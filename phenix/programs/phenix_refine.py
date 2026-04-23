"""
  phenix.refine
"""
from __future__ import division, print_function

from iotbx.cli_parser import CCTBXParser, run_program
from phenix.program_template import ProgramTemplate
import mmtbx.utils
import iotbx.extract_xtal_data
import phenix.refinement
import phenix.refinement.driver
import sys
import iotbx.phil # XXX
from phenix.refinement.fsr import io
from libtbx.utils import multi_out
from libtbx.utils import Sorry
from six.moves import cStringIO as StringIO
import libtbx.load_env, os
from copy import deepcopy
from libtbx.utils import null_out

import phenix.refinement.fmodels


def _inject_user_bulk_solvent(fmodel, params, log):
  """Read user-supplied bulk solvent structure factors from an MTZ file and
  inject them into fmodel, replacing the internally computed flat mask."""
  if params.file_name is None:
    return
  import iotbx.mtz
  import math
  mtz_obj = iotbx.mtz.object(file_name=params.file_name)
  ma_list = mtz_obj.as_miller_arrays()
  amp_array = None
  phi_array = None
  for ma in ma_list:
    if params.amplitudes_label in ma.info().labels:
      amp_array = ma
    if params.phases_label in ma.info().labels:
      phi_array = ma
  if amp_array is None:
    raise Sorry(
      "Column '%s' not found in %s" % (params.amplitudes_label, params.file_name))
  if phi_array is None:
    raise Sorry(
      "Column '%s' not found in %s" % (params.phases_label, params.file_name))
  f_complex = amp_array.phase_transfer(phase_source=phi_array, deg=True)
  f_complex = f_complex.map_to_asu()
  # Re-express on exactly f_obs's index set (same indices, same canonical order).
  # Uses match_indices so the result satisfies f_calc.indices().all_eq(fm.indices()).
  # Any f_obs reflections absent from the user map get zero contribution.
  from cctbx import miller as cctbx_miller
  from cctbx.array_family import flex
  f_obs = fmodel.f_obs()
  zero_data = flex.complex_double(f_obs.indices().size(), 0+0j)
  match = cctbx_miller.match_indices(f_obs.indices(), f_complex.indices())
  for i_fobs, i_user in match.pairs():
    zero_data[i_fobs] = f_complex.data()[i_user]
  f_mask = f_obs.array(data=zero_data)
  fmodel.set_user_f_masks([f_mask])
  print("Using user-supplied bulk solvent map from: %s" % params.file_name,
    file=log)
  print("  amplitude column: %s   phase column: %s" % (
    params.amplitudes_label, params.phases_label), file=log)


# logic for handling the default behavior for
#
#   phenix.refine <model file> <data file>
#
# without needing to provide DataManager parameters, but enable the
# setting of the data type (e.g. x_ray, neutron, electron) with the
# refinement.main.scattering_table parameter

def customize_and_process_single(parser):
  """
  Customizes and processes the input for single (non-joint) refinement.

  Applies default map settings, handles reference models, and sets data types
  based on scattering table if only one model and data file are provided.

  Args:
    parser (CCTBXParser): The command-line argument parser.
  """
  working_extract = parser.working_phil.extract()
  dm = parser.data_manager
  # Use pre-defined map settings; user can overwrite.
  emaps_params = working_extract.refinement.electron_density_maps
  if emaps_params.apply_default_maps is not False:
    if len(emaps_params.map_coefficients) == 0 \
      and len(emaps_params.map) == 0:
      emaps_params.apply_default_maps = True
      print('''\
No user-defined map coefficients or files defined; will use default map
outputs instead.''', file=parser.logger)
    elif emaps_params.apply_default_maps :
      print('''Will apply parameters for default map types.''', file=parser.logger)
    if emaps_params.apply_default_maps:
      maps_par = os.path.abspath(
        os.path.join(__file__, '..', '..', 'refinement', 'customizations', 'maps.params'))
      assert(os.path.isfile(maps_par))
      maps_par_phil = iotbx.phil.parse(file_name=maps_par)
      parser.working_phil = parser.working_phil.fetch(sources=[maps_par_phil])
  # Fetch DEN params
  if("den" in working_extract.refinement.refine.strategy):
    den_par = os.path.abspath(
      os.path.join(__file__, '..', '..', 'refinement', 'customizations', 'den.params'))
    assert(os.path.isfile(den_par))
    den_par_phil = iotbx.phil.parse(file_name=den_par)
    parser.working_phil = parser.working_phil.fetch(sources=[den_par_phil])
  # Deal with reference model
  for fn in working_extract.refinement.reference_model.file:
    dm.set_model_type(filename=fn, model_type=['reference'])
  # check that there are exactly one model and one data file
  if len(dm.get_model_names()) == 1 and len(dm.get_miller_array_names()) == 1:
    data_type = working_extract.refinement.main.scattering_table
    if data_type not in ['electron', 'neutron']:
      data_type = 'x_ray'
    if len(dm.get_model_type()) == 1:
      dm.set_model_type(model_type=[data_type])
    for label in dm.get_miller_array_labels():
      dm.set_miller_array_type(label=label, array_type=[data_type])

def is_joint_refinement(parser):
  """
  Checks if the input criteria for joint refinement are met.

  Criteria: one X-ray model, one neutron model, X-ray data, and neutron data.

  Args:
    parser (CCTBXParser): The command-line argument parser.

  Returns:
    bool: True if joint refinement criteria are met, False otherwise.
  """
  # check for joint refinement criteria
  #   1 xray model and 1 neutron model
  #   xray and neutron data available
  dm = parser.data_manager
  if len(dm.get_model_names()) == 2 or len(dm.get_model_type()) == 2:
    # check models first (1 xray, 1 neutron)
    model_types = set()
    for filename in dm.get_model_names():
      for model_type in dm.get_model_type(filename):
        model_types.add(model_type)
    # then check data
    if 'x_ray' in model_types and 'neutron' in model_types:
      xray_available = False
      neutron_available = False
      for filename in dm.get_miller_array_names():
        data_types = dm.get_miller_array_types(filename).values()
        # convert list of lists into a set
        # [['x_ray'], ['x_ray', 'neutron']] becomes {'x_ray', 'neutron'}
        data_types = {dt for dt_list in data_types for dt in dt_list}
        xray_available = xray_available or 'x_ray' in data_types
        neutron_available = neutron_available or 'neutron' in data_types
        if xray_available and neutron_available:
          return True
  return False

def customize_and_process_joint(parser):
  """
  Customizes and processes the input for joint X-ray/neutron refinement.

  Applies map and refinement customizations specific to joint refinement,
  sets a new master PHIL, and handles single parameters if provided.

  Args:
    parser (CCTBXParser): The command-line argument parser.
  """
  # store original changed parameters
  single_extract = None
  if len(parser.unused_phil) == 0:  # joint data provided, but only common parameters provided
    single_extract = parser.master_phil.fetch(parser.working_phil).extract()

  # Apply map customizations
  maps_par = os.path.abspath(
    os.path.join(__file__, '..', '..', 'refinement', 'customizations',
    'maps_joint.params'))
  assert(os.path.isfile(maps_par))
  maps_par_phil = iotbx.phil.parse(file_name=maps_par)
  # Apply general refinement customizations
  ref_par = os.path.abspath(
    os.path.join(__file__, '..', '..', 'refinement', 'customizations',
    'ref_joint.params'))
  assert(os.path.isfile(maps_par))
  ref_par_phil = iotbx.phil.parse(file_name=ref_par)
  #
  joint_master_phil_str = Program.create_joint_master_phil_str(Program)
  # Set new master phil
  parser.master_phil = iotbx.phil.parse(
    joint_master_phil_str, process_includes=True)
  required_output_phil = iotbx.phil.parse(ProgramTemplate.output_phil_str)
  parser.master_phil.adopt_scope(required_output_phil)
  # Fetch all customizations
  parser.master_phil = parser.master_phil.fetch(
    sources=[maps_par_phil, ref_par_phil])

  # reset parser for joint refinement PHIL
  parser.unused_phil_raises_sorry = True
  parser.unused_phil = []

  # copy original changed parameters to new joint parameters if needed
  # retain neutron scattering table
  if single_extract is not None:
    joint_extract = parser.master_phil.fetch().extract()
    joint_extract.xray.refinement = deepcopy(single_extract.refinement)
    joint_extract.neutron.refinement = deepcopy(single_extract.refinement)
    joint_extract.neutron.refinement.main.scattering_table = 'neutron'
    joint_extract.xray.output = single_extract.output
    joint_extract.neutron.output = single_extract.output
    joint_extract.output = single_extract.output
    parser.working_phil = parser.master_phil.format(python_object=joint_extract)
  # or reparse joint input parameters
  else:
    parser.working_phil = None
    parser.process_phil(parser.namespace.phil)

def custom_process_arguments(parser):
  """
  Processes command-line arguments, handling both single and joint refinement.

  First processes arguments for single refinement. Then, if joint refinement
  criteria are met, switches to joint refinement mode and re-processes.

  Args:
    parser (CCTBXParser): The command-line argument parser.
  """
  # Single (not joint XN) refinement
  customize_and_process_single(parser)
  single_unused_phil = parser.unused_phil
  # Joint XN refinement: modify phil if joint refinement criteria satisfied
  if(is_joint_refinement(parser)):
    print('Switching to joint x-ray/neutron refinement mode', file=parser.logger)
    print('-'*parser.text_width, file=parser.logger)
    print()
    customize_and_process_joint(parser)
  else:
    parser.unused_phil_raises_sorry = True
    parser.unused_phil = single_unused_phil
  parser.raise_Sorry_for_unused_phil()


# programs directory
class Program(ProgramTemplate):
  """
  Program for phenix.refine.
  """

  datatypes=['model', 'phil', 'miller_array', 'restraint']

  description = '''
  phenix.refine
'''

  master_phil_str = """
%s
"""%phenix.refinement.master_params().as_str(attributes_level=9)

  # xray parameters for joint refinement
  xray_master_phil_str = """
  xray {
    %s
  }
  """ % master_phil_str

  # neutron parameters for joint refinement
  neutron_master_phil_str = """
  neutron {
    %s
  }
  """ % master_phil_str

  # common parameters for joint refinement
  joint_master_phil_str = """
  joint {
    dummy_parameter = 0
      .type = int
      .help = common parameter for joint refinements
  }
  """

  show_data_manager_scope_by_default = True

  @staticmethod
  def create_joint_master_phil_str(Program):
    '''
    Class function to construct the joint master PHIL from class variables
    '''
    # xray
    xray_master_phil = iotbx.phil.parse(
      Program.xray_master_phil_str, process_includes=True)
    xray_output_phil = iotbx.phil.parse(
      'xray { %s }' % ProgramTemplate.output_phil_str)
    xray_master_phil.adopt_scope(xray_output_phil)

    # neutron
    neutron_master_phil = iotbx.phil.parse(
      Program.neutron_master_phil_str, process_includes=True)
    neutron_output_phil = iotbx.phil.parse(
      'neutron { %s }' % ProgramTemplate.output_phil_str)
    neutron_master_phil.adopt_scope(neutron_output_phil)

    # complete set of joint parameters
    joint_master_phil_str = xray_master_phil.as_str(attributes_level=9) \
                          + neutron_master_phil.as_str(attributes_level=9) \
                          + Program.joint_master_phil_str

    return joint_master_phil_str

  def custom_init(self):
    """
    Performs custom initialization tasks.

    Checks for PHENIX_OVERWRITE_ALL and joint refinement mode.
    Handles model file names and output prefixes.
    Sets up logging.
    """
    # check for PHENIX_OVERWRITE_ALL
    super(Program, self).custom_init()

    # check if joint refinement mode is set
    self.joint_xn = False
    self.single_master_phil = iotbx.phil.parse(self.master_phil_str, process_includes=True)
    self.single_master_phil.adopt_scope(iotbx.phil.parse(
      ProgramTemplate.output_phil_str, process_includes=True))
    if hasattr(self.params, 'neutron'):
      self.joint_xn = True

      # apply PHENIX_OVERWRITE_ALL
      self.params.xray.output.overwrite = self.params.output.overwrite
      self.params.neutron.output.overwrite = self.params.output.overwrite

    params = self.params
    if self.joint_xn:  # why is this needed for joint refinement?
      params = self.params.xray

    # XXX backwards compatibility, remove later! XXX
    model_file_names = []
    for model_type in ["x_ray", "neutron", "electron"]:
      for fn in self.data_manager.get_model_names(model_type=model_type):
        model_file_names.append(fn)
    # XXX backwards compatibility, remove later! XXX
    if(self.params.output.serial==0): self.params.output.serial=1
    # XXX
    #
    if(len(model_file_names)==0):
      raise Sorry("An atomic model is required.")
    #
    prefix_from_file = io.get_prefix(
      file_name  = model_file_names[0],
      add_suffix = True)#self.params.refinement.output.prefix is not None)
    io.determine_output_prefix_and_serial(
      params           = self.params,  # if joint_xn, the main output scope needs to be updated
      prefix_from_file = prefix_from_file)
#XXX    # XXX backwards compatibility, remove later! XXX
#XXX    self.params.output.prefix = params.refinement.output.prefix
#XXX    self.params.output.serial = params.refinement.output.serial
#XXX    if self.joint_xn:
#XXX      self.params.xray.output.prefix = params.refinement.output.prefix
#XXX      self.params.xray.output.serial = params.refinement.output.serial
#XXX      self.params.neutron.refinement.output.prefix = params.refinement.output.prefix
#XXX      self.params.neutron.refinement.output.serial = params.refinement.output.serial
#XXX      self.params.neutron.output.prefix = params.refinement.output.prefix
#XXX      self.params.neutron.output.serial = params.refinement.output.serial

    # XXX
    if self.joint_xn:
      self.params.xray.output.prefix = self.params.output.prefix
      self.params.xray.output.serial = self.params.output.serial
      self.params.neutron.output.prefix = self.params.output.prefix
      self.params.neutron.output.serial = self.params.output.serial
      osf = self.params.xray.output.serial_format
      prefix = self.params.xray.output.prefix
      serial = self.params.xray.output.serial
      overwrite = self.params.xray.output.overwrite
    else:
      osf = self.params.output.serial_format
      prefix = self.params.output.prefix
      serial = self.params.output.serial
      overwrite = self.params.output.overwrite
    file_name = io.output_file_name_manager(
      prefix        = prefix,
      serial        = serial,
      overwrite     = overwrite,
      serial_format = osf
        ).get_file_name(ext=".log")
    self.log = self.logger
    if(isinstance(self.logger, multi_out)):
      log_fname = file_name
      log_fname = self.data_manager._update_default_output_filename(log_fname)
      log_f = open(log_fname, 'w')
      self.logger.replace_stringio(
        old_label       = "parser_log",
        new_label       = "log_file",
        new_file_object = log_f)
      self.log = self.logger
      # NEED FIX: Never overwrite stderr without restoring it.
      # Upstream callers will be srewed.
      sys.stderr = self.log

  def validate(self):
    """
    Validates the program parameters.

    Currently, it checks if any models are provided.
    """
    # get models
    params = self.params
    if self.joint_xn:
      params = self.params.xray
    #
    def _get_models(model_type):
      result = []
      for filename in self.data_manager.get_model_names(model_type=model_type):
        if(filename in params.refinement.reference_model.file): continue
        result.append(self.data_manager.get_model(filename, model_type))
      return result
    #
    self.xray_models     = _get_models(model_type='x_ray')
    self.neutron_models  = _get_models(model_type='neutron')
    self.electron_models = _get_models(model_type='electron')
    #
    all_models = self.xray_models+self.neutron_models+self.electron_models
    self.model = None
    if(len(all_models)==1):
      self.model = all_models[0]

    # check models
    if(len(all_models)==0):
      raise Sorry("An atomic model is required.")
    if(not(len(self.xray_models) in [0,1] and
           len(self.neutron_models) in [0,1])):
      raise Sorry("Wrong number of models of each type supplied.")

    self.fmodel_xray, self.fmodel_neutron, self.fmodel = None, None, None
    if(not self.joint_xn):
      #self.data_manager.update_all_defaults(self.params.refinement.main.scattering_table)

      rp = self.params.refinement
      sf_accuracy_params = rp.structure_factors_and_gradients_accuracy
      self.fmodel = self.data_manager.get_fmodel(
        crystal_symmetry           = rp.crystal_symmetry,
        #experimental_phases_params = rp.input.experimental_phases,
        mask_params                = rp.mask,
        sf_accuracy_params         = sf_accuracy_params,
        scattering_table           = rp.main.scattering_table,
        free_r_flags_scope = 'miller_array.labels.name')
      self.fmodel.update(mask_params = rp.mask)
      _inject_user_bulk_solvent(
        fmodel = self.fmodel,
        params = rp.input.bulk_solvent_map,
        log    = self.log)
    else:
      rp = self.params.xray.refinement
      sf_accuracy_params = rp.structure_factors_and_gradients_accuracy
      self.fmodel_xray = self.data_manager.get_fmodel(
        crystal_symmetry           = rp.crystal_symmetry,
        sf_accuracy_params         = sf_accuracy_params,
        #experimental_phases_params = rp.input.experimental_phases,
        mask_params                = rp.mask,
        scattering_table           = rp.main.scattering_table,
        free_r_flags_scope = 'xray_data')
      rp = self.params.neutron.refinement
      sf_accuracy_params = rp.structure_factors_and_gradients_accuracy
      self.fmodel_neutron = self.data_manager.get_fmodel(
        crystal_symmetry           = rp.crystal_symmetry,
        #experimental_phases_params = rp.input.experimental_phases,
        mask_params                = rp.mask,
        sf_accuracy_params         = sf_accuracy_params,
        scattering_table           = rp.main.scattering_table,
        free_r_flags_scope = 'neutron_data')

  def _run_xn_joint(self):
    params_x = self.params.xray
    params_n = self.params.neutron
    #
    params_x.output.prefix            += "_xray"
    params_n.output.prefix            += "_neutron"
    #
    return phenix.refinement.driver.run0_joint(
      master_params = self.single_master_phil,
      params_x        = params_x,
      params_n        = params_n,
      dm_working_phil = self.data_manager.export_phil_scope(),
      model_x         = self.xray_models[0],
      model_n         = self.neutron_models[0],
      log             = self.log,
      fmodel_x        = self.fmodel_xray,
      fmodel_n        = self.fmodel_neutron,
      overwrite       = self.params.output.overwrite,
      dry_run         = self.params.xray.refinement.dry_run)

  def _run_single(self):
    return phenix.refinement.driver.run0(
      master_params   = self.master_phil,
      params          = self.params,
      dm_working_phil = self.data_manager.export_phil_scope(),
      model           = self.model,
      log             = self.log,
      fmodel          = self.fmodel,
      overwrite       = self.params.output.overwrite,
      dry_run         = self.params.refinement.dry_run)

  def run(self):
    if(self.joint_xn): self.refined = self._run_xn_joint()
    else:              self.refined = self._run_single()
    if 0:
      print(self.refined.fmodel.show_call_stats_sorted_and_rounded())
      print()
      print(self.refined.model.show_call_stats_sorted_and_rounded())
      print()
      print(self.refined.macro_cycle.show_call_stats_sorted_and_rounded())

  def get_results(self):
    return self.refined


# subclass for CCTBXParser to enable displaying of joint refinement parameters
class PhenixRefineParser(CCTBXParser):

  def add_default_options(self):
    super(PhenixRefineParser, self).add_default_options()

    # --show-joint-defaults by itself is set to 0
    # --show-joint-defaults=n sets it to n and it can only be {0, 1, 2, 3}
    self.add_argument(
      '--show-joint-defaults', '--show_joint_defaults',
      nargs='?', const=0, type=int, choices=range(0,4),
      help='show default parameters with expert level (default=0)'
    )

    # # move flag to below --show-defaults
    # new_flag = self.__dict__['_actions'].pop(-1)
    # for i, _action in enumerate(self.__dict__['_actions']):
    #   if _action.dest == 'show_defaults':
    #     self.__dict__['_actions'].insert(i+1, new_flag)
    #     break

    # new_flag_0 = self.__dict__['_actions'].pop(-1)
    # new_flag_1 = self.__dict__['_actions'].pop(-1)
    # new_dict = {}
    # for _action in self.__dict__['_option_string_actions'].keys():
    #   if _action == '--show_defaults':
    #     new_dict[_action] = self.__dict__['_option_string_actions'][_action]
    #     new_dict['--show-joint-defaults'] = new_flag_0
    #     new_dict['--show_joint_defaults'] = new_flag_1
    #   else:
    #     new_dict[_action] = self.__dict__['_option_string_actions'][_action]
    # self.__dict__['_option_string_actions'] = new_dict

    # for _action in self.__dict__['_option_string_actions'].keys():
    #   print(_action)

    # print(dir(self))

  def parse_args(self, args):
    '''
    The first few steps from CCTBXParser is copied to handle some basic options
    '''
    # default behavior with no arguments
    if len(args) == 0:
      self.print_help()
      self.exit()

    # parse arguments
    if sys.version_info >= (3, 7):
      # https://docs.python.org/3/library/argparse.html#argparse.ArgumentParser.parse_intermixed_args
      # https://bugs.python.org/issue9338
      # https://bugs.python.org/issue15112
      self.namespace = super(CCTBXParser, self).parse_intermixed_args(args)
    else:
      self.namespace = super(CCTBXParser, self).parse_args(args)

    # process command-line options
    if self.namespace.attributes_level is not None:
      if self.namespace.show_defaults is None and self.namespace.show_joint_defaults is None:
        self.error('--attributes-level requires --show-defaults or --show-joint-defaults to be set')
    if self.namespace.show_joint_defaults is not None:
      if self.namespace.attributes_level is None:
        self.namespace.attributes_level = 0
      if self.program_class.show_data_manager_scope_by_default:
        self.data_manager.master_phil.show(
          expert_level=self.namespace.show_defaults,
          attributes_level=self.namespace.attributes_level,
          out=self.logger)
      self.master_phil = iotbx.phil.parse(Program.create_joint_master_phil_str(Program),
                                          process_includes=True)
      self.master_phil.show(expert_level=self.namespace.show_joint_defaults,
                            attributes_level=self.namespace.attributes_level,
                            out=self.logger)
      self.exit()

    self.namespace = super(PhenixRefineParser, self).parse_args(args)

    return self.namespace

def run_phenix_refine(args=None, json=False, logger=None, hide_parsing_output=False):
  '''
  Convenience function for running phenix.refine. Due to the custom behavior for
  handling joint x-ray/neutron refinements, there is a custom parser and extra
  parsing code.

  Parameters
  ----------
  args : list, optional
      list of command-line arguments
  json : bool, optional
      switch for JSON output, by default False
  logger : multi_out, optional
      custom logger, by default None
  hide_parsing_output : bool, optional
      switch for command-line parsing output, by default False
  '''
  from iotbx.cli_parser import run_program
  run_program(program_class=Program,
              parser_class=PhenixRefineParser,
              custom_process_arguments=custom_process_arguments,
              unused_phil_raises_sorry=False,
              args=args,
              json=json,
              logger=logger,
              hide_parsing_output=hide_parsing_output)

if (__name__ == '__main__'):
  run_phenix_refine()
