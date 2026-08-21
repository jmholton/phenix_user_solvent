#!/usr/bin/env phenix.python
"""
combine_solvent: fold a squish solvent INCREMENT into the current bulk-solvent map to
produce ONE updated user Fpart for patched phenix.refine.

Patched phenix supports exactly one user Fpart, so (unlike refmac's scpart, which keeps
each increment as a separately-scaled partial) the old solvent and the new increment must
be pre-combined into a single Fpart/PHIpart column. Everything is in electron units.

  old solvent (electron units)     :  K_MASK * FMASK          (from the phenix f_model)
  residual the increment explains  :  fofc = FOBS - FMODEL    (model phases; electron units)
  squish increment ("the divot")   :  Fpart,PHIpart in refme_Fpart.mtz (electron units;
                                       generate with fofc_prescale=1 and the loB>=0 guard)

The divot is re-scaled by a single k fit to fofc over the WORK set (the phenix-side analog
of refmac's scpart):   k = Re(<fofc, incr*>) / <|incr|^2>.

The increment is at the DATA (FOBS) scale; FMASK is at phenix's PRE-overall scale, so the
increment is divided by the exact per-reflection K_ISO*K_ANISO to land on FMASK's scale
(then phenix's fresh k_mask starts ~1):

  new_FMASK = K_MASK*FMASK  +  k * increment / (K_ISO*K_ANISO)      (complex, per reflection)

Output: an MTZ with Fpart,PHIpart (|new_FMASK|, phase) to hand phenix as
  refinement.input.bulk_solvent_map.file_name=... amplitudes_label=Fpart phases_label=PHIpart

Usage:
  phenix.python combine_solvent.py  fmodel.mtz  refme_Fpart.mtz  [out.mtz]
     [k=auto|<float>] [fobs=FOBS fmodel=FMODEL fmask=FMASK kmask=K_MASK
      kiso=K_ISOTROPIC kaniso=K_ANISOTROPIC free=R_FREE_FLAGS incr=Fpart]
"""

import sys
import iotbx.mtz
from cctbx import miller as cctbx_miller
from cctbx.array_family import flex


def run(fmodel_mtz, incr_mtz, out_mtz=None, k='auto',
        fobs_l='FOBS', fmodel_l='FMODEL', fmask_l='FMASK', kmask_l='K_MASK',
        kiso_l='K_ISOTROPIC', kaniso_l='K_ANISOTROPIC', free_l='R_FREE_FLAGS',
        incr_l='Fpart'):

    def load(mtz):
        o = iotbx.mtz.object(mtz)
        by = {}
        for ma in o.as_miller_arrays():
            for lab in ma.info().labels:
                by.setdefault(lab, ma)
        return by

    fm = load(fmodel_mtz)
    inc = load(incr_mtz)

    def col(by, label, src):
        if label not in by:
            raise RuntimeError('column %r not found in %s (have: %s)'
                               % (label, src, sorted(by)))
        return by[label].map_to_asu()

    fo = col(fm, fobs_l, fmodel_mtz)
    ref = fo.indices()
    n = ref.size()

    def a_real(by, label, src):
        m = col(by, label, src)
        out = flex.double(n, 0.0)
        for i, j in cctbx_miller.match_indices(ref, m.indices()).pairs():
            out[i] = m.data()[j]
        return out

    def a_cplx(by, label, src):
        m = col(by, label, src)
        out = flex.complex_double(n, 0j)
        for i, j in cctbx_miller.match_indices(ref, m.indices()).pairs():
            out[i] = m.data()[j]
        return out

    def a_int(by, label, src):
        m = col(by, label, src)
        out = flex.int(n, -1)
        for i, j in cctbx_miller.match_indices(ref, m.indices()).pairs():
            out[i] = int(round(float(m.data()[j])))
        return out

    fo_amp = fo.data()
    fmodel = a_cplx(fm, fmodel_l, fmodel_mtz)     # complex, model phases
    fmask  = a_cplx(fm, fmask_l, fmodel_mtz)      # complex solvent mask
    kmask  = a_real(fm, kmask_l, fmodel_mtz)
    kiso   = a_real(fm, kiso_l, fmodel_mtz)
    kaniso = a_real(fm, kaniso_l, fmodel_mtz)
    free   = a_int(fm, free_l, fmodel_mtz)        # phenix binary: free=1, work=0
    incr   = a_cplx(inc, incr_l, incr_mtz)        # squish increment (electron units)

    # fofc = (FOBS - |FMODEL|) with model phase  (electron units)
    fofc = flex.complex_double(n, 0j)
    for i in range(n):
        fm_i = fmodel[i]
        amp = abs(fm_i)
        if amp > 1e-9:
            ph = fm_i / amp                       # unit phasor e^{i phi_model}
        else:
            ph = 1+0j
        fofc[i] = (fo_amp[i] - amp) * ph

    # fit single divot scale k over the WORK set:  k = Re<fofc,incr*> / <|incr|^2>
    if str(k) == 'auto':
        num = 0.0
        den = 0.0
        for i in range(n):
            if free[i] == 0:                      # work only
                num += (fofc[i] * incr[i].conjugate()).real
                den += abs(incr[i]) ** 2
        k_fit = (num / den) if den > 0 else 0.0
    else:
        k_fit = float(k)

    # new_FMASK = K_MASK*FMASK + k*incr/(K_ISO*K_ANISO)   (pre-overall scale)
    new_fmask = flex.complex_double(n, 0j)
    for i in range(n):
        kk = kiso[i] * kaniso[i]
        if kk < 1e-6:
            kk = 1.0
        new_fmask[i] = kmask[i] * fmask[i] + k_fit * incr[i] / kk

    fpart_amp = flex.abs(new_fmask)
    # write Fpart/PHIpart via a complex miller array (root 'Fpart' -> Fpart, PHIFpart);
    # rename the phase to PHIpart so it matches the established example convention.
    ds = fo.set().array(data=new_fmask).as_mtz_dataset(column_root_label='Fpart')
    mtz = ds.mtz_object()
    for col_obj in mtz.columns():
        if col_obj.label() == 'PHIFpart':
            col_obj.set_label('PHIpart')

    if out_mtz is None:
        out_mtz = 'solvent_Fpart_updated.mtz'
    mtz.write(out_mtz)

    # report
    old_amp = flex.double([kmask[i] * abs(fmask[i]) for i in range(n)])
    print('combine_solvent: %s + %s -> %s' % (fmodel_mtz, incr_mtz, out_mtz))
    print('  reflections: %d   divot scale k = %.4f (%s)'
          % (n, k_fit, 'fit to fofc, work set' if str(k) == 'auto' else 'user'))
    print('  mean |K_MASK*FMASK| (old solvent) = %.2f' % flex.mean(old_amp))
    print('  mean |Fpart|        (new solvent) = %.2f' % flex.mean(fpart_amp))
    print('  columns: Fpart PHIpart')
    print('  phenix: refinement.input.bulk_solvent_map.file_name=%s '
          'amplitudes_label=Fpart phases_label=PHIpart' % out_mtz)


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print(__doc__); sys.exit(1)
    kw = {}
    pos = []
    for a in sys.argv[1:]:
        if '=' in a and not a.endswith('.mtz'):
            key, val = a.split('=', 1)
            m = {'k': 'k', 'fobs': 'fobs_l', 'fmodel': 'fmodel_l', 'fmask': 'fmask_l',
                 'kmask': 'kmask_l', 'kiso': 'kiso_l', 'kaniso': 'kaniso_l',
                 'free': 'free_l', 'incr': 'incr_l'}.get(key)
            if m:
                kw[m] = val
            else:
                print('unknown option', a); sys.exit(1)
        else:
            pos.append(a)
    run(pos[0], pos[1], pos[2] if len(pos) > 2 else None, **kw)
