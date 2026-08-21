#!/usr/bin/env phenix.python
"""
phenix2squish: adapt a phenix.refine _f_model.mtz so squish_solvent_runme.com can
consume it (see the phenix-input mode added to that script).

squish was written for refmac output. The only real differences for a phenix
_f_model.mtz are (a) the free-flag convention and (b) column names. This tool writes
a clean MTZ that keeps phenix's own column names (squish's phenix mode knows them) and
fixes the free flag:

  FOBS, SIGFOBS          passed through (squish FP, SIGFP)
  FreeR_flag  = 1 - R_FREE_FLAGS   phenix binary free=1  ->  CCP4/refmac free=0
                                   (so `fft ... FREE=FreeR_flag` excludes the test set)
  FMODEL, PHIFMODEL      passed through -- the total scaled model
                         K_ISO*K_ANISO*(FCALC + K_MASK*FMASK); squish subtracts this to
                         form Fo-Fc and adds it back to rebuild the solvent-removed Fobs
  FMASK, PHIFMASK        passed through -- the (user) bulk solvent mask, unscaled
  K_MASK                 passed through -- per-bin solvent scale
  Fpart, PHIFpart = K_MASK*FMASK, PHIFMASK   the CURRENT physical solvent as a partial
                         (needed only for the one-Fpart hand-back: new solvent =
                          this + squish's increment, combined into a single column,
                          because patched phenix supports exactly one user Fpart)

Everything is re-expressed on FOBS's index set (match_indices), so all columns share
one HKL list. Fobs is NOT wilsonified here (squish needs the true residual Fo-FMODEL);
use wilsonify.py instead when the goal is refmac scale-transfer.

Usage:
  phenix.python phenix2squish.py  fmodel.mtz  [out.mtz]
      [fobs=FOBS free=R_FREE_FLAGS fmodel=FMODEL fmask=FMASK kmask=K_MASK]
Then:
  squish_solvent_runme.com  model.pdb  out.mtz     # squish auto-detects phenix columns
"""

import sys
import iotbx.mtz
from cctbx import miller as cctbx_miller
from cctbx.array_family import flex


def run(fmodel_mtz, out_mtz=None,
        fobs_l='FOBS', free_l='R_FREE_FLAGS',
        fmodel_l='FMODEL', fmask_l='FMASK', kmask_l='K_MASK'):

    obj = iotbx.mtz.object(fmodel_mtz)
    by_label = {}
    for ma in obj.as_miller_arrays():
        for lab in ma.info().labels:
            by_label.setdefault(lab, ma)

    def col(label, required=True):
        if label not in by_label:
            if required:
                raise RuntimeError('column %r not found in %s (have: %s)'
                                   % (label, fmodel_mtz, sorted(by_label)))
            return None
        return by_label[label].map_to_asu()

    fo = col(fobs_l)
    if fo.sigmas() is None:
        raise RuntimeError('%s has no sigmas (need an F,Q column)' % fobs_l)
    ref = fo.indices()
    n = ref.size()

    def aligned_real(m):
        out = flex.double(n, 0.0)
        for i, j in cctbx_miller.match_indices(ref, m.indices()).pairs():
            out[i] = m.data()[j]
        return out

    def aligned_complex(m):
        out = flex.complex_double(n, 0j)
        for i, j in cctbx_miller.match_indices(ref, m.indices()).pairs():
            out[i] = m.data()[j]
        return out

    def aligned_int(m):
        out = flex.int(n, 0)
        for i, j in cctbx_miller.match_indices(ref, m.indices()).pairs():
            out[i] = int(round(float(m.data()[j])))
        return out

    # --- total scaled model FMODEL (complex) ---
    fmodel = aligned_complex(col(fmodel_l))
    # --- bulk solvent mask FMASK (complex) and per-bin scale K_MASK (real) ---
    fmask = aligned_complex(col(fmask_l))
    kmask = aligned_real(col(kmask_l))
    # --- current physical solvent partial: Fpart = K_MASK * FMASK (complex) ---
    fpart = flex.complex_double([kmask[i] * fmask[i] for i in range(n)])

    # --- free flag: flip phenix binary (free=1) to CCP4 convention (free=0) ---
    free = aligned_int(col(free_l))
    freer = flex.int([1 - free[i] for i in range(n)])   # 1->0 (free), 0->1 (work)
    n_free = (freer == 0).count(True)

    # --- assemble the output MTZ (all arrays share FOBS's index set) ---
    ds = fo.as_mtz_dataset(column_root_label='FOBS', column_types='FQ')
    ds.add_miller_array(fo.set().array(data=freer),
                        column_root_label='FreeR_flag')
    ds.add_miller_array(fo.set().array(data=fmodel),
                        column_root_label='FMODEL')      # -> FMODEL, PHIFMODEL
    ds.add_miller_array(fo.set().array(data=fmask),
                        column_root_label='FMASK')        # -> FMASK, PHIFMASK
    ds.add_miller_array(fo.set().array(data=kmask),
                        column_root_label='K_MASK')
    ds.add_miller_array(fo.set().array(data=fpart),
                        column_root_label='Fpart')        # -> Fpart, PHIFpart

    if out_mtz is None:
        base = fmodel_mtz[:-4] if fmodel_mtz.endswith('.mtz') else fmodel_mtz
        out_mtz = base + '_squishin.mtz'
    ds.mtz_object().write(out_mtz)

    print('phenix2squish %s -> %s' % (fmodel_mtz, out_mtz))
    print('  %d reflections, %d free (FreeR_flag==0, CCP4 convention)' % (n, n_free))
    print('  columns: FOBS SIGFOBS FreeR_flag FMODEL PHIFMODEL '
          'FMASK PHIFMASK K_MASK Fpart PHIFpart')
    print('  squish:  squish_solvent_runme.com model.pdb %s' % out_mtz)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(1)
    kw = {}
    pos = []
    for a in sys.argv[1:]:
        if a.startswith('fobs='):    kw['fobs_l'] = a.split('=', 1)[1]
        elif a.startswith('free='):  kw['free_l'] = a.split('=', 1)[1]
        elif a.startswith('fmodel='):kw['fmodel_l'] = a.split('=', 1)[1]
        elif a.startswith('fmask='): kw['fmask_l'] = a.split('=', 1)[1]
        elif a.startswith('kmask='): kw['kmask_l'] = a.split('=', 1)[1]
        else: pos.append(a)
    run(pos[0], pos[1] if len(pos) > 1 else None, **kw)
