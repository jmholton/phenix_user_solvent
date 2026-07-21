#!/usr/bin/env phenix.python
"""
Wilsonify: transfer phenix's per-resolution-bin scaling into a refmac-ready MTZ.

Divides Fobs by phenix's EXACT per-reflection K_ISOTROPIC * K_ANISOTROPIC (its full
scaling envelope), so that refmac's flat overall scaling reproduces phenix's per-bin
scaling. Also packages the physical bulk solvent as Fpart = K_MASK * FMASK.

Result: refmac (SOLVENT NO, SHANNON 4-6, no SCPART) reproduces phenix's R-free to the
digit. Called "Wilsonification" because dividing out K_ISO straightens the Wilson plot
of Fobs, leaving a flat envelope refmac's single k*exp(-B s^2/4) can match.

Unlike Fpart = FMODEL - FCALC_R (a fixed-geometry validator that packs a non-physical
correction into the partial), Wilsonification keeps the solvent PHYSICAL and transfers
the scaling through the DATA -- so it also works as a refinement target (NCYC>0).

Input columns are read from a phenix.refine _f_model.mtz
(output.export_final_f_model=True): FOBS/SIGFOBS, R_FREE_FLAGS, K_ISOTROPIC,
K_ANISOTROPIC, K_MASK, and FMASK.

Usage:
    phenix.python wilsonify.py  phenix_f_model.mtz  [out.mtz]
        [fobs=FOBS] [free=R_FREE_FLAGS] [kiso=K_ISOTROPIC] [kaniso=K_ANISOTROPIC]
        [kmask=K_MASK] [fmask=FMASK] [free_in=1] [free_out=0] [no_solvent]

Then run refmac (feed it the CIF, not a lossy PDB, for a multi-conformer ensemble):
    refmac5 XYZIN model.cif HKLIN out.mtz HKLOUT refined.mtz << eof
    LABIN FP=FP SIGFP=SIGFP FREE=FreeR_flag FPART1=Fpart PHIP1=PHIFpart
    SOLVENT NO
    SHANNON 6
    NCYC 0            # or >0 to refine
    eof

NOTE (refinement): the K_ISO/K_ANISO are fixed at the input model's geometry. For
NCYC>0 they slowly go stale as coordinates move; refresh by re-exporting a phenix
f_model and re-Wilsonifying every few macrocycles.
"""

import sys, math
import iotbx.mtz
from cctbx import miller as cctbx_miller
from cctbx.array_family import flex


def run(fmodel_mtz, out_mtz=None,
        fobs_l='FOBS', free_l='R_FREE_FLAGS',
        kiso_l='K_ISOTROPIC', kaniso_l='K_ANISOTROPIC',
        kmask_l='K_MASK', fmask_l='FMASK',
        free_in=1, free_out=0, want_solvent=True):

    obj = iotbx.mtz.object(fmodel_mtz)
    by_label = {}
    for ma in obj.as_miller_arrays():
        for lab in ma.info().labels:
            by_label.setdefault(lab, ma)

    def col(label):
        if label not in by_label:
            raise RuntimeError('column %r not found in %s (have: %s)'
                               % (label, fmodel_mtz, sorted(by_label)))
        return by_label[label].map_to_asu()

    fo = col(fobs_l)
    if fo.sigmas() is None:
        raise RuntimeError('%s has no sigmas (need an F,Q column)' % fobs_l)
    ref = fo.indices()
    uc = fo.unit_cell()

    def aligned_real(label):
        m = col(label)
        out = flex.double(ref.size(), float('nan'))
        for i, j in cctbx_miller.match_indices(ref, m.indices()).pairs():
            out[i] = m.data()[j]
        return out

    def aligned_complex(label):
        m = col(label)
        out = flex.complex_double(ref.size(), 0j)
        for i, j in cctbx_miller.match_indices(ref, m.indices()).pairs():
            out[i] = m.data()[j]
        return out

    kiso = aligned_real(kiso_l)
    kaniso = aligned_real(kaniso_l)
    kika = flex.double([kiso[i] * kaniso[i] for i in range(ref.size())])

    # --- Wilsonify Fobs and its sigma by the exact per-reflection envelope ---
    fo_d = fo.data()
    fo_s = fo.sigmas()
    fo_w = flex.double(ref.size())
    sg_w = flex.double(ref.size())
    for i in range(ref.size()):
        k = kika[i]
        if k <= 0 or k != k:            # guard bad/absent scale -> leave unscaled
            k = 1.0
        fo_w[i] = fo_d[i] / k
        sg_w[i] = fo_s[i] / k
    fo_wil = fo.customized_copy(data=fo_w, sigmas=sg_w)

    # --- free flags -> output convention (refmac test set = free_out, default 0) ---
    fr = col(free_l)
    frr = fr.data()
    frv = flex.int(ref.size(), -999)
    for i, j in cctbx_miller.match_indices(ref, fr.indices()).pairs():
        frv[i] = int(round(float(frr[j])))
    work_val = 1 if free_out == 0 else 0
    fr_out = flex.int(ref.size())
    for i in range(ref.size()):
        fr_out[i] = free_out if frv[i] == free_in else work_val

    # --- assemble output ---
    ds = fo_wil.as_mtz_dataset(column_root_label='FP', column_types='FQ')
    ds.add_miller_array(fo.set().array(data=fr_out), column_root_label='FreeR_flag')

    n_free = (fr_out == free_out).count(True)
    if want_solvent:
        kmask = aligned_real(kmask_l)
        fmask = aligned_complex(fmask_l)
        fpart = flex.complex_double([kmask[i] * fmask[i] for i in range(ref.size())])
        ds.add_miller_array(fo.set().array(data=fpart),
                            column_root_label='Fpart')

    if out_mtz is None:
        base = fmodel_mtz[:-4] if fmodel_mtz.endswith('.mtz') else fmodel_mtz
        out_mtz = base + '_wilsonified.mtz'
    ds.mtz_object().write(out_mtz)

    # --- report: Wilson-plot straightening + the LABIN line ---
    ss = flex.double([uc.d_star_sq(h) for h in ref])
    import numpy as np
    o = np.argsort(list(ss)); n = ref.size()
    fo_o = np.array(list(fo_d)); fo_ww = np.array(list(fo_w))
    print('Wilsonified %s -> %s  (%d reflns, %d free[=%d])'
          % (fmodel_mtz, out_mtz, n, n_free, free_out))
    print('  K_ISO*K_ANISO range %.3f - %.3f' % (min(kika), max(kika)))
    print('  Wilson plot <Fobs> per shell (straightening):')
    print('    %-13s %9s %9s' % ('res(A)', 'orig', 'wilsonified'))
    for b in range(6):
        s = o[b*n//6:(b+1)*n//6]
        d1 = 1.0/math.sqrt(ss[int(s[0])]); d2 = 1.0/math.sqrt(ss[int(s[-1])])
        print('    %5.2f-%5.2f  %9.2f %9.2f'
              % (d1, d2, fo_o[s].mean(), fo_ww[s].mean()))
    print('  refmac LABIN: FP=FP SIGFP=SIGFP FREE=FreeR_flag'
          + (' FPART1=Fpart PHIP1=PHIFpart' if want_solvent else ''))
    print('  refmac keywords: SOLVENT NO ; SHANNON 6 ; (no SCPART)')


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(1)
    kw = dict(free_in=1, free_out=0, want_solvent=True)
    pos = []
    for a in sys.argv[1:]:
        if a == 'no_solvent': kw['want_solvent'] = False
        elif a.startswith('fobs='): kw['fobs_l'] = a.split('=', 1)[1]
        elif a.startswith('free='): kw['free_l'] = a.split('=', 1)[1]
        elif a.startswith('kiso='): kw['kiso_l'] = a.split('=', 1)[1]
        elif a.startswith('kaniso='): kw['kaniso_l'] = a.split('=', 1)[1]
        elif a.startswith('kmask='): kw['kmask_l'] = a.split('=', 1)[1]
        elif a.startswith('fmask='): kw['fmask_l'] = a.split('=', 1)[1]
        elif a.startswith('free_in='): kw['free_in'] = int(a.split('=', 1)[1])
        elif a.startswith('free_out='): kw['free_out'] = int(a.split('=', 1)[1])
        else: pos.append(a)
    fmodel_mtz = pos[0]
    out_mtz = pos[1] if len(pos) > 1 else None
    run(fmodel_mtz, out_mtz, **kw)
