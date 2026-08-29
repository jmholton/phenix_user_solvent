#!/usr/bin/env phenix.python
"""
Package phenix bulk solvent for a refmac run.

Extrapolates K_MASK from F_obs positions to all ASU reflections using per-bin
means (same method as extrapolate_kmask.py), applies it to the original
complete Fpart/PHIpart, then merges into the main reflection MTZ.

    K_MASK_complete(h) = bin_mean(K_MASK)   extrapolated
                       = K_MASK(h) exactly  where F_obs exists

    Fpart_out(h) = K_MASK_complete(h) * Fpart_original(h)

Usage:
    phenix.python package_for_refmac.py  main.mtz  solvent.mtz  fmodel.mtz  [output.mtz] \
                  [amp=Fpart] [phi=PHIpart] [kmask=K_MASK] [n_bins=20]

Defaults: amp=Fpart  phi=PHIpart  kmask=K_MASK  n_bins=20
Output MTZ adds:  Fpart (type F)  PHIFpart (type P)

Refmac run:
    refmac5 hklin output.mtz xyzin model.pdb xyzout refined.pdb hklout out.mtz << EOF
    LABIN FP=FP SIGFP=SIGFP FREE=FreeR_flag FPART1=Fpart PHIPART1=PHIFpart
    LABOUT FP=FP SIGFP=SIGFP FREE=FreeR_flag
    NCYC 10
    END
    EOF
"""

import sys
import iotbx.mtz
from cctbx import miller as cctbx_miller
from cctbx.array_family import flex


def _extrapolate_kmask(kmask_ma, n_bins=20):
    """Extrapolate K_MASK (at F_obs positions) to a complete ASU set.

    Bins the existing K_MASK, computes per-bin means, assigns those means to
    every reflection in the complete ASU at the same d_min, then restores
    original values exactly where F_obs existed.  Identical algorithm to
    extrapolate_kmask.py.
    """
    cs    = kmask_ma.crystal_symmetry()
    d_min = kmask_ma.d_min()

    complete = cctbx_miller.build_set(
        crystal_symmetry=cs, anomalous_flag=False, d_min=d_min)

    binner = kmask_ma.setup_binner(n_bins=n_bins)

    bin_ranges = []
    bin_means  = []
    print("\n  bin   d_max    d_min    n_refl  mean_K_MASK")
    print("  " + "-" * 46)
    for i_bin in binner.range_used():
        sel = binner.selection(i_bin)
        if sel.count(True) < 1:
            continue
        d_range = binner.bin_d_range(i_bin)
        km      = flex.mean(kmask_ma.data().select(sel))
        bin_ranges.append(d_range)
        bin_means.append(km)
        print("  %3d  %7.3f  %7.3f  %7d  %9.4f" % (
            i_bin, d_range[0], d_range[1], sel.count(True), km))

    ds_complete = complete.d_star_sq().data()
    kmask_out   = flex.double(complete.indices().size(), 0.0)
    assigned    = flex.bool(complete.indices().size(), False)

    for (d_max, d_min_b), km in zip(bin_ranges, bin_means):
        ds_lo = 1.0 / (d_max   * d_max)
        ds_hi = 1.0 / (d_min_b * d_min_b)
        sel   = (ds_complete >= ds_lo) & (ds_complete <= ds_hi)
        kmask_out = kmask_out.set_selected(sel, km)
        assigned  = assigned.set_selected(sel, True)

    unassigned = ~assigned
    if unassigned.count(True) > 0:
        ds_lo_first = 1.0 / (bin_ranges[0][0]  ** 2)
        ds_hi_last  = 1.0 / (bin_ranges[-1][1] ** 2)
        kmask_out = kmask_out.set_selected(
            unassigned & (ds_complete < ds_lo_first), bin_means[0])
        kmask_out = kmask_out.set_selected(
            unassigned & (ds_complete > ds_hi_last),  bin_means[-1])
        n_clamped = (unassigned & (ds_complete < ds_lo_first)).count(True) + \
                    (unassigned & (ds_complete > ds_hi_last)).count(True)
        if n_clamped:
            print("  Warning: %d reflections clamped to endpoint bin" % n_clamped)

    # Restore original K_MASK exactly at F_obs positions.
    match = cctbx_miller.match_indices(complete.indices(), kmask_ma.indices())
    for i_complete, i_orig in match.pairs():
        kmask_out[i_complete] = kmask_ma.data()[i_orig]

    return complete.array(data=kmask_out)   # complete ASU, extrapolated K_MASK


def run(main_mtz, solvent_mtz, fmodel_mtz, output_file=None,
        amp_label='Fpart', phi_label='PHIpart', kmask_label='K_MASK', n_bins=20):

    # --- read original (complete) Fpart/PHIpart ---
    sol_obj = iotbx.mtz.object(file_name=solvent_mtz)
    amp_array = phi_array = None
    for ma in sol_obj.as_miller_arrays():
        labels = ma.info().labels
        if amp_label in labels:
            amp_array = ma
        if phi_label in labels:
            phi_array = ma
    if amp_array is None:
        raise RuntimeError("Column '%s' not found in %s" % (amp_label, solvent_mtz))
    if phi_array is None:
        raise RuntimeError("Column '%s' not found in %s" % (phi_label, solvent_mtz))
    f_complex = amp_array.phase_transfer(phase_source=phi_array, deg=True).map_to_asu()
    print("Original solvent (complete): %d reflections" % f_complex.indices().size())

    # --- read K_MASK from _f_model.mtz and extrapolate to complete ASU ---
    fm_obj = iotbx.mtz.object(file_name=fmodel_mtz)
    kmask_ma = None
    for ma in fm_obj.as_miller_arrays():
        if kmask_label in ma.info().labels and ma.is_real_array():
            kmask_ma = ma
            break
    if kmask_ma is None:
        raise RuntimeError("Column '%s' not found in %s" % (kmask_label, fmodel_mtz))
    kmask_ma = kmask_ma.map_to_asu()
    print("K_MASK at F_obs positions   : %d" % kmask_ma.indices().size())

    print("\nExtrapolating K_MASK to complete ASU...")
    kmask_complete = _extrapolate_kmask(kmask_ma, n_bins=n_bins)
    print("K_MASK complete ASU         : %d" % kmask_complete.indices().size())

    # --- apply extrapolated K_MASK to original Fpart ---
    # Match kmask_complete onto f_complex (both should be the same complete ASU).
    match_km = cctbx_miller.match_indices(f_complex.indices(), kmask_complete.indices())
    scaled_data = f_complex.data().deep_copy()
    for i_fc, i_km in match_km.pairs():
        scaled_data[i_fc] = kmask_complete.data()[i_km] * f_complex.data()[i_fc]
    n_applied = match_km.pairs().size()
    n_missing = f_complex.indices().size() - n_applied
    print("\nK_MASK applied: %d  (unmatched, left at zero: %d)" % (n_applied, n_missing))
    f_scaled = f_complex.array(data=scaled_data)

    # --- read main MTZ and build output ---
    main_obj = iotbx.mtz.object(file_name=main_mtz)
    main_arrays = main_obj.as_miller_arrays()
    if not main_arrays:
        raise RuntimeError("No miller arrays found in %s" % main_mtz)

    ref = main_arrays[0]
    n   = ref.indices().size()

    complex_out = flex.complex_double(n, 0+0j)
    match = cctbx_miller.match_indices(ref.indices(), f_scaled.indices())
    for i_ref, i_sol in match.pairs():
        complex_out[i_ref] = f_scaled.data()[i_sol]
    print("Matched to main MTZ: %d / %d" % (match.pairs().size(), n))

    # Exclude any existing Fpart/PHIFpart columns — they will be replaced.
    out_labels = {'Fpart', 'PHIFpart', 'PHIpart'}
    filtered = [ma for ma in main_arrays
                if not any(l in out_labels for l in ma.info().labels)]
    if not filtered:
        raise RuntimeError("All columns filtered out — nothing left from main MTZ")

    mtz_dataset = filtered[0].as_mtz_dataset(
        column_root_label=filtered[0].info().labels[0])
    for ma in filtered[1:]:
        mtz_dataset.add_miller_array(
            miller_array=ma,
            column_root_label=ma.info().labels[0])

    fpart = ref.array(data=complex_out)
    mtz_dataset.add_miller_array(miller_array=fpart, column_root_label='Fpart')

    if output_file is None:
        base = main_mtz[:-4] if main_mtz.endswith('.mtz') else main_mtz
        output_file = base + '_with_fpart.mtz'

    mtz_dataset.mtz_object().write(file_name=output_file)
    print("Wrote: %s  (%d reflections)" % (output_file, n))
    print("\nRefmac LABIN card:")
    print("  LABIN FP=FP SIGFP=SIGFP FREE=FreeR_flag FPART1=Fpart PHIPART1=PHIFpart")


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit(1)
    main_mtz    = sys.argv[1]
    solvent_mtz = sys.argv[2]
    fmodel_mtz  = sys.argv[3]
    output_file = None
    amp_label   = 'Fpart'
    phi_label   = 'PHIpart'
    kmask_label = 'K_MASK'
    n_bins      = 20
    for arg in sys.argv[4:]:
        if arg.startswith('amp='):
            amp_label = arg[4:]
        elif arg.startswith('phi='):
            phi_label = arg[4:]
        elif arg.startswith('kmask='):
            kmask_label = arg[6:]
        elif arg.startswith('n_bins='):
            n_bins = int(arg[7:])
        elif arg.endswith('.mtz'):
            output_file = arg
    run(main_mtz, solvent_mtz, fmodel_mtz, output_file=output_file,
        amp_label=amp_label, phi_label=phi_label,
        kmask_label=kmask_label, n_bins=n_bins)
