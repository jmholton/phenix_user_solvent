#!/usr/bin/env phenix.python
"""
Extrapolate K_MASK from a phenix.refine _f_model.mtz to all ASU reflections.

phenix.refine writes K_MASK only where F_obs exists. This script fits a
per-bin step function to the existing K_MASK values and assigns the
appropriate bin mean to every reflection in the complete ASU. Original
K_MASK values are preserved exactly at f_obs positions.

Usage:
    phenix.python extrapolate_kmask.py <f_model.mtz> [n_bins] [output.mtz]

n_bins defaults to 20.
Output MTZ column: K_MASK
"""

import sys
import iotbx.mtz
from cctbx import miller as cctbx_miller
from cctbx.array_family import flex


def run(mtz_file, n_bins=20, output_file=None):
    # --- read K_MASK from _f_model.mtz ---
    mtz_obj = iotbx.mtz.object(file_name=mtz_file)
    ma_list = mtz_obj.as_miller_arrays()

    kmask = None
    for ma in ma_list:
        if 'K_MASK' in ma.info().labels and ma.is_real_array():
            kmask = ma
            break

    if kmask is None:
        raise RuntimeError("K_MASK real array not found in %s" % mtz_file)

    cs    = kmask.crystal_symmetry()
    d_min = kmask.d_min()

    # --- build complete ASU at same resolution ---
    complete = cctbx_miller.build_set(
        crystal_symmetry=cs, anomalous_flag=False, d_min=d_min)

    n_obs   = kmask.indices().size()
    n_total = complete.indices().size()
    print("\nK_MASK at f_obs positions : %d" % n_obs)
    print("Complete ASU (d >= %.3f A): %d" % (d_min, n_total))

    # --- bin the existing K_MASK ---
    binner = kmask.setup_binner(n_bins=n_bins)

    bin_ranges = []   # (d_max, d_min) per bin, low-res → high-res
    bin_means  = []

    print("\n  bin   d_max    d_min    n_refl  mean_K_MASK")
    print("  " + "-" * 46)
    for i_bin in binner.range_used():
        sel = binner.selection(i_bin)
        n   = sel.count(True)
        if n < 1:
            continue
        d_range = binner.bin_d_range(i_bin)
        km      = flex.mean(kmask.data().select(sel))
        bin_ranges.append(d_range)
        bin_means.append(km)
        print("  %3d  %7.3f  %7.3f  %7d  %9.4f" % (
            i_bin, d_range[0], d_range[1], n, km))

    # --- assign extrapolated K_MASK to every complete-set reflection ---
    # Work in d_star_sq (= 1/d^2) for vectorised flex selection.
    ds_complete = complete.d_star_sq().data()
    kmask_out   = flex.double(n_total, 0.0)
    assigned    = flex.bool(n_total, False)

    for (d_max, d_min_b), km in zip(bin_ranges, bin_means):
        # d_max  = low-res edge  → small d_star_sq
        # d_min_b = high-res edge → large d_star_sq
        ds_lo = 1.0 / (d_max   * d_max)
        ds_hi = 1.0 / (d_min_b * d_min_b)
        sel   = (ds_complete >= ds_lo) & (ds_complete <= ds_hi)
        kmask_out = kmask_out.set_selected(sel, km)
        assigned  = assigned.set_selected(sel, True)

    # Clamp any reflections that fell outside all bin ranges to the
    # nearest endpoint value (happens when complete set extends slightly
    # beyond the bin boundaries due to floating-point rounding).
    unassigned = ~assigned
    if unassigned.count(True) > 0:
        ds_lo_first = 1.0 / (bin_ranges[0][0]  ** 2)
        ds_hi_last  = 1.0 / (bin_ranges[-1][1] ** 2)
        sel_below = unassigned & (ds_complete < ds_lo_first)
        sel_above = unassigned & (ds_complete > ds_hi_last)
        kmask_out = kmask_out.set_selected(sel_below, bin_means[0])
        kmask_out = kmask_out.set_selected(sel_above, bin_means[-1])
        n_clamped = sel_below.count(True) + sel_above.count(True)
        if n_clamped:
            print("  Warning: %d reflections outside all bins — clamped to "
                  "endpoint value" % n_clamped)

    # --- override with original values where F_obs existed ---
    match  = cctbx_miller.match_indices(complete.indices(), kmask.indices())
    n_orig = match.pairs().size()
    for i_complete, i_orig in match.pairs():
        kmask_out[i_complete] = kmask.data()[i_orig]

    n_filled = n_total - n_orig
    print("\n  %d at f_obs positions (original K_MASK preserved)" % n_orig)
    print("  %d filled by bin extrapolation" % n_filled)

    # --- write output ---
    if output_file is None:
        base        = mtz_file[:-4] if mtz_file.endswith('.mtz') else mtz_file
        output_file = base + '_kmask_complete.mtz'

    out_array   = complete.array(data=kmask_out)
    mtz_dataset = out_array.as_mtz_dataset(column_root_label='K_MASK')
    mtz_dataset.mtz_object().write(file_name=output_file)
    print("  wrote %d reflections → %s" % (n_total, output_file))

    return out_array


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    mtz_file    = sys.argv[1]
    n_bins      = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    output_file = sys.argv[3] if len(sys.argv) > 3 else None
    run(mtz_file, n_bins=n_bins, output_file=output_file)
