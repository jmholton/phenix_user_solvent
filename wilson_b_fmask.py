#!/usr/bin/env phenix.python
"""
Compute the Wilson B factor of the bulk solvent contribution k_mask * F_mask
from a phenix _f_model.mtz file.

Usage:
    phenix.python wilson_b_fmask.py <f_model.mtz> [n_bins]

n_bins defaults to 20.
"""
import sys
import math
import iotbx.mtz
from cctbx import miller as cctbx_miller
from cctbx.array_family import flex
import numpy as np


def run(mtz_file, n_bins=20):
    mtz_obj = iotbx.mtz.object(file_name=mtz_file)

    # Read as miller arrays so F+PHI pairs are automatically combined
    ma_list = mtz_obj.as_miller_arrays()

    fmask = None
    kmask = None
    for ma in ma_list:
        labels = ma.info().labels
        if "FMASK" in labels and ma.is_complex_array():
            fmask = ma
        if "K_MASK" in labels and ma.is_real_array():
            kmask = ma

    if fmask is None:
        raise RuntimeError("FMASK complex array not found in %s" % mtz_file)
    if kmask is None:
        raise RuntimeError("K_MASK array not found in %s" % mtz_file)

    # Align K_MASK to FMASK index set (should already match, but be safe)
    kmask = kmask.common_set(fmask)
    fmask = fmask.common_set(kmask)

    # Bulk solvent amplitudes: |k_mask * F_mask|
    fbulk_amp = flex.abs(fmask.data()) * kmask.data()

    # d_star_sq = (2 sin θ / λ)² = 1/d²  (cctbx convention)
    d_star_sq = fmask.d_star_sq().data()

    # Bin by resolution and compute mean intensity per bin
    binner = fmask.setup_binner(n_bins=n_bins)

    s2_vals  = []   # mean d_star_sq per bin
    logI_vals = []  # log(mean |F_bulk|²) per bin

    print("\n  bin   d_min(A)  d_max(A)   <|kF|>    <|kF|²>   log(<|kF|²>)")
    print("  " + "-"*62)

    for i_bin in binner.range_used():
        sel = binner.selection(i_bin)
        if sel.count(True) < 3:
            continue
        fb   = fbulk_amp.select(sel)
        s2   = d_star_sq.select(sel)
        mean_I = flex.mean(fb * fb)
        if mean_I <= 0:
            continue
        mean_amp = flex.mean(fb)
        mean_s2  = flex.mean(s2)
        d_range = binner.bin_d_range(i_bin)
        s2_vals.append(mean_s2)
        logI_vals.append(math.log(mean_I))
        print("  %3d   %7.3f   %7.3f   %8.4f  %9.4f   %8.4f" % (
            i_bin, d_range[0], d_range[1], mean_amp, mean_I, math.log(mean_I)))

    if len(s2_vals) < 3:
        raise RuntimeError("Too few resolution bins with data to fit Wilson B")

    # Wilson plot: log(<I>) = log(C) - (B/2) * d_star_sq
    # slope = -B/2  →  B = -2 * slope
    # (because I ~ exp(-2*B*s²) and s² = d_star_sq/4,
    #  so log I = C - B/2 * d_star_sq)
    s2_arr   = np.array(s2_vals)
    logI_arr = np.array(logI_vals)
    slope, intercept = np.polyfit(s2_arr, logI_arr, 1)

    B_wilson = -2.0 * slope

    print("\n  Linear fit: log(<|kF|²>) = %.4f - %.4f * d_star_sq" % (
        intercept, -slope))
    print("  Wilson B of k_mask * F_mask: %.2f A²" % B_wilson)
    print("  (negative B = rising with resolution; unusual but possible for")
    print("   bulk solvent that has structure at medium resolution)\n")
    return B_wilson


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    mtz_file = sys.argv[1]
    n_bins   = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    run(mtz_file, n_bins=n_bins)
