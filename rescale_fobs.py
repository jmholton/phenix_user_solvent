#!/usr/bin/env phenix.python
"""
Rescale Fobs by the inverse of phenix K_ISOTROPIC * K_ANISOTROPIC.

Extrapolates both scale arrays from F_obs positions to all reflections:
  K_ISOTROPIC  : per-bin means  (same algorithm as extrapolate_kmask.py)
  K_ANISOTROPIC: fitted 6-parameter u_star tensor
    ln(K_ANISO) = -2pi^2 (U11 h^2 + U22 k^2 + U33 l^2
                          + 2U12 hk + 2U13 hl + 2U23 kl)

Output MTZ replaces FP -> FP/(K_ISO*K_ANISO) and SIGFP likewise.
All other columns (Fpart, FREE, etc.) are passed through unchanged.

Usage:
    phenix.python rescale_fobs.py  main.mtz  phenix_f_model.mtz  [output.mtz]
                  [fp=FP] [sigfp=SIGFP] [iso=K_ISOTROPIC] [aniso=K_ANISOTROPIC]
                  [n_bins=20]
"""

import sys, math
import numpy as np
import iotbx.mtz
from cctbx import miller as cctbx_miller
from cctbx.array_family import flex

TWO_PI_SQ = 2.0 * math.pi ** 2


# ---------------------------------------------------------------------------
def _extrapolate_kiso(kiso_ma, n_bins=20):
    """Per-bin mean extrapolation to complete ASU; restores originals in place."""
    cs    = kiso_ma.crystal_symmetry()
    d_min = kiso_ma.d_min()
    complete = cctbx_miller.build_set(cs, anomalous_flag=False, d_min=d_min)

    binner     = kiso_ma.setup_binner(n_bins=n_bins)
    bin_ranges = []
    bin_means  = []
    print("  K_ISO bin means:")
    for i_bin in binner.range_used():
        sel = binner.selection(i_bin)
        if sel.count(True) < 1:
            continue
        d_range = binner.bin_d_range(i_bin)
        km      = flex.mean(kiso_ma.data().select(sel))
        bin_ranges.append(d_range)
        bin_means.append(km)
        print("    %6.3f - %6.3f A   n=%5d   mean=%.5f" % (
            d_range[0], d_range[1], sel.count(True), km))

    ds         = complete.d_star_sq().data()
    kiso_out   = flex.double(complete.indices().size(), 0.0)
    assigned   = flex.bool(complete.indices().size(), False)

    for (d_max, d_min_b), km in zip(bin_ranges, bin_means):
        ds_lo = 1.0 / (d_max   ** 2)
        ds_hi = 1.0 / (d_min_b ** 2)
        sel   = (ds >= ds_lo) & (ds <= ds_hi)
        kiso_out = kiso_out.set_selected(sel, km)
        assigned = assigned.set_selected(sel, True)

    unassigned = ~assigned
    if unassigned.count(True) > 0:
        ds_lo_first = 1.0 / bin_ranges[0][0]  ** 2
        ds_hi_last  = 1.0 / bin_ranges[-1][1] ** 2
        kiso_out = kiso_out.set_selected(
            unassigned & (ds < ds_lo_first), bin_means[0])
        kiso_out = kiso_out.set_selected(
            unassigned & (ds > ds_hi_last),  bin_means[-1])

    match = cctbx_miller.match_indices(complete.indices(), kiso_ma.indices())
    for i_c, i_o in match.pairs():
        kiso_out[i_c] = kiso_ma.data()[i_o]

    return complete.array(data=kiso_out)


# ---------------------------------------------------------------------------
def _fit_aniso_tensor(kaniso_ma):
    """Fit intercept + 6 u_star components by linear least squares from K_ANISO values.

    Returns (c, u_star) where c = ln(k_overall) and u_star = [U11..U23].
    The full model is: ln(K_ANISO) = c + (-2pi^2) * h^T U* h
    """
    hkl  = kaniso_ma.map_to_asu().indices()
    lnK  = np.log(np.array(list(kaniso_ma.data()), dtype=float))

    n    = len(lnK)
    A    = np.zeros((n, 7), dtype=float)
    A[:, 0] = 1.0  # intercept — absorbs ln(k_overall)
    for i, (h, k, l) in enumerate(hkl):
        A[i, 1:] = [-TWO_PI_SQ * h*h, -TWO_PI_SQ * k*k, -TWO_PI_SQ * l*l,
                    -TWO_PI_SQ * 2*h*k, -TWO_PI_SQ * 2*h*l, -TWO_PI_SQ * 2*k*l]

    params, _, _, _ = np.linalg.lstsq(A, lnK, rcond=None)
    c      = params[0]
    u_star = params[1:]
    return c, u_star


def _apply_aniso_tensor(c, u_star, miller_set):
    """Evaluate exp(c - 2pi^2 h^T U* h) for every reflection in miller_set."""
    U11, U22, U33, U12, U13, U23 = u_star
    kaniso = flex.double(miller_set.indices().size())
    for i, (h, k, l) in enumerate(miller_set.indices()):
        exp_val = c - TWO_PI_SQ * (U11*h*h + U22*k*k + U33*l*l +
                                    2*U12*h*k + 2*U13*h*l + 2*U23*k*l)
        kaniso[i] = math.exp(exp_val)
    return miller_set.array(data=kaniso)


# ---------------------------------------------------------------------------
def run(main_mtz, phenix_mtz, output_mtz=None,
        fp_label='FP', sigfp_label='SIGFP',
        iso_label='K_ISOTROPIC', aniso_label='K_ANISOTROPIC',
        n_bins=20):

    # --- read phenix scale columns ---
    ph = iotbx.mtz.object(file_name=phenix_mtz)
    kiso_ma = kaniso_ma = None
    for ma in ph.as_miller_arrays():
        labels = ma.info().labels
        if iso_label   in labels and ma.is_real_array(): kiso_ma   = ma
        if aniso_label in labels and ma.is_real_array(): kaniso_ma = ma
    if kiso_ma   is None: raise RuntimeError('%s not found in %s' % (iso_label,   phenix_mtz))
    if kaniso_ma is None: raise RuntimeError('%s not found in %s' % (aniso_label, phenix_mtz))

    kiso_ma   = kiso_ma.map_to_asu()
    kaniso_ma = kaniso_ma.map_to_asu()
    print("K_ISO:   %d reflections  range %.4f - %.4f" % (
        kiso_ma.indices().size(),
        flex.min(kiso_ma.data()), flex.max(kiso_ma.data())))
    print("K_ANISO: %d reflections  range %.4f - %.4f" % (
        kaniso_ma.indices().size(),
        flex.min(kaniso_ma.data()), flex.max(kaniso_ma.data())))

    # --- extrapolate K_ISO ---
    print("\nExtrapolating K_ISO to complete ASU...")
    kiso_complete = _extrapolate_kiso(kiso_ma, n_bins=n_bins)

    # --- fit and extrapolate K_ANISO ---
    print("\nFitting u_star tensor to K_ANISO (%d values)..." % kaniso_ma.indices().size())
    c, u_star = _fit_aniso_tensor(kaniso_ma)
    print("  fitted k_overall = exp(c) = %.5f  (phenix reports k_overall in log)" % math.exp(c))
    print("  U_star (frac): U11=%+.6f U22=%+.6f U33=%+.6f" % (u_star[0], u_star[1], u_star[2]))
    print("                 U12=%+.6f U13=%+.6f U23=%+.6f" % (u_star[3], u_star[4], u_star[5]))
    kaniso_complete = _apply_aniso_tensor(c, u_star, kiso_complete)

    # verify fit residual
    match_v = cctbx_miller.match_indices(kaniso_ma.indices(), kaniso_complete.indices())
    diffs = flex.double()
    for i_orig, i_fit in match_v.pairs():
        diffs.append(abs(kaniso_ma.data()[i_orig] - kaniso_complete.data()[i_fit]))
    print("  fit residual |K_ANISO_orig - K_ANISO_fit|: mean=%.6f  max=%.6f" % (
        flex.mean(diffs), flex.max(diffs)))

    # combined scale at complete ASU
    combined = kiso_complete.array(data=kiso_complete.data() * kaniso_complete.data())
    print("\nCombined scale (K_ISO*K_ANISO) range: %.4f - %.4f  mean: %.4f" % (
        flex.min(combined.data()), flex.max(combined.data()), flex.mean(combined.data())))

    # --- read and rescale main MTZ ---
    main_obj    = iotbx.mtz.object(file_name=main_mtz)
    main_arrays = main_obj.as_miller_arrays()

    out_arrays = []
    for ma in main_arrays:
        labels = ma.info().labels
        is_fp    = fp_label    in labels and ma.is_real_array()
        is_sigfp = sigfp_label in labels and ma.is_real_array()
        if is_fp or is_sigfp:
            ma_asu = ma.map_to_asu()
            match  = cctbx_miller.match_indices(ma_asu.indices(), combined.indices())
            data_new   = ma_asu.data().deep_copy()
            sigmas_new = ma_asu.sigmas().deep_copy() if ma_asu.sigmas() is not None else None
            n_scaled = 0
            for i_ma, i_sc in match.pairs():
                scale = combined.data()[i_sc]
                data_new[i_ma] = ma_asu.data()[i_ma] / scale
                if sigmas_new is not None:
                    sigmas_new[i_ma] = ma_asu.sigmas()[i_ma] / scale
                n_scaled += 1
            sig_note = " (sigmas also rescaled)" if sigmas_new is not None else ""
            print("Rescaled %-8s: %d / %d reflections divided by K_ISO*K_ANISO%s" % (
                labels[0], n_scaled, ma_asu.indices().size(), sig_note))
            out_arrays.append(ma_asu.array(data=data_new, sigmas=sigmas_new))
        else:
            out_arrays.append(ma)

    # --- write output ---
    mtz_dataset = out_arrays[0].as_mtz_dataset(
        column_root_label=main_arrays[0].info().labels[0])
    for ma, ma_orig in zip(out_arrays[1:], main_arrays[1:]):
        mtz_dataset.add_miller_array(
            miller_array=ma,
            column_root_label=ma_orig.info().labels[0])

    if output_mtz is None:
        base       = main_mtz[:-4] if main_mtz.endswith('.mtz') else main_mtz
        output_mtz = base + '_rescaled.mtz'

    mtz_dataset.mtz_object().write(file_name=output_mtz)
    print("\nWrote: %s" % output_mtz)
    print("Refmac LABIN: FP=FP SIGFP=SIGFP FREE=FreeR_flag FPART1=Fpart PHIPART1=PHIFpart")
    print("Add:  SCAL MLSC FIXB B 0.0    (or let refmac fit; should converge near B=0)")


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    main_mtz   = sys.argv[1]
    phenix_mtz = sys.argv[2]
    output_mtz = None
    fp_label    = 'FP'
    sigfp_label = 'SIGFP'
    iso_label   = 'K_ISOTROPIC'
    aniso_label = 'K_ANISOTROPIC'
    n_bins      = 20
    for arg in sys.argv[3:]:
        if   arg.startswith('fp='):    fp_label    = arg[3:]
        elif arg.startswith('sigfp='): sigfp_label = arg[6:]
        elif arg.startswith('iso='):   iso_label   = arg[4:]
        elif arg.startswith('aniso='): aniso_label = arg[6:]
        elif arg.startswith('n_bins='): n_bins     = int(arg[7:])
        elif arg.endswith('.mtz'):     output_mtz  = arg
    run(main_mtz, phenix_mtz, output_mtz=output_mtz,
        fp_label=fp_label, sigfp_label=sigfp_label,
        iso_label=iso_label, aniso_label=aniso_label,
        n_bins=n_bins)
