#!/usr/bin/env phenix.python
"""
sigma_power: per-resolution-shell noise power Pnn = <sigF^2>, binned by stol^2, for
squish's data-driven resolution filters (resfilter=wiener|snr).

Wiener filter:  H(s) = Pss/(Pss+Pnn) = 1 - Pnn/Ptot,  Ptot = <|Fdiff|^2> (squish's
wilson_<map>.txt), Pnn = <sigF^2> (this tool). The experimental sigma is the noise floor
of the difference map, so shells where the difference is all noise (Ptot ~ Pnn) get H->0
and shells with real signal (Ptot >> Pnn) get H->1 -- a self-tuning band-pass.

Bins exactly like squish's Wilson plot: bin = round(stol^2 / binsize), stol^2 = d*^2/4,
free reflections excluded (CCP4 convention free = 0), origin bin dropped. Output columns:
  stol^2   <sigF^2>
matched to wilson_<map>.txt by bin index (bs) on the squish side.

Usage:
  phenix.python sigma_power.py  data.mtz  binsize  out.txt  [free=FreeR_flag] [sig=SIGFOBS]
"""

import sys
import iotbx.mtz
from cctbx.array_family import flex


def run(mtz, binsize, out, free_l='FreeR_flag', sig_l='SIGFOBS'):
    binsize = float(binsize)
    o = iotbx.mtz.object(mtz)
    by = {}
    for ma in o.as_miller_arrays():
        for lab in ma.info().labels:
            by.setdefault(lab, ma)

    # the F,Q array carrying sigmas (SIGFOBS lives on the FOBS amplitude array)
    fo = None
    for lab, ma in by.items():
        if ma.sigmas() is not None and (sig_l in ma.info().labels or lab == 'FOBS'):
            fo = ma.map_to_asu(); break
    if fo is None or fo.sigmas() is None:
        raise RuntimeError('no amplitude array with sigmas (%s) in %s' % (sig_l, mtz))

    free = by.get(free_l)
    free_by_h = {}
    if free is not None:
        fr = free.map_to_asu()
        free_by_h = {h: int(round(v)) for h, v in zip(fr.indices(), fr.data())}

    uc = fo.unit_cell()
    sums = {}
    counts = {}
    for h, sig in zip(fo.indices(), fo.sigmas()):
        if free_by_h.get(h, 1) == 0:        # exclude the free set (CCP4 free = 0)
            continue
        stol2 = uc.d_star_sq(h) / 4.0
        b = int(round(stol2 / binsize))
        if b <= 0:                          # drop origin bin (matches Wilson binning)
            continue
        sums[b] = sums.get(b, 0.0) + sig * sig
        counts[b] = counts.get(b, 0) + 1

    with open(out, 'w') as f:
        for b in sorted(sums):
            f.write('%g %g\n' % (b * binsize, sums[b] / counts[b]))
    print('sigma_power: %s -> %s  (%d shells, binsize=%g, free excluded)'
          % (mtz, out, len(sums), binsize))


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print(__doc__); sys.exit(1)
    kw = {}
    pos = []
    for a in sys.argv[1:]:
        if a.startswith('free='): kw['free_l'] = a.split('=', 1)[1]
        elif a.startswith('sig='): kw['sig_l'] = a.split('=', 1)[1]
        else: pos.append(a)
    run(pos[0], pos[1], pos[2], **kw)
