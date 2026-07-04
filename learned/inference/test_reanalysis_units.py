#!/usr/bin/env python
# SPDX-License-Identifier: MIT
"""
Unit tests for the DATA-INDEPENDENT parts of the new reanalysis count path (no real caches needed):
  1. inpolygon ROI rasterization (square polygon -> mask + in-ROI count) + 'full' sentinel
  2. loc_sub aggregation (variable blocks-per-rung + pedestal subtract + floor-at-0)
  3. lf fit (a known log-log slope recovers beta; n<3 -> NaN; x<=0/y<=0 filtered)
  4. datasets.py block-tag expansion (apr17 24 rung tags; jun23 76 all tags) + key/label helpers

Run: python test_reanalysis_units.py
"""
import os, sys, traceback
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import datasets as DS
import count_locrate as CL

_results = []


def check(name, cond, detail=""):
    _results.append((name, bool(cond), detail))
    print("  [%s] %s%s" % ("PASS" if cond else "FAIL", name, (" -- " + detail) if detail else ""))


# ---------------------------------------------------------------- 1. ROI rasterization
def test_rasterize():
    print("\n== 1. inpolygon ROI rasterization ==")
    xGrid = np.arange(0.0, 4.0001, 0.25)        # 17 cols, indices 0..16
    zGrid = np.arange(8.0, 14.0001, 0.25)       # 25 rows, indices 0..24
    poly = np.array([[1, 10], [3, 10], [3, 12], [1, 12]], float)   # x in [1,3], z in [10,12]
    mask = CL.rasterize_roi(poly, xGrid, zGrid)
    check("mask shape == [nZ x nX]", mask.shape == (len(zGrid), len(xGrid)), str(mask.shape))

    # interior pixel x=2 (col 8), z=11 (row 12) -> inside
    check("interior pixel True", mask[12, 8], "mask[12,8]=%s" % mask[12, 8])
    # exterior pixel x=0.5 (col 2), z=9 (row 4) -> outside
    check("exterior pixel False", not mask[4, 2], "mask[4,2]=%s" % mask[4, 2])

    # synthetic detections: 4 strictly inside, 3 strictly outside
    inside = [(12, 8), (10, 6), (14, 10), (12, 6)]      # z 11,10.5,11.5,11 ; x 2,1.5,2.5,1.5
    outside = [(4, 2), (24, 16), (0, 0)]
    rows = np.array([r for r, c in inside + outside])
    cols = np.array([c for r, c in inside + outside])
    # verify each "inside" pixel is actually inside the mask we built (geometry sanity)
    inside_ok = all(mask[r, c] for r, c in inside)
    outside_ok = all(not mask[r, c] for r, c in outside)
    check("inside pts land in mask", inside_ok)
    check("outside pts land outside mask", outside_ok)
    cnt = CL.count_in_roi(rows, cols, mask)
    check("count_in_roi == #inside (4)", cnt == 4, "got %d" % cnt)

    # 'full' sentinel: poly None -> all true; count == # in-bounds detections (7)
    full = CL.rasterize_roi(None, xGrid, zGrid)
    check("'full' sentinel all-true", full.all() and full.shape == (25, 17))
    check("count in 'full' == 7", CL.count_in_roi(rows, cols, full) == 7,
          "got %d" % CL.count_in_roi(rows, cols, full))

    # loc_rate divides by nFrames
    lr = CL.loc_rate(rows, cols, mask, nFrames=2)
    check("loc_rate == count/nFrames", np.isclose(lr, 4 / 2.0), "got %g" % lr)

    # empty detections -> 0
    check("empty detections -> 0 count", CL.count_in_roi(np.array([]), np.array([]), mask) == 0)


# ---------------------------------------------------------------- 2. loc_sub aggregation
def test_aggregate():
    print("\n== 2. loc_sub aggregation (variable blocks/rung + pedestal + floor) ==")
    # 3 rungs, variable blocks [2,3,1]; one NaN block to exercise omitnan
    rung_rates = [
        [10.0, 12.0],            # mean 11
        [20.0, np.nan, 24.0],    # mean 22 (omitnan)
        [3.0],                   # mean 3
    ]
    ped = 5.0
    loc_mean, loc_sub = CL.aggregate_loc_sub(rung_rates, ped)
    check("loc_mean omitnan", np.allclose(loc_mean, [11.0, 22.0, 3.0]), str(loc_mean))
    # loc_sub = max(mean-ped,0): [6,17,0] (rung2 3-5=-2 -> floored to 0)
    check("loc_sub subtract+floor", np.allclose(loc_sub, [6.0, 17.0, 0.0]), str(loc_sub))

    # pedestal helper: mean over bgFlow rates, omitnan; empty -> 0
    check("pedestal mean omitnan", np.isclose(CL.pedestal([4.0, np.nan, 6.0]), 5.0))
    check("pedestal empty -> 0", CL.pedestal([]) == 0.0)
    check("pedestal all-nan -> 0", CL.pedestal([np.nan, np.nan]) == 0.0)

    # whole-rung NaN -> loc_mean NaN, loc_sub max(nan-ped,0) -> nan stays (not a crash)
    lm, lsub = CL.aggregate_loc_sub([[np.nan, np.nan]], 1.0)
    check("all-nan rung -> loc_mean NaN", np.isnan(lm[0]))


# ---------------------------------------------------------------- 3. lf fit
def test_lf_fit():
    print("\n== 3. lf fit (log-log slope recovery) ==")
    x = np.array([1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0])
    beta_true = 1.5
    y = 3.0 * x ** beta_true                       # exact power law
    fit = CL.lf_fit(x, y, n_boot=200, seed=1)
    check("beta recovers slope 1.5", np.isclose(fit["beta"], beta_true, atol=1e-9),
          "beta=%.10f" % fit["beta"])
    check("R2 ~ 1", fit["R2"] > 1 - 1e-12, "R2=%.12f" % fit["R2"])
    check("n == 8", fit["n"] == 8)
    check("ci brackets beta", fit["ci"][0] <= fit["beta"] <= fit["ci"][1],
          "ci=%s" % fit["ci"])

    # n<3 -> NaN
    fitn = CL.lf_fit([1.0, 2.0], [1.0, 4.0])
    check("n<3 -> NaN beta", np.isnan(fitn["beta"]) and fitn["n"] == 2)

    # x<=0 / y<=0 filtered (drop the bad points, keep the rest)
    xf = np.array([0.0, 1.0, 2.0, 4.0, 8.0])
    yf = np.array([5.0, 3.0, 3.0 * 2 ** 1.5, 3.0 * 4 ** 1.5, 3.0 * 8 ** 1.5])
    fitf = CL.lf_fit(xf, yf)
    check("x<=0 dropped (n==4)", fitf["n"] == 4, "n=%d" % fitf["n"])
    check("slope still 1.5 after filter", np.isclose(fitf["beta"], 1.5, atol=1e-9),
          "beta=%.10f" % fitf["beta"])


# ---------------------------------------------------------------- 4. block-tag expansion
def test_block_tags():
    print("\n== 4. datasets.py block-tag expansion ==")
    a = DS.get_dataset("apr17")
    j = DS.get_dataset("jun23")

    a_rung = DS.rung_block_tags(a)
    a_all = DS.all_block_tags(a)
    j_rung = DS.rung_block_tags(j)
    j_all = DS.all_block_tags(j)

    check("apr17 rung tags == 24", len(a_rung) == 24, "got %d" % len(a_rung))
    check("apr17 first/last == C1b1/C8b3", a_rung[0] == "C1b1" and a_rung[-1] == "C8b3",
          "%s..%s" % (a_rung[0], a_rung[-1]))
    check("apr17 lowercase 'b' format", all("b" in t and "B" not in t for t in a_rung))
    check("apr17 all (rung+bg) == 27", len(a_all) == 27, "got %d" % len(a_all))
    check("apr17 conc len == 8", len(a["conc"]) == 8)

    check("jun23 rung tags == 68", len(j_rung) == 68, "got %d" % len(j_rung))
    check("jun23 ALL tags == 76", len(j_all) == 76, "got %d" % len(j_all))
    check("jun23 first rung == L1B1", j_rung[0] == "L1B1", j_rung[0])
    check("jun23 M1B1 present (uppercase B)", "M1B1" in j_rung)
    check("jun23 blocksPerRung sums to 68", sum(j["blocksPerRung"]) == 68)
    check("jun23 conc len == 15", len(j["conc"]) == 15)
    check("jun23 bgFlow has 6, bgStatic 2", len(j["bgFlow"]) == 6 and len(j["bgStatic"]) == 2)

    # helpers
    check("window_label C3-C7", DS.window_label(a, [3, 7]) == "C3-C7", DS.window_label(a, [3, 7]))
    check("window_label M1-M5", DS.window_label(j, [6, 10]) == "M1-M5", DS.window_label(j, [6, 10]))
    key = CL.make_beta_key("PI", "combinedTube", "locRate", "C3-C7")
    check("beta key PI_combinedTube_locRate_C3_C7", key == "PI_combinedTube_locRate_C3_C7", key)
    check("pol_cache_path basename", os.path.basename(DS.pol_cache_path(a, "C3b2")) == "C3b2_POL.mat",
          os.path.basename(DS.pol_cache_path(a, "C3b2")))
    check("rung_index_of_tag C8b3 -> 7", DS.rung_index_of_tag(a, "C8b3") == 7)
    check("rung_index_of_tag Bg1 -> None", DS.rung_index_of_tag(a, "Bg1") is None)


def main():
    print("=" * 70)
    print("  Phase 3 reanalysis unit tests (data-independent)")
    print("=" * 70)
    for fn in (test_rasterize, test_aggregate, test_lf_fit, test_block_tags):
        try:
            fn()
        except Exception:
            _results.append((fn.__name__, False, "EXCEPTION"))
            print("  [FAIL] %s raised:" % fn.__name__)
            traceback.print_exc()
    n_pass = sum(1 for _, ok, _ in _results if ok)
    n_tot = len(_results)
    print("\n" + "=" * 70)
    print("  RESULT: %d/%d checks passed%s" % (n_pass, n_tot, "" if n_pass == n_tot else "  <-- FAILURES"))
    print("=" * 70)
    return 0 if n_pass == n_tot else 1


if __name__ == "__main__":
    sys.exit(main())
