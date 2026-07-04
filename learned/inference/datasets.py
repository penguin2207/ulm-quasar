#!/usr/bin/env python
# SPDX-License-Identifier: MIT
"""
Phase 3 reanalysis dataset config (FROZEN spec, REANALYSIS_INPUTS_RESPONSE.md Section 5).

ONE struct per dataset (apr17, jun23) carrying every value the new in-tube-localization-rate
count path needs: paths, polarity-cache naming, block-tag format, rung labels + concentration
axis (verbatim, in order), variable blocks-per-rung, Bg tags, fit windows, headline domain/ROI,
SVD cutoff/seed, derived domains, FISTA lambda.

NO tracking, NO QC anywhere in this pipeline -- the count is `detections-in-ROI / nFrames`
(Bg-pedestal-subtracted). See RESPONSE.md Sections 2, 4, 6, 7.

Block-tag construction (Section 5):
    rung block tag = blockFmt % (rungLabel, blockIdx)
      apr17 blockFmt '%sb%d' (LOWERCASE b) -> C3b2
      jun23 blockFmt '%sB%d' (UPPERCASE B) -> M3B4
    Bg tags (bgFlow + bgStatic) are used verbatim.
    Polarity cache file = '<tag>_POL.mat' in <polCacheDir>.
"""

import os
from paths import OUTPUT_ROOT, RAW_ROOT  # machine-local roots (learned_config.py / env); see learned/README.md

# domain order is FROZEN (matches readout.domains, RESPONSE.md Section 6a)
DOMAINS = ["PI", "fundamental", "singlepol"]

# ROI order is FROZEN (matches roi.names / readout.roiNames, RESPONSE.md Section 4b/6a)
ROI_NAMES = ["full", "combinedTube", "tubeL", "tubeR", "background"]

# SVD clutter filter (RESPONSE.md Sections 2, 5): cutoff 8, seed 12345.
SVD_CUTOFF = 8
SVD_SEED = 12345

# FISTA sparse-deconvolution regularizer (RESPONSE.md Section 5).
FISTA_LAMBDA = 0.10


DATASETS = {
    # ------------------------------------------------------------------ apr17
    "apr17": dict(
        name="apr17",
        rawRoot=os.path.join(RAW_ROOT, "04-17-raw-data"),
        outDir=os.path.join(OUTPUT_ROOT, "apr17"),
        polCacheDir=os.path.join(OUTPUT_ROOT, "apr17", "pol_cache"),
        cacheFmt="%s_POL.mat",            # %s = userLabel (block tag or Bg tag)
        blockFmt="%sb%d",                 # lowercase b  -> C1b1 .. C8b3
        roiFile=os.path.join(OUTPUT_ROOT, "apr17", "roi_polys_apr17.mat"),
        readoutFile=os.path.join(OUTPUT_ROOT, "apr17", "readout_apr17.mat"),
        betasFile=os.path.join(OUTPUT_ROOT, "apr17", "betas_apr17.mat"),
        nRung=8,
        rungLabels=["C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"],
        # MB/mL, in order (RESPONSE.md Section 5)
        conc=[2.9e5, 4.4e5, 6.6e5, 9.9e5, 1.5e6, 2.2e6, 3.3e6, 5.0e6],
        blocksPerRung=[3, 3, 3, 3, 3, 3, 3, 3],     # 8*3 = 24 rung block tags
        bgFlow=["Bg1", "Bg2"],                       # sweep + pedestal (subtracted)
        bgStatic=["C5NoFlow"],                       # reported, NOT subtracted
        nTubes=2,
        fitWindows=[[3, 7], [3, 8], [1, 8]],         # 1-based rung-index ranges
        headlineWindow=[3, 7],                       # C3-C7
        headlineDomain="PI",
        headlineROI="combinedTube",
        domains=list(DOMAINS),
        svdCutoff=SVD_CUTOFF,
        svdSeed=SVD_SEED,
        fistaLambda=FISTA_LAMBDA,
        # soft consistency check only (RESPONSE.md Section 5)
        target=dict(roi="tubeL", locBeta=1.19, ampBeta=0.94, tol=0.25, window="C3-C7"),
    ),
    # ------------------------------------------------------------------ jun23
    "jun23": dict(
        name="jun23",
        rawRoot=os.path.join(RAW_ROOT, "06-23-raw-data"),
        outDir=os.path.join(OUTPUT_ROOT, "jun23"),
        polCacheDir=os.path.join(OUTPUT_ROOT, "jun23", "pol_cache"),
        cacheFmt="%s_POL.mat",
        blockFmt="%sB%d",                 # UPPERCASE B  -> M1B1 ..
        roiFile=os.path.join(OUTPUT_ROOT, "jun23", "roi_polys_jun23.mat"),
        readoutFile=os.path.join(OUTPUT_ROOT, "jun23", "readout_jun23.mat"),
        betasFile=os.path.join(OUTPUT_ROOT, "jun23", "betas_jun23.mat"),
        nRung=15,
        rungLabels=["L1", "L2", "L3", "L4", "L5",
                    "M1", "M2", "M3", "M4", "M5",
                    "U1", "U2", "U3", "U4", "U5"],
        # LOCKED axis (RESPONSE.md Section 5): L1-U3 nominal, U4/U5 direct Countess reads
        conc=[2.5e5, 3.1e5, 3.9e5, 4.8e5, 6.0e5,
              7.0e5, 1.15e6, 1.9e6, 3.0e6, 5.0e6,
              6.5e6, 9.5e6, 1.4e7, 2.81e7, 3.68e7],
        blocksPerRung=[4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 4, 4, 4, 3, 3],  # sum = 68 rung blocks
        bgFlow=["BGTF1", "BGTF2", "BGTF3", "BGTF4", "BGTF2B1", "BGTF2B2"],  # all 6
        bgStatic=["M2SB1", "M2SB2"],
        nTubes=2,
        fitWindows=[[6, 10], [1, 5], [11, 15], [1, 15]],
        headlineWindow=[6, 10],                      # M1-M5 (cross-cal)
        headlineDomain="PI",
        headlineROI=["tubeL", "tubeR"],              # per-tube for jun23
        domains=list(DOMAINS),
        svdCutoff=SVD_CUTOFF,
        svdSeed=SVD_SEED,
        fistaLambda=FISTA_LAMBDA,
        target=None,                                  # consistency report only
    ),
}


# --------------------------------------------------------------------------- helpers
def get_dataset(name):
    """Return the frozen config dict for a dataset name ('apr17' | 'jun23')."""
    if name not in DATASETS:
        raise KeyError("unknown dataset %r (have %s)" % (name, list(DATASETS)))
    return DATASETS[name]


def rung_blocks(cfg):
    """List-of-lists: rung_blocks(cfg)[r] = ordered block tags for rung r (0-based r).

    Expands the rung+block table (RESPONSE.md Section 5) honoring the variable
    blocks-per-rung. Block index is 1-based in the tag (blockFmt %d)."""
    fmt = cfg["blockFmt"]
    out = []
    for label, nb in zip(cfg["rungLabels"], cfg["blocksPerRung"]):
        out.append([fmt % (label, b) for b in range(1, nb + 1)])
    return out


def rung_block_tags(cfg):
    """Flat, rung-major ordered list of every RUNG block tag (no Bg tags).

    apr17 -> 24 (8 rungs x 3 blocks); jun23 -> 68 (variable blocks-per-rung)."""
    return [t for rung in rung_blocks(cfg) for t in rung]


def all_block_tags(cfg, include_bg_static=True):
    """Full ordered block-tag list the sweep must beamform/infer:
    rung blocks + bgFlow + (optionally) bgStatic.

    apr17 -> 24 + 2 (+1) = 27 ; jun23 -> 68 + 6 (+2) = 76.
    (bgFlow blocks are required for the pedestal; bgStatic is reported-only.)"""
    tags = list(rung_block_tags(cfg))
    tags += list(cfg["bgFlow"])
    if include_bg_static:
        tags += list(cfg["bgStatic"])
    return tags


def pol_cache_path(cfg, tag):
    """Absolute path to a block's polarity cache '<tag>_POL.mat' (Section 3/5)."""
    import os
    return os.path.join(cfg["polCacheDir"], cfg["cacheFmt"] % tag)


def window_label(cfg, window):
    """Human window name from a 1-based [lo hi] rung-index range, e.g. [3,7] -> 'C3-C7'.
    Matches the LABEL-based window names in betas_<ds>.mat (RESPONSE.md Section 6b)."""
    lo, hi = window
    return "%s-%s" % (cfg["rungLabels"][lo - 1], cfg["rungLabels"][hi - 1])


def rung_index_of_tag(cfg, tag):
    """0-based rung index that a rung block tag belongs to, or None for Bg tags."""
    for ri, rung in enumerate(rung_blocks(cfg)):
        if tag in rung:
            return ri
    return None


if __name__ == "__main__":
    for nm, cfg in DATASETS.items():
        rb = rung_block_tags(cfg)
        ab = all_block_tags(cfg)
        print("%-6s nRung=%-2d rungBlocks=%-3d allBlocks(+bg)=%-3d  conc[0..-1]=%g..%g" % (
            nm, cfg["nRung"], len(rb), len(ab), cfg["conc"][0], cfg["conc"][-1]))
        print("        first/last rung tag: %s .. %s   bgFlow=%s  bgStatic=%s" % (
            rb[0], rb[-1], cfg["bgFlow"], cfg["bgStatic"]))
