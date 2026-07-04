#!/usr/bin/env python
# SPDX-License-Identifier: MIT
"""
Phase 3: build the FINAL (stopgap) training PSF artifact from the apr17 validation.

  >>> STOPGAP PSF <<<  See phase3/UPGRADE_LATER.md. Decision (Eli, 2026-06-18): sanity-check +
  reuse. The apr17 single-bubble measurement (measure_psf.py, C1-C6) VALIDATED the analytical
  model laterally (measured FWHM_lat 0.232 mm, concentration-independent, vs model 0.245 mm @20mm
  / 0.220 mm @ the tube depth ~18mm). Axial is at the ~2px sampling floor (unmeasurable on this
  grid) -> keep the model's lambda. So the training PSF = analytical model with the lateral FWHM
  lightly refined to the measured value. UPGRADE to a measured LINEAR PSF from the followup later.

Emits phase3/artifacts/psf_training.npz:
  sigma_x_mm, sigma_z_mm  (anisotropic Gaussian; x=lateral, z=axial)  -- the SOURCE OF TRUTH
  fwhm_lat_mm, fwhm_ax_mm
  kernel_native           rendered kernel on the native grid (dx=dz=0.0346875 mm)
  meta (json)             full provenance + the loud upgrade note
The synthetic generator / deconvolution dictionary RENDER the PSF from (sigma_x_mm, sigma_z_mm)
at their own working resolution; the stored native kernel is for convenience / QC only.
"""
import os, json
import numpy as np

DX_MM = DZ_MM = 0.0346875
FWHM_K = 2.0 * np.sqrt(2.0 * np.log(2.0))
HERE = os.path.dirname(__file__)
ART = os.path.normpath(os.path.join(HERE, "..", "artifacts"))

# --- measured (apr17 C1-C6, measure_psf.py pooled) ---
MEAS_FWHM_LAT_MM = 0.2323          # pooled median, conc-independent (C1-C4 cleaner subset ~0.228)
MEAS_FWHM_LAT_IQR_MM = 0.105
# --- analytical model (psf.mat / build_sushi_psf.m) ---
MODEL_FWHM_LAT_MM = 0.24470898     # @ representative depth 20 mm
MODEL_FWHM_AX_MM = 0.069375        # = lambda (c=1480, 21.333 MHz)

# --- FINAL training PSF: measured lateral, model axial (axial undersampled) ---
FWHM_LAT_MM = MEAS_FWHM_LAT_MM
FWHM_AX_MM = MODEL_FWHM_AX_MM
SIGMA_X_MM = FWHM_LAT_MM / FWHM_K   # lateral
SIGMA_Z_MM = FWHM_AX_MM / FWHM_K    # axial


def render_kernel(sigma_z_mm, sigma_x_mm, dz_mm, dx_mm, nsig=4.0):
    sz, sx = sigma_z_mm / dz_mm, sigma_x_mm / dx_mm
    hz, hx = int(np.ceil(nsig * sz)), int(np.ceil(nsig * sx))
    zz, xx = np.mgrid[-hz:hz + 1, -hx:hx + 1]
    k = np.exp(-(zz ** 2) / (2 * sz ** 2) - (xx ** 2) / (2 * sx ** 2))
    return (k / k.sum()).astype(np.float32)


def main():
    kernel = render_kernel(SIGMA_Z_MM, SIGMA_X_MM, DZ_MM, DX_MM)
    meta = dict(
        provenance="STOPGAP — apr17 sanity-check + reuse (Eli decision 2026-06-18). "
                   "Lateral FWHM refined to apr17 single-bubble measurement; axial kept at model lambda "
                   "(axial FWHM ~2px is at the grid sampling floor, not measurable here).",
        upgrade_later="REPLACE with a MEASURED LINEAR PSF from the followup acquisition (raw .vada, "
                      "both PI polarities, bottom rungs C01/C02 + Bg controls). Then regenerate the "
                      "synthetic training set and RETRAIN both nets. See phase3/UPGRADE_LATER.md.",
        domain="PI-sum (nonlinear/contrast), post-SVD envelope — matches what the localizers see at "
               "inference on apr17. NOT a linear-imaging PSF.",
        grid=dict(dx_mm=DX_MM, dz_mm=DZ_MM, nZ=404, nX=193),
        measured=dict(fwhm_lat_mm=MEAS_FWHM_LAT_MM, fwhm_lat_iqr_mm=MEAS_FWHM_LAT_IQR_MM,
                      conc_independent=True, blocks="C1b1..C6b1 x512 frames, tube ROI",
                      n_fits=2572, note="C1-C4 clean ~0.228mm; C5-C6 broaden ~0.241mm (overlap onset)"),
        model=dict(fwhm_lat_mm_at20mm=MODEL_FWHM_LAT_MM, fwhm_ax_mm=MODEL_FWHM_AX_MM,
                   note="diffraction-limited separable Gaussian, build_sushi_psf.m, depth 20mm"),
        final=dict(fwhm_lat_mm=FWHM_LAT_MM, fwhm_ax_mm=FWHM_AX_MM,
                   sigma_x_mm=SIGMA_X_MM, sigma_z_mm=SIGMA_Z_MM,
                   sigma_x_px=SIGMA_X_MM / DX_MM, sigma_z_px=SIGMA_Z_MM / DZ_MM,
                   axis="x=lateral(nX,193), z=axial(nZ,404)"),
    )
    out = os.path.join(ART, "psf_training.npz")
    np.savez(out, sigma_x_mm=SIGMA_X_MM, sigma_z_mm=SIGMA_Z_MM,
             fwhm_lat_mm=FWHM_LAT_MM, fwhm_ax_mm=FWHM_AX_MM,
             kernel_native=kernel, meta=json.dumps(meta))
    print("[psf] FINAL training PSF (STOPGAP):")
    print(f"      lateral: FWHM {FWHM_LAT_MM:.4f} mm  sigma {SIGMA_X_MM:.4f} mm ({SIGMA_X_MM/DX_MM:.3f} px)")
    print(f"      axial  : FWHM {FWHM_AX_MM:.4f} mm  sigma {SIGMA_Z_MM:.4f} mm ({SIGMA_Z_MM/DZ_MM:.3f} px)  [undersampled; =model lambda]")
    print(f"      native kernel shape {kernel.shape}, sum {kernel.sum():.4f}")
    print(f"[psf] saved -> {out}")

    # --- QC figure (internal validation, NOT a thesis figure per feedback_resolvability_framing) ---
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        mf = os.path.join(ART, "psf_measure_C1b1_C2b1_C3b1_C4b1_C5b1_C6b1.npz")
        summ = json.loads(np.load(mf, allow_pickle=True)["summary"].item())
        tags = [s["tag"] for s in summ]; fl = [s["fwhm_lat_mm"] for s in summ]
        fig, (a0, a1) = plt.subplots(1, 2, figsize=(10, 4))
        a0.bar(tags, fl, color="#0B5563"); a0.axhline(FWHM_LAT_MM, color="#1C1C1C", ls="-", label=f"final {FWHM_LAT_MM:.3f}")
        a0.axhline(MODEL_FWHM_LAT_MM, color="#b00", ls="--", label=f"model@20mm {MODEL_FWHM_LAT_MM:.3f}")
        a0.set_ylabel("lateral FWHM (mm)"); a0.set_title("apr17 single-bubble lateral PSF vs C (conc-independent)")
        a0.set_ylim(0, 0.32); a0.legend(fontsize=8)
        a1.imshow(kernel, aspect="auto", cmap="magma",
                  extent=[-kernel.shape[1]//2*DX_MM, kernel.shape[1]//2*DX_MM,
                          kernel.shape[0]//2*DZ_MM, -kernel.shape[0]//2*DZ_MM])
        a1.set_title("final training PSF kernel (native grid)"); a1.set_xlabel("lateral (mm)"); a1.set_ylabel("axial (mm)")
        fig.suptitle("STOPGAP PSF — upgrade to measured linear PSF from followup (UPGRADE_LATER.md)", fontsize=9, color="#b00")
        fig.tight_layout(); figp = os.path.join(ART, "psf_validation_QC.png")
        fig.savefig(figp, dpi=130); print(f"[psf] QC figure -> {figp}")
    except Exception as e:
        print(f"[psf] (QC figure skipped: {e})")


if __name__ == "__main__":
    main()
