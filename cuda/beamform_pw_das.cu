/*
 * beamform_pw_das.cu
 * ------------------
 * Production MEX-CUDA plane-wave delay-and-sum beamformer for ULM.
 *
 * Two-stage GPU pipeline:
 *   Stage 1: Batch cuFFT Hilbert transform (real RF -> complex IQ)
 *   Stage 2: DAS with Catmull-Rom cubic interpolation + f-number apodization
 *
 * MATLAB signature (R2018a+ interleaved complex):
 *   bf_iq = beamform_pw_das(rf_data, elem_pos_x, grid_x, grid_z, ...
 *                           tx_angle, c, fs, t0, fnum)
 *
 * Inputs:
 *   rf_data     [nSamples x nChannels x nFrames] real single (raw RF)
 *   elem_pos_x  [nChannels x 1] single, element x-positions [m]
 *   grid_x      [nX x 1] single, lateral pixel positions [m]
 *   grid_z      [nZ x 1] single, axial pixel positions [m]
 *   tx_angle    scalar double, steering angle [radians]
 *   c           scalar double, speed of sound [m/s]
 *   fs          scalar double, sampling frequency [Hz]
 *   t0          scalar double, first sample time [s]
 *   fnum        scalar double, f-number for apodization (<=0 to disable)
 *
 * Output:
 *   bf_iq       [nZ x nX x nFrames] complex single (beamformed IQ)
 *
 * Compile:
 *   mexcuda -R2018a NVCCFLAGS="$NVCCFLAGS --allow-unsupported-compiler" beamform_pw_das.cu -lcufft
 */

/* Bypass MSVC STL version check when CUDA toolkit is older than what
 * the installed VS headers expect (e.g. MATLAB's bundled CUDA 12.2
 * with VS 2022 17.12+).  Must appear before any STL includes. */
#define _ALLOW_COMPILER_AND_STL_VERSION_MISMATCH

#include "mex.h"
#include "gpu/mxGPUArray.h"
#include <cuda_runtime.h>
#include <cufft.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

/* ====================================================================
 *  Error handling -- goto-cleanup pattern for leak-free error exits.
 *  mexErrMsgIdAndTxt does a longjmp, so we must free GPU resources
 *  BEFORE calling it.  We store the message, jump to cleanup, free
 *  everything, then report.
 * ==================================================================== */

#define GPU_CHECK(call) do { \
    cudaError_t _e = (call); \
    if (_e != cudaSuccess) { \
        snprintf(error_msg, sizeof(error_msg), \
            "CUDA error at line %d: %s", __LINE__, \
            cudaGetErrorString(_e)); \
        goto cleanup; \
    } \
} while(0)

#define FFT_CHECK(call) do { \
    cufftResult _e = (call); \
    if (_e != CUFFT_SUCCESS) { \
        snprintf(error_msg, sizeof(error_msg), \
            "cuFFT error at line %d: code %d", __LINE__, (int)_e); \
        goto cleanup; \
    } \
} while(0)

#define KERNEL_CHECK() do { \
    cudaError_t _e = cudaGetLastError(); \
    if (_e != cudaSuccess) { \
        snprintf(error_msg, sizeof(error_msg), \
            "Kernel launch error at line %d: %s", __LINE__, \
            cudaGetErrorString(_e)); \
        goto cleanup; \
    } \
} while(0)

/* ====================================================================
 *  Constants
 * ==================================================================== */

/* Max elements in constant memory: 64 KB / 4 bytes = 16384 floats.
 * UHF29x has 256 elements, L38xp has 128.  Plenty of headroom. */
#define MAX_ELEM 16384

__constant__ float d_elem_x[MAX_ELEM];

/* ====================================================================
 *  Device helper: Catmull-Rom cubic interpolation
 * ====================================================================
 * Given fractional sample index `idx` into complex array `data` of
 * length `n`, returns interpolated complex value.  Uses 4-point
 * support window [i0-1, i0, i0+1, i0+2].  Returns zero if the
 * support window falls outside the array.
 *
 * Catmull-Rom preserves PSF shape better than linear interpolation,
 * which directly impacts ULM localization precision.
 */

__device__ __forceinline__
void cubic_interp(const float2 *data, int n, float idx,
                  float *out_re, float *out_im)
{
    int i0 = __float2int_rd(idx);   /* floor via hardware intrinsic */
    float t = idx - (float)i0;

    if (i0 < 1 || i0 + 2 >= n) {
        *out_re = 0.0f;
        *out_im = 0.0f;
        return;
    }

    float t2 = t * t;
    float t3 = t2 * t;

    float w0 = -0.5f*t3 +       t2 - 0.5f*t;
    float w1 =  1.5f*t3 - 2.5f*t2           + 1.0f;
    float w2 = -1.5f*t3 + 2.0f*t2 + 0.5f*t;
    float w3 =  0.5f*t3 - 0.5f*t2;

    float2 s0 = data[i0 - 1];
    float2 s1 = data[i0];
    float2 s2 = data[i0 + 1];
    float2 s3 = data[i0 + 2];

    *out_re = w0*s0.x + w1*s1.x + w2*s2.x + w3*s3.x;
    *out_im = w0*s0.y + w1*s1.y + w2*s2.y + w3*s3.y;
}

/* ====================================================================
 *  Kernel: real RF -> complex  (set imaginary part to zero)
 * ==================================================================== */

__global__ void kernel_real_to_complex(const float *rf, float2 *iq, int n)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    iq[i].x = rf[i];
    iq[i].y = 0.0f;
}

/* ====================================================================
 *  Kernel: Hilbert multiplier (frequency domain)
 * ====================================================================
 * Applies the analytic-signal filter in-place and folds in 1/N
 * normalization (cuFFT's inverse does not normalize).
 *
 *   DC:         multiply by 1/N
 *   Positive:   multiply by 2/N
 *   Nyquist:    multiply by 1/N  (even N only)
 *   Negative:   multiply by 0
 */

__global__ void kernel_hilbert(float2 *spectrum, int nSamples, int nChannels)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int total = nSamples * nChannels;
    if (i >= total) return;

    int k = i % nSamples;          /* frequency bin within channel */
    int half = nSamples / 2;       /* integer division */
    float inv_n = 1.0f / (float)nSamples;

    float m;
    if (k == 0) {
        m = inv_n;                                      /* DC */
    } else if (nSamples % 2 == 0) {
        /* Even N */
        if      (k <  half) m = 2.0f * inv_n;          /* positive freq */
        else if (k == half) m = inv_n;                  /* Nyquist */
        else                m = 0.0f;                   /* negative freq */
    } else {
        /* Odd N */
        if (k <= half)      m = 2.0f * inv_n;          /* positive freq */
        else                m = 0.0f;                   /* negative freq */
    }

    spectrum[i].x *= m;
    spectrum[i].y *= m;
}

/* ====================================================================
 *  Kernel: plane-wave DAS beamforming
 * ====================================================================
 * One thread per pixel (iz, ix).  Loops over receive channels,
 * computes round-trip delay, interpolates channel IQ, accumulates.
 *
 * Thread mapping:
 *   blockIdx.x * blockDim.x + threadIdx.x  ->  iz (axial)
 *   blockIdx.y * blockDim.y + threadIdx.y  ->  ix (lateral)
 */

__global__ void kernel_das(
    const float2 *ch_iq,        /* [nSamples x nCh] interleaved complex */
    float2       *bf_out,       /* [nZ x nX] interleaved complex */
    const float  *grid_x,       /* [nX] lateral positions [m] */
    const float  *grid_z,       /* [nZ] axial positions [m] */
    int nSamples, int nCh, int nX, int nZ,
    float sin_a, float cos_a,
    float c, float fs, float t0, float fnum)
{
    int iz = blockIdx.x * blockDim.x + threadIdx.x;
    int ix = blockIdx.y * blockDim.y + threadIdx.y;
    if (iz >= nZ || ix >= nX) return;

    float px = grid_x[ix];
    float pz = grid_z[iz];

    /* TX delay: plane-wave projection onto propagation direction */
    float t_tx = (pz * cos_a + px * sin_a) / c;

    /* F-number half-aperture (precompute outside channel loop) */
    float half_ap = (fnum > 0.0f) ? (pz / (2.0f * fnum)) : 1e30f;

    float sum_re = 0.0f;
    float sum_im = 0.0f;

    for (int ch = 0; ch < nCh; ch++) {
        float dx = px - d_elem_x[ch];

        /* F-number apodization: skip elements beyond the aperture */
        if (fabsf(dx) > half_ap) continue;

        /* RX delay: one-way distance from pixel to element */
        float dist = sqrtf(dx * dx + pz * pz);
        float t_rx = dist / c;

        /* Round-trip delay -> fractional sample index (0-based) */
        float sample_idx = (t_tx + t_rx - t0) * fs;

        /* Cubic interpolation into this channel's IQ data */
        float val_re, val_im;
        cubic_interp(ch_iq + ch * nSamples, nSamples,
                     sample_idx, &val_re, &val_im);

        sum_re += val_re;
        sum_im += val_im;
    }

    /* Column-major output to match MATLAB layout */
    int out_idx = iz + ix * nZ;
    bf_out[out_idx].x = sum_re;
    bf_out[out_idx].y = sum_im;
}

/* ====================================================================
 *  MEX gateway
 * ==================================================================== */

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* ------------------------------------------------------------ */
    /*  ALL variable declarations up front (C89 style).              */
    /*  Required because goto-cleanup can't jump over C++ inits.     */
    /* ------------------------------------------------------------ */

    char error_msg[512];
    error_msg[0] = '\0';

    /* GPU resource pointers -- NULL-init for safe cleanup.
     * cudaFree(NULL) is a documented no-op. */
    float  *d_rf  = NULL;
    float2 *d_iq  = NULL;
    float  *d_gx  = NULL;
    float  *d_gz  = NULL;
    float2 *d_bf  = NULL;
    cufftHandle plan = 0;
    int plan_ok = 0;

    /* Parsed from inputs (assigned after validation) */
    const mwSize *dims = NULL;
    int ndim      = 0;
    int nSamples  = 0;
    int nChannels = 0;
    int nFrames   = 0;
    const float *h_elem = NULL;
    int nX = 0, nZ = 0;
    const float *h_gx = NULL;
    const float *h_gz = NULL;
    float tx_angle = 0, c_val = 0, fs_val = 0, t0_val = 0, fnum = 0;
    float sin_a = 0, cos_a = 0;

    /* Sizes computed after parsing */
    size_t frame_n  = 0;
    size_t rf_bytes = 0;
    size_t iq_bytes = 0;
    size_t bf_n     = 0;
    size_t bf_bytes = 0;

    /* Output and launch config */
    mxComplexSingle *h_out = NULL;
    int thr1d = 256;
    int blk1d = 0;
    dim3 das_thr(16, 16);
    dim3 das_blk(1, 1);
    const float *h_rf = NULL;
    int n_fft = 0;
    mwSize out_dims[3];
    int out_ndim = 0;
    int f = 0;

    /* ------------------------------------------------------------ */
    /*  Input validation (safe to use mexErrMsgIdAndTxt directly --  */
    /*  no GPU resources allocated yet)                              */
    /* ------------------------------------------------------------ */

    if (nrhs != 9)
        mexErrMsgIdAndTxt("beamform:nrhs",
            "Nine inputs required: rf_data, elem_pos_x, grid_x, "
            "grid_z, tx_angle, c, fs, t0, fnum.");
    if (nlhs > 1)
        mexErrMsgIdAndTxt("beamform:nlhs",
            "One output required.");

    /* rf_data: must be real single (Hilbert is done internally) */
    if (!mxIsSingle(prhs[0]) || mxIsComplex(prhs[0]))
        mexErrMsgIdAndTxt("beamform:input",
            "rf_data must be real single (Hilbert transform is applied internally).");

    /* ------------------------------------------------------------ */
    /*  Parse dimensions                                             */
    /* ------------------------------------------------------------ */

    dims      = mxGetDimensions(prhs[0]);
    ndim      = (int)mxGetNumberOfDimensions(prhs[0]);
    nSamples  = (int)dims[0];
    nChannels = (ndim > 1) ? (int)dims[1] : 1;
    nFrames   = (ndim > 2) ? (int)dims[2] : 1;

    if (nSamples < 4)
        mexErrMsgIdAndTxt("beamform:input",
            "Need >= 4 samples for cubic interpolation (got %d).", nSamples);
    if (nChannels > MAX_ELEM)
        mexErrMsgIdAndTxt("beamform:input",
            "nChannels (%d) exceeds constant-memory limit (%d).",
            nChannels, MAX_ELEM);

    /* elem_pos_x */
    if (!mxIsSingle(prhs[1]))
        mexErrMsgIdAndTxt("beamform:input",
            "elem_pos_x must be single.");
    if ((int)mxGetNumberOfElements(prhs[1]) != nChannels)
        mexErrMsgIdAndTxt("beamform:input",
            "elem_pos_x has %d elements but rf_data has %d channels.",
            (int)mxGetNumberOfElements(prhs[1]), nChannels);
    h_elem = mxGetSingles(prhs[1]);

    /* grid_x, grid_z */
    if (!mxIsSingle(prhs[2]))
        mexErrMsgIdAndTxt("beamform:input", "grid_x must be single.");
    if (!mxIsSingle(prhs[3]))
        mexErrMsgIdAndTxt("beamform:input", "grid_z must be single.");
    nX = (int)mxGetNumberOfElements(prhs[2]);
    nZ = (int)mxGetNumberOfElements(prhs[3]);
    if (nX < 1 || nZ < 1)
        mexErrMsgIdAndTxt("beamform:input",
            "Grid dimensions must be >= 1 (got nX=%d, nZ=%d).", nX, nZ);
    h_gx = mxGetSingles(prhs[2]);
    h_gz = mxGetSingles(prhs[3]);

    /* Scalar parameters */
    tx_angle = (float)mxGetScalar(prhs[4]);
    c_val    = (float)mxGetScalar(prhs[5]);
    fs_val   = (float)mxGetScalar(prhs[6]);
    t0_val   = (float)mxGetScalar(prhs[7]);
    fnum     = (float)mxGetScalar(prhs[8]);

    if (c_val <= 0.0f)
        mexErrMsgIdAndTxt("beamform:input",
            "Speed of sound must be positive (got %.2f).", (double)c_val);
    if (fs_val <= 0.0f)
        mexErrMsgIdAndTxt("beamform:input",
            "Sampling frequency must be positive (got %.2f).", (double)fs_val);

    sin_a = sinf(tx_angle);
    cos_a = cosf(tx_angle);

    /* ------------------------------------------------------------ */
    /*  Initialize GPU                                               */
    /* ------------------------------------------------------------ */
    mxInitGPU();

    /* ------------------------------------------------------------ */
    /*  Compute sizes                                                */
    /* ------------------------------------------------------------ */

    frame_n  = (size_t)nSamples * nChannels;
    rf_bytes = frame_n * sizeof(float);
    iq_bytes = frame_n * sizeof(float2);
    bf_n     = (size_t)nZ * nX;
    bf_bytes = bf_n * sizeof(float2);

    /* ------------------------------------------------------------ */
    /*  GPU memory allocation (goto-cleanup from here on)            */
    /* ------------------------------------------------------------ */

    GPU_CHECK(cudaMemcpyToSymbol(d_elem_x, h_elem,
                                 nChannels * sizeof(float)));
    GPU_CHECK(cudaMalloc(&d_rf, rf_bytes));
    GPU_CHECK(cudaMalloc(&d_iq, iq_bytes));
    GPU_CHECK(cudaMalloc(&d_gx, nX * sizeof(float)));
    GPU_CHECK(cudaMalloc(&d_gz, nZ * sizeof(float)));
    GPU_CHECK(cudaMalloc(&d_bf, bf_bytes));

    GPU_CHECK(cudaMemcpy(d_gx, h_gx, nX * sizeof(float),
                         cudaMemcpyHostToDevice));
    GPU_CHECK(cudaMemcpy(d_gz, h_gz, nZ * sizeof(float),
                         cudaMemcpyHostToDevice));

    /* ------------------------------------------------------------ */
    /*  cuFFT plan: batch C2C along sample dimension                 */
    /* ------------------------------------------------------------ */
    n_fft = nSamples;
    FFT_CHECK(cufftPlanMany(&plan, 1, &n_fft,
        NULL, 1, nSamples,      /* inembed, istride, idist */
        NULL, 1, nSamples,      /* onembed, ostride, odist */
        CUFFT_C2C, nChannels));
    plan_ok = 1;

    /* ------------------------------------------------------------ */
    /*  Allocate MATLAB output array                                 */
    /* ------------------------------------------------------------ */
    out_dims[0] = (mwSize)nZ;
    out_dims[1] = (mwSize)nX;
    out_dims[2] = (mwSize)nFrames;
    out_ndim = (nFrames > 1) ? 3 : 2;
    plhs[0] = mxCreateNumericArray(out_ndim, out_dims,
                                   mxSINGLE_CLASS, mxCOMPLEX);
    h_out = mxGetComplexSingles(plhs[0]);

    /* ------------------------------------------------------------ */
    /*  Kernel launch configurations                                 */
    /* ------------------------------------------------------------ */

    blk1d = ((int)frame_n + thr1d - 1) / thr1d;

    das_blk.x = (nZ + (int)das_thr.x - 1) / (int)das_thr.x;
    das_blk.y = (nX + (int)das_thr.y - 1) / (int)das_thr.y;

    /* ------------------------------------------------------------ */
    /*  Host input pointer                                           */
    /* ------------------------------------------------------------ */
    h_rf = mxGetSingles(prhs[0]);

    /* ------------------------------------------------------------ */
    /*  Frame processing loop                                        */
    /*                                                                */
    /*  Each iteration:                                               */
    /*    1. H->D copy of real RF                                     */
    /*    2. Real -> complex kernel                                   */
    /*    3. Forward FFT (in-place)                                   */
    /*    4. Hilbert multiplier (in-place, includes 1/N)              */
    /*    5. Inverse FFT (in-place) -> analytic signal                */
    /*    6. DAS kernel -> beamformed IQ                              */
    /*    7. D->H copy of result                                      */
    /*                                                                */
    /*  Per-frame GPU memory: O(nSamples*nChannels + nZ*nX) -- the    */
    /*  same buffers are reused across frames.                        */
    /* ------------------------------------------------------------ */

    for (f = 0; f < nFrames; f++) {

        /* 1. Copy real RF frame to GPU */
        GPU_CHECK(cudaMemcpy(d_rf, h_rf + (size_t)f * frame_n,
                             rf_bytes, cudaMemcpyHostToDevice));

        /* 2. Convert real -> complex (imaginary = 0) */
        kernel_real_to_complex<<<blk1d, thr1d>>>(d_rf, d_iq, (int)frame_n);
        KERNEL_CHECK();

        /* 3. Forward C2C FFT (in-place) */
        FFT_CHECK(cufftExecC2C(plan,
            (cufftComplex *)d_iq, (cufftComplex *)d_iq, CUFFT_FORWARD));

        /* 4. Apply Hilbert multiplier (in-place) */
        kernel_hilbert<<<blk1d, thr1d>>>(d_iq, nSamples, nChannels);
        KERNEL_CHECK();

        /* 5. Inverse C2C FFT (in-place) -> analytic signal */
        FFT_CHECK(cufftExecC2C(plan,
            (cufftComplex *)d_iq, (cufftComplex *)d_iq, CUFFT_INVERSE));

        /* 6. DAS beamforming with f-number apodization */
        kernel_das<<<das_blk, das_thr>>>(
            d_iq, d_bf, d_gx, d_gz,
            nSamples, nChannels, nX, nZ,
            sin_a, cos_a, c_val, fs_val, t0_val, fnum);
        KERNEL_CHECK();

        /* 7. Copy beamformed frame to host output.
         * cudaMemcpy is synchronous on the default stream, so this
         * implicitly waits for all preceding kernels to finish. */
        GPU_CHECK(cudaMemcpy(
            (void *)(h_out + (size_t)f * bf_n),
            d_bf, bf_bytes,
            cudaMemcpyDeviceToHost));
    }

    /* Final sync to catch any lingering async errors */
    GPU_CHECK(cudaDeviceSynchronize());

    /* ------------------------------------------------------------ */
    /*  Cleanup -- always reached, even on error (via goto)          */
    /* ------------------------------------------------------------ */
cleanup:
    if (plan_ok) cufftDestroy(plan);
    cudaFree(d_rf);     /* cudaFree(NULL) is a no-op */
    cudaFree(d_iq);
    cudaFree(d_gx);
    cudaFree(d_gz);
    cudaFree(d_bf);

    /* Report error AFTER freeing GPU resources */
    if (error_msg[0] != '\0')
        mexErrMsgIdAndTxt("beamform:cuda", "%s", error_msg);
}
