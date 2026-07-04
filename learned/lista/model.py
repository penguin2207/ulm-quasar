#!/usr/bin/env python
# SPDX-License-Identifier: MIT
"""
Deep-unfolded ULM (LISTA) — van Sloun, Cohen, Eldar 2020 (Proc IEEE) Sec IV-C2, eq (17):

    x_{k+1} = T( W1 y + W2 x_k ) ,   init W1 = A^T,  W2 = I - A^T A
    T = smooth sigmoid-based soft-threshold = Zhang 2001 eq (6) [ref 103]   (see LOG_IMPLEMENTATION_DETAILS.md)

  10 layers, 5x5 conv kernels, per-layer trainable W1_k, W2_k, threshold t_k. ~506 parameters.
  NAMING: van Sloun's per-layer "lambda_k" is the THRESHOLD t_k here; Zhang's lambda is the SMOOTHNESS
  (a fixed hyperparameter `lam_smooth`). Single-channel conv unrolling at the SR grid; input pre-upsampled.

  CAVEAT (our system): lateral PSF ~53 SR-px >> van Sloun's; a 5x5/10-layer receptive field is
  borderline laterally. We keep the spec by default (n_layers/ksize configurable); under-resolution
  laterally would itself be a finding (consistent with the "less robust" caveat, review Fig 7d).
"""
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F


def smooth_soft_threshold(x, t, lam):
    """Zhang 2001 eq (6). t=threshold (>=0), lam=smoothness (>0). -> standard soft-threshold as lam->0."""
    return x + 0.5 * (torch.sqrt((x - t) ** 2 + lam) - torch.sqrt((x + t) ** 2 + lam))


class DeepUnfoldedULM(nn.Module):
    def __init__(self, n_layers=10, ksize=5, sr=8, lam_smooth=1e-2, nonneg=True):
        super().__init__()
        self.sr, self.lam, self.nonneg = sr, lam_smooth, nonneg
        pad = ksize // 2
        self.W1 = nn.ModuleList([nn.Conv2d(1, 1, ksize, padding=pad, bias=False) for _ in range(n_layers)])
        self.W2 = nn.ModuleList([nn.Conv2d(1, 1, ksize, padding=pad, bias=False) for _ in range(n_layers)])
        self.t = nn.ParameterList([nn.Parameter(torch.tensor(0.005)) for _ in range(n_layers)])
        self._init(ksize)

    def _init(self, ksize):
        c = ksize // 2
        for w in self.W1:                       # W1 ~ A^T : centered positive (stronger init so signal
            nn.init.zeros_(w.weight); w.weight.data[0, 0, c, c] = 0.3   # propagates through the 10 layers
        for w in self.W2:                       # W2 ~ I - A^T A : start at identity (delta)
            nn.init.zeros_(w.weight); w.weight.data[0, 0, c, c] = 1.0

    def forward(self, y):
        yu = F.interpolate(y, scale_factor=self.sr, mode="nearest")   # pre-upsample to SR grid
        x = torch.zeros_like(yu)
        for k in range(len(self.t)):
            x = smooth_soft_threshold(self.W1[k](yu) + self.W2[k](x), self.t[k].clamp_min(0.0), self.lam)
        return F.relu(x) if self.nonneg else x


if __name__ == "__main__":
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    net = DeepUnfoldedULM().to(dev)
    n = sum(p.numel() for p in net.parameters())
    x = torch.rand(2, 1, 64, 64, device=dev)
    y = net(x)
    print(f"[lista] params = {n} (target ~506); in {tuple(x.shape)} -> out {tuple(y.shape)} on {dev}")
    print(f"[lista] thresholds init {[round(float(p),3) for p in net.t]}")
