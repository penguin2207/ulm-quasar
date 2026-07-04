#!/usr/bin/env python
# SPDX-License-Identifier: MIT
"""
Deep-ULM — van Sloun et al. 2021 (TMI 40(3):829-839) Sec II.B.

Fully-convolutional encoder-decoder (U-Net-based), leaky-ReLU everywhere, batchnorm before
activations, x8 net upsampling, single-channel LINEAR output (a super-res intensity map).
  - Encoder: 3 blocks, each = two 3x3 conv (+BN+leakyReLU) then 2x2 max-pool.
  - Latent : two 3x3 conv (+BN+leakyReLU) + dropout(0.5).
  - Decoder: 3 blocks, each upsamples x4 (ConvTranspose stride2 + nearest x2) with 5x5 convs;
             final 5x5 conv -> 1 channel, linear. Net: input dims x8 (64 -> 512).
Encoder downsamples x8 and the decoder upsamples x64, so this is the asymmetric Deep-STORM/SegNet-style
encoder-decoder (NO symmetric U-Net skip connections -- the x8 output size precludes them); the U-Net
reference [25] is for the encoder-decoder idea. Exact channel widths are not in the paper text; `base`
is tuned to land near the paper's ~700k params while fitting 8 GB. Loss (eq 2) lives in the trainer.
"""
import torch
import torch.nn as nn
import torch.nn.functional as F


def cba(ci, co, k=3):
    return nn.Sequential(nn.Conv2d(ci, co, k, padding=k // 2),
                         nn.BatchNorm2d(co), nn.LeakyReLU(0.1, inplace=True))


class EncBlock(nn.Module):
    def __init__(self, ci, co):
        super().__init__()
        self.body = nn.Sequential(cba(ci, co), cba(co, co))
        self.pool = nn.MaxPool2d(2)

    def forward(self, x):
        return self.pool(self.body(x))


class DecBlock(nn.Module):
    """Upsample x4 = ConvTranspose(stride2) [x2] + nearest [x2], with 5x5 convs."""
    def __init__(self, ci, co):
        super().__init__()
        self.up = nn.ConvTranspose2d(ci, co, 5, stride=2, padding=2, output_padding=1)
        self.bn = nn.BatchNorm2d(co)
        self.act = nn.LeakyReLU(0.1, inplace=True)
        self.conv = cba(co, co, k=5)

    def forward(self, x):
        x = self.act(self.bn(self.up(x)))                       # x2
        x = F.interpolate(x, scale_factor=2, mode="nearest")    # x2 -> x4 total
        return self.conv(x)


class DeepULM(nn.Module):
    def __init__(self, base=28, sr=8, dropout=0.5):
        super().__init__()
        assert sr == 8
        self.enc1 = EncBlock(1, base)
        self.enc2 = EncBlock(base, 2 * base)
        self.enc3 = EncBlock(2 * base, 4 * base)
        self.latent = nn.Sequential(cba(4 * base, 4 * base), cba(4 * base, 4 * base), nn.Dropout2d(dropout))
        self.dec1 = DecBlock(4 * base, 2 * base)
        self.dec2 = DecBlock(2 * base, base)
        self.dec3 = DecBlock(base, base // 2)
        self.head = nn.Conv2d(base // 2, 1, 5, padding=2)        # linear output (Deep-STORM)

    def forward(self, x):
        x = self.enc1(x); x = self.enc2(x); x = self.enc3(x)     # /8
        x = self.latent(x)
        x = self.dec1(x); x = self.dec2(x); x = self.dec3(x)     # x64 -> input x8
        return self.head(x)


if __name__ == "__main__":
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    net = DeepULM().to(dev)
    n = sum(p.numel() for p in net.parameters())
    x = torch.rand(2, 1, 64, 64, device=dev)
    with torch.no_grad():
        y = net(x)
    print(f"[deepulm] params = {n:,} (target ~700k); in {tuple(x.shape)} -> out {tuple(y.shape)} on {dev}")
