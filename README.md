# E101 Homework 9 – Part 1 of 3
**Due Date:** Friday, November 14th, 2025  
**Course:** E101 – Signals and Systems  
**Assignment:** Pulse Train Spectrum (FFT Analysis)

---

## Overview
This assignment is a **pre-lab design and analysis** for the final project. You’ll analyze a sampled periodic pulse train signal using Fourier methods. Work **with your final project partner** — only one submission per team is required. **Teams of three are not allowed.**

You’ll use MATLAB to:
1. Inspect the time-domain data.
2. Estimate basic parameters of the pulse train.
3. Compute and interpret FFT spectra.
4. Compare results with analytical DTFS predictions.

Before starting, **read `FFT-Reference.pdf`** to review the relevant Fourier analysis concepts.

---

## Files Provided
- `HW9.mat` → contains the sampled signal `x[n]`
- `fdomain.m` → function that computes `[X,f] = fdomain(x,Fs)`  
  Returns the FFT result and corresponding frequency vector (Hz)

