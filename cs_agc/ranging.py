"""IFFT ranging back-end applied after stitching (common to every method)."""
import numpy as np

from .sim_model import C, DF

N_FFT = 1024
N_BINS = 79                      # CS channel index 0..78 on a 1 MHz grid
FIRST_PATH_DB = 6.0              # leading-edge threshold below the strongest peak
D_MAX = 60.0                     # search window [m]


def to_grid(y, ch):
    """(A, N) samples in hop order -> (A, 79) 1 MHz grid, zeros where unmeasured."""
    X = np.zeros((y.shape[0], N_BINS), complex)
    X[:, ch] = y
    return X


def delay_profile(X):
    """Non-coherent (over antennas) IFFT power profile; bin -> distance."""
    w = np.hanning(N_BINS + 2)[1:-1]
    P = np.abs(np.fft.ifft(X * w, N_FFT, axis=1)) ** 2
    prof = P.sum(0)
    dist = np.arange(N_FFT) * C / (2 * N_FFT * DF)      # round-trip -> one-way
    return dist, prof


def estimate_distance(y, ch):
    dist, prof = delay_profile(to_grid(y, ch))
    n_max = int(np.searchsorted(dist, D_MAX))
    p = prof[:n_max]
    thr = p.max() * 10 ** (-FIRST_PATH_DB / 10)
    # earliest local maximum above threshold (leading edge / first path)
    k = int(np.argmax(p))
    for i in range(1, n_max - 1):
        if p[i] >= thr and p[i] >= p[i - 1] and p[i] >= p[i + 1]:
            k = i
            break
    # parabolic interpolation on the log profile
    if 0 < k < n_max - 1:
        a, b, c = np.log(p[k - 1:k + 2] + 1e-30)
        den = a - 2 * b + c
        k = k + (0.5 * (a - c) / den if den != 0 else 0.0)
    return k * C / (2 * N_FFT * DF)
