"""Doppler / AGC-phase compensation algorithms for multi-subevent CS.

All functions take the raw measurement (y in hop order) and return the
compensated IQ in the same order, ready to be put on the 1 MHz grid and fed
into the IFFT ranging back-end.

Signal model (step n, subevent k, antenna a), tau_n = t_n - t_ref:
    y[a,n] = H_a(f_n) * exp(-j*4*pi*f_n*v*tau_n/c) * exp(j*psi_k) + noise
psi_k is the unknown AGC phase of subevent k (initiator + reflector).
"""
from dataclasses import dataclass

import numpy as np

from .sim_model import C, DF

F_C = 2440e6
K_DOP = 4 * np.pi / C            # round-trip phase per (Hz * m)


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def adjacent_pairs(ch, orders=(1,)):
    """Index pairs (n, m) with ch[m] = ch[n] + o for o in orders.

    Pairs that would straddle the 23..25 hole simply do not exist, so every
    pair has an exactly known frequency spacing o * 1 MHz.
    """
    pos = {int(c): i for i, c in enumerate(ch)}
    n_idx, m_idx, o_idx = [], [], []
    for o in orders:
        for c, i in pos.items():
            j = pos.get(c + o)
            if j is not None:
                n_idx.append(i), m_idx.append(j), o_idx.append(o)
    return np.array(n_idx, int), np.array(m_idx, int), np.array(o_idx, int)


def parabolic_peak(grid, metric):
    k = int(np.argmax(metric))
    if 0 < k < len(grid) - 1:
        a, b, c = metric[k - 1:k + 2]
        den = a - 2 * b + c
        if den < 0:
            return grid[k] + 0.5 * (a - c) / den * (grid[1] - grid[0])
    return grid[k]


def doppler_phasor(f, tau, v):
    """exp(+j 4 pi f v tau / c): removes Doppler AND range migration exactly."""
    return np.exp(1j * K_DOP * f * v * tau)


# --------------------------------------------------------------------------- #
# reference / baselines
# --------------------------------------------------------------------------- #
def genie(meas):
    """Oracle compensation with the true v and true AGC phases."""
    tau = meas.t - meas.t_ref
    th = meas.theta
    agc = th[meas.se] if th.ndim == 1 else th[:, meas.se]
    return meas.y * doppler_phasor(meas.f, tau, meas.v) * np.exp(-1j * agc)


def raw(meas):
    return meas.y


def legacy(meas, v_grid, sel=None):
    """Current single-subevent scheme (as described): phase of IQ only,
    adjacent-channel pairs sorted by frequency, NDFT over the time difference
    with omega = 4 pi f_c v / c, then de-rotate IQ by omega * t.
    Applied blindly across subevents it ignores the AGC jumps."""
    sel = np.arange(meas.y.shape[1]) if sel is None else sel
    y, f, t, ch = meas.y[:, sel], meas.f[sel], meas.t[sel], meas.ch[sel]
    n, m, _ = adjacent_pairs(ch, (1,))
    u = np.exp(1j * np.angle(y))
    z = u[:, m] * np.conj(u[:, n])                        # (A, P)
    dt = t[m] - t[n]
    w = K_DOP * F_C * v_grid                              # omega grid
    E = np.exp(1j * np.outer(w, dt))                      # (G, P)
    metric = np.abs(E @ z.T).sum(1)
    v_hat = parabolic_peak(v_grid, metric)
    y_c = y * np.exp(1j * K_DOP * F_C * v_hat * (t - meas.t_ref))
    return y_c, sel, v_hat


# --------------------------------------------------------------------------- #
# proposed: Joint Doppler & Phase Stitching (JDPS)
# --------------------------------------------------------------------------- #
@dataclass
class JdpsCfg:
    v_max: float = 10.0
    v_step: float = 0.25
    v_guard: float = 1.0         # search beyond +-v_max so true |v| ~ v_max is not an edge hit
    min_v_score: float = 5.0     # velocity accepted if noise-normalised peak score >= this
    orders: tuple = (1, 2)       # adjacent-channel spacings used (MHz)
    n_power_iter: int = 8
    n_refine_iter: int = 3
    per_ant_theta: bool = False  # estimate an AGC phase per antenna path

    @property
    def v_grid(self):
        n = int(round((self.v_max + self.v_guard) / self.v_step))
        return np.arange(-n, n + 1) * self.v_step


def _group_index(se, n, m, o, n_se, orders):
    """Group id for (subevent of lower-freq tone, subevent of higher, order)."""
    oi = np.searchsorted(np.array(orders), o)
    return (oi * n_se + se[n]) * n_se + se[m]


def phase_sync(Z, n_se, orders, n_power_iter=8, n_refine_iter=3):
    """Estimate per-subevent phases from group sums.

    Z : (A, n_orders, K, K) complex, Z[a, o, k, j] = sum over pairs whose
        lower-frequency tone is in subevent k and higher in subevent j.
        Model: Z[a,o,k,j] ~ |.| * exp(j(c[a,o] + psi_j - psi_k)).
    Returns psi (K,) with psi[0] = 0 and c (A, n_orders).
    """
    K = n_se
    diag = np.einsum("aokk->ao", Z)
    c = np.angle(diag)                                   # intercept from intra-subevent pairs
    psi = np.zeros(K)
    if K == 1:
        return psi, c
    for it in range(n_refine_iter):
        # K x K Hermitian "relative phase" matrix, principal eigenvector (power iteration)
        Q = np.einsum("aokj,ao->kj", Z, np.exp(-1j * c))
        Hm = Q + Q.conj().T                              # Hm[k,j] ~ exp(j(psi_j - psi_k))
        # zero diagonal would give eigenvalues +-lambda (K=2 oscillates);
        # shift so the wanted (positive) eigenvector dominates
        np.fill_diagonal(Hm, np.abs(Hm).sum(1).max())
        x = np.exp(-1j * psi)
        for _ in range(n_power_iter):
            x = Hm @ x                                   # Hm = conj(u) u^T -> principal vector conj(u)
            x /= np.abs(x) + 1e-30
        psi = -np.angle(x)
        psi = np.angle(np.exp(1j * (psi - psi[0])))
        # intercept update with all groups
        rot = np.exp(-1j * (psi[None, :] - psi[:, None]))  # exp(-j(psi_j - psi_k)) [k,j]
        c = np.angle(np.einsum("aokj,kj->ao", Z, rot))
    return psi, c


def velocity_score(metric_peak, z, g, n_groups):
    """Noise-normalised peak score (K independent).

    For pure noise each |group sum| is Rayleigh with mean sqrt(pi/4*S_g) and
    variance (1-pi/4)*S_g, S_g = sum |z|^2 of the group. score = (peak - mu) / sigma
    stays below ~4.5 for noise whatever the number of segments.
    Returns None when there is no signal energy at all.
    """
    S = np.stack([np.bincount(g, weights=np.abs(za) ** 2, minlength=n_groups) for za in z])
    sigma = np.sqrt(((1 - np.pi / 4) * S).sum())
    if sigma == 0:
        return None
    mu = np.sqrt(np.pi / 4 * S).sum()
    return (metric_peak - mu) / sigma


def jdps(meas, cfg: JdpsCfg = JdpsCfg()):
    """Returns (compensated y, info).

    info: v_hat (raw search result), v_used (applied: v_hat if valid else 0),
    v_valid, v_score, v_at_edge, psi, metric, v_grid, status ('ok' / 'no_signal').
    With status 'no_signal' (all-zero input, e.g. remote IQ missing) y is returned unchanged.
    """
    y, f, ch, se = meas.y, meas.f, meas.ch, meas.se
    tau = meas.t - meas.t_ref
    K = int(se.max()) + 1
    A = y.shape[0]
    n, m, o = adjacent_pairs(ch, cfg.orders)
    g = _group_index(se, n, m, o, K, cfg.orders)
    n_groups = len(cfg.orders) * K * K

    # ---- stage 1: velocity search (group-wise coherent, non-coherent across groups/antennas)
    z = y[:, m] * np.conj(y[:, n])                                     # (A, P) amplitude-weighted
    alpha = K_DOP * (f[m] * tau[m] - f[n] * tau[n])                   # exact per-pair Doppler rate
    se_mid = np.array([tau[se == k].mean() for k in range(K)])
    alpha_red = alpha - K_DOP * F_C * (se_mid[se[m]] - se_mid[se[n]])  # drop group-common part
    v_grid = cfg.v_grid
    E = np.exp(1j * np.outer(v_grid, alpha_red))                       # (G, P)
    onehot = np.zeros((len(g), n_groups))
    onehot[np.arange(len(g)), g] = 1.0
    metric = np.zeros(len(v_grid))
    for a in range(A):
        metric += np.abs(E @ (z[a][:, None] * onehot)).sum(1)
    v_hat = parabolic_peak(v_grid, metric)

    # ---- reliability: no energy -> give up; weak peak or peak on the grid edge -> v = 0
    score = velocity_score(metric.max(), z, g, n_groups) if len(g) else None
    if score is None:
        K0 = np.zeros((A, K)) if cfg.per_ant_theta else np.zeros(K)
        return y.copy(), dict(status="no_signal", v_hat=np.nan, v_used=0.0, v_valid=False,
                              v_score=0.0, v_at_edge=False, psi=K0, metric=metric, v_grid=v_grid)
    k_peak = int(np.argmax(metric))
    at_edge = k_peak in (0, len(v_grid) - 1)
    v_valid = bool(score >= cfg.min_v_score and not at_edge)
    v_used = v_hat if v_valid else 0.0

    # ---- stage 2: Doppler + range-migration compensation with v_used
    y1 = y * doppler_phasor(f, tau, v_used)

    # ---- stage 3: AGC (+ inter-subevent Doppler residual) phase estimation
    z1 = y1[:, m] * np.conj(y1[:, n])
    Zg = np.zeros((A, n_groups), complex)
    for a in range(A):
        Zg[a] = np.bincount(g, weights=z1[a].real, minlength=n_groups) \
            + 1j * np.bincount(g, weights=z1[a].imag, minlength=n_groups)
    Zg = Zg.reshape(A, len(cfg.orders), K, K)
    if cfg.per_ant_theta:
        psi = np.stack([phase_sync(Zg[a:a + 1], K, cfg.orders,
                                   cfg.n_power_iter, cfg.n_refine_iter)[0] for a in range(A)])
        y2 = y1 * np.exp(-1j * psi[:, se])
    else:
        psi, _ = phase_sync(Zg, K, cfg.orders, cfg.n_power_iter, cfg.n_refine_iter)
        y2 = y1 * np.exp(-1j * psi[se])[None, :]
    return y2, dict(status="ok", v_hat=v_hat, v_used=v_used, v_valid=v_valid, v_score=float(score),
                    v_at_edge=at_edge, psi=psi, metric=metric, v_grid=v_grid)
