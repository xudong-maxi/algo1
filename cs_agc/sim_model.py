"""BLE Channel Sounding (Mode-2 PBR) multi-subevent measurement model.

Generates the per-step round-trip IQ (initiator tone x reflector tone) for a
CS procedure whose 72 channels are randomly hopped and split over K subevents.
Each subevent carries an unknown constant phase caused by the AGC gear chosen
on both sides; the target moves with radial velocity v (Doppler + range
migration); the channel is a sparse multipath model.
"""
from dataclasses import dataclass, field

import numpy as np

C = 299_792_458.0
F0 = 2402e6                      # CS channel 0
DF = 1e6                         # channel spacing
CS_CHANNELS = np.array([ch for ch in range(2, 77) if ch not in (23, 24, 25)])  # 72 ch


@dataclass
class TimingCfg:
    t_meas: float = 565e-6           # measurement duration of one step
    t_gap: float = 150e-6            # gap between steps
    se_gap: float = 40e-3            # idle gap between end of one subevent and start of the next

    @property
    def t_step(self):
        return self.t_meas + self.t_gap


@dataclass
class ChannelCfg:
    """Sparse multipath channel (one-way)."""
    n_nlos: int = 4                  # number of NLOS paths
    k_factor_db: float = 3.0         # LOS power / total NLOS power
    mean_excess_m: float = 4.0       # mean excess path length (exponential)
    decay_m: float = 6.0             # NLOS power ~ exp(-excess / decay_m)
    per_path_velocity: bool = False  # NLOS radial velocity = v*cos(beta)


# Named multipath scenarios used by the simulation
SCENARIOS = {
    "LOS":      ChannelCfg(n_nlos=2, k_factor_db=10.0, mean_excess_m=3.0),
    "MP_mod":   ChannelCfg(n_nlos=4, k_factor_db=3.0, mean_excess_m=4.0),
    "MP_severe": ChannelCfg(n_nlos=6, k_factor_db=0.0, mean_excess_m=3.0),
}


@dataclass
class SimCfg:
    n_se: int = 3                    # number of subevents
    n_ant: int = 4                   # antenna paths
    snr_db: float = 20.0             # per one-way tone SNR
    d_range: tuple = (1.0, 20.0)
    v_range: tuple = (-10.0, 10.0)
    agc_shared_across_ant: bool = True
    timing: TimingCfg = field(default_factory=TimingCfg)
    channel: ChannelCfg = field(default_factory=lambda: SCENARIOS["MP_mod"])


@dataclass
class Measurement:
    y: np.ndarray        # (A, N) complex, PBR product IQ, in hop (time) order
    f: np.ndarray        # (N,) Hz
    ch: np.ndarray       # (N,) channel index
    t: np.ndarray        # (N,) s, step centre time (absolute)
    se: np.ndarray       # (N,) subevent id
    t_ref: float         # reference time the distance is reported at
    # ground truth
    d0: float            # LOS distance at t_ref
    v: float             # LOS radial velocity
    theta: np.ndarray    # (K,) or (A, K) AGC phase per subevent
    h_clean: np.ndarray  # (A, N) noiseless static round-trip CFR at t_ref (no AGC)


def split_steps(n_ch, n_se):
    """Near-equal step counts per subevent, e.g. 72 over 3 -> 24/24/24, 72 over 5 -> 15/15/14/14/14."""
    return [len(a) for a in np.array_split(np.arange(n_ch), n_se)]


def draw_paths(rng, cc: ChannelCfg, d0, v):
    excess = np.concatenate([[0.0], np.sort(rng.exponential(cc.mean_excess_m, cc.n_nlos))])
    p_nlos = np.exp(-excess[1:] / cc.decay_m)
    p_nlos = p_nlos / p_nlos.sum() * 10 ** (-cc.k_factor_db / 10) if cc.n_nlos else p_nlos
    amp = np.sqrt(np.concatenate([[1.0], p_nlos]))
    if cc.per_path_velocity:
        vel = v * np.cos(np.concatenate([[0.0], rng.uniform(0, np.pi, cc.n_nlos)]))
    else:
        vel = np.full(cc.n_nlos + 1, v)
    return d0 + excess, amp, vel


def generate(rng, cfg: SimCfg, d0=None, v=None) -> Measurement:
    tm, cc = cfg.timing, cfg.channel
    d0 = rng.uniform(*cfg.d_range) if d0 is None else d0
    v = rng.uniform(*cfg.v_range) if v is None else v

    # random hopping of the 72 channels, split over the subevents
    ch = rng.permutation(CS_CHANNELS)
    counts = split_steps(len(ch), cfg.n_se)
    se = np.repeat(np.arange(cfg.n_se), counts)
    idx_in_se = np.concatenate([np.arange(c) for c in counts])
    # subevent duration = c steps (last step has no trailing inter-step gap), then se_gap idle
    se_dur = np.array(counts) * tm.t_step - tm.t_gap
    se_start = np.concatenate([[0.0], np.cumsum(se_dur + tm.se_gap)[:-1]])
    t = se_start[se] + idx_in_se * tm.t_step + tm.t_meas / 2
    t_ref = t.mean()
    f = F0 + ch * DF

    d_p, amp, vel = draw_paths(rng, cc, d0, v)
    A = cfg.n_ant
    ph = rng.uniform(0, 2 * np.pi, (A, len(d_p)))          # per-antenna path phases
    g = amp[None, :] * np.exp(1j * ph)                        # (A, P)
    # one-way CFR at the time of each step: (A, N)
    dist = d_p[None, :] + vel[None, :] * (t - t_ref)[:, None]          # (N, P)
    h = np.einsum("ap,np->an", g, np.exp(-2j * np.pi * f[:, None] * dist / C))
    dist0 = np.broadcast_to(d_p, (len(f), len(d_p)))
    h0 = np.einsum("ap,np->an", g, np.exp(-2j * np.pi * f[:, None] * dist0 / C))

    # AGC phase: each side contributes an unknown constant per subevent
    if cfg.agc_shared_across_ant:
        theta = rng.uniform(0, 2 * np.pi, (2, cfg.n_se))
        th_i, th_r = theta[0][se][None, :], theta[1][se][None, :]
        theta_tot = np.angle(np.exp(1j * theta.sum(0)))
    else:
        theta = rng.uniform(0, 2 * np.pi, (2, A, cfg.n_se))
        th_i, th_r = theta[0][:, se], theta[1][:, se]
        theta_tot = np.angle(np.exp(1j * theta.sum(0)))

    p_sig = np.mean(np.abs(h) ** 2)
    sig_n = np.sqrt(p_sig / 10 ** (cfg.snr_db / 10) / 2)
    n_i = sig_n * (rng.standard_normal(h.shape) + 1j * rng.standard_normal(h.shape))
    n_r = sig_n * (rng.standard_normal(h.shape) + 1j * rng.standard_normal(h.shape))
    y = (h * np.exp(1j * th_i) + n_i) * (h * np.exp(1j * th_r) + n_r)

    return Measurement(y=y, f=f, ch=ch, t=t, se=se, t_ref=t_ref,
                       d0=d0, v=v, theta=theta_tot, h_clean=h0 ** 2)
