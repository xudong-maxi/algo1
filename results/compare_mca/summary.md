## Baseline: K=3, MP_mod, SNR 20 dB, v~U[-10,10]

| method | P50 err [m] | P90 err [m] | >1 m [%] | P90 dev vs genie [m] | dev>0.5 m [%] | P90 phase RMS [rad] | P90 peak loss [dB] | P90 v err [m/s] |
|---|---|---|---|---|---|---|---|---|
| Genie (true v, true AGC) | 0.50 | 1.24 | 16.50 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| Existing motion_correct_alg | 0.49 | 1.34 | 17.50 | 0.21 | 0.50 | 0.81 | 2.43 | 0.71 |
| JDPS, PAIR_ORDER=1 | 0.51 | 1.24 | 16.25 | 0.00 | 0.00 | 0.02 | 0.00 | 0.03 |
| JDPS, PAIR_ORDER=2 (default) | 0.51 | 1.24 | 16.25 | 0.00 | 0.00 | 0.02 | 0.00 | 0.02 |

## sweep K: P90 range error [m] / P90 dev vs genie [m] / P90 v err [m/s]

| K | Genie (true v, true AGC) | Existing motion_correct_alg | JDPS, PAIR_ORDER=1 | JDPS, PAIR_ORDER=2 (default) |
|---|---|---|---|---|
| 1 | 1.27 / 0.00 / 0.00 | 1.29 / 0.03 / 0.11 | 1.27 / 0.00 / 0.01 | 1.27 / 0.00 / 0.01 |
| 2 | 1.23 / 0.00 / 0.00 | 1.21 / 0.09 / 0.45 | 1.23 / 0.00 / 0.02 | 1.23 / 0.00 / 0.02 |
| 3 | 1.24 / 0.00 / 0.00 | 1.23 / 0.22 / 0.75 | 1.24 / 0.00 / 0.03 | 1.23 / 0.00 / 0.02 |
| 4 | 1.17 / 0.00 / 0.00 | 1.27 / 0.34 / 1.03 | 1.16 / 0.01 / 0.04 | 1.16 / 0.01 / 0.03 |

## sweep speed: P90 range error [m] / P90 dev vs genie [m] / P90 v err [m/s]

| speed | Genie (true v, true AGC) | Existing motion_correct_alg | JDPS, PAIR_ORDER=1 | JDPS, PAIR_ORDER=2 (default) |
|---|---|---|---|---|
| 0 | 1.29 / 0.00 / 0.00 | 1.29 / 0.20 / 0.80 | 1.29 / 0.01 / 0.03 | 1.29 / 0.00 / 0.02 |
| 2 | 1.23 / 0.00 / 0.00 | 1.26 / 0.20 / 0.80 | 1.23 / 0.00 / 0.03 | 1.23 / 0.00 / 0.02 |
| 5 | 1.25 / 0.00 / 0.00 | 1.24 / 0.21 / 0.80 | 1.25 / 0.00 / 0.03 | 1.25 / 0.00 / 0.02 |
| 8 | 1.15 / 0.00 / 0.00 | 1.25 / 0.22 / 0.80 | 1.14 / 0.00 / 0.03 | 1.14 / 0.00 / 0.02 |
| 10 | 1.14 / 0.00 / 0.00 | 1.23 / 0.27 / 0.80 | 1.14 / 0.00 / 0.03 | 1.14 / 0.00 / 0.02 |

## sweep snr: P90 range error [m] / P90 dev vs genie [m] / P90 v err [m/s]

| snr | Genie (true v, true AGC) | Existing motion_correct_alg | JDPS, PAIR_ORDER=1 | JDPS, PAIR_ORDER=2 (default) |
|---|---|---|---|---|
| 0 | 1.31 / 0.00 / 0.00 | 12.43 / 12.56 / 15.96 | 13.49 / 13.89 / 8.62 | 13.23 / 13.55 / 8.33 |
| 5 | 1.24 / 0.00 / 0.00 | 1.26 / 0.26 / 11.02 | 1.24 / 0.03 / 0.17 | 1.24 / 0.03 / 0.15 |
| 10 | 1.24 / 0.00 / 0.00 | 1.26 / 0.21 / 0.92 | 1.24 / 0.02 / 0.09 | 1.24 / 0.01 / 0.08 |
| 20 | 1.15 / 0.00 / 0.00 | 1.22 / 0.22 / 0.77 | 1.15 / 0.00 / 0.03 | 1.15 / 0.00 / 0.02 |
| 30 | 1.15 / 0.00 / 0.00 | 1.27 / 0.22 / 0.71 | 1.15 / 0.00 / 0.01 | 1.15 / 0.00 / 0.01 |

## scenarios (K=3 unless noted, SNR 20 dB): P90 range error / P90 dev vs genie / P90 v err

| scenario | Genie (true v, true AGC) | Existing motion_correct_alg | JDPS, PAIR_ORDER=1 | JDPS, PAIR_ORDER=2 (default) |
|---|---|---|---|---|
| LOS | 0.37 / 0.00 / 0.00 | 0.43 / 0.20 / 0.70 | 0.37 / 0.00 / 0.02 | 0.37 / 0.00 / 0.02 |
| MP_severe | 1.55 / 0.00 / 0.00 | 1.59 / 0.25 / 0.73 | 1.55 / 0.01 / 0.03 | 1.55 / 0.00 / 0.02 |
| K=1 (single subevent) | 1.24 / 0.00 / 0.00 | 1.23 / 0.03 / 0.12 | 1.24 / 0.00 / 0.01 | 1.24 / 0.00 / 0.01 |

## Sample log proc 7121 (uncalibrated distance)

| method | distance [m] | speed [m/s] |
|---|---|---|
| raw | 1.472 | - |
| mca | 1.451 | -0.00 / -0.00 / -0.60 / -0.00 |
| jdps_o1 | 1.378 | +0.10 (score 5.8) |
| jdps_o2 | 1.387 | +0.22 (score 9.5) |
