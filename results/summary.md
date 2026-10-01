## Baseline: K=3, MP_mod, SNR 20 dB, v~U[-10,10] m/s

| method | P50 err [m] | P90 err [m] | P95 err [m] | >1 m [%] | P90 dev vs genie [m] | dev>0.5 m [%] | P90 phase RMS [rad] | P90 peak loss [dB] | P90 v err [m/s] |
|---|---|---|---|---|---|---|---|---|---|
| Genie (true v, true AGC) | 0.47 | 1.29 | 1.62 | 18.25 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| No compensation | 5.22 | 14.60 | 16.73 | 82.00 | 14.83 | 85.50 | 1.78 | 10.48 | - |
| True Doppler, AGC not handled | 0.77 | 8.53 | 12.42 | 39.50 | 8.72 | 31.50 | 1.66 | 9.53 | 0.00 |
| Legacy, 1st subevent only | 0.58 | 1.58 | 2.48 | 22.25 | 0.55 | 17.25 | - | - | 0.14 |
| Legacy per subevent + stitch | 0.82 | 7.97 | 12.12 | 42.50 | 8.46 | 41.25 | 1.66 | 9.81 | - |
| Legacy over whole procedure | 0.47 | 1.35 | 1.70 | 19.25 | 0.21 | 1.00 | 0.78 | 2.51 | 0.78 |
| Proposed JDPS | 0.47 | 1.29 | 1.62 | 18.25 | 0.00 | 0.00 | 0.02 | 0.00 | 0.02 |

## speed_sweep: P90 range error [m] / P90 deviation from genie [m]

| |v| [m/s] | Genie (true v, true AGC) | Legacy per subevent + stitch | Legacy over whole procedure | Proposed JDPS |
|---|---|---|---|---|
| 0 | 1.33 / 0.00 | 9.51 / 9.97 | 1.38 / 0.21 | 1.33 / 0.00 |
| 2 | 1.25 / 0.00 | 9.55 / 9.89 | 1.24 / 0.21 | 1.25 / 0.00 |
| 5 | 1.26 / 0.00 | 9.13 / 9.35 | 1.31 / 0.22 | 1.26 / 0.00 |
| 8 | 1.26 / 0.00 | 10.23 / 10.56 | 1.32 / 0.23 | 1.26 / 0.00 |
| 10 | 1.18 / 0.00 | 10.18 / 10.60 | 1.29 / 0.25 | 1.19 / 0.00 |

## snr_sweep: P90 range error [m] / P90 deviation from genie [m]

| SNR [dB] | Genie (true v, true AGC) | Legacy per subevent + stitch | Legacy over whole procedure | Proposed JDPS |
|---|---|---|---|---|
| 5 | 1.10 / 0.00 | 11.80 / 12.24 | 1.12 / 0.24 | 1.11 / 0.03 |
| 10 | 1.12 / 0.00 | 9.94 / 10.30 | 1.17 / 0.23 | 1.12 / 0.01 |
| 15 | 1.12 / 0.00 | 9.11 / 9.67 | 1.15 / 0.22 | 1.12 / 0.01 |
| 20 | 1.14 / 0.00 | 10.42 / 10.15 | 1.20 / 0.22 | 1.14 / 0.00 |
| 25 | 1.19 / 0.00 | 7.67 / 8.00 | 1.25 / 0.21 | 1.19 / 0.00 |
| 30 | 1.19 / 0.00 | 7.87 / 8.60 | 1.27 / 0.20 | 1.19 / 0.00 |

## K_sweep: P90 range error [m] / P90 deviation from genie [m]

| K | Genie (true v, true AGC) | Legacy per subevent + stitch | Legacy over whole procedure | Proposed JDPS |
|---|---|---|---|---|
| 1 | 1.24 / 0.00 | 1.25 / 0.01 | 1.25 / 0.01 | 1.24 / 0.00 |
| 2 | 1.29 / 0.00 | 7.56 / 8.04 | 1.30 / 0.09 | 1.29 / 0.00 |
| 3 | 1.35 / 0.00 | 9.55 / 10.12 | 1.30 / 0.23 | 1.35 / 0.00 |
| 4 | 1.25 / 0.00 | 11.80 / 11.95 | 1.38 / 0.39 | 1.24 / 0.01 |
| 5 | 1.24 / 0.00 | 10.83 / 11.23 | 1.65 / 0.84 | 1.24 / 0.01 |

## Scenarios (K=3, SNR 20 dB): P90 range error [m] / P90 dev vs genie [m] / P90 v err [m/s]

| scenario | Genie (true v, true AGC) | Legacy per subevent + stitch | Legacy over whole procedure | Proposed JDPS |
|---|---|---|---|---|
| LOS | 0.37 / 0.00 / 0.00 | 8.14 / 8.26 / - | 0.43 / 0.20 / 0.72 | 0.37 / 0.00 / 0.02 |
| MP_mod | 1.28 / 0.00 / 0.00 | 9.68 / 10.29 / - | 1.28 / 0.21 / 0.75 | 1.28 / 0.00 / 0.03 |
| MP_severe | 1.53 / 0.00 / 0.00 | 9.27 / 9.77 / - | 1.57 / 0.24 / 0.74 | 1.53 / 0.00 / 0.02 |
| MP_mod + per-path velocity | 0.51 / 0.00 / 0.00 | 11.35 / 11.04 / - | 1.24 / 1.20 / 3.03 | 1.15 / 1.06 / 4.27 |
| MP_mod + 100 ms SE gap | 1.18 / 0.00 / 0.00 | 7.64 / 8.25 / - | 1.48 / 0.77 / 1.03 | 1.18 / 0.00 / 0.02 |
| MP_mod + AGC per antenna | 1.16 / 0.00 / 0.00 | 5.65 / 5.92 / - | 1.23 / 0.46 / 0.61 | 1.18 / 0.27 / 0.03 (per-ant option: 1.16 / 0.00 / 0.03) |

## Complexity knobs (K=3, SNR 10 dB)

| config | P90 err [m] | P90 dev vs genie [m] | P90 phase RMS [rad] | P90 v err [m/s] | est. time 4 ant [ms] |
|---|---|---|---|---|---|
| step0.25_orders12 | 1.19 | 0.01 | 0.06 | 0.08 | 8.56 |
| step0.5_orders12 | 1.19 | 0.01 | 0.06 | 0.08 | 5.03 |
| step1.0_orders12 | 1.20 | 0.01 | 0.06 | 0.09 | 3.26 |
| step0.5_orders1 | 1.19 | 0.02 | 0.07 | 0.09 | 2.64 |

## Estimated MCU load (default config, 4 antennas, 128 MHz)

| K | pairs | v grid | groups | cycles | time [ms] |
|---|---|---|---|---|---|
| 1 | 138 | 89 | 2 | 978224 | 7.64 |
| 2 | 138 | 89 | 8 | 1022096 | 7.99 |
| 3 | 138 | 89 | 18 | 1095216 | 8.56 |
| 5 | 138 | 89 | 50 | 1329200 | 10.38 |

Breakdown for K=3:

| stage | cycles |
|---|---|
| pass1: pair products | 4416 |
| pass1: rotator init (2 sincos / pair) | 71760 |
| pass1: rotator recursion | 393024 |
| pass1: group accumulation | 393024 |
| pass1: |Z| per group | 128160 |
| pass2: Doppler rot (2 sincos / pair) + pair products | 79488 |
| phase sync (K x K power iteration) | 3456 |
| pass3: Doppler + AGC de-rotation | 21888 |
