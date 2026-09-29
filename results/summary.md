## Baseline: K=3, MP_mod, SNR 20 dB, v~U[-10,10] m/s

| method | P50 err [m] | P90 err [m] | P95 err [m] | >1 m [%] | P90 dev vs genie [m] | dev>0.5 m [%] | P90 phase RMS [rad] | P90 peak loss [dB] | P90 v err [m/s] |
|---|---|---|---|---|---|---|---|---|---|
| Genie (true v, true AGC) | 0.48 | 1.26 | 1.61 | 18.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| No compensation | 5.47 | 15.34 | 17.47 | 80.50 | 15.87 | 85.75 | 1.79 | 10.43 | - |
| True Doppler, AGC not handled | 0.76 | 8.63 | 12.16 | 39.50 | 9.22 | 31.50 | 1.66 | 9.51 | 0.00 |
| Legacy, 1st subevent only | 0.54 | 1.56 | 2.41 | 21.50 | 0.42 | 6.25 | - | - | 0.13 |
| Legacy per subevent + stitch | 0.76 | 8.14 | 11.70 | 38.75 | 8.17 | 36.25 | 1.65 | 9.64 | - |
| Legacy over whole procedure | 0.49 | 1.30 | 1.69 | 19.00 | 0.20 | 0.50 | 0.74 | 2.36 | 0.71 |
| Proposed JDPS | 0.48 | 1.28 | 1.62 | 18.25 | 0.00 | 0.25 | 0.02 | 0.00 | 0.02 |

## speed_sweep: P90 range error [m] / P90 deviation from genie [m]

| |v| [m/s] | Genie (true v, true AGC) | Legacy per subevent + stitch | Legacy over whole procedure | Proposed JDPS |
|---|---|---|---|---|
| 0 | 1.33 / 0.00 | 9.83 / 10.05 | 1.36 / 0.19 | 1.33 / 0.00 |
| 2 | 1.22 / 0.00 | 10.47 / 10.47 | 1.24 / 0.20 | 1.23 / 0.00 |
| 5 | 1.28 / 0.00 | 10.73 / 11.26 | 1.33 / 0.21 | 1.28 / 0.00 |
| 8 | 1.26 / 0.00 | 9.01 / 9.44 | 1.33 / 0.20 | 1.26 / 0.00 |
| 10 | 1.17 / 0.00 | 10.34 / 10.62 | 1.24 / 0.25 | 1.17 / 0.00 |

## snr_sweep: P90 range error [m] / P90 deviation from genie [m]

| SNR [dB] | Genie (true v, true AGC) | Legacy per subevent + stitch | Legacy over whole procedure | Proposed JDPS |
|---|---|---|---|---|
| 5 | 1.12 / 0.00 | 12.16 / 12.66 | 1.20 / 0.22 | 1.12 / 0.03 |
| 10 | 1.13 / 0.00 | 10.20 / 10.73 | 1.17 / 0.22 | 1.13 / 0.02 |
| 15 | 1.11 / 0.00 | 8.60 / 9.58 | 1.14 / 0.20 | 1.12 / 0.01 |
| 20 | 1.13 / 0.00 | 7.35 / 8.06 | 1.18 / 0.18 | 1.14 / 0.00 |
| 25 | 1.18 / 0.00 | 7.29 / 7.97 | 1.27 / 0.18 | 1.18 / 0.00 |
| 30 | 1.19 / 0.00 | 7.67 / 8.40 | 1.29 / 0.18 | 1.19 / 0.00 |

## K_sweep: P90 range error [m] / P90 deviation from genie [m]

| K | Genie (true v, true AGC) | Legacy per subevent + stitch | Legacy over whole procedure | Proposed JDPS |
|---|---|---|---|---|
| 1 | 1.24 / 0.00 | 1.25 / 0.01 | 1.25 / 0.01 | 1.24 / 0.00 |
| 2 | 1.29 / 0.00 | 6.79 / 7.49 | 1.29 / 0.13 | 1.29 / 0.00 |
| 3 | 1.33 / 0.00 | 9.26 / 9.88 | 1.28 / 0.19 | 1.34 / 0.00 |
| 4 | 1.27 / 0.00 | 11.89 / 12.52 | 1.36 / 0.32 | 1.27 / 0.01 |
| 5 | 1.24 / 0.00 | 10.81 / 11.35 | 1.51 / 0.56 | 1.24 / 0.01 |

## Scenarios (K=3, SNR 20 dB): P90 range error [m] / P90 dev vs genie [m] / P90 v err [m/s]

| scenario | Genie (true v, true AGC) | Legacy per subevent + stitch | Legacy over whole procedure | Proposed JDPS |
|---|---|---|---|---|
| LOS | 0.37 / 0.00 / 0.00 | 9.88 / 10.04 / - | 0.42 / 0.18 / 0.69 | 0.37 / 0.00 / 0.02 |
| MP_mod | 1.27 / 0.00 / 0.00 | 9.44 / 9.59 / - | 1.27 / 0.19 / 0.71 | 1.27 / 0.00 / 0.02 |
| MP_severe | 1.50 / 0.00 / 0.00 | 9.35 / 10.14 / - | 1.56 / 0.24 / 0.67 | 1.50 / 0.00 / 0.03 |
| MP_mod + per-path velocity | 0.47 / 0.00 / 0.00 | 10.53 / 10.45 / - | 1.05 / 1.00 / 2.87 | 1.04 / 0.94 / 4.19 |
| MP_mod + 100 ms SE interval | 1.15 / 0.00 / 0.00 | 7.41 / 7.97 / - | 1.45 / 0.49 / 0.63 | 1.15 / 0.00 / 0.02 |
| MP_mod + AGC per antenna | 1.15 / 0.00 / 0.00 | 4.13 / 4.77 / - | 1.18 / 0.44 / 0.66 | 1.19 / 0.23 / 0.03 (per-ant option: 1.16 / 0.01 / 0.03) |

## Complexity knobs (K=3, SNR 10 dB)

| config | P90 err [m] | P90 dev vs genie [m] | P90 phase RMS [rad] | P90 v err [m/s] | est. time 4 ant [ms] |
|---|---|---|---|---|---|
| step0.25_orders12 | 1.18 | 0.01 | 0.06 | 0.07 | 4.75 |
| step0.5_orders12 | 1.18 | 0.01 | 0.06 | 0.08 | 2.57 |
| step1.0_orders12 | 1.18 | 0.02 | 0.07 | 0.10 | 1.48 |
| step0.5_orders1 | 1.18 | 0.02 | 0.07 | 0.09 | 1.34 |

## Estimated MCU load (default config, 4 antennas, 128 MHz)

| K | pairs | v grid | groups | cycles | time [ms] |
|---|---|---|---|---|---|
| 1 | 138 | 81 | 2 | 500640 | 3.91 |
| 2 | 138 | 81 | 8 | 540732 | 4.22 |
| 3 | 138 | 81 | 18 | 607512 | 4.75 |
| 5 | 138 | 81 | 50 | 821136 | 6.42 |

Breakdown for K=3:

| stage | cycles |
|---|---|
| pair products (y_m*conj(y_n)) | 4416 |
| pair phase-rate + start/step phasors | 17940 |
| search: phasor recursion | 89424 |
| search: group accumulation | 357696 |
| search: |Z| per group | 116640 |
| Doppler+migration compensation | 6624 |
| re-pair + group sums | 8832 |
| phase sync (K x K power iteration) | 3456 |
| AGC de-rotation | 2484 |
