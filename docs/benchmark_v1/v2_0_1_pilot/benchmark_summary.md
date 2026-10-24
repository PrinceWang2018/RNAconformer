# Benchmark summary — run `v2_0_1_pilot`

Engine rnafold3d-v2-foldtree 2.0.1; 8 trajectories and up to 5 models per target and mode; seed 42; ViennaRNA 2.7.2.

C4' RMSD in Å over residues observed in the reference; **top-1** = rank-1 model, **best** = best of the exported models, **cand** = best of all valid trajectory models (sampling ceiling); **null** = median C4' RMSD of 30 random arrangements with the correct (oracle) stems; INF-WC on canonical pairs; ViennaRNA F1 = MFE base-pair F1 vs. the reference nested pairs. A dash means no valid model.

## dev split — sequence only (ViennaRNA 2D)

| target | nt | family | valid traj | top-1 | best | cand | null | eRMSD best | INF top-1 | Vienna F1 | time s |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 2KOC_A | 14 | cUUCGg tetraloop hairpin | 8/8 | 1.8 | 1.8 | 1.8 | 3.3 | 1.01 | 1.00 | 1.00 | 17 |
| 480D_A | 27 | sarcin-ricin domain | 5/8 | 6.2 | 4.4 | 4.4 | 7.2 | 1.40 | 0.77 | 0.77 | 57 |
| 1KXK_A | 70 | group II intron D5/D6 | 8/8 | 10.7 | 10.5 | 9.8 | 9.7 | 1.33 | 0.78 | 0.81 | 115 |
| 1Y26_A | 71 | adenine riboswitch | 7/8 | 13.0 | 13.0 | 13.0 | 15.1 | 1.45 | 0.81 | 0.98 | 118 |
| 1EHZ_A | 76 | tRNA-Phe | 5/8 | 18.7 | 10.8 | 10.8 | 20.6 | 1.33 | 0.50 | 0.98 | 134 |
| 1GID_A | 158 | group I intron P4-P6 domain | 3/8 | 31.3 | 24.3 | 24.3 | 36.8 | 1.55 | 0.79 | 0.83 | 348 |

Targets with ≥1 valid model: 6/6. Median top-1 11.9 Å, median best-of-K 10.7 Å (null median 12.4 Å over the same targets). Top-1 better than null: 5/6; best-of-K better than null: 5/6. Best-of-K ≤ 5 Å: 2/6; ≤ 10 Å: 2/6.

## dev split — oracle 2D (from reference; upper bound)

| target | nt | family | valid traj | top-1 | best | cand | null | eRMSD best | INF top-1 | Vienna F1 | time s |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 2KOC_A | 14 | cUUCGg tetraloop hairpin | 8/8 | 1.8 | 1.8 | 1.8 | 3.3 | 1.01 | 1.00 | 1.00 | 18 |
| 480D_A | 27 | sarcin-ricin domain | 7/8 | 3.5 | 3.5 | 3.5 | 7.2 | 1.35 | 1.00 | 0.77 | 48 |
| 1KXK_A | 70 | group II intron D5/D6 | 7/8 | 12.5 | 10.9 | 10.2 | 9.7 | 1.17 | 0.98 | 0.81 | 118 |
| 1Y26_A | 71 | adenine riboswitch | 7/8 | 15.2 | 15.2 | 15.2 | 15.1 | 1.41 | 0.94 | 0.98 | 131 |
| 1EHZ_A | 76 | tRNA-Phe | 5/8 | 12.4 | 12.4 | 12.4 | 20.6 | 1.34 | 0.95 | 0.98 | 131 |
| 1GID_A | 158 | group I intron P4-P6 domain | 7/8 | 26.5 | 23.4 | 23.4 | 36.8 | 1.43 | 0.97 | 0.83 | 413 |

Targets with ≥1 valid model: 6/6. Median top-1 12.5 Å, median best-of-K 11.7 Å (null median 12.4 Å over the same targets). Top-1 better than null: 4/6; best-of-K better than null: 4/6. Best-of-K ≤ 5 Å: 2/6; ≤ 10 Å: 2/6.

## test split — sequence only (ViennaRNA 2D)

| target | nt | family | valid traj | top-1 | best | cand | null | eRMSD best | INF top-1 | Vienna F1 | time s |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1ANR_A | 29 | HIV-1 TAR | 8/8 | 4.1 | 3.4 | 3.4 | 7.9 | 0.90 | 0.95 | 0.95 | 56 |
| 1RNK_A | 34 | MMTV frameshift pseudoknot | 8/8 | 23.9 | 15.9 | 15.9 | 18.4 | 1.51 | 0.57 | 0.00 | 58 |
| 1F1T_A | 38 | malachite green aptamer | 8/8 | 8.0 | 6.8 | 6.2 | 6.8 | 1.33 | 1.00 | 1.00 | 74 |
| 5NWQ_A | 41 | guanidine-III riboswitch | 8/8 | 19.3 | 17.4 | 17.4 | 19.7 | 1.89 | 0.57 | 0.67 | 70 |
| 4ENC_A | 52 | fluoride riboswitch | 7/8 | 18.9 | 15.3 | 15.3 | 17.7 | 1.82 | 0.37 | 0.43 | 95 |
| 3E5C_A | 53 | SAM-III riboswitch | 7/8 | 17.8 | 12.2 | 12.2 | 14.1 | 1.19 | 0.69 | 0.69 | 99 |
| 4XW7_A | 64 | ZTP riboswitch | 7/8 | 15.7 | 10.5 | 10.5 | 13.7 | 1.39 | 0.76 | 0.76 | 108 |
| 5KPY_A | 71 | 5-hydroxytryptophan aptamer | 8/8 | 21.0 | 19.6 | 19.6 | 23.9 | 1.58 | 0.80 | 0.84 | 123 |
| 5TPY_A | 71 | Zika xrRNA | 7/8 | 21.5 | 17.6 | 17.6 | 16.5 | 1.51 | 0.59 | 0.65 | 107 |
| 2OIU_A | 71 | L1 ligase ribozyme | 7/8 | 19.2 | 13.3 | 6.8 | 18.0 | 1.17 | 0.91 | 0.87 | 128 |
| 2GDI_A | 80 | TPP riboswitch | 7/8 | 21.6 | 18.9 | 18.4 | 17.1 | 1.64 | 0.79 | 0.79 | 179 |
| 3MXH_B | 92 | c-di-GMP-I riboswitch | 7/8 | 26.4 | 25.4 | 25.4 | 22.0 | 1.97 | 0.33 | 0.33 | 209 |
| 4L81_A | 96 | SAM-I/IV riboswitch | 7/8 | 26.2 | 21.9 | 17.3 | 20.6 | 1.46 | 0.82 | 0.93 | 198 |
| 4QLM_A | 125 | ydaO (c-di-AMP) riboswitch | 5/8 | 27.9 | 19.9 | 19.9 | 29.7 | 1.52 | 0.80 | 0.70 | 272 |
| 3V7E_C | 126 | SAM-I riboswitch (YbxF-bound) | 6/8 | 24.0 | 22.7 | 22.5 | 19.1 | 1.23 | 0.65 | 0.94 | 215 |
| 1U9S_A | 161 | RNase P specificity domain | 5/8 | 27.9 | 23.2 | 23.2 | 27.9 | 1.65 | 0.77 | 0.81 | 393 |
| 3D0U_A | 161 | lysine riboswitch | 1/8 | 29.0 | 29.0 | 29.0 | 29.3 | 1.62 | 0.68 | 0.74 | 412 |
| 2QBZ_A | 161 | M-box riboswitch | 1/8 | 21.8 | 21.8 | 21.8 | 25.7 | 1.47 | 0.91 | 0.90 | 376 |
| 4GXY_A | 172 | AdoCbl riboswitch | 5/8 | 32.6 | 23.7 | 23.7 | 29.9 | 1.63 | 0.59 | 0.59 | 390 |
| 4P9R_A | 192 | lariat-capping ribozyme | 4/8 | 28.2 | 28.2 | 28.2 | 29.9 | 1.83 | 0.44 | 0.47 | 446 |

Targets with ≥1 valid model: 20/20. Median top-1 21.7 Å, median best-of-K 19.3 Å (null median 19.4 Å over the same targets). Top-1 better than null: 8/20; best-of-K better than null: 14/20. Best-of-K ≤ 5 Å: 1/20; ≤ 10 Å: 2/20.

## test split — oracle 2D (from reference; upper bound)

| target | nt | family | valid traj | top-1 | best | cand | null | eRMSD best | INF top-1 | Vienna F1 | time s |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1ANR_A | 29 | HIV-1 TAR | 8/8 | 5.5 | 4.2 | 4.1 | 7.9 | 0.96 | 1.00 | 0.95 | 65 |
| 1RNK_A | 34 | MMTV frameshift pseudoknot | 8/8 | 17.6 | 14.3 | 14.3 | 18.4 | 1.45 | 0.74 | 0.00 | 72 |
| 1F1T_A | 38 | malachite green aptamer | 8/8 | 8.0 | 6.8 | 6.2 | 6.8 | 1.34 | 1.00 | 1.00 | 80 |
| 5NWQ_A | 41 | guanidine-III riboswitch | 8/8 | 22.9 | 21.6 | 21.4 | 19.7 | 1.71 | 0.80 | 0.67 | 85 |
| 4ENC_A | 52 | fluoride riboswitch | 8/8 | 18.0 | 17.0 | 15.4 | 17.7 | 1.38 | 0.86 | 0.43 | 98 |
| 3E5C_A | 53 | SAM-III riboswitch | 8/8 | 13.1 | 12.2 | 12.2 | 14.1 | 1.19 | 1.00 | 0.69 | 110 |
| 4XW7_A | 64 | ZTP riboswitch | 7/8 | 12.0 | 8.0 | 8.0 | 13.7 | 1.06 | 1.00 | 0.76 | 109 |
| 5KPY_A | 71 | 5-hydroxytryptophan aptamer | 8/8 | 24.5 | 14.3 | 14.3 | 23.9 | 1.42 | 0.95 | 0.84 | 134 |
| 5TPY_A | 71 | Zika xrRNA | 7/8 | 21.9 | 14.0 | 14.0 | 16.5 | 1.35 | 0.87 | 0.65 | 123 |
| 2OIU_A | 71 | L1 ligase ribozyme | 6/8 | 15.7 | 13.4 | 10.1 | 18.0 | 1.18 | 0.98 | 0.87 | 141 |
| 2GDI_A | 80 | TPP riboswitch | 4/8 | 16.1 | 16.1 | 16.1 | 17.1 | 1.48 | 1.00 | 0.79 | 168 |
| 3MXH_B | 92 | c-di-GMP-I riboswitch | 7/8 | 19.3 | 14.4 | 14.4 | 22.0 | 1.47 | 0.98 | 0.33 | 233 |
| 4L81_A | 96 | SAM-I/IV riboswitch | 5/8 | 31.4 | 19.3 | 19.3 | 20.6 | 1.32 | 0.93 | 0.93 | 217 |
| 4QLM_A | 125 | ydaO (c-di-AMP) riboswitch | 4/8 | 21.7 | 21.7 | 21.7 | 29.7 | 1.39 | 0.88 | 0.70 | 352 |
| 3V7E_C | 126 | SAM-I riboswitch (YbxF-bound) | 3/8 | 24.0 | 13.7 | 13.7 | 19.1 | 1.28 | 0.96 | 0.94 | 307 |
| 1U9S_A | 161 | RNase P specificity domain | 2/8 | 28.4 | 28.4 | 28.4 | 27.9 | 1.51 | 1.00 | 0.81 | 442 |
| 3D0U_A | 161 | lysine riboswitch | 3/8 | 19.4 | 19.4 | 19.4 | 29.3 | 1.36 | 0.93 | 0.74 | 386 |
| 2QBZ_A | 161 | M-box riboswitch | 2/8 | 25.4 | 25.4 | 25.4 | 25.7 | 1.41 | 0.98 | 0.90 | 367 |
| 4GXY_A | 172 | AdoCbl riboswitch | 2/8 | 28.0 | 28.0 | 28.0 | 29.9 | 1.52 | 0.96 | 0.59 | 437 |
| 4P9R_A | 192 | lariat-capping ribozyme | 0/8 | — | — | — | 29.9 | — | — | 0.47 | 506 |

Targets with ≥1 valid model: 19/20. Median top-1 19.4 Å, median best-of-K 14.4 Å (null median 19.1 Å over the same targets). Top-1 better than null: 11/19; best-of-K better than null: 16/19. Best-of-K ≤ 5 Å: 1/20; ≤ 10 Å: 3/20.

## Aggregate (all targets; failures counted as failures)

```json
{
 "independent": {
  "targets": 26,
  "yield": 1.0,
  "median_top1_c4_rmsd": 21.23,
  "median_best_c4_rmsd": 17.48,
  "median_top1_ermsd": 1.59,
  "median_best_ermsd": 1.49,
  "median_top1_inf_wc": 0.77,
  "median_best_inf_wc": 0.8,
  "success_best_c4_le_5A": 0.115,
  "success_best_c4_le_10A": 0.154,
  "success_top1_c4_le_10A": 0.154,
  "success_best_ermsd_le_1": 0.038,
  "median_wall_time_s": 125.4,
  "by_stratum": {
   "s100_149": {
    "targets": 2,
    "with_models": 2,
    "median_best_c4_rmsd": 21.33,
    "median_top1_c4_rmsd": 25.92
   },
   "s150_200": {
    "targets": 6,
    "with_models": 6,
    "median_best_c4_rmsd": 23.99,
    "median_top1_c4_rmsd": 28.62
   },
   "s30_49": {
    "targets": 3,
    "with_models": 3,
    "median_best_c4_rmsd": 15.9,
    "median_top1_c4_rmsd": 19.27
   },
   "s50_99": {
    "targets": 12,
    "with_models": 12,
    "median_best_c4_rmsd": 14.28,
    "median_top1_c4_rmsd": 19.05
   },
   "smoke_lt30": {
    "targets": 3,
    "with_models": 3,
    "median_best_c4_rmsd": 3.41,
    "median_top1_c4_rmsd": 4.1
   }
  }
 },
 "oracle": {
  "targets": 26,
  "yield": 0.962,
  "median_top1_c4_rmsd": 17.96,
  "median_best_c4_rmsd": 14.29,
  "median_top1_ermsd": 1.42,
  "median_best_ermsd": 1.36,
  "median_top1_inf_wc": 0.97,
  "median_best_inf_wc": 0.97,
  "success_best_c4_le_5A": 0.115,
  "success_best_c4_le_10A": 0.192,
  "success_top1_c4_le_10A": 0.154,
  "success_best_ermsd_le_1": 0.038,
  "median_wall_time_s": 132.7,
  "by_stratum": {
   "s100_149": {
    "targets": 2,
    "with_models": 2,
    "median_best_c4_rmsd": 17.71,
    "median_top1_c4_rmsd": 22.85
   },
   "s150_200": {
    "targets": 6,
    "with_models": 5,
    "median_best_c4_rmsd": 25.44,
    "median_top1_c4_rmsd": 26.48
   },
   "s30_49": {
    "targets": 3,
    "with_models": 3,
    "median_best_c4_rmsd": 14.28,
    "median_top1_c4_rmsd": 17.62
   },
   "s50_99": {
    "targets": 12,
    "with_models": 12,
    "median_best_c4_rmsd": 14.17,
    "median_top1_c4_rmsd": 15.94
   },
   "smoke_lt30": {
    "targets": 3,
    "with_models": 3,
    "median_best_c4_rmsd": 3.51,
    "median_top1_c4_rmsd": 3.51
   }
  }
 }
}
```
