# Benchmark summary — run `v2_0_pilot`

Engine rnafold3d-v2-foldtree 2.0.0; 8 trajectories and up to 5 models per target and mode; seed 42; ViennaRNA 2.7.2.

C4' RMSD in Å over residues observed in the reference; **top-1** = rank-1 model, **best** = best of the exported models, **cand** = best of all valid trajectory models (sampling ceiling); **null** = median C4' RMSD of 30 random arrangements with the correct (oracle) stems; INF-WC on canonical pairs; ViennaRNA F1 = MFE base-pair F1 vs. the reference nested pairs. A dash means no valid model.

## dev split — sequence only (ViennaRNA 2D)

| target | nt | family | valid traj | top-1 | best | cand | null | eRMSD best | INF top-1 | Vienna F1 | time s |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 2KOC_A | 14 | cUUCGg tetraloop hairpin | 8/8 | 1.8 | 1.8 | 1.8 | 3.3 | 1.01 | 1.00 | 1.00 | 12 |
| 480D_A | 27 | sarcin-ricin domain | 5/8 | 6.2 | 4.4 | 4.4 | 7.2 | 1.40 | 0.77 | 0.77 | 44 |
| 1KXK_A | 70 | group II intron D5/D6 | 6/8 | 11.4 | 10.4 | 10.4 | 9.7 | 1.31 | 0.78 | 0.81 | 92 |
| 1Y26_A | 71 | adenine riboswitch | 8/8 | 13.0 | 13.0 | 13.0 | 15.1 | 1.42 | 0.81 | 0.98 | 92 |
| 1EHZ_A | 76 | tRNA-Phe | 4/8 | 14.0 | 10.8 | 10.8 | 20.6 | 1.33 | 0.93 | 0.98 | 107 |
| 1GID_A | 158 | group I intron P4-P6 domain | 3/8 | 28.6 | 23.1 | 23.1 | 36.8 | 1.56 | 0.80 | 0.83 | 287 |

Targets with ≥1 valid model: 6/6. Median top-1 12.2 Å, median best-of-K 10.6 Å (null median 12.4 Å over the same targets). Top-1 better than null: 5/6; best-of-K better than null: 5/6. Best-of-K ≤ 5 Å: 2/6; ≤ 10 Å: 2/6.

## dev split — oracle 2D (from reference; upper bound)

| target | nt | family | valid traj | top-1 | best | cand | null | eRMSD best | INF top-1 | Vienna F1 | time s |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 2KOC_A | 14 | cUUCGg tetraloop hairpin | 8/8 | 1.8 | 1.8 | 1.8 | 3.3 | 1.01 | 1.00 | 1.00 | 13 |
| 480D_A | 27 | sarcin-ricin domain | 7/8 | 3.5 | 3.5 | 3.5 | 7.2 | 1.35 | 1.00 | 0.77 | 36 |
| 1KXK_A | 70 | group II intron D5/D6 | 7/8 | 10.8 | 10.8 | 9.9 | 9.7 | 1.18 | 0.98 | 0.81 | 92 |
| 1Y26_A | 71 | adenine riboswitch | 7/8 | 15.7 | 13.5 | 13.5 | 15.1 | 1.43 | 0.94 | 0.98 | 102 |
| 1EHZ_A | 76 | tRNA-Phe | 5/8 | 12.4 | 12.4 | 12.4 | 20.6 | 1.34 | 0.95 | 0.98 | 101 |
| 1GID_A | 158 | group I intron P4-P6 domain | 7/8 | 26.5 | 23.4 | 23.4 | 36.8 | 1.43 | 0.97 | 0.83 | 335 |

Targets with ≥1 valid model: 6/6. Median top-1 11.6 Å, median best-of-K 11.6 Å (null median 12.4 Å over the same targets). Top-1 better than null: 4/6; best-of-K better than null: 5/6. Best-of-K ≤ 5 Å: 2/6; ≤ 10 Å: 2/6.

## test split — sequence only (ViennaRNA 2D)

| target | nt | family | valid traj | top-1 | best | cand | null | eRMSD best | INF top-1 | Vienna F1 | time s |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1ANR_A | 29 | HIV-1 TAR | 8/8 | 4.1 | 3.4 | 3.4 | 7.9 | 0.90 | 0.95 | 0.95 | 40 |
| 1RNK_A | 34 | MMTV frameshift pseudoknot | 8/8 | 23.9 | 15.9 | 15.9 | 18.4 | 1.51 | 0.57 | 0.00 | 49 |
| 1F1T_A | 38 | malachite green aptamer | 8/8 | 8.0 | 6.8 | 6.2 | 6.8 | 1.33 | 1.00 | 1.00 | 60 |
| 5NWQ_A | 41 | guanidine-III riboswitch | 8/8 | 19.3 | 17.4 | 17.4 | 19.7 | 1.89 | 0.57 | 0.67 | 56 |
| 4ENC_A | 52 | fluoride riboswitch | 7/8 | 18.9 | 15.3 | 15.3 | 17.7 | 1.82 | 0.37 | 0.43 | 76 |
| 3E5C_A | 53 | SAM-III riboswitch | 7/8 | 17.8 | 12.2 | 12.2 | 14.1 | 1.19 | 0.69 | 0.69 | 80 |
| 4XW7_A | 64 | ZTP riboswitch | 6/8 | 14.2 | 7.4 | 7.4 | 13.7 | 1.39 | 0.76 | 0.76 | 89 |
| 5KPY_A | 71 | 5-hydroxytryptophan aptamer | 8/8 | 21.0 | 19.6 | 19.6 | 23.9 | 1.58 | 0.80 | 0.84 | 99 |
| 5TPY_A | 71 | Zika xrRNA | 6/8 | 20.5 | 20.5 | 18.8 | 16.5 | 1.53 | 0.76 | 0.65 | 80 |
| 2OIU_A | 71 | L1 ligase ribozyme | 7/8 | 19.2 | 13.3 | 6.8 | 18.0 | 1.17 | 0.91 | 0.87 | 105 |
| 2GDI_A | 80 | TPP riboswitch | 7/8 | 21.6 | 18.9 | 18.9 | 17.1 | 1.64 | 0.79 | 0.79 | 146 |
| 3MXH_B | 92 | c-di-GMP-I riboswitch | 0/8 | — | — | — | 22.0 | — | — | 0.33 | 168 |
| 4L81_A | 96 | SAM-I/IV riboswitch | 7/8 | 26.2 | 21.9 | 17.3 | 20.6 | 1.46 | 0.82 | 0.93 | 154 |
| 4QLM_A | 125 | ydaO (c-di-AMP) riboswitch | 7/8 | 30.6 | 18.5 | 18.5 | 29.7 | 1.49 | 0.68 | 0.70 | 267 |
| 3V7E_C | 126 | SAM-I riboswitch (YbxF-bound) | 8/8 | 21.7 | 16.6 | 16.6 | 19.1 | 1.23 | 0.65 | 0.94 | 202 |
| 1U9S_A | 161 | RNase P specificity domain | 2/8 | 28.1 | 28.1 | 28.1 | 27.9 | 1.66 | 0.77 | 0.81 | 293 |
| 3D0U_A | 161 | lysine riboswitch | 5/8 | 25.0 | 25.0 | 25.0 | 29.3 | 1.58 | 0.73 | 0.74 | 316 |
| 2QBZ_A | 161 | M-box riboswitch | 4/8 | 17.6 | 17.6 | 17.6 | 25.7 | 1.44 | 0.91 | 0.90 | 269 |
| 4GXY_A | 172 | AdoCbl riboswitch | 1/8 | 32.2 | 32.2 | 32.2 | 29.9 | 1.64 | 0.60 | 0.59 | 268 |
| 4P9R_A | 192 | lariat-capping ribozyme | 4/8 | 32.1 | 29.2 | 29.2 | 29.9 | 1.83 | 0.44 | 0.47 | 365 |

Targets with ≥1 valid model: 19/20. Median top-1 21.0 Å, median best-of-K 17.6 Å (null median 19.1 Å over the same targets). Top-1 better than null: 5/19; best-of-K better than null: 13/19. Best-of-K ≤ 5 Å: 1/20; ≤ 10 Å: 3/20.

## test split — oracle 2D (from reference; upper bound)

| target | nt | family | valid traj | top-1 | best | cand | null | eRMSD best | INF top-1 | Vienna F1 | time s |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1ANR_A | 29 | HIV-1 TAR | 8/8 | 5.5 | 4.2 | 4.1 | 7.9 | 0.96 | 1.00 | 0.95 | 52 |
| 1RNK_A | 34 | MMTV frameshift pseudoknot | 8/8 | 14.0 | 11.7 | 11.7 | 18.4 | 1.42 | 0.74 | 0.00 | 58 |
| 1F1T_A | 38 | malachite green aptamer | 8/8 | 8.0 | 6.8 | 6.2 | 6.8 | 1.34 | 1.00 | 1.00 | 61 |
| 5NWQ_A | 41 | guanidine-III riboswitch | 8/8 | 22.9 | 21.6 | 21.4 | 19.7 | 1.71 | 0.80 | 0.67 | 67 |
| 4ENC_A | 52 | fluoride riboswitch | 8/8 | 18.0 | 17.0 | 15.4 | 17.7 | 1.38 | 0.86 | 0.43 | 75 |
| 3E5C_A | 53 | SAM-III riboswitch | 8/8 | 13.1 | 12.2 | 12.2 | 14.1 | 1.19 | 1.00 | 0.69 | 86 |
| 4XW7_A | 64 | ZTP riboswitch | 8/8 | 14.7 | 8.3 | 7.9 | 13.7 | 1.04 | 1.00 | 0.76 | 88 |
| 5KPY_A | 71 | 5-hydroxytryptophan aptamer | 8/8 | 24.5 | 14.3 | 14.3 | 23.9 | 1.42 | 0.95 | 0.84 | 115 |
| 5TPY_A | 71 | Zika xrRNA | 7/8 | 21.9 | 14.0 | 14.0 | 16.5 | 1.35 | 0.87 | 0.65 | 95 |
| 2OIU_A | 71 | L1 ligase ribozyme | 6/8 | 15.7 | 13.4 | 10.1 | 18.0 | 1.18 | 0.98 | 0.87 | 111 |
| 2GDI_A | 80 | TPP riboswitch | 4/8 | 16.1 | 16.1 | 16.1 | 17.1 | 1.48 | 1.00 | 0.79 | 134 |
| 3MXH_B | 92 | c-di-GMP-I riboswitch | 0/8 | — | — | — | 22.0 | — | — | 0.33 | 189 |
| 4L81_A | 96 | SAM-I/IV riboswitch | 5/8 | 31.4 | 19.3 | 19.3 | 20.6 | 1.32 | 0.93 | 0.93 | 170 |
| 4QLM_A | 125 | ydaO (c-di-AMP) riboswitch | 2/8 | 27.1 | 22.1 | 22.1 | 29.7 | 1.42 | 0.88 | 0.70 | 280 |
| 3V7E_C | 126 | SAM-I riboswitch (YbxF-bound) | 4/8 | 18.6 | 18.6 | 18.6 | 19.1 | 1.25 | 0.96 | 0.94 | 247 |
| 1U9S_A | 161 | RNase P specificity domain | 4/8 | 33.1 | 30.9 | 30.9 | 27.9 | 1.49 | 1.00 | 0.81 | 332 |
| 3D0U_A | 161 | lysine riboswitch | 0/8 | — | — | — | 29.3 | — | — | 0.74 | 313 |
| 2QBZ_A | 161 | M-box riboswitch | 2/8 | 25.4 | 25.4 | 25.4 | 25.7 | 1.41 | 0.98 | 0.90 | 274 |
| 4GXY_A | 172 | AdoCbl riboswitch | 0/8 | — | — | — | 29.9 | — | — | 0.59 | 335 |
| 4P9R_A | 192 | lariat-capping ribozyme | 1/8 | 28.6 | 28.6 | 28.6 | 29.9 | 1.55 | 0.91 | 0.47 | 368 |

Targets with ≥1 valid model: 17/20. Median top-1 18.6 Å, median best-of-K 16.1 Å (null median 18.4 Å over the same targets). Top-1 better than null: 9/17; best-of-K better than null: 14/17. Best-of-K ≤ 5 Å: 1/20; ≤ 10 Å: 3/20.

## Aggregate (all targets; failures counted as failures)

```json
{
 "independent": {
  "targets": 26,
  "yield": 0.962,
  "median_top1_c4_rmsd": 19.27,
  "median_best_c4_rmsd": 16.56,
  "median_top1_ermsd": 1.55,
  "median_best_ermsd": 1.46,
  "median_top1_inf_wc": 0.77,
  "median_best_inf_wc": 0.81,
  "success_best_c4_le_5A": 0.115,
  "success_best_c4_le_10A": 0.192,
  "success_top1_c4_le_10A": 0.154,
  "success_best_ermsd_le_1": 0.038,
  "median_wall_time_s": 102.1,
  "by_stratum": {
   "s100_149": {
    "targets": 2,
    "with_models": 2,
    "median_best_c4_rmsd": 17.54,
    "median_top1_c4_rmsd": 26.16
   },
   "s150_200": {
    "targets": 6,
    "with_models": 6,
    "median_best_c4_rmsd": 26.54,
    "median_top1_c4_rmsd": 28.32
   },
   "s30_49": {
    "targets": 3,
    "with_models": 3,
    "median_best_c4_rmsd": 15.9,
    "median_top1_c4_rmsd": 19.27
   },
   "s50_99": {
    "targets": 12,
    "with_models": 11,
    "median_best_c4_rmsd": 13.26,
    "median_top1_c4_rmsd": 18.9
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
  "yield": 0.885,
  "median_top1_c4_rmsd": 16.14,
  "median_best_c4_rmsd": 14.05,
  "median_top1_ermsd": 1.42,
  "median_best_ermsd": 1.35,
  "median_top1_inf_wc": 0.97,
  "median_best_inf_wc": 0.97,
  "success_best_c4_le_5A": 0.115,
  "success_best_c4_le_10A": 0.192,
  "success_top1_c4_le_10A": 0.154,
  "success_best_ermsd_le_1": 0.038,
  "median_wall_time_s": 106.5,
  "by_stratum": {
   "s100_149": {
    "targets": 2,
    "with_models": 2,
    "median_best_c4_rmsd": 20.36,
    "median_top1_c4_rmsd": 22.86
   },
   "s150_200": {
    "targets": 6,
    "with_models": 4,
    "median_best_c4_rmsd": 27.01,
    "median_top1_c4_rmsd": 27.53
   },
   "s30_49": {
    "targets": 3,
    "with_models": 3,
    "median_best_c4_rmsd": 11.67,
    "median_top1_c4_rmsd": 14.03
   },
   "s50_99": {
    "targets": 12,
    "with_models": 11,
    "median_best_c4_rmsd": 13.48,
    "median_top1_c4_rmsd": 15.74
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
