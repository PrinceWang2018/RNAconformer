# RCV demo data

Opened by *Try the demo ensembles* on the start page or `rcv --demo`; bundled with the installers.

| Item | Content | Source |
|---|---|---|
| `HIV1_TAR_NMR_1ANR.pdb` | 20 NMR models of the HIV-1 TAR RNA (29 nt), heavy atoms, one file with MODEL records | PDB entry 1ANR (Aboul-ela, Karn & Varani, *Nucleic Acids Res.* 1996); PDB data are CC0 |
| `HIV1_TAR_prediction/` | rnafold3d v2.0.1 prediction from the 1ANR sequence (`HIV1_TAR.fasta`): 24 trajectories, 10 conformers over two ViennaRNA 2D hypotheses, seed 7 | generated with `python -m rnafold3d predict --fasta HIV1_TAR.fasta --models 10 --trajectories 24 --seed 7` |
| `1EHZ_tRNA_prediction/` | rnafold3d v2.0.1 sequence-only prediction of yeast tRNA-Phe (76 nt), 5 conformers | copy of the benchmark pilot run `benchmark/runs/v2_0_1_pilot/independent/1EHZ_A` |

Predicted models are physics-inspired predictions, not experimental structures (see docs/V2_STATUS.md for
their accuracy).
