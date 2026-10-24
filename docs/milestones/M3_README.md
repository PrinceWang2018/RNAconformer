# M3 — Experimental Atomistic Reconstruction and PDB Export

This milestone introduces a **fail-closed, template-fitted atomistic export adapter**, not a chemically validated RNA reconstruction algorithm.

## Scope

- **M3.1 (experimental):** fit chemically complete, user-provided canonical RNA residue PDB templates to P/C4'/N1-or-N9 guide coordinates using a proper rigid Kabsch rotation.
- **M3.2 (diagnostic only):** inspect inter-residue O3'–P linkage distances, nonlocal heavy-atom clashes, and per-residue anchor fitting residuals. **No molecular mechanics refinement is implemented.**
- **M3.3 (partial):** emit PDB files only when basic checks pass, plus a JSON manifest and explicit rejection reasons. mmCIF and complete stereochemical validation are pending.

## Why external templates?

The M2 guide sites are synthetic and do not determine ribose pucker, phosphate torsions, or chemically correct atom positions. Independent rigid fitting of full-residue templates preserves each residue's internal geometry, but **cannot guarantee covalent continuity across residues**. This module therefore rejects most arbitrary M2 outputs instead of silently creating invalid full-atom structures. It is a scaffold for a future chain-aware reconstruction/refinement backend (e.g., an independently licensed RNA reconstruction tool and validated force field).

## Usage

Supply one chemically complete, correctly named single-residue PDB template for each base present in the sequence. The templates must contain P, C4', N9 (A/G) or N1 (C/U), O3', O5', C3', C5', C1' and enough additional heavy atoms. Input is an M2 JSON result.

```bash
python -m rnafold3d.m3_cli --input ensemble.json \
  --template-a templates/A.pdb --template-c templates/C.pdb \
  --template-g templates/G.pdb --template-u templates/U.pdb \
  --outdir atomistic_out
```

`--allow-invalid-debug` is only for inspecting rejected geometries; such files are labeled invalid in `manifest.json` and **must not** be treated as valid structures.

## Remaining before claiming M3 completion

1. Validated residue templates with proper atom names and stereochemistry.
2. Chain-aware inverse kinematics or restrained atomistic minimization that closes all O3'–P bonds without destroying base-pair geometry.
3. Complete bond/angle/torsion, sugar-pucker, chirality and clash validation.
4. mmCIF writer and standards-compliant ensemble provenance.
5. Benchmarking against experimental RNA structures and independent validators.

## M3.2 update — constrained rigid-residue chain refinement

An **experimental** `refinement.py` module performs chain-aware rigid-body
least-squares optimization. It minimizes deviations from a target O3'–P distance
(1.61 Å), two approximate linkage angles, original P/C4'/base-N anchors, and
nonlocal heavy-atom clashes. Each nucleotide remains a rigid unit, so internal
geometry and chirality of supplied templates are preserved. The optimizer does
**not** optimize sugar puckers or the alpha/beta/gamma/delta/epsilon/zeta/chi
backbone torsions; linkage angle targets are approximate and are not validated
force-field parameters. Its numerical success flag is **not** chemical validity.

Use `--refine-rigid-chain` with `python -m rnafold3d.m3_cli` to opt in.
The manifest stores before/after validation and optimizer diagnostics. Export
still fails closed if basic O3'–P distances or nonlocal clashes are invalid.

The finite-difference implementation is intentionally limited to 40 residues;
**200-nt atomistic refinement is not yet supported**. For production-scale RNA,
a validated force field and torsion-space refinement backend are required.

Test with `PYTHONPATH=src python -m unittest discover -s tests -v`.

## M3.2b: experimental torsion-space refinement

New module: `src/rnafold3d/torsion_refinement.py`.

Run the existing atomistic CLI with `--refine-torsions` (optionally after
`--refine-rigid-chain`). The torsion module requires complete, canonical RNA
residue templates and currently limits optimization to 20 nt. It rotates
acyclic P-O5', O5'-C5', C3'-O3', and O3'-P(next) bonds; optional chi rotation
is available through the Python API. Cyclic sugar bonds and sugar pucker are
not optimized. The `alpha` and `zeta` dihedral definitions share the
O3'(previous)-P bond axis, so that axis is optimized once, not twice.

Example:

```bash
python -m rnafold3d.m3_cli --input ensemble.json \
  --template-a templates/A.pdb --template-c templates/C.pdb \
  --template-g templates/G.pdb --template-u templates/U.pdb \
  --refine-torsions --outdir refined_output
```

The optimizer reports objective before/after, selected torsion axes,
maximum rotation, linkage diagnostics, and its numerical convergence flag.
**Convergence is not evidence of stereochemical validity.** Output remains
blocked when basic atomistic checks fail, unless the explicit debugging
export option is set. In particular, this is not full backbone refinement:
no sugar pucker, ring closure optimization, explicit solvent, validated force
field, or independent RNA stereochemical validation is included.

For a practical 50-200 nt pipeline, this experimental local optimizer needs
replacement or substantial extension before production use.
