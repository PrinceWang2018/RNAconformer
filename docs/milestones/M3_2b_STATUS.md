# M3.2b Status — RNA Backbone Torsion-Space Refinement

Implemented an experimental acyclic-bond torsion optimizer and wired it to
`rnafold3d.m3_cli` with `--refine-torsions`.

## Included
- Graph-based downstream atom selection; ring bonds cannot rotate.
- Beta, gamma, epsilon and zeta bond-axis rotations; optional chi.
- Rotation limits, guide anchors, inter-residue linkage and angle restraints,
  and nonlocal clash penalties.
- Before/after diagnostic report, rejection on failed basic atomistic checks.
- Unit tests for graph rotor selection, length invariance, optimizer objective,
  malformed inputs and single-residue handling.

## Not included
- Sugar-pucker sampling, delta torsion, complete alpha/beta/gamma/delta/
  epsilon/zeta/chi RNA internal-coordinate representation.
- Full atomistic force field or validated RNA stereochemistry.
- 50-200 nt production optimization (current cap: 20 nt).
- Evidence of success on experimental RNA structures.

This is a prototype and does not guarantee exportable PDB structures.
