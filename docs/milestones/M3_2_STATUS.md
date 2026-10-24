# M3.2 — Chain-Aware Stereochemical Refinement: Implementation Status

## Implemented

- A joint, rigid-residue, restrained least-squares optimizer (`refinement.py`).
- O3'–P distance restraint and approximate C3'–O3'–P / O3'–P–O5' angular restraints.
- Guide anchor restraints to reduce loss of predicted overall conformation.
- Nonlocal heavy-atom soft-clash penalties.
- Preservation of the supplied nucleotide templates' internal geometry and chirality.
- Before/after linkage and clash diagnostics, optimizer status and objective.
- Optional CLI integration with fail-closed export behavior.
- Four new unit tests; 73 tests total passing.

## Not yet implemented / scientific caveats

- No internal backbone torsion, sugar pucker, or phosphate geometry optimization.
- No independently validated RNA force-field minimization.
- No full stereochemical validation or experimental structure benchmark.
- No assurance that optimizer convergence implies a valid RNA structure.
- Finite-difference prototype limited to 40 residues; 200 nt is not supported.
- Approximate restraint constants must not be interpreted as calibrated force-field parameters.

## Next implementation step

Implement torsion-space chain closure with fixed covalent geometry and explicit
sugar pucker states, followed by a validated force-field minimization and
independent stereochemical checks. Only then attempt 50–200 nt benchmarks.
