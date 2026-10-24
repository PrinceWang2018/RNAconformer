# M3.3 — Independent Stereochemical Validation

## Implemented

- `src/rnafold3d/stereochemistry.py`: independent atomistic diagnostics.
- Checks: atom completeness, finite coordinates, selected intranucleotide bond distances and angles, phosphodiester linkage, backbone torsion measurement, ribose ring torsion measurement, near-planar sugar stereocenter, and nonlocal heavy atom clashes.
- `export_ensemble()` includes stereochemistry results in the manifest and gates the default PDB export on both the previous basic checks and this new validator.
- Debug export remains explicitly possible with `--allow-invalid-debug` and is marked invalid.

## Important limitations

- Reference intervals are provisional broad screening bounds, **not** calibrated against high-quality RNA structures.
- Sugar-pucker phase/amplitude and canonical conformer classification are not yet implemented.
- Signed sugar volume is a degeneracy check, not a reference-verified absolute stereochemistry assignment.
- Not all atomistic bonds and angles are covered; adjacent-residue clash exclusions are coarse.
- A passing result is **not** a scientifically validated RNA structure.
- No experimental-structure benchmarking was performed.

## Tests

Run `PYTHONPATH=src python -m unittest discover -s tests -q`.

Observed: 87 tests passed, including 5 new stereochemistry unit tests.
