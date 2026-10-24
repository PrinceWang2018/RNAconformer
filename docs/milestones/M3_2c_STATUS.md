# M3.2c — Sugar Pucker & Backbone Geometry (experimental increment)

Added `sugar_pucker.py`: five-membered ribose ring pseudorotation diagnostics and restrained ring coordinate optimization. Added optional CLI `--refine-sugar-pucker`, five tests, and per-residue reports in the manifest.

**Important:** This is NOT complete backbone geometry. The ring optimizer holds substituents fixed and may worsen substituent bond geometry; it must not be treated as a validated RNA all-atom repair. Full α/β/γ/δ/ε/ζ/χ coupling, sugar stereocenters, phosphate ester geometry, explicit force-field minimization, and independent stereochemical validation remain outstanding.

No pretrained model is used; atomistic reconstruction still requires externally supplied nucleotide templates.
