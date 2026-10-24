"""M2.1 experimental coarse-grained RNA guide-site energy engine.

Units: coordinates in Angstrom; energies are arbitrary model units, NOT kcal/mol.
P/C4/N sites from M1 are synthetic guides, not chemically calibrated atoms.
This module is for testing energy behavior, not ranking real RNA conformations.
"""
from dataclasses import dataclass, asdict
import math
import numpy as np

@dataclass(frozen=True)
class EnergyConfig:
    bond_k: float = 1.0
    bond_target: float = 6.0
    pair_k: float = 1.0
    pair_target: float = 6.0
    orientation_k: float = 0.5
    stack_strength: float = 1.0
    stack_target: float = 4.0
    stack_sigma: float = 2.0
    steric_k: float = 2.0
    steric_cutoff: float = 4.0
    electrostatic_strength: float = 1.0
    screening_length: float = 10.0
    min_separation: float = 1e-6

    def __post_init__(self):
        for name, val in asdict(self).items():
            if not math.isfinite(val) or val <= 0:
                raise ValueError(f'{name} must be finite and positive')

@dataclass(frozen=True)
class EnergyBreakdown:
    bonded: float
    pairing: float
    orientation: float
    stacking: float
    steric: float
    electrostatic: float
    total: float
    n_residues: int
    n_pairs: int
    n_stacks: int

    def as_dict(self):
        return asdict(self)


def _extract(sites, sequence):
    if not sites:
        raise ValueError('At least one nucleotide site is required')
    mapping = {}
    for s in sites:
        i = int(s.residue)
        if i < 0 or i >= len(sequence) or i in mapping:
            raise ValueError('Out-of-range or duplicate residue index')
        xyz = np.asarray([s.P, s.C4, s.N], dtype=float)
        if xyz.shape != (3, 3) or not np.isfinite(xyz).all():
            raise ValueError('Invalid site coordinates')
        mapping[i] = xyz
    return mapping


def score_sites(model, sites, config=EnergyConfig(), *, require_complete=True):
    """Score P/C4/N guide sites; missing residues optionally allowed.

    Pairing is a harmonic anchor-distance restraint to the *input secondary structure*.
    Stacking uses a directional base-anchor proxy, not real base-plane normals.
    Excludes sequence neighbors and specified paired residues from nonbonded sterics.
    """
    xyz = _extract(sites, model.sequence)
    if require_complete and set(xyz) != set(range(len(model.sequence))):
        raise ValueError('Complete residue coverage required')
    ids = sorted(xyz)
    c4 = {i: xyz[i][1] for i in ids}
    p = {i: xyz[i][0] for i in ids}
    n = {i: xyz[i][2] for i in ids}
    unit = {}
    for i in ids:
        v = n[i] - c4[i]
        length = np.linalg.norm(v)
        if length <= config.min_separation:
            raise ValueError('Degenerate C4-to-base direction')
        unit[i] = v / length

    bonded = 0.0
    for i in ids:
        if i + 1 in xyz:
            d = np.linalg.norm(c4[i] - c4[i+1])
            bonded += config.bond_k * (d - config.bond_target)**2

    pairs = [tuple(pair) for pair in model.pairs if pair[0] in xyz and pair[1] in xyz]
    pair_set = {tuple(sorted(pair)) for pair in model.pairs}
    pairing = orientation = 0.0
    for i, j in pairs:
        d = np.linalg.norm(n[i] - n[j])
        pairing += config.pair_k * (d - config.pair_target)**2
        # Facing base anchors should point towards each other.
        orientation += config.orientation_k * (1.0 + float(np.dot(unit[i], unit[j])))**2

    # Stacking is applied to successive residues on each strand of an input stem.
    # Gaussian attraction with a parallel-direction preference.
    stacks = set()
    for stem in model.stems:
        for (i, j), (ip, jp) in zip(stem.pairs, stem.pairs[1:]):
            for a, b in ((i, ip), (j, jp)):
                if a in xyz and b in xyz:
                    stacks.add(tuple(sorted((a, b))))
    stacking = 0.0
    for i, j in stacks:
        d = np.linalg.norm(n[i] - n[j])
        direction = max(0.0, float(np.dot(unit[i], unit[j])))
        stacking -= config.stack_strength * math.exp(-0.5*((d-config.stack_target)/config.stack_sigma)**2) * direction

    steric = electrostatic = 0.0
    for k, i in enumerate(ids):
        for j in ids[k+1:]:
            if abs(i-j) <= 2 or (i,j) in pair_set:
                continue
            d = max(float(np.linalg.norm(c4[i]-c4[j])), config.min_separation)
            steric += config.steric_k * max(0.0, config.steric_cutoff-d)**2
            # screened P-P repulsion; excluded for local sequence neighbors
            dp = max(float(np.linalg.norm(p[i]-p[j])), config.min_separation)
            electrostatic += config.electrostatic_strength * math.exp(-dp/config.screening_length)/dp

    total = bonded+pairing+orientation+stacking+steric+electrostatic
    return EnergyBreakdown(*(float(x) for x in (bonded,pairing,orientation,stacking,steric,electrostatic,total)),len(ids),len(pairs),len(stacks))
