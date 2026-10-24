r"""Example analysis plugin for RCV (RNA Conformation Viewer).

Copy this file to the `plugins` folder of the RCV user data folder
(Windows %APPDATA%\RCV\plugins, macOS ~/Library/Application Support/RCV/plugins,
Linux ~/.local/share/rcv/plugins) and restart RCV. Its metrics then appear in the
energy panel (group "Plugins") for every ensemble opened afterwards, and in the CSV export.

A Python plugin defines PLUGIN (id, name, metrics) and compute(sequence, residues) -> {key: value}.
`residues` is a list (one entry per nucleotide, 5'->3') of {atom_name: numpy xyz} after superposition.
"""
import numpy as np

PLUGIN = {
    'id': 'contacts',
    'name': 'Tertiary contacts',
    'metrics': [
        {'key': 'tertiary_contacts', 'label': 'Tertiary contacts', 'unit': '', 'lower_is_better': False,
         'description': 'Residue pairs more than 3 apart in sequence with any heavy atoms within 4.0 Å, '
                        'excluding Watson-Crick partners (example plugin).'},
        {'key': 'end_to_end', 'label': "5'-3' end distance", 'unit': 'Å', 'lower_is_better': None,
         'description': "Distance between the C4' atoms of the first and last nucleotide (example plugin)."},
    ],
}

CUTOFF = 4.0


def compute(sequence, residues):
    centres, coords = [], []
    for r in residues:
        xyz = np.array(list(r.values())) if r else np.zeros((0, 3))
        coords.append(xyz)
        centres.append(xyz.mean(axis=0) if len(xyz) else None)
    n = len(residues)
    contacts = 0
    for i in range(n):
        if centres[i] is None:
            continue
        for j in range(i + 4, n):
            if centres[j] is None or np.linalg.norm(centres[i] - centres[j]) > 18.0:
                continue
            d = np.linalg.norm(coords[i][:, None, :] - coords[j][None, :, :], axis=2)
            if d.min() < CUTOFF and not _wc_partner(sequence, residues, i, j):
                contacts += 1
    first = residues[0].get("C4'") if residues and residues[0] else None
    last = residues[-1].get("C4'") if residues and residues[-1] else None
    end = float(np.linalg.norm(first - last)) if first is not None and last is not None else None
    return {'tertiary_contacts': contacts, 'end_to_end': None if end is None else round(end, 2)}


def _wc_partner(seq, residues, i, j):
    pair = ''.join(sorted(seq[i] + seq[j]))
    if pair not in ('AU', 'CG', 'GU'):
        return False
    a, b = residues[i], residues[j]
    for x, y in (('N1', 'N3'), ('N3', 'N1'), ('N1', 'O2'), ('O2', 'N1'), ('O6', 'N3'), ('N3', 'O6')):
        if x in a and y in b and np.linalg.norm(a[x] - b[y]) < 3.3:
            return True
    return False
