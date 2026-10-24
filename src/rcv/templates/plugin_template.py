r"""__NAME__ — an RCV analysis plugin (template written by Plugins ▸ New plugin).

Python plugins add metrics that are computed when an ensemble is opened (on this computer, with numpy).
Edit this file and restart RCV. Guide: Help ▸ Plugin developer guide (docs/PLUGINS.md).
"""
import numpy as np

PLUGIN = {
    'id': '__ID__',
    'name': '__NAME__',
    'metrics': [
        {'key': 'radius', 'label': 'Max. distance from centre', 'unit': 'Å', 'lower_is_better': None,
         'description': 'Largest distance of a C4\' atom from the centroid (template plugin).'},
    ],
}


def compute(sequence, residues):
    """sequence: str (5'->3'); residues: list of {atom_name: numpy xyz} after superposition.
    Return {metric_key: number or None}."""
    c4 = np.array([r["C4'"] for r in residues if "C4'" in r])
    if len(c4) < 2:
        return {'radius': None}
    return {'radius': round(float(np.linalg.norm(c4 - c4.mean(axis=0), axis=1).max()), 2)}
