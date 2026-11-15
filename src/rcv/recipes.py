"""Recipes: named parameter sets for a whole analysis, saved and shared as a small file (.rcvrecipe, JSON).

A recipe has a name, a description, optional prediction settings (used when a prediction is started with the
recipe chosen on the start page) and steps: console commands (or plain-language lines) run on the ensemble, in
order — conditions, colours, tabs, views, comparisons. A lab can keep its way of looking at an RNA in one file
and pass it around. User recipes live in <user data>/recipes/; the built-in ones are below.

    {"rcv_recipe": 1, "name": "Riboswitch report", "description": "...", "author": "",
     "prediction": {"models": 10, "sampling": "thorough", "templates": "auto"},
     "steps": ["conditions standard", "color switching", "tab switch"]}
"""
from __future__ import annotations

import json
import re
from pathlib import Path

from . import paths as P

FORMAT = 1
EXT = '.rcvrecipe'
SAMPLING = ('fast', 'standard', 'thorough')
MAX_STEPS = 200

BUILTIN = [
    {'id': 'builtin-riboswitch-report', 'name': 'Riboswitch report',
     'description': 'For a two-state switch: thorough sampling, then the states and their populations under standard '
                    'conditions, the residues that switch, the mutations that shift the balance, and a saved overview '
                    'view.',
     'prediction': {'models': 10, 'sampling': 'thorough', 'templates': 'auto'},
     'steps': ['conditions standard', 'order rank', 'color switching', 'show ensemble', 'select switching',
               'tab switch', 'zoom', 'view save overview',
               'echo Switch tab: states, populations and mutations; Pathway tab: run "pathway first, last" for the barrier.']},
    {'id': 'builtin-pull-and-compare', 'name': 'Pull and compare',
     'description': 'Compare the first and the last conformation: superpose them, list the residues that move most, '
                    'colour by displacement and play the morph; then drag a residue with the Pull tool.',
     'prediction': {},
     'steps': ['conf first', 'conf last', 'align first, last', 'color displacement', 'show labels, moved > 4',
               'tab energy', 'morph first, last, 2',
               'echo Press G and drag a residue to pull the RNA into another conformation.']},
    {'id': 'builtin-flexibility', 'name': 'Flexibility overview',
     'description': 'Where the ensemble is flexible: colour by RMSF, flexibility putty, the most flexible residues '
                    'labelled, Details tab.',
     'prediction': {},
     'steps': ['color flexibility', 'show putty', 'label rmsf > 3', 'tab details', 'zoom']},
    {'id': 'builtin-publication-figure', 'name': 'Publication figure',
     'description': 'A clean figure: illustrative style on white, helices coloured, every 10th residue labelled, the '
                    'molecule oriented along its long axis. Then Export > Image or Figure panels.',
     'prediction': {},
     'steps': ['style illustrative', 'bg white', 'rep cartoon', 'color helix', 'set labels, every10', 'orient',
               'echo Ready: Export > Image (Ctrl+E) or Export > Figure panels.']},
]


def recipes_dir():
    return P.user_data_dir() / 'recipes'


def slug(name):
    s = re.sub(r'[^A-Za-z0-9]+', '-', str(name)).strip('-').lower()
    return s[:60] or 'recipe'


def validate(obj):
    """A clean recipe dict from user input; ValueError for anything that is not a recipe."""
    if not isinstance(obj, dict):
        raise ValueError('a recipe is a JSON object')
    if obj.get('rcv_recipe', FORMAT) != FORMAT:
        raise ValueError(f'recipe format {obj.get("rcv_recipe")} is newer than this RCV understands')
    name = str(obj.get('name') or '').strip()
    if not name:
        raise ValueError('the recipe has no name')
    steps = obj.get('steps') or []
    if isinstance(steps, str):
        steps = steps.splitlines()
    if not isinstance(steps, list):
        raise ValueError('steps must be a list of commands')
    steps = [str(s).strip()[:300] for s in steps if str(s).strip()][:MAX_STEPS]
    pred_in = obj.get('prediction') or {}
    if not isinstance(pred_in, dict):
        raise ValueError('prediction must be an object')
    pred = {}
    if pred_in.get('models') not in (None, ''):
        pred['models'] = max(1, min(50, int(pred_in['models'])))
    if pred_in.get('sampling'):
        if pred_in['sampling'] not in SAMPLING:
            raise ValueError(f'sampling must be one of {", ".join(SAMPLING)}')
        pred['sampling'] = pred_in['sampling']
    if pred_in.get('seed') not in (None, ''):
        pred['seed'] = int(pred_in['seed'])
    if pred_in.get('templates'):
        pred['templates'] = 'off' if str(pred_in['templates']).lower() in ('off', 'false', 'no') else 'auto'
    if not steps and not pred:
        raise ValueError('the recipe has neither steps nor prediction settings')
    return {'rcv_recipe': FORMAT, 'name': name[:80], 'description': str(obj.get('description') or '').strip()[:1000],
            'author': str(obj.get('author') or '').strip()[:120], 'prediction': pred, 'steps': steps}


def list_recipes():
    out = [{**r, 'rcv_recipe': FORMAT, 'author': 'RCV', 'builtin': True} for r in BUILTIN]
    d = recipes_dir()
    if d.is_dir():
        for f in sorted(d.glob('*' + EXT)):
            try:
                r = validate(json.loads(f.read_text(encoding='utf-8')))
            except (OSError, ValueError, TypeError) as exc:
                out.append({'id': f.stem, 'name': f.stem, 'error': str(exc), 'path': str(f), 'builtin': False, 'steps': []})
                continue
            out.append({**r, 'id': f.stem, 'path': str(f), 'builtin': False})
    return out


def get(rid_or_name):
    key = str(rid_or_name).strip().lower()
    for r in list_recipes():
        if r['id'].lower() == key or r['name'].lower() == key:
            return r
    raise KeyError(rid_or_name)


def save(obj, rid=None):
    """Save a user recipe (a new file, or replace rid); returns it with its id and path."""
    r = validate(obj)
    if rid and str(rid).startswith('builtin-'):
        rid = None                                          # a changed built-in recipe is saved as the user's own
    d = recipes_dir()
    d.mkdir(parents=True, exist_ok=True)
    if rid and re.fullmatch(r'[A-Za-z0-9_-]+', str(rid)):
        f = d / f'{rid}{EXT}'
    else:
        base, k = slug(r['name']), 1
        f = d / f'{base}{EXT}'
        while f.exists():
            k += 1
            f = d / f'{base}-{k}{EXT}'
    f.write_text(json.dumps(r, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    return {**r, 'id': f.stem, 'path': str(f), 'builtin': False}


def import_text(text):
    try:
        obj = json.loads(text)
    except ValueError as exc:
        raise ValueError(f'not a recipe file (JSON): {exc}') from None
    return save(obj)


def delete(rid):
    if str(rid).startswith('builtin-'):
        raise ValueError('built-in recipes cannot be deleted')
    if not re.fullmatch(r'[A-Za-z0-9_-]+', str(rid)):
        raise ValueError('bad recipe id')
    f = recipes_dir() / f'{rid}{EXT}'
    if not f.is_file():
        raise KeyError(rid)
    f.unlink()
    return {'deleted': rid}
