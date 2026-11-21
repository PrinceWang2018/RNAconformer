"""End-to-end v2 prediction: sequence -> 2D hypotheses -> 3D ensemble of validated all-atom models.

    rnafold3d predict --fasta input.fa [--dbn input.dbn] --models 5 --seed 42 --out run_001

Pipeline per trajectory (independent, seeded by (seed, trajectory index)):
  1. FoldModel for the trajectory's 2D hypothesis (rigid A-form stems + fold tree)
  2. random torsions from the calibrated loop-conformer library, simulated annealing
  3. CCD closure of junction and hairpin loops (several starts, least clash / torsion strain)
  4. glycosidic (chi) clean-up of loop residues
  5. independent all-atom validation (checks.validate_model)
Valid models are ranked by 3D score + 2D free-energy offset and selected greedily for
structural diversity (C4' RMSD). Invalid models are never exported as predictions.
"""
from __future__ import annotations

import csv
import datetime
import hashlib
import json
import math
import os
import sys
import platform
import time
from dataclasses import dataclass, asdict, field
from pathlib import Path

import numpy as np

from . import __version__ as PACKAGE_VERSION
from .annotate import check_dotbracket, dotbracket_to_pairs, nested_subset, pairs_to_dotbracket
from .checks import validate_model
from .fold import FoldModel
from .loops import TorsionPrior, close_linker, final_relief, finish_closure
from .search import Scorer, ScoreWeights, AnnealConfig, anneal, fit_restraints, relax_modified
from .structure_io import write_pdb, write_mmcif

ENGINE = 'rnafold3d-v2-foldtree'
ENGINE_VERSION = '2.5.0'
KCAL_TO_KT_37C = 1.0 / (0.0019872 * 310.15)


@dataclass
class PredictConfig:
    n_models: int = 5
    trajectories: int = 0            # 0 -> max(8, 4 * n_models)
    seed: int = 42
    max_structures: int = 3          # 2D hypotheses when predicting from sequence
    mfe_share: float = 0.6           # fraction of trajectories for the MFE hypothesis
    temperature: float = 37.0
    closure_attempts: int = 24
    cluster_rmsd: float = 0.0        # 0 -> automatic, 0.15 * Rg0 clamped to [2, 6] A
    short_loop_max: int = 4          # closure loops up to this length are explicit from the start
    pseudoknot_restraints: bool = False  # EXPERIMENTAL: crossing pairs of a user 2D -> distance restraints
    anneal_steps: int = 0            # 0 -> automatic
    workers: int = 1
    export_rejected: bool = False
    save_candidates: bool = False    # also write every valid trajectory model to candidates/
    weights: ScoreWeights = field(default_factory=ScoreWeights)
    # guided exploration (RCV "pull"): extra 2D hypotheses tried first, distance restraints
    # [{'i','j' (0-based residues), 'a','b' (atom names), 'target', 'tol'}] applied as a bias during
    # sampling, then removed for an unbiased relaxation of `relax_steps` annealing steps
    seed_structures: tuple = ()
    pull_restraints: tuple = ()
    relax_steps: int = 1500
    # template-based modelling from the built-in library of experimental RNA structures (templates.py)
    templates: str = 'auto'          # 'auto': use significant hits; 'off': physics only
    template_cutoff: str = ''        # only templates released before this date (YYYY-MM-DD), for blind tests
    template_exclude: tuple = ()     # PDB ids that are never used as templates
    max_templates: int = 3
    template_share: float = 0.7      # fraction of trajectories for template hypotheses when hits exist
    motifs: str = 'auto'             # junction geometry from non-homologous library structures ('off': none)
    template_structure: bool = True  # also accept remote homologs confirmed by their secondary structure
    template_fit_attempts: int = 4   # fits onto the template atoms: after a fit that fails validation, the residues
                                     # with the problems are released and the fit repeated; models whose fit was
                                     # rejected rank after the others (1: one fit, order by score, as up to 2.4.0)
    domains: str = 'off'             # 'auto': chains >= 240 nt fold their domains separately, then assemble them
                                     # (experimental: worse than the whole-chain search in benchmark/long_rna.py)
    pseudoknots: str = 'off'         # 'auto': 2D hypothesis with a crossing helix (pseudoknots.py) when it beats the
                                     # MFE (no gain on top of fragments in benchmark/dev_templates.py, slower)
    fragments: str = 'auto'          # motif fragments (hairpin loops, junction strands) grafted as torsions ('off': none)
    # experimental and evolutionary information for the 2D hypotheses (sequence-only predictions)
    reactivities: tuple = ()         # SHAPE-like reactivity per nucleotide (None = missing), see ss.read_reactivities
    alignment: tuple = ()            # aligned homologous sequences, query first, see ss.read_alignment
    # modified nucleotides [{'pos' (1-based), 'code'}] (modifications.py); a sequence with [code] tokens also works
    modifications: tuple = ()
    # Z-RNA (zrna.py): 'auto' adds a hypothesis with one helix segment built left-handed when its estimated
    # A -> Z free energy is below zrna.Z_PROPOSE_MAX under these conditions (extra trajectories, about +25 % time)
    zrna: str = 'off'
    zrna_salt: float = 0.15          # Na+ (mol/L) for the Z-RNA estimate
    zrna_binder: bool = False        # a Z-binding protein is present (ADAR1 Zalpha, ZBP1)


def normalize_sequence(seq, allow_t=True):
    # position numbers are dropped: sequences copied from GenBank or a paper read "1 ggcagatctg 11 agcctgggag"
    s = ''.join(ch for ch in seq if not ch.isspace() and not ch.isdigit()).upper()
    if allow_t:
        s = s.replace('T', 'U')
    bad = sorted(set(s) - set('ACGU'))
    if bad:
        raise ValueError(f'sequence contains unsupported symbols {bad}; only A/C/G/U (T->U) are modelled')
    if not s:
        raise ValueError('empty sequence')
    return s


def read_fasta(path):
    name, seq = None, []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if name is not None and seq:
                break
            name = line[1:].strip() or 'query'
            continue
        seq.append(line)
    return name or 'query', ''.join(seq)


def read_dbn(path, n):
    lines = [l.strip() for l in Path(path).read_text().splitlines() if l.strip() and not l.startswith('>')]
    for l in reversed(lines):
        tok = l.split()[0]
        if len(tok) == n and set(tok) <= set('.()[]{}<>'):
            return tok
    raise ValueError(f'no dot-bracket line of length {n} found in {path}')


# ---------------------------------------------------------------------------- trajectory worker
def _c4(residues):
    return np.array([r["C4'"] for r in residues])


def _seed_hypotheses(seq, seeds, temperature, hyps):
    """Put user/seed dot-brackets first (deduplicated against `hyps`), renumbering ids ss1, ss2, ..."""
    out, seen = [], set()
    for k, dbn in enumerate(seeds):
        if len(dbn) != len(seq):
            raise ValueError('seed structure length differs from the sequence length')
        pairs = [(i - 1, j - 1) for i, j in dotbracket_to_pairs(dbn.replace('[', '.').replace(']', '.')
                                                                  .replace('{', '.').replace('}', '.'))]
        pairs = [(i, j) for i, j in pairs if seq[i] + seq[j] in ('AU', 'UA', 'GC', 'CG', 'GU', 'UG')]
        nested, _ = nested_subset(pairs)
        d = pairs_to_dotbracket(len(seq), nested, offset=0)
        if d in seen:
            continue
        seen.add(d)
        e = None
        try:
            import RNA
            md = RNA.md()
            md.temperature = temperature
            e = round(RNA.fold_compound(seq, md).eval_structure(d), 2)
        except Exception:
            pass
        out.append({'id': '', 'dbn': d, 'pairs': nested, 'energy_kcal_mol': e, 'source': 'seed'})
    out += [h for h in hyps if h['dbn'] not in seen]
    for k, h in enumerate(out, 1):
        h['id'] = f'ss{k}'
    return out


def _r3(v):
    return None if v is None else round(float(v), 3)


def _round_stats(st):
    return None if st is None else {k: (round(v, 3) if isinstance(v, float) else v) for k, v in st.items()}


def run_trajectory(job):
    """Worker: one independent search trajectory. Returns a JSON-friendly candidate dict.

    Every BLAS/LAPACK library that the worker may use (numpy's and scipy's bundled OpenBLAS,
    loaded lazily by scipy.optimize/scipy.linalg) is imported first and then limited to one
    thread; otherwise N workers x M BLAS threads oversubscribe the CPU by an order of magnitude.
    """
    import scipy.linalg       # noqa: F401  (load all BLAS libraries before limiting them)
    import scipy.optimize     # noqa: F401
    try:
        from threadpoolctl import threadpool_limits
    except ImportError:
        return _run_trajectory(job)
    with threadpool_limits(limits=1):
        return _run_trajectory(job)


def _worker_init():
    """Start-up of a trajectory worker: stop as soon as the engine process is gone. Without it, the workers of a
    cancelled or killed run finish their trajectory (minutes on every core for a long RNA) before they notice."""
    import multiprocessing as mp
    import signal
    import threading
    try:
        signal.signal(signal.SIGTERM, signal.SIG_DFL)     # not the engine's handler (inherited with fork)
        signal.signal(signal.SIGINT, signal.SIG_IGN)      # Ctrl+C is handled by the engine process, which ends the pool
    except (ValueError, OSError):
        pass
    parent = mp.parent_process()
    if parent is None:
        return
    ppid = os.getppid()

    def watch():
        # the sentinel covers Windows (and spawn); the parent id covers POSIX, where an orphan gets a new parent
        while parent.is_alive() and (os.name == 'nt' or os.getppid() == ppid):
            time.sleep(1.0)
        os._exit(1)
    threading.Thread(target=watch, daemon=True, name='engine-watch').start()


class _StopOnSignal:
    """While the worker pool runs, a request to stop (SIGTERM; SIGBREAK on Windows) ends the run through the normal
    exit path, so the pool terminates its workers. Only in the main thread (signal handlers cannot be set elsewhere);
    the previous handlers are put back afterwards."""

    def __enter__(self):
        import signal
        import threading
        self.saved = {}
        if threading.current_thread() is not threading.main_thread():
            return self

        def stop(signum, _frame):
            raise SystemExit(128 + signum)
        for name in ('SIGTERM', 'SIGBREAK'):
            sig = getattr(signal, name, None)
            if sig is None:
                continue
            try:
                self.saved[sig] = signal.signal(sig, stop)
            except (ValueError, OSError):
                pass
        return self

    def __exit__(self, *exc):
        import signal
        for sig, handler in self.saved.items():
            try:
                signal.signal(sig, handler if handler is not None else signal.SIG_DFL)
            except (ValueError, OSError, TypeError):
                pass
        return False


LOCAL_ISSUES = {'bond_angle', 'bond_length', 'chain_break', 'severe_clash'}
# fit onto the template atoms (PredictConfig.template_fit_attempts): when the fitted model fails validation, the
# residues with the problems (and their neighbours, one more on each side per attempt) are released from the
# template and the fit is repeated; given up when more than this fraction of the chain would be released
TEMPLATE_FIT_MAX_RELEASED = 0.35


def _release_for_refit(issues, attempt, released, n):
    """Residues (0-based) left out of the next fit onto the template: those already released plus the residues
    of the validation issues with attempt + 1 neighbours on each side. None when nothing new would be released
    or more than TEMPLATE_FIT_MAX_RELEASED of the chain would be: the fit is then given up."""
    bad = {int(i['residue']) - 1 for i in issues}
    grow = {r + d for r in bad for d in range(-1 - attempt, 2 + attempt) if 0 <= r + d < n} - set(released)
    if not grow or len(set(released) | grow) > TEMPLATE_FIT_MAX_RELEASED * n:
        return None
    return set(released) | grow


def _order_template_models(by_template, refit):
    """Rank keys of the template models, {template rank: [trajectory results]}. Within a template the models are
    ordered by the 3D score. A model whose fit onto the template atoms was rejected kept the shape of the
    torsion-space search, which follows the template only loosely: with refit (template_fit_attempts > 1) such
    models come after all models with an accepted fit, of every template."""
    unfit = lambda r: int(refit and (r.get('template_fit') or {}).get('accepted') is False)
    for t_rank, rs in by_template.items():
        for pos, r in enumerate(sorted(rs, key=lambda r: (unfit(r), r['score_3d']))):
            r['rank_key'] = (0, unfit(r), pos, t_rank, r['score_3d'])


def _written_atom_names(interactions, mods):
    """The interaction report of the search with the atom names the files carry (a residue written as a PDB chemical
    component has that component's atom names, modifications.COMPONENTS)."""
    from .modifications import written_names
    names = written_names(mods)
    if not names:
        return interactions
    out = []
    for row in interactions:
        row = dict(row)
        own = names.get(row.get('residue'), {})
        hb = []
        for h in row.get('hbonds') or ():
            res, _, atom = str(h.get('partner', '')).partition(':')
            other = names.get(int(res), {}) if res.isdigit() else {}
            hb.append({**h, 'atom': own.get(h.get('atom'), h.get('atom')), 'partner': f'{res}:{other.get(atom, atom)}'})
        row['hbonds'] = hb
        out.append(row)
    return out


def _run_trajectory(job):
    t0 = time.time()
    seq, hyp, k, seed, cfg = job['sequence'], job['hypothesis'], job['index'], job['seed'], job['config']
    weights = ScoreWeights(**cfg['weights'])
    rng = np.random.default_rng([seed, k])
    out = {'trajectory': k, 'ss_id': hyp['id'], 'seed': [seed, k]}
    if hyp.get('motifs'):
        out['junction_motifs'] = hyp['motifs']
    try:
        model = FoldModel(seq, hyp['pairs'], explicit_loops=True, modifications=cfg.get('modifications'),
                          mod_pairs=hyp.get('mod_pairs') or (), z_pairs=hyp.get('z_pairs') or ())
        prior = TorsionPrior(model.geom)
        scorer = Scorer(model, weights, prior, restraints=hyp.get('restraints', ()),
                        pull_restraints=cfg.get('pull_restraints') or (),
                        distance_restraints=hyp.get('distance_restraints'))
        # stage 1: global arrangement; short closure loops explicit, long loops and hairpins implicit
        short = {L.id for L in model.loop_linkers if L.kind == 'closure' and L.n <= cfg.get('short_loop_max', 4)}
        # pseudoknot helix (pseudoknots.py): loops holding one of its strands are explicit from the start, so that the
        # fit and the search can bring the two strands together; its residues start from A-form torsions
        pk_res = {int(r) for p in (hyp.get('restraints') or ()) for r in p} if not hyp.get('torsions') else set()
        short |= {L.id for L in model.loop_linkers if any(r in pk_res for r in L.loop)}
        model.set_active_loops(short)
        model.sample_torsions(rng)
        frag_res = set()
        if hyp.get('fragment_torsions'):         # motif fragments from library structures (fragments.py)
            model.apply_torsions(hyp['fragment_torsions'])
            frag_res = {int(r) for r in hyp['fragment_torsions']}
        if pk_res:
            af = model.geom.aform
            model.apply_torsions({r: [af['alpha'], af['beta'], af['gamma'], 80.0, af['epsilon'], af['zeta'], af['chi']]
                                  for r in pk_res})
        guided = bool(hyp.get('torsions'))
        if guided:
            model.apply_torsions(hyp['torsions'])
            scorer.dr_scale = float(cfg.get('template_weight', 5.0))
        model.build()
        if scorer.dr is not None:              # template / junction-motif restraints: least-squares fit first
            scorer.clash_scale = 0.3
            out['fit'] = _round_stats(fit_restraints(model, scorer, iters=60 if guided else 30))
        # a template that covers (nearly) everything fixes the fold: small moves only; a partial template
        # leaves the arrangement of the other parts open: broader search
        close_fit = guided and float(hyp.get('template_coverage') or 1.0) >= 0.85
        if close_fit:
            a1 = AnnealConfig(steps=cfg['anneal_steps'], t_start=0.8, sigma_start=10.0, sigma_end=3.0,
                              p_residue_jump=0.02, p_junction_aform=0.0)
        elif guided:
            a1 = AnnealConfig(steps=cfg['anneal_steps'], t_start=3.0, sigma_start=25.0, sigma_end=4.0,
                              p_residue_jump=0.1, p_junction_aform=0.0)
        else:
            a1 = AnnealConfig(steps=cfg['anneal_steps'])
        stats = anneal(model, scorer, rng, a1)
        out['search'] = {k2: (round(v, 3) if isinstance(v, float) else v) for k2, v in stats.items()}
        # stage 2: grow + CCD-close the remaining loops
        closures = []
        order = sorted(model.loop_linkers, key=lambda L: (L.kind != 'closure', L.n))
        for L in order:
            if L.id in model.active_loops:
                continue
            grafted = scorer.dr is not None and L.n >= 1 and all(r in frag_res for r in L.loop)
            if (guided or grafted) and scorer.dr is not None:   # keep the template / fragment loop shape
                closures.append(close_linker(model, L, rng, prior, attempts=cfg['closure_attempts'], from_current=4,
                                             extra=lambda: scorer.dr_scale * scorer.w.template
                                             * scorer.distance_restraint_energy(), restraints=scorer.dr))
            else:
                closures.append(close_linker(model, L, rng, prior, attempts=cfg['closure_attempts']))
            model.set_active_loops(model.active_loops | {L.id})
        if not all(c.get('closed') for c in closures):
            out['closures'] = closures
            out['status'] = 'rejected'
            out['reasons'] = ['loop_closure_failed']
            out['time_s'] = round(time.time() - t0, 2)
            return out
        # stage 3: joint annealing with every loop explicit (closure-preserving moves)
        model.set_loops_active(True)
        if scorer.dr is not None:
            scorer.clash_scale = 2.0
            out['fit2'] = _round_stats(fit_restraints(model, scorer, iters=40 if guided else 20))
        ref = anneal(model, scorer, rng, AnnealConfig(steps=min(8000, 2000 + 12 * len(model.active_dofs)),
                                                      t_start=0.6 if close_fit else 1.5, t_end=0.1,
                                                      sigma_start=8.0 if close_fit else 15.0,
                                                      sigma_end=3.0 if close_fit else 4.0,
                                                      clash_scale_start=2.0, clash_scale_end=25.0,
                                                      p_residue_jump=0.01 if close_fit else 0.05, p_junction_aform=0.0))
        out['refine'] = {k2: (round(v, 3) if isinstance(v, float) else v) for k2, v in ref.items()}
        if scorer.pull:
            # stage 3b: remove the pull bias and relax into the nearest stable state of the unbiased score
            out['pull'] = {'rms_biased_A': _r3(scorer.pull_deviation())}
            scorer.pull_scale = 0.0
            anneal(model, scorer, rng, AnnealConfig(steps=int(cfg.get('relax_steps') or 1500), t_start=0.5, t_end=0.05,
                                                     sigma_start=6.0, sigma_end=2.0, clash_scale_start=25.0,
                                                     clash_scale_end=25.0, p_residue_jump=0.02, p_junction_aform=0.0))
            out['pull']['rms_relaxed_A'] = _r3(scorer.pull_deviation())
        # stage 3c: modified nucleotides — local relaxation with their atoms in the contact scoring
        if model.riders or model.pref_chi:
            out['modification_relax'] = relax_modified(model, scorer, rng)
        # stage 4: exact closure from the refined state
        closures = []
        for L in order:
            keep_shape = scorer.dr is not None and L.n >= 1 and (any(r in pk_res for r in L.loop)
                                                                 or all(r in frag_res for r in L.loop))
            fb = {'from_current': 4, 'restraints': scorer.dr,
                  'extra': lambda: scorer.dr_scale * scorer.w.template * scorer.distance_restraint_energy()} \
                if keep_shape else {}
            closures.append(finish_closure(model, L, rng, prior, **fb))
        out['closures'] = closures
        if not all(c.get('closed') for c in closures):
            out['status'] = 'rejected'
            out['reasons'] = ['loop_closure_failed']
            out['time_s'] = round(time.time() - t0, 2)
            return out
        out['relief'] = final_relief(model, rng, prior)
        if scorer.modint is not None:
            out['modification_interactions'] = scorer.modint.report(scorer.base_centroids())
        residues = model.residues_atoms()
        ideal = residues                     # covalent geometry of the torsion model (ideal bonds and angles)
        val = validate_model(seq, residues, hyp['pairs'], model.geom)
        ta = hyp.get('template_atoms')
        if ta and ta.get('atoms'):
            # stage 5: Cartesian refinement onto the template atoms (keeps covalent geometry, removes overlaps)
            from .cartesian import refine_cartesian, template_targets
            targets = template_targets(seq, residues, {int(r): v for r, v in ta['atoms'].items()}, ta['names'],
                                       group={int(r): g for r, g in (ta.get('group') or {}).items()})
            from .cartesian import regularize
            first, released = None, set()
            for attempt in range(max(1, int(cfg.get('template_fit_attempts', 4)))):
                use = {k2: v for k2, v in targets.items() if k2[0] not in released} if released else targets
                refined, fit = refine_cartesian(seq, residues, use, geom=model.geom)
                if fit is None:
                    break
                # polish: ideal covalent geometry of the torsion model, shape of the refined model
                refined, pol = regularize(seq, refined, model.geom, reference=residues)
                fit['polish_moved_rms_A'] = pol['atoms_moved_rms_A']
                val2 = validate_model(seq, refined, hyp['pairs'], model.geom)
                fit['accepted'] = bool(val2['passed'] or (not val['passed'] and len(val2['issues']) <= len(val['issues'])))
                if released:
                    fit['released_residues'] = len(released)
                first = first or (refined, fit, val2)
                if val2['passed']:
                    first = (refined, fit, val2)
                    break
                # the template asks for a geometry the chain cannot take at these residues (typically next to a gap
                # of the alignment): leave them and their neighbours to the covalent terms and fit the rest again
                wider = _release_for_refit(val2['issues'], attempt, released, len(seq))
                if wider is None:
                    break
                released = wider
            if first is not None:
                refined, fit, val2 = first      # the first fit that passes validation, else the fit onto all targets
                if fit['accepted']:
                    residues, val = refined, val2
                out['template_fit'] = fit
        if not val['passed'] and set(val['issue_counts']) <= LOCAL_ISSUES:
            # local strain only (closure junction angles, overlaps): Cartesian regularisation keeps the shape
            from .cartesian import regularize
            fixed, reg = regularize(seq, residues, model.geom, reference=ideal)
            val2 = validate_model(seq, fixed, hyp['pairs'], model.geom)
            reg['issues_before'] = dict(val['issue_counts'])
            reg['passed'] = bool(val2['passed'])
            if val2['passed']:
                residues, val = fixed, val2
            out['regularized'] = reg
        if hyp.get('restraints'):
            rv = validate_model(seq, residues, hyp['restraints'], model.geom, max_report=0)['pair_recovery']
            val['restraint_pair_recovery'] = rv
        detail = scorer(detail=True)
        if scorer.dr is not None:
            out['template_restraints'] = scorer.distance_restraint_stats()
        loop_tor = float(sum(c.get('torsion_energy', 0.0) for c in closures))
        score3d = detail['total'] + loop_tor + 10.0 * val['severe_clashes']
        out.update({
            'status': 'valid' if val['passed'] else 'rejected',
            'reasons': sorted(val['issue_counts']) if not val['passed'] else [],
            'score_terms': {k2: round(v, 3) if isinstance(v, float) else v for k2, v in detail.items()},
            'loop_torsion_energy': round(loop_tor, 3),
            'score_3d': round(score3d, 3),
            'validation': val,
            'residues': [{a: [round(float(c), 3) for c in xyz] for a, xyz in r.items()} for r in residues],
            'model_summary': model.summary(),
        })
        if model.modifications:
            from .modifications import conformation_report
            out['modification_conformations'] = conformation_report(seq, residues, model.modifications)
    except Exception as exc:          # never let one trajectory kill the run; record it
        out['status'] = 'error'
        out['reasons'] = [f'{type(exc).__name__}: {exc}']
    out['time_s'] = round(time.time() - t0, 2)
    return out


# ---------------------------------------------------------------------------- selection
def kabsch_rmsd(a, b):
    a = a - a.mean(axis=0)
    b = b - b.mean(axis=0)
    u, _, vt = np.linalg.svd(a.T @ b)
    d = np.sign(np.linalg.det(u @ vt))
    r = u @ np.diag([1.0, 1.0, d]) @ vt
    return float(np.sqrt(np.mean(np.sum((a @ r - b) ** 2, axis=1))))


def select_diverse(cands, n_models, threshold):
    order = sorted(range(len(cands)), key=lambda i: cands[i].get('rank_key', (1, 0, cands[i]['total_score'])))
    reps = []
    for i in order:
        x = cands[i]['_c4']
        d = [kabsch_rmsd(x, cands[j]['_c4']) for j in reps]
        if not d or min(d) >= threshold:
            if len(reps) < n_models:
                reps.append(i)
    # cluster membership: nearest representative
    for i, c in enumerate(cands):
        d = [kabsch_rmsd(c['_c4'], cands[j]['_c4']) for j in reps]
        k = int(np.argmin(d))
        c['cluster'] = k + 1
        c['rmsd_to_representative'] = round(d[k], 3)
    return reps


# ---------------------------------------------------------------------------- templates
def _template_hypotheses(seq, config, say, lib, scan, meta):
    """Search the built-in structure library; one 2D/3D hypothesis per significant hit."""
    from .templates import search, hypothesis, rank_templates, segment_hits
    meta = {**meta, 'mode': config.templates, 'cutoff': config.template_cutoff or None,
            'excluded': list(config.template_exclude), 'library': lib.version if lib else None}
    if lib is None:
        say('template library not installed; physics-only prediction')
        return [], [], meta
    t0 = time.time()
    hits = search(seq, lib, cutoff=config.template_cutoff or None, exclude=config.template_exclude,
                  max_hits=max(8, config.max_templates), scan=scan, structure=config.template_structure)
    n_found = len(hits)
    hits = rank_templates(seq, hits, lib, config.max_templates)
    meta['search_time_s'] = round(time.time() - t0 + meta.get('scan_time_s', 0.0), 2)
    say(f'template search: {n_found} significant hit(s) in {meta["search_time_s"]} s'
        + (f' (released before {config.template_cutoff})' if config.template_cutoff else ''))
    hyps, out = [], []
    for rank, h in enumerate(hits):
        # the best template may leave a part uncovered (an insert in a scaffold): add templates for that part
        extra = segment_hits(seq, h, lib, cutoff=config.template_cutoff or None,
                             exclude=config.template_exclude) if rank == 0 else []
        hy = hypothesis(seq, h, lib, config.temperature, extra=extra)
        label = '+'.join(f"{x['pdb']}_{x['chain']}" for x in (h, *extra))
        say(f"  template {h['pdb']}_{h['chain']} E={h['evalue']:.2g} identity {h['identity']:.0%} coverage {h['coverage']:.0%} "
            f"pairs {h['pairs_compatible']}/{h['pairs_mapped']} [{h.get('found_by', 'sequence')}] '{h['title'][:60]}'")
        for x in extra:
            say(f"  template {x['pdb']}_{x['chain']} E={x['evalue']:.2g} identity {x['identity']:.0%} coverage "
                f"{len(x['aligned']) / len(seq):.0%} pairs {x['pairs_compatible']}/{x['pairs_mapped']} '{x['title'][:60]}'")
        info = {k: h[k] for k in ('pdb', 'chain', 'model', 'release', 'method', 'resolution', 'title', 'evalue',
                                  'bits', 'identity', 'coverage', 'consensus', 'rank_score', 'pairs_mapped',
                                  'pairs_compatible', 'found_by')}
        if extra:
            info['parts'] = [{k: x[k] for k in ('pdb', 'chain', 'model', 'release', 'method', 'resolution', 'title',
                                                'evalue', 'identity', 'segment')} for x in extra]
            info['coverage'] = round(len({a for x in (h, *extra) for a, _ in x['aligned']}) / len(seq), 3)
        info['label'] = label
        info['aligned_residues'] = sum(len(x['aligned']) for x in (h, *extra))
        info['distance_restraints'] = len(hy['distance_restraints']['i'])
        out.append(info)
        e = hy['energy_kcal_mol']
        hyps.append({'id': '', 'dbn': hy['dbn'], 'pairs': hy['pairs'], 'energy_kcal_mol': e,
                     'source': f'template {label}', 'restraints': hy['restraints'],
                     'pseudoknot_pairs': [(i + 1, j + 1) for i, j in hy['restraints']],
                     'distance_restraints': hy['distance_restraints'], 'torsions': hy['torsions'],
                     'template_atoms': hy['template_atoms'],
                     'template': {**info, 'rank': rank}})
    return out, hyps, meta


def _pseudoknot_hypothesis(seq, hyps, config, ss_meta):
    """Put a 2D with one crossing helix first if it is more stable than the nested MFE (pseudoknots.py)."""
    from . import pseudoknots as PK
    reac = list(config.reactivities) or None
    try:
        ph = PK.hypothesis(seq, config.temperature, reactivities=reac)
    except Exception as exc:                      # the pseudoknot search must never break a prediction
        ss_meta['pseudoknots'] = {'error': f'{type(exc).__name__}: {exc}'}
        return hyps
    ss_meta['pseudoknots'] = {'penalty_kcal_mol': PK.PENALTY, 'found': ph is not None}
    if ph is None:
        return hyps
    if reac:                                      # ranking energy includes the pseudo-energies, as for the others
        import RNA
        md = RNA.md()
        md.temperature = config.temperature
        nested = ''.join(c if c in '().' else '.' for c in ph['dbn'])
        ph['rank_energy_kcal_mol'] = ph['energy_kcal_mol']
        ph['energy_kcal_mol'] = round(RNA.fold_compound(seq, md).eval_structure(nested)
                                      + ph['pseudoknot']['helix_energy_kcal_mol'] + PK.PENALTY, 2)
        ph['source'] += '_probing'
    ss_meta['pseudoknots'].update({'pairs': ph['pseudoknot_pairs'], **ph['pseudoknot']})
    out = [ph] + hyps
    if len(out) > max(config.max_structures, 2):
        out = out[:max(config.max_structures, 2)]
    for k, h in enumerate(out, 1):
        h['id'] = f'ss{k}'
    return out


def _domains():
    from . import domains
    return domains


def _domain_hypothesis(seq, hyps, config, out, say, meta):
    """Fold the domains of the best whole-chain 2D separately; their models become a multi-part template."""
    D = _domains()
    from .annotate import nested_subset
    base = next((h for h in hyps if not h.get('template')), None)
    if base is None:
        return None
    t0 = time.time()
    kept, models, report = D.fold_domains(seq, base['pairs'], config, out, predict, say)
    meta.update({'domains': report, 'based_on': base['id'], 'time_s': round(time.time() - t0, 1)})
    if not kept:
        say('domains: none could be folded; whole-chain search only')
        return None
    asm = D.assembly(seq, kept, models)
    pairs, _crossing = nested_subset(sorted(set(asm['pairs']) | set(D.outer_pairs(len(seq), base['pairs'], kept))))
    dbn = pairs_to_dotbracket(len(seq), pairs, offset=0)
    e = None
    try:
        import RNA
        md = RNA.md()
        md.temperature = config.temperature
        e = round(RNA.fold_compound(seq, md).eval_structure(dbn), 2)
    except Exception:
        pass
    n_tmpl = sum(1 for h in hyps if h.get('template'))
    label = 'domains ' + ', '.join(f'{s + 1}-{e2}' for s, e2 in kept)
    meta['covered'] = round(asm['covered'], 3)
    say(f'{label}: {asm["covered"]:.0%} of the chain folded as domains in {meta["time_s"]} s')
    return {'id': '', 'dbn': dbn, 'pairs': pairs, 'energy_kcal_mol': e, 'source': 'domain assembly',
            'restraints': [], 'pseudoknot_pairs': [], 'distance_restraints': asm['distance_restraints'],
            'torsions': asm['torsions'], 'template_atoms': asm['template_atoms'],
            'template': {'label': label, 'rank': n_tmpl, 'coverage': round(asm['covered'], 3), 'domains': True}}


def _allocate(hyps, n_traj, config):
    """Trajectories per hypothesis: templates first (best hit most), then MFE and suboptimal 2D."""
    alloc = [0] * len(hyps)
    t_idx = [k for k, h in enumerate(hyps) if h.get('template')]
    p_idx = [k for k, h in enumerate(hyps) if not h.get('template')]
    n_t = 0
    if t_idx:
        n_t = n_traj if not p_idx else min(n_traj - 1, max(1, int(round(config.template_share * n_traj))))
        w = np.array([0.6 ** r for r in range(len(t_idx))])
        share = np.floor(w / w.sum() * n_t).astype(int)
        share[0] += n_t - share.sum()
        for k, c in zip(t_idx, share):
            alloc[k] = int(c)
    rest = n_traj - n_t
    if p_idx and rest > 0:
        if len(p_idx) == 1:
            alloc[p_idx[0]] = rest
        else:
            alloc[p_idx[0]] = max(1, int(round(config.mfe_share * rest)))
            for k in range(rest - alloc[p_idx[0]]):
                alloc[p_idx[1 + k % (len(p_idx) - 1)]] += 1
    return alloc


# ---------------------------------------------------------------------------- main entry
def predict(sequence, out_dir, name='query', dbn=None, config: PredictConfig = PredictConfig(), log=print):
    t_start = time.time()
    from . import modifications as modmod
    mods = list(config.modifications or ())
    if modmod.has_notation(sequence):
        sequence, parsed, problems = modmod.parse(sequence)
        if problems:
            raise ValueError('sequence: ' + '; '.join(problems))
        mods = parsed + mods
    seq = normalize_sequence(sequence)
    n = len(seq)
    mods = modmod.normalize(mods, seq)
    out = Path(out_dir)
    for sub in ('input', 'conformations', 'metadata', 'logs'):
        (out / sub).mkdir(parents=True, exist_ok=True)
    logf = (out / 'logs' / 'run.log').open('w', encoding='utf-8')

    def say(msg):
        line = f'[{time.strftime("%H:%M:%S")}] {msg}'
        logf.write(line + '\n')
        logf.flush()
        if log:
            log(line)
    warnings = []
    if n > 200:
        warnings.append('sequence longer than the 200-nt design range; results are untested at this size')
    # -------------------------------------------------------------- 2D hypotheses
    if dbn is not None:
        problem = check_dotbracket(dbn, n)
        if problem:
            raise ValueError(problem)
        pairs = [(i - 1, j - 1) for i, j in dotbracket_to_pairs(dbn)]
        canon = {'AU', 'UA', 'GC', 'CG', 'GU', 'UG'}
        if mods:                      # pairs only a modification allows stay (as restraints); blocked bases are freed
            _rigid, _special, _dropped = modmod.split_pairs(seq, pairs, mods)
            blocked = modmod.blocked_positions(mods)
            freed = [(i, j) for i, j in _dropped if i in blocked or j in blocked]
            if freed:
                warnings.append(f'{len(freed)} pairs of the user 2D involve a nucleotide whose modification blocks '
                                f'Watson-Crick pairing; they were dropped: {[(i + 1, j + 1) for i, j in freed][:10]}')
                pairs = [p for p in pairs if p not in set(freed)]
            pairs = [p for p in pairs if p not in set(_special)]
        noncanon = [(i, j) for i, j in pairs if seq[i] + seq[j] not in canon]
        if noncanon:
            pairs = [p for p in pairs if p not in set(noncanon)]
            warnings.append(f'{len(noncanon)} non-canonical pairs in the user 2D were dropped (only A-U, G-C, G-U '
                            f'helices are modelled): {[(i + 1, j + 1, seq[i] + seq[j]) for i, j in noncanon][:10]}')
        nested, removed = nested_subset(pairs)
        if removed:
            warnings.append(f'{len(removed)} crossing (pseudoknot) pairs are not part of the nested fold tree; '
                            + ('they are applied as base-pair distance restraints' if config.pseudoknot_restraints
                               else 'they were dropped'))
        hyps = [{'id': 'ss1', 'dbn': pairs_to_dotbracket(n, nested, offset=0), 'pairs': nested,
                 'energy_kcal_mol': None, 'source': 'user', 'restraints': removed if config.pseudoknot_restraints else [],
                 'pseudoknot_pairs': [(i + 1, j + 1) for i, j in removed],
                 'pseudoknot_treatment': 'distance restraints' if config.pseudoknot_restraints else 'dropped'}]
        ss_meta = {'program': 'user-supplied'}
        try:
            import RNA
            md = RNA.md()
            md.temperature = config.temperature
            fc_e = modmod.fold_compound(seq, mods, md)[0] if mods else RNA.fold_compound(seq, md)
            ev_dbn = pairs_to_dotbracket(n, sorted(set(nested) | set(_special)), offset=0) if mods else hyps[0]['dbn']
            hyps[0]['energy_kcal_mol'] = round(fc_e.eval_structure(ev_dbn), 2)
        except Exception:
            pass
    else:
        from .ss import predict_structures
        try:
            if config.reactivities and len(config.reactivities) != n:
                raise ValueError(f'probing data: {len(config.reactivities)} values for {n} nucleotides')
            hyps, ss_meta = predict_structures(seq, config.max_structures, config.temperature,
                                               reactivities=list(config.reactivities) or None,
                                               alignment=list(config.alignment) or None,
                                               modifications=mods or None)
        except RuntimeError:
            if not config.seed_structures:
                raise
            hyps, ss_meta = [], {'program': 'seed structures only (ViennaRNA unavailable)'}
        if config.pseudoknots != 'off' and hyps:
            hyps = _pseudoknot_hypothesis(seq, hyps, config, ss_meta)
        if config.seed_structures:
            hyps = _seed_hypotheses(seq, config.seed_structures, config.temperature, hyps)
    # -------------------------------------------------------------- templates (experimental structures)
    template_hits, template_meta = [], {'mode': config.templates}
    lib, scan = None, None
    guided_run = bool(config.seed_structures or config.pull_restraints)
    if (config.templates != 'off' or config.motifs != 'off') and not guided_run:
        from .templates import load_library, sw_scan
        lib = load_library()
        if lib is not None:
            t0 = time.time()
            scan = sw_scan(seq, lib)
            template_meta['scan_time_s'] = round(time.time() - t0, 2)
    if config.templates != 'off' and dbn is None and not guided_run:
        template_hits, tmpl_hyps, template_meta = _template_hypotheses(seq, config, say, lib, scan, template_meta)
        if tmpl_hyps:
            hyps = tmpl_hyps + hyps
            for k2, h in enumerate(hyps, 1):
                h['id'] = f'ss{k2}'
    # -------------------------------------------------------------- long chains: domains folded separately
    domain_meta = {'mode': config.domains}
    if config.domains != 'off' and dbn is None and not guided_run and n >= _domains().MIN_LENGTH and hyps:
        cover = max([float((h.get('template') or {}).get('coverage') or 0.0) for h in hyps] + [0.0])
        if cover >= 0.85:
            domain_meta['skipped'] = 'a template covers the chain'
        else:
            dh = _domain_hypothesis(seq, hyps, config, out, say, domain_meta)
            if dh is not None:
                hyps = [h for h in hyps if h.get('template')] + [dh] + [h for h in hyps if not h.get('template')]
                for k2, h in enumerate(hyps, 1):
                    h['id'] = f'ss{k2}'
    # -------------------------------------------------------------- modified nucleotides: pairs in 3D
    mod_meta = None
    if mods:
        if dbn is not None and _special:
            hyps[0]['pairs'] = sorted(set(hyps[0]['pairs']) | set(_special))
        n_special = n_dropped = 0
        for h in hyps:
            rigid, special, dropped = modmod.split_pairs(seq, h['pairs'], mods)
            h['pairs'] = rigid
            if special:
                h['modification_pairs'] = [(i + 1, j + 1) for i, j in special]
                h['mod_restraints'] = modmod.special_pair_restraints(seq, rigid, special)
            n_special += len(special)
            n_dropped += len([p for p in dropped if p[0] in modmod.blocked_positions(mods) or p[1] in modmod.blocked_positions(mods)])
        mod_meta = {'sites': [{**m, 'name': modmod.BY_CODE[m['code']]['name'], 'ccd': modmod.BY_CODE[m['code']]['ccd'],
                               'parent': modmod.BY_CODE[m['code']]['parent']} for m in mods],
                    'sequence': modmod.format_sequence(seq, mods), 'summary': modmod.summary(mods),
                    'pairs_as_restraints': n_special, 'pairs_dropped_blocked': n_dropped,
                    'method': '2D: ViennaRNA soft constraints (measured parameters where ViennaRNA has them, otherwise '
                              'Turner 2004 stacks of a reference pair plus an estimated correction) and hard constraints '
                              'for blocked Watson-Crick edges; 3D: modifications inside the search — sugar pucker '
                              'and syn/anti preferences in the conformers drawn for non-helical residues and in the '
                              'torsion prior, modification atoms placed from their parent groups and scored for '
                              'contacts at every step, hydrogen bonds, the pseudouridine water bridge and the stacking '
                              'change of the modification scored in the search (modinteract.py), a local relaxation '
                              'around the modified residues, and pairs only '
                              'a modification allows held by A-form distance restraints; helical residues keep the '
                              "A-form C3'-endo sugar",
                    'preferences': {'pucker': {str(r + 1): {'pucker': "C2'-endo" if v[0] == 'C2endo' else "C3'-endo",
                                                            'weight': v[1]}
                                               for r, v in modmod.sampling_preferences(mods)['pucker'].items()},
                                    'glycosidic': {str(r + 1): {'mode': v[0], 'weight': v[1], 'preference_kT': v[2]}
                                                   for r, v in modmod.sampling_preferences(mods)['chi'].items()}}}
        say(f'modified nucleotides: {mod_meta["summary"]}')
    # -------------------------------------------------------------- Z-RNA hypothesis (left-handed helix segment)
    zrna_meta = {'mode': config.zrna}
    if config.zrna != 'off' and hyps and not guided_run:
        from . import zrna as zmod
        zrna_meta.update({'salt_M': config.zrna_salt, 'binder': config.zrna_binder,
                          'propose_below_kcal_mol': zmod.Z_PROPOSE_MAX})
        found = []
        for h in hyps:
            for c in zmod.candidates(seq, h['pairs'], config.zrna_salt, mods, config.zrna_binder):
                found.append((c['dG_kcal_mol'], h['id'], c, h))
        found.sort(key=lambda x: x[0])
        zrna_meta['candidates'] = [{'hypothesis': hid, 'pairs': [(i + 1, j + 1) for i, j in c['pairs']],
                                    'dG_kcal_mol': dg, 'z_fraction_estimate': c['z_fraction_estimate']}
                                   for dg, hid, c, _ in found[:5]]
        if found and found[0][0] <= zmod.Z_PROPOSE_MAX:
            dg, hid, c, h = found[0]
            zh = zmod.hypothesis(seq, h, c)
            zh['id'] = f'ss{len(hyps) + 1}'
            hyps.append(zh)
            zrna_meta['hypothesis'] = zh['id']
            say(f"Z-RNA: {zh['id']} has pairs {zh['zrna']['pairs'][0]}..{zh['zrna']['pairs'][-1]} left-handed "
                f"(estimated dG(A->Z) {dg:+.1f} kcal/mol)")
        elif found:
            zrna_meta['note'] = (f'best segment dG(A->Z) {found[0][0]:+.1f} kcal/mol is above '
                                 f'{zmod.Z_PROPOSE_MAX} kcal/mol: no Z-RNA hypothesis')
        else:
            zrna_meta['note'] = 'no helix with >= 4 alternating purine-pyrimidine Watson-Crick pairs'
    say(f'{name}: {n} nt; {len(hyps)} secondary-structure hypotheses ({ss_meta.get("program")})')
    for h in hyps:
        say(f"  {h['id']} {h['dbn']} {h.get('energy_kcal_mol')} kcal/mol [{h['source']}]")
    # -------------------------------------------------------------- trajectories
    n_traj = config.trajectories or max(8, 4 * config.n_models)
    z_ids = [k2 for k2, h in enumerate(hyps) if h.get('z_pairs')]
    alloc = _allocate([h for h in hyps if not h.get('z_pairs')], n_traj, config)
    for k2 in z_ids:                                     # extra trajectories: the other hypotheses keep theirs
        alloc.insert(k2, max(2, int(round(0.25 * n_traj))))
    cfg_dict = asdict(config)
    cfg_dict['modifications'] = mods
    if config.alignment:                      # recorded in input/ and in the 2D method; not copied into every job
        cfg_dict['alignment'] = [f'{len(config.alignment)} aligned sequences']
    motif_mask, motif_cache, motif_meta = None, {}, {'mode': config.motifs, 'fragments': config.fragments}
    if config.motifs != 'off' and lib is not None:
        from .templates import homologs
        motif_mask = lib.allowed(config.template_cutoff or None, config.template_exclude) & ~homologs(seq, lib, scan)
        motif_meta.update({'library_chains': int(motif_mask.sum()), 'homologs_excluded': int(homologs(seq, lib, scan).sum())})
    jobs = []
    k = 0
    for h, cnt in zip(hyps, alloc):
        for _ in range(cnt):
            hyp_job = {'id': h['id'], 'pairs': h['pairs'], 'restraints': h.get('restraints', []),
                       'distance_restraints': h.get('distance_restraints'), 'torsions': h.get('torsions'),
                       'template_atoms': h.get('template_atoms'),
                       'template_coverage': (h.get('template') or {}).get('coverage'),
                       'mod_pairs': [(i - 1, j - 1) for i, j in h.get('modification_pairs', ())],
                       'z_pairs': h.get('z_pairs') or []}
            if motif_mask is not None and not h.get('template'):
                from .motifs import junction_restraints
                tors = {} if config.fragments != 'off' else None
                dr, used = junction_restraints(seq, h['pairs'], lib, motif_mask, np.random.default_rng([config.seed, k, 7]),
                                               cache=motif_cache, torsions=tors)
                if config.fragments != 'off':
                    from .fragments import hairpin_fragments, merge_restraints
                    hdr, htors, _loops, hused = hairpin_fragments(seq, h['pairs'], lib, motif_mask,
                                                                  np.random.default_rng([config.seed, k, 11]),
                                                                  cache=motif_cache)
                    dr = merge_restraints(dr, hdr)
                    tors.update(htors)
                    used = used + hused
                    if tors:
                        hyp_job['fragment_torsions'] = tors
                if h.get('restraints') and h.get('source', '').startswith('vienna_pseudoknot'):
                    from .pseudoknots import helix_restraints
                    from .fragments import merge_restraints
                    dr = merge_restraints(dr, helix_restraints(seq, h['restraints']))
                if dr['i']:
                    hyp_job['distance_restraints'] = dr
                    hyp_job['motifs'] = used
            if h.get('mod_restraints'):
                from .fragments import merge_restraints
                hyp_job['distance_restraints'] = merge_restraints(hyp_job.get('distance_restraints'), h['mod_restraints'])
            jobs.append({'sequence': seq, 'hypothesis': hyp_job, 'index': k, 'seed': config.seed, 'config': cfg_dict})
            k += 1
    say(f'running {len(jobs)} trajectories on {config.workers} worker(s)')
    if config.workers > 1:
        import multiprocessing as mp
        for var in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
            os.environ.setdefault(var, '1')   # inherited by workers; libraries loaded later honour it
        # fork is fast and safe on Linux; macOS app bundles and Windows need spawn
        method = os.environ.get('RNAFOLD3D_MP_START') or ('fork' if sys.platform.startswith('linux') else 'spawn')
        ctx = mp.get_context(method)
        results = []
        # leaving the block for any reason (finished, error, Ctrl+C, a stop signal) terminates the workers
        with _StopOnSignal(), ctx.Pool(config.workers, initializer=_worker_init) as pool:
            for r in pool.imap_unordered(run_trajectory, jobs):
                results.append(r)
                say(f'  progress {len(results)}/{len(jobs)} trajectories finished')
        results.sort(key=lambda r: r['trajectory'])      # completion order must not affect the result
    else:
        results = []
        for j in jobs:
            results.append(run_trajectory(j))
            say(f'  progress {len(results)}/{len(jobs)} trajectories finished')
    # -------------------------------------------------------------- ranking & selection
    # 2D energy used for ranking: with probing data the energy that includes the pseudo-energies
    rank_e = lambda h: h.get('rank_energy_kcal_mol', h.get('energy_kcal_mol'))
    e_mfe = min((rank_e(h) for h in hyps if rank_e(h) is not None), default=None)
    hyp_by_id = {h['id']: h for h in hyps}
    valid = []
    for r in results:
        say(f"  traj {r['trajectory']:3d} [{r['ss_id']}] {r['status']:8s} score3d={r.get('score_3d')} "
            f"{'' if r['status'] == 'valid' else r.get('reasons')} {r.get('time_s')}s")
        if r['status'] != 'valid':
            continue
        h = hyp_by_id[r['ss_id']]
        dg = rank_e(h)
        if h.get('template'):
            r['ss_offset_kT'] = 0.0          # 2D taken from an experimental structure
            r['rank_key'] = (0, h['template']['rank'], r['score_3d'])
        else:
            r['ss_offset_kT'] = round((dg - e_mfe) * KCAL_TO_KT_37C, 3) if (dg is not None and e_mfe is not None) else 0.0
        # junction geometry from other structures only guides the sampling: its restraint energy is not part of the
        # score that ranks models (otherwise models of hypotheses without junctions would always rank first)
        guide = float((r.get('score_terms') or {}).get('template') or 0.0) if r.get('junction_motifs') else 0.0
        if guide:
            r['guide_energy_excluded'] = round(guide, 3)
        r['total_score'] = r['score_3d'] - guide + r['ss_offset_kT']
        r.setdefault('rank_key', (1, 0, r['total_score']))
        r['_residues'] = [{a: np.array(xyz) for a, xyz in res.items()} for res in r['residues']]
        r['_c4'] = _c4(r['_residues'])
        valid.append(r)
    # Z-RNA hypotheses: the 3D score is calibrated on right-handed helices and scores a left-handed stem tens of
    # units worse, whatever its free energy. The estimated dG(A->Z) is already in the 2D energy of the hypothesis
    # (ss_offset_kT). What the 3D score adds for the helix form alone is taken out by setting the best Z model level
    # with the best model of the hypothesis it was derived from: the Z model then ranks first exactly when the
    # estimate favours Z (dG < 0) and next to the first models when Z is a minor form, instead of always last.
    for h in hyps:
        zs = [r for r in valid if r['ss_id'] == h['id']] if h.get('z_pairs') else []
        ref = [r for r in valid if r['ss_id'] == h.get('parent') and r['rank_key'][0] == 1] or \
              [r for r in valid if r['rank_key'][0] == 1 and not hyp_by_id[r['ss_id']].get('z_pairs')]
        if zs and ref:
            shift = (min(r['total_score'] - r['ss_offset_kT'] for r in zs)
                     - min(r['total_score'] - r['ss_offset_kT'] for r in ref))
            for r in zs:
                r['helix_form_score_excluded'] = round(shift, 3)
                r['total_score'] = r['total_score'] - shift
                r['rank_key'] = (1, 0, r['total_score'])
    # template models: the best model of every template first (templates in rank order), then the second-best of
    # each, so that the exported models cover the different templates instead of one template several times
    by_template = {}
    for r in valid:
        if r['rank_key'][0] == 0:
            by_template.setdefault(r['rank_key'][1], []).append(r)
    _order_template_models(by_template, config.template_fit_attempts > 1)
    rg0 = config.weights.rg_coefficient * n ** (1 / 3)
    threshold = config.cluster_rmsd or float(np.clip(0.15 * rg0, 2.0, 6.0))
    reps = select_diverse(valid, config.n_models, threshold) if valid else []
    if len(reps) < config.n_models:
        warnings.append(f'only {len(reps)} distinct valid conformers (requested {config.n_models}); '
                        f'{len(valid)} of {len(results)} trajectories produced valid models')
    # -------------------------------------------------------------- write
    (out / 'input' / 'sequence.fasta').write_text(f'>{name}\n{seq}\n')
    if mods:
        (out / 'input' / 'sequence_modified.txt').write_text(f'>{name}\n{mod_meta["sequence"]}\n')
    with (out / 'input' / 'secondary_structures.dbn').open('w', encoding='utf-8') as f:
        for h in hyps:
            f.write(f">{name}|{h['id']}|{h['source']}|{h.get('energy_kcal_mol')}\n{seq}\n{h['dbn']}\n")
    conformers = []
    for rank, i in enumerate(reps, 1):
        r = valid[i]
        cid = f'conf_{rank:03d}'
        h = hyp_by_id[r['ss_id']]
        remarks = [f'RNAFOLD3D {ENGINE} {ENGINE_VERSION} PREDICTED MODEL {cid} RANK {rank}',
                   f'SECONDARY STRUCTURE {h["id"]} ({h["source"]}) {h.get("energy_kcal_mol")} KCAL/MOL',
                   f'MODEL SCORE {r["total_score"]:.2f} (ARBITRARY UNITS, NOT A FREE ENERGY)',
                   'PHYSICS-INSPIRED PREDICTION; NOT AN EXPERIMENTAL STRUCTURE']
        residues_out, comps, mod_report = r['_residues'], None, None
        if mods:
            residues_out, comps, mod_report = modmod.decorate(seq, r['_residues'], mods,
                                                              list(h['pairs']) + [(i - 1, j - 1) for i, j in h.get('modification_pairs', ())],
                                                              sampled=True)
            mod_report['conformations'] = r.get('modification_conformations') or []
            if r.get('modification_relax'):
                mod_report['relaxation'] = r['modification_relax']
            if r.get('modification_interactions'):
                mod_report['interactions'] = _written_atom_names(r['modification_interactions'], mods)
            remarks.append(f'MODIFIED NUCLEOTIDES {mod_meta["summary"]}'[:69])
        from .zrna import detect as _zdetect
        z_regions = _zdetect(seq, r['_residues'], h['pairs'])
        if z_regions:
            remarks.append(f"Z-RNA (LEFT-HANDED) PAIRS {z_regions[0]['pairs'][0]}..{z_regions[0]['pairs'][-1]}"[:69])
        write_pdb(out / 'conformations' / f'{cid}.pdb', residues_out, seq, remarks=remarks, comps=comps)
        write_mmcif(out / 'conformations' / f'{cid}.cif', residues_out, seq, entry_id=cid.upper(),
                    title=f'{name} {cid} predicted by {ENGINE} {ENGINE_VERSION}', comps=comps)
        cluster_size = sum(1 for c in valid if c['cluster'] == rank)
        conformers.append({
            'id': cid, 'rank': rank, 'files': {'pdb': f'conformations/{cid}.pdb', 'mmcif': f'conformations/{cid}.cif'},
            'secondary_structure_id': h['id'], 'secondary_structure': h['dbn'],
            **({'template': h['template']['label'], 'template_restraints': r.get('template_restraints'),
                'template_fit': r.get('template_fit')} if h.get('template') else {}),
            **({'junction_motifs': r['junction_motifs'], 'motif_restraints': r.get('template_restraints'),
                'guide_energy_excluded': r.get('guide_energy_excluded', 0.0)} if r.get('junction_motifs') else {}),
            'ss_free_energy_kcal_mol': h.get('energy_kcal_mol'),
            **({'helix_form_score_excluded': r['helix_form_score_excluded'],
                'z_fraction_estimate': h['zrna'].get('z_fraction_estimate')} if 'helix_form_score_excluded' in r else {}),
            'total_score': round(r['total_score'], 3), 'score_3d': r['score_3d'], 'ss_offset_kT': r['ss_offset_kT'],
            'score_terms': r['score_terms'], 'trajectory': r['trajectory'], 'seed': r['seed'],
            'cluster_id': rank, 'cluster_size': cluster_size,
            'validation': {k2: r['validation'].get(k2) for k2 in ('passed', 'severe_clashes',
                           'heavy_atom_overlaps_per_1000_atoms', 'bond_outliers_4sd', 'angle_outliers_4sd',
                           'sugar_puckers', 'pair_recovery', 'restraint_pair_recovery')},
            'loop_closures': r['closures'],
            **({'pull': r['pull']} if r.get('pull') else {}),
            **({'modifications_3d': mod_report} if mod_report else {}),
            'zrna_regions': z_regions,
        })
    if config.save_candidates:
        cdir = out / 'candidates'
        cdir.mkdir(exist_ok=True)
        for r in valid:
            write_pdb(cdir / f"traj_{r['trajectory']:03d}.pdb", r['_residues'], seq,
                      remarks=[f"CANDIDATE TRAJECTORY {r['trajectory']} SS {r['ss_id']} SCORE {r['total_score']:.2f}"])
    scores_path = out / 'metadata' / 'scores.csv'
    with scores_path.open('w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['trajectory', 'ss_id', 'status', 'selected_as', 'total_score', 'score_3d', 'ss_offset_kT',
                    'clash', 'stacking', 'compaction', 'electro', 'torsion', 'rg', 'severe_clashes',
                    'cluster', 'rmsd_to_representative', 'time_s', 'reasons'])
        sel = {valid[i]['trajectory']: f'conf_{k:03d}' for k, i in enumerate(reps, 1)}
        for r in results:
            st = r.get('score_terms', {})
            w.writerow([r['trajectory'], r['ss_id'], r['status'], sel.get(r['trajectory'], ''), r.get('total_score', ''),
                        r.get('score_3d', ''), r.get('ss_offset_kT', ''), st.get('clash', ''), st.get('stacking', ''),
                        st.get('compaction', ''), st.get('electro', ''), st.get('torsion', ''), st.get('rg', ''),
                        r.get('validation', {}).get('severe_clashes', ''), r.get('cluster', ''),
                        r.get('rmsd_to_representative', ''), r.get('time_s', ''), ';'.join(r.get('reasons', []))])
    clusters = {'method': 'greedy score-ordered selection; C4\' Kabsch RMSD threshold', 'threshold_angstrom': round(threshold, 3),
                'clusters': [{'id': k, 'representative': f'conf_{k:03d}',
                              'members': [c['trajectory'] for c in valid if c['cluster'] == k]} for k in range(1, len(reps) + 1)],
                'note': 'Cluster sizes are sampling statistics, not thermodynamic populations.'}
    (out / 'metadata' / 'clusters.json').write_text(json.dumps(clusters, indent=2))
    validation = {r['trajectory']: {'status': r['status'], 'reasons': r.get('reasons', []),
                                    'validation': r.get('validation'), 'closures': r.get('closures')} for r in results}
    (out / 'metadata' / 'validation.json').write_text(json.dumps(validation, indent=1, default=str))
    from .build import load_geometry
    geom = load_geometry()
    manifest = {
        'schema': 'rnafold3d.ensemble/1',
        'created': datetime.datetime.now().isoformat(timespec='seconds'),
        'software': {'name': 'rnafold3d', 'package_version': PACKAGE_VERSION, 'engine': ENGINE,
                     'engine_version': ENGINE_VERSION, 'geometry_parameters': geom.version,
                     'python': platform.python_version()},
        'input': {'name': name, 'sequence': seq, 'length': n, 'sha256': hashlib.sha256(seq.encode()).hexdigest(),
                  'secondary_structure_mode': 'user' if dbn is not None else 'predicted_from_sequence',
                  'secondary_structure_method': ss_meta,
                  'secondary_structures': [{k2: v for k2, v in h.items() if k2 not in ('pairs', 'restraints', 'distance_restraints', 'torsions', 'mod_restraints')}
                                           for h in hyps],
                  'templates': {**template_meta, 'hits': template_hits}, 'junction_motifs': motif_meta, 'domains': domain_meta,
                  **({'modifications': mod_meta} if mod_meta else {}), 'zrna': zrna_meta},
        'parameters': {**{k2: v for k2, v in cfg_dict.items() if k2 != 'weights'}, 'weights': cfg_dict['weights'],
                       'trajectories_run': len(jobs), 'trajectory_allocation': dict(zip([h['id'] for h in hyps], alloc)),
                       'cluster_rmsd_threshold': round(threshold, 3), 'ionic_conditions': 'implicit; Debye length 10 A '
                       '(~0.1 M monovalent); no explicit Mg2+'},
        'chain': {'id': 'A', 'numbering': '1..N matches input sequence positions', 'atoms': 'heavy atoms only; '
                  "5'-terminal residue has no phosphate"},
        'conformers': conformers,
        'statistics': {'trajectories': len(results), 'valid': len(valid),
                       'rejected': sum(r['status'] == 'rejected' for r in results),
                       'errors': sum(r['status'] == 'error' for r in results),
                       'wall_time_s': round(time.time() - t_start, 1)},
        'warnings': warnings,
        'notes': ['Scores are dimensionless model scores (kT-like units), not free energies; ranks are not populations.',
                  'Secondary-structure energies are ViennaRNA nearest-neighbour free energies (kcal/mol).',
                  ('Models built on a template follow an experimental structure of a related RNA (see input.templates); '
                   'the others are physics-inspired predictions. No deep learning is used.' if template_hits else
                   'Models are physics-inspired predictions without deep learning or templates; '
                   'accuracy must be judged from the benchmark report.')],
    }
    (out / 'metadata' / 'manifest.json').write_text(json.dumps(manifest, indent=2))
    (out / 'config.json').write_text(json.dumps(cfg_dict, indent=2))
    say(f'{len(valid)}/{len(results)} valid models; exported {len(reps)} conformers to {out} '
        f'in {time.time() - t_start:.1f}s')
    for wmsg in warnings:
        say('WARNING: ' + wmsg)
    logf.close()
    return manifest
