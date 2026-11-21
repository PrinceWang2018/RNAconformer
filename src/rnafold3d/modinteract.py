"""Interactions made by the modified nucleotides themselves, scored in the 3D search (0.7.6).

The search has no general hydrogen-bond or loop-stacking term (helices come from the 2D, whose free energy already
holds their pairing and stacking, including the measured parameters of m6A, pseudouridine, inosine ...). This module
adds only what is specific to a modification, so that nothing is counted twice:

* hydrogen bonds of the polar groups a modification brings (side-chain carbonyls, hydroxyls, amines, carboxylates of
  t6A, Q, mnm5U, cmo5U, ac4C ...; thio groups as weak acceptors) and of the donors it creates on the parent ring
  (N1-H of pseudouridine, N7-H of 8-oxo purines), with any atom of another nucleotide that has the matching role;
  the roles of the parent atoms follow the modification (a methyl on G N1 or U N3 removes the donor, two methyls on
  an amino group too, a methyl on a ring N acceptor removes it);
* the water bridge of pseudouridine: N1-H ... water ... phosphate of its own or the preceding nucleotide (Arnez and
  Steitz 1994; absent in m1Ψ, whose N1 carries the methyl);
* the change in base stacking for residues outside helices: stronger for Ψ, m1Ψ, m5C, m6A (unpaired), i6A/t6A and
  their relatives, wybutosine, 5-halo-pyrimidines; weaker for dihydrouridine (non-aromatic ring).

Geometry is reduced to distances (hydrogen bonds; no angle terms) and to centroid distance, parallel planes and
overlap along the normal (stacking). Energies are in the kT-like units of the search; the strengths are estimates.
"""
from __future__ import annotations

import numpy as np

from .build import BASE_ATOMS

W_HB = 0.8            # per hydrogen bond of a modification atom
W_WATER = 0.6         # pseudouridine N1-H water bridge
HB_NO = (2.6, 3.1, 3.5)      # full reward between the first two distances, none beyond the third
HB_S = (3.0, 3.5, 3.9)
WB = (4.3, 4.4, 5.4, 6.0)    # N1 ... OP distance window of a water-mediated contact

# stacking change for residues outside helices (kT; negative: stacks less than the parent)
STACKING = {'Y': 0.5, 'Ym': 0.5, 'm1Y': 0.6, 'm3Y': 0.4, 'm5C': 0.4, 'm5Cm': 0.4, 'm5U': 0.3, 'm5Um': 0.3,
            's2U': 0.3, 's2Um': 0.3, 'm5s2U': 0.4, 'm6A': 0.4, 'm6Am': 0.4, 'i6A': 0.8, 'ms2i6A': 1.0, 't6A': 0.8,
            'ms2t6A': 0.9, 'm6t6A': 0.8, 'yW': 1.2, 'imG': 1.0, 'Q': 0.4, '5BrU': 0.5, '5IU': 0.6, '5BrC': 0.5,
            '2AP': 0.3, 'DAP': 0.3, 'D': -0.8}
STACK_D = (3.0, 3.3, 4.1, 4.7)

PARENT_ACC = {'A': {'N1', 'N3', 'N7'}, 'G': {'O6', 'N3', 'N7'}, 'C': {'O2', 'N3'}, 'U': {'O2', 'O4'}}
PARENT_DON = {'A': {'N6'}, 'G': {'N1', 'N2'}, 'C': {'N4'}, 'U': {'N3'}}
BACKBONE_ACC = {'OP1', 'OP2', "O3'", "O5'", "O4'", "O2'"}
BACKBONE_DON = {"O2'"}
AMINO = {'N6', 'N2', 'N4'}
# rider atoms: N donors unless listed as acceptors; O acceptors, also donors when hydroxyl; S weak acceptors
RIDER_N_ACCEPTOR = {('k2C', 'N2')}
HYDROXYL = {('t6A', 'O13'), ('ms2t6A', 'O13'), ('m6t6A', 'O13'), ('hm5C', 'O5B'), ('ho5C', 'O5A'),
            ('ho5U', 'O5A'), ('Q', 'O15'), ('Q', 'O16')}
RING = ('C2', 'C4', 'C6')


def _ramp(d, lo, hi, far):
    """1 between lo and hi, linear to 0 at far (and below lo as well, the clash term takes over)."""
    return np.clip(np.minimum((d - (lo - 0.3)) / 0.3, (far - d) / (far - hi)), 0.0, 1.0)


class ModInteractions:
    """Interaction energy of the modifications of a FoldModel (riders and parent atoms)."""

    def __init__(self, model):
        from .modifications import BY_CODE, PSI_MAP
        m = self.m = model
        I = model.idx
        n = model.n
        seq = model.seq
        by_res = {}
        for x in model.modifications:
            e = BY_CODE.get(x['code'])
            if e is not None and e['kind'] != 'backbone':
                by_res[x['pos'] - 1] = e
        # roles of every atom (parent tables with the modification's changes, then riders)
        acc, don = {}, {}
        for r, b in enumerate(seq):
            for nm in BACKBONE_ACC:
                acc[(r, nm)] = True
            for nm in BACKBONE_DON:
                don[(r, nm)] = True
            for nm in PARENT_ACC[b]:
                acc[(r, nm)] = True
            for nm in PARENT_DON[b]:
                don[(r, nm)] = True
        special = []                                  # (atom key, label): donors a modification creates
        self.water = []                               # (donor key, residue)
        for r, e in by_res.items():
            b = seq[r]
            methyl_on = {}
            for op in e['atoms']:
                if op[0] in ('z', 'ring') and str(op[1]).startswith('C'):
                    methyl_on[op[2]] = methyl_on.get(op[2], 0) + 1
            psi = ('psi',) in e['atoms']
            if psi:
                # parent names of the isomer: the new N1 sits where U has C5; ring O2/O4 trade places (both acceptors)
                inv = {v: k for k, v in PSI_MAP.items()}
                methyl_parent = {inv.get(k, k): v for k, v in methyl_on.items()}
                if not methyl_parent.get('C5'):
                    special.append(((r, 'C5'), 'N1'))
                    don[(r, 'C5')] = True
                    self.water.append(((r, 'C5'), r))
                if methyl_parent.get('N3'):
                    don.pop((r, 'N3'), None)
            else:
                for nm, k in methyl_on.items():
                    if nm in AMINO:
                        if k >= 2:
                            don.pop((r, nm), None)
                    elif nm in PARENT_DON[b]:
                        don.pop((r, nm), None)
                    if nm in PARENT_ACC[b] and nm not in AMINO:
                        acc.pop((r, nm), None)
            if e['code'] in ('8oxoG', '8oxoA'):
                acc.pop((r, 'N7'), None)
                don[(r, 'N7')] = True
                special.append(((r, 'N7'), 'N7'))
        q = []                                        # query atoms: (index, donor, acceptor, weight, sulfur, res, label)
        for rd in model.riders:
            code = by_res.get(rd['res'], {}).get('code') if rd['group'] == 'base' else None
            el = model.element[rd['k']]
            if el not in ('N', 'O', 'S'):
                continue
            if rd['group'] == 'phosphate':
                continue                              # phosphorothioate S: the phosphate's role is unchanged
            if el == 'S':
                d_, a_, w_, s_ = False, True, 0.5, True
            elif el == 'N':
                a_ = (code, rd['name']) in RIDER_N_ACCEPTOR
                d_, w_, s_ = not a_, 1.0, False
            else:
                d_, a_, w_, s_ = (code, rd['name']) in HYDROXYL, True, 1.0, False
            q.append((rd['k'], d_, a_, w_, s_, rd['res'], rd['name']))
            acc[(rd['res'], rd['name'])] = a_
            don[(rd['res'], rd['name'])] = d_
        for key, label in special:
            k = I.get(key)
            if k is not None:
                q.append((k, True, False, 1.0, False, key[0], label))
        ghost = model.ghost
        part = sorted({I[k] for k in set(acc) | set(don) if k in I and (acc.get(k) or don.get(k))})
        part = [k for k in part if not ghost[k]]
        self.p_idx = np.array(part, int)
        self.p_res = model.atom_res[self.p_idx] if len(part) else np.zeros(0, int)
        keys = {v: k for k, v in I.items()}
        self.p_acc = np.array([bool(acc.get(keys[k])) for k in part], bool)
        self.p_don = np.array([bool(don.get(keys[k])) for k in part], bool)
        self.q = q
        if q:
            self.q_idx = np.array([t[0] for t in q], int)
            self.q_don = np.array([t[1] for t in q], bool)
            self.q_acc = np.array([t[2] for t in q], bool)
            self.q_w = np.array([t[3] for t in q])
            self.q_s = np.array([t[4] for t in q], bool)
            self.q_res = np.array([t[5] for t in q], int)
            # compatible donor/acceptor pairs of a query atom and a partner of another residue
            self.compat = ((self.q_don[:, None] & self.p_acc[None, :]) | (self.q_acc[:, None] & self.p_don[None, :])) \
                & (self.q_res[:, None] != self.p_res[None, :])
        # water bridges: N1-H of pseudouridine to OP1/OP2 of its own and the preceding phosphate
        self.wb = []
        for key, r in self.water:
            ops = [I[(rr, nm)] for rr in (r, r - 1) if rr >= 1 for nm in ('OP1', 'OP2') if (rr, nm) in I]
            if key in I and ops:
                self.wb.append((I[key], np.array(ops, int), r))
        # stacking change for residues outside helices
        self.ring = np.array([[I[(r, a)] for a in RING] for r in range(n)], int)
        self.stack = [(r, STACKING[e['code']]) for r, e in by_res.items()
                      if e['code'] in STACKING and model.res_body[r] < 0]
        self.active = bool(self.q or self.wb or self.stack)

    # ------------------------------------------------------------------ terms
    def _hbonds(self):
        x, pres = self.m.xyz, self.m.present
        qok = pres[self.q_idx]
        pok = pres[self.p_idx]
        if not qok.any() or not pok.any():
            return np.zeros(len(self.q)), np.full(len(self.q), -1), np.zeros(len(self.q))
        D = np.linalg.norm(x[self.q_idx][:, None, :] - x[self.p_idx][None, :, :], axis=2)
        g = np.where(self.q_s[:, None], _ramp(D, *HB_S), _ramp(D, *HB_NO))
        g = np.where(self.compat & qok[:, None] & pok[None, :], g, 0.0)
        best = g.argmax(axis=1)
        return g[np.arange(len(self.q)), best], best, D[np.arange(len(self.q)), best]

    def _water(self):
        x, pres = self.m.xyz, self.m.present
        out = []
        for k, ops, r in self.wb:
            ok = pres[ops]
            if not pres[k] or not ok.any():
                out.append((0.0, None, None))
                continue
            d = np.linalg.norm(x[ops[ok]] - x[k], axis=1)
            lo0, lo, hi, far = WB
            g = np.clip(np.minimum((d - lo0) / (lo - lo0), (far - d) / (far - hi)), 0.0, 1.0)
            j = int(np.argmax(g))
            out.append((float(g[j]), int(ops[ok][j]), float(d[j])))
        return out

    def _stacking(self, bc):
        x, pres = self.m.xyz, self.m.present
        rok = pres[self.ring].all(axis=1)
        out = []
        if not self.stack:
            return out
        a, b, c = x[self.ring[:, 0]], x[self.ring[:, 1]], x[self.ring[:, 2]]
        nrm = np.cross(b - a, c - a)
        nrm /= np.linalg.norm(nrm, axis=1)[:, None] + 1e-12
        for r, s in self.stack:
            if not rok[r]:
                out.append((r, s, 0.0, []))
                continue
            v = bc - bc[r]
            d = np.linalg.norm(v, axis=1)
            cand = np.where(rok & (d < STACK_D[3]) & (np.arange(len(d)) != r))[0]
            if not len(cand):
                out.append((r, s, 0.0, []))
                continue
            lo0, lo, hi, far = STACK_D
            dc = d[cand]
            gd = np.clip(np.minimum((dc - lo0) / (lo - lo0), (far - dc) / (far - hi)), 0.0, 1.0)
            par = np.clip((np.abs(nrm[cand] @ nrm[r]) - 0.65) / 0.25, 0.0, 1.0)
            over = np.clip((np.abs(v[cand] @ nrm[r]) / dc - 0.6) / 0.25, 0.0, 1.0)
            f = gd * par * over
            tot = float(f.sum())
            partners = [(int(j), round(float(fj), 2)) for j, fj in zip(cand, f) if fj > 0]
            out.append((r, s, min(2.0, tot), partners))
        return out

    def energy(self, bc):
        """Interaction energy (negative: favourable) for the current coordinates; bc: base centroids."""
        if not self.active:
            return 0.0
        self.m.update_riders()
        e = 0.0
        if self.q:
            g, _, _ = self._hbonds()
            e -= W_HB * float(np.dot(g, self.q_w))
        for g, _, _ in self._water():
            e -= W_WATER * g
        for _r, s, tot, _p in self._stacking(bc):
            e -= s * tot
        return e

    def report(self, bc):
        """Per modified residue: hydrogen bonds of its atoms, water bridge, stacking (for the manifest)."""
        m = self.m
        rows = {}

        def row(r):
            return rows.setdefault(r, {'residue': r + 1, 'hbonds': [], 'water_bridge': None, 'stacking': None})
        m.update_riders()
        if self.q:
            g, best, dist = self._hbonds()
            for t, (gi, bi, di) in enumerate(zip(g, best, dist)):
                if gi > 0.5:
                    k = self.p_idx[bi]
                    row(self.q[t][5])['hbonds'].append({'atom': self.q[t][6], 'partner': f'{int(m.atom_res[k]) + 1}:'
                                                        f'{m.atom_name[k]}', 'distance': round(float(di), 2)})
        for (k, _ops, r), (g, j, d) in zip(self.wb, self._water()):
            row(r)['water_bridge'] = ({'phosphate_oxygen': f'{int(m.atom_res[j]) + 1}:{m.atom_name[j]}',
                                       'distance': round(d, 2)} if g > 0.5 else False)
        for r, s, tot, partners in self._stacking(bc):
            row(r)['stacking'] = {'change_kT': s, 'stacked_faces': round(float(tot), 2),
                                  'partners': [j + 1 for j, f in partners if f > 0.3]}
        return [rows[r] for r in sorted(rows)]
