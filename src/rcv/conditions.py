"""Solution conditions for the secondary-structure analyses (Switch and Pathway tabs).

One place builds the ViennaRNA energy model for all of them:
  * temperature (°C);
  * monovalent salt: Na+ concentration in mol/L (ViennaRNA salt correction; 1.021 M is the standard state of the
    Turner parameters);
  * ligand binding to an aptamer that forms a hairpin or interior loop (ViennaRNA hairpin/interior-loop motif soft
    constraint). With a dissociation constant Kd and a ligand concentration c, every structure that forms the
    binding-competent motif gains the binding free energy  -RT ln(1 + c/Kd)  (two-state binding, ligand in excess);
  * chemical probing data (SHAPE or DMS reactivities) as pseudo-energies after Deigan et al. 2009. These are not
    physical energies: they are used for data-informed populations only, never for barriers or rates.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

R_KCAL = 0.0019872036            # kcal/(mol K)
STANDARD_SALT = 1.021            # mol/L, ViennaRNA default (Turner 2004 standard conditions)

# Aptamers that bind in a hairpin or interior loop (the shapes the ViennaRNA motif constraint supports).
LIGAND_PRESETS = {
    'theophylline': {
        'name': 'Theophylline aptamer', 'motif_seq': 'GAUACCAG&CCCUUGGCAGC', 'motif_struct': '(...((((&)...)))...)',
        'kd_uM': 0.32,
        'reference': 'Zimmermann et al. 1997 (Kd 0.32 µM); motif as in the ViennaRNA documentation (Hofacker lab)'},
}


@dataclass
class Conditions:
    temperature: float = 37.0
    sodium: float = STANDARD_SALT
    ligand: dict | None = None                  # motif_seq, motif_struct, kd_uM (or dg), conc_uM, name
    probing: dict | None = None                 # kind ('SHAPE'|'DMS'), values [float|None per nucleotide], m, b
    notes: list = field(default_factory=list)

    @classmethod
    def from_json(cls, d, n=None):
        d = d or {}
        t = float(d.get('temperature', 37.0))
        if not -10.0 <= t <= 100.0:
            raise ValueError('temperature must be between -10 and 100 °C')
        na = float(d.get('sodium', STANDARD_SALT) or STANDARD_SALT)
        if not 0.005 <= na <= 3.0:
            raise ValueError('Na+ concentration must be between 0.005 and 3 mol/L')
        lig = d.get('ligand') or None
        if lig:
            preset = LIGAND_PRESETS.get(str(lig.get('preset') or ''))
            lig = {**(preset or {}), **{k: v for k, v in lig.items() if v not in (None, '')}}
            if not lig.get('motif_seq') or not lig.get('motif_struct'):
                raise ValueError('a ligand needs an aptamer motif (sequence and structure)')
            lig['motif_seq'] = str(lig['motif_seq']).upper().replace('T', 'U').replace(' ', '')
            lig['motif_struct'] = str(lig['motif_struct']).replace(' ', '')
            if len(lig['motif_seq']) != len(lig['motif_struct']) or lig['motif_seq'].count('&') != lig['motif_struct'].count('&') or lig['motif_seq'].count('&') > 1:
                raise ValueError('the aptamer motif sequence and structure must have the same length (one "&" for an interior loop)')
            if set(lig['motif_seq']) - set('ACGU&') or set(lig['motif_struct']) - set('().&'):
                raise ValueError('the aptamer motif may only contain A, C, G, U, "&" and the structure "(", ")", ".", "&"')
            lig['conc_uM'] = max(0.0, float(lig.get('conc_uM', 0.0)))
            if lig.get('dg') is None and not lig.get('kd_uM'):
                raise ValueError('give the ligand Kd (µM) or a binding free energy (kcal/mol)')
        probe = d.get('probing') or None
        if probe:
            vals = probe.get('values') or []
            if n is not None and len(vals) != n:
                raise ValueError(f'probing data has {len(vals)} values for {n} nucleotides')
            kind = str(probe.get('kind') or 'SHAPE').upper()
            if kind not in ('SHAPE', 'DMS'):
                raise ValueError('probing data must be SHAPE or DMS')
            probe = {'kind': kind, 'values': [None if v is None or not math.isfinite(float(v)) or float(v) < -500 else float(v) for v in vals],
                     'm': float(probe.get('m', 1.8)), 'b': float(probe.get('b', -0.6))}
        return cls(temperature=t, sodium=na, ligand=lig, probing=probe)

    def to_json(self):
        return {'temperature': self.temperature, 'sodium': self.sodium, 'ligand': self.ligand,
                'probing': {k: v for k, v in self.probing.items() if k != 'values'} if self.probing else None}

    @property
    def kT(self):
        return R_KCAL * (self.temperature + 273.15)

    def model(self, uniq_ml=False):
        import RNA
        md = RNA.md()
        md.temperature = float(self.temperature)
        md.salt = float(self.sodium)
        if uniq_ml:
            md.uniq_ML = 1
        return md

    def binding_energy(self):
        """Free energy (kcal/mol) that a binding-competent structure gains from the ligand at this concentration."""
        lig = self.ligand
        if not lig:
            return 0.0
        if lig.get('dg') is not None and not lig.get('kd_uM'):
            return float(lig['dg']) if lig.get('conc_uM', 0) > 0 else 0.0
        c, kd = float(lig.get('conc_uM', 0.0)), float(lig['kd_uM'])
        return -self.kT * math.log1p(c / kd) if c > 0 else 0.0

    def fold_compound(self, seq, eval_only=False, uniq_ml=False, probing=False):
        """ViennaRNA fold compound for seq under these conditions. probing=True adds the reactivity pseudo-energies
        (populations only)."""
        import RNA
        md = self.model(uniq_ml=uniq_ml)
        fc = RNA.fold_compound(seq, md, RNA.OPTION_EVAL_ONLY) if eval_only else RNA.fold_compound(seq, md)
        dg = self.binding_energy()
        if self.ligand and dg < 0:
            fc.sc_add_hi_motif(self.ligand['motif_seq'], self.ligand['motif_struct'], float(dg))
        if probing and self.probing:
            vals = self.probing['values']
            if len(vals) == len(seq):
                data = [-999.0] + [(-999.0 if v is None else max(0.0, v)) for v in vals]
                if self.probing['kind'] == 'DMS':          # DMS reports on A and C only
                    data = [-999.0] + [(-999.0 if v is None or seq[i] not in 'AC' else max(0.0, v)) for i, v in enumerate(vals)]
                fc.sc_add_SHAPE_deigan(data, self.probing['m'], self.probing['b'])
        return fc

    # ViennaRNA's eval_structure() leaves the motif soft constraint out (fold and partition function include it), so
    # structure energies add the binding energy themselves when the structure forms the motif.
    def motif_sites(self, seq):
        """Places where the aptamer motif fits seq: list of (positions, motif pairs as position pairs, unpaired)."""
        lig = self.ligand
        if not lig:
            return []
        segs, sts = lig['motif_seq'].split('&'), lig['motif_struct'].split('&')
        flat = ''.join(sts)
        stack, mpairs = [], []
        for i, ch in enumerate(flat):
            if ch == '(':
                stack.append(i)
            elif ch == ')' and stack:
                mpairs.append((stack.pop(), i))
        unp = [i for i, ch in enumerate(flat) if ch == '.']
        n, out = len(seq), []
        starts0 = [p for p in range(n - len(segs[0]) + 1) if seq[p:p + len(segs[0])] == segs[0]]
        for p in starts0:
            if len(segs) == 1:
                pos_lists = [list(range(p, p + len(segs[0])))]
            else:
                pos_lists = [list(range(p, p + len(segs[0]))) + list(range(q, q + len(segs[1])))
                             for q in range(p + len(segs[0]) + 3, n - len(segs[1]) + 1) if seq[q:q + len(segs[1])] == segs[1]]
            for pos in pos_lists:
                out.append(([(pos[a], pos[b]) for a, b in mpairs], [pos[u] for u in unp]))
        return out

    def motif_bonus(self, dbn, sites):
        """Binding energy if dbn forms the motif at one of the sites, else 0."""
        if not sites:
            return 0.0
        pt = [-1] * len(dbn)
        stack = []
        for i, ch in enumerate(dbn):
            if ch == '(':
                stack.append(i)
            elif ch == ')' and stack:
                j = stack.pop()
                pt[i], pt[j] = j, i
        for pairs, unp in sites:
            if all(pt[a] == b for a, b in pairs) and all(pt[u] < 0 for u in unp):
                return self.binding_energy()
        return 0.0

    def energy(self, fc, seq, dbn, sites=None):
        """Free energy of dbn under these conditions (Turner + salt + ligand binding [+ probing if fc has it])."""
        e = fc.eval_structure(dbn)
        if self.ligand:
            e += self.motif_bonus(dbn, self.motif_sites(seq) if sites is None else sites)
        return e

    def describe(self):
        parts = [f'{self.temperature:g} °C', f'{self.sodium:g} M Na+' + (' (standard)' if abs(self.sodium - STANDARD_SALT) < 1e-6 else '')]
        if self.ligand:
            lig = self.ligand
            nm = lig.get('name') or 'ligand'
            if lig.get('kd_uM'):
                parts.append(f"{nm} {lig.get('conc_uM', 0):g} µM (Kd {lig['kd_uM']:g} µM, ΔG bind {self.binding_energy():.2f} kcal/mol)")
            else:
                parts.append(f"{nm} (ΔG bind {self.binding_energy():.2f} kcal/mol)")
        if self.probing:
            parts.append(f"{self.probing['kind']} pseudo-energies (Deigan 2009, m {self.probing['m']:g}, b {self.probing['b']:g})")
        return ', '.join(parts)

    def key(self):
        """Hashable summary for caches."""
        lig = self.ligand or {}
        pr = self.probing or {}
        return (round(self.temperature, 3), round(self.sodium, 5), lig.get('motif_seq'), lig.get('motif_struct'),
                lig.get('kd_uM'), lig.get('dg'), lig.get('conc_uM'), pr.get('kind'), tuple(pr.get('values') or ()), pr.get('m'), pr.get('b'))
