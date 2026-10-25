// Methods paragraphs and references for the analyses of RCV 0.7 that were run on the open ensemble (pathway,
// phi-value scan, ensemble map, probing data, confidence model) - built from the stored results, so they quote
// the exact settings. Appended to Export ▸ Methods text & citation.
import { fmt } from './util.js';
import { conditionsOf, conditionsText } from './conditions.js';

export const REFERENCES = {
  vienna: 'Lorenz R, Bernhart SH, Höner zu Siederdissen C, et al. ViennaRNA Package 2.0. Algorithms Mol Biol. 2011;6:26.',
  turner: 'Turner DH, Mathews DH. NNDB: the nearest neighbor parameter database for predicting stability of nucleic acid secondary structure. Nucleic Acids Res. 2010;38:D280–D282.',
  findpath: 'Flamm C, Hofacker IL, Maurer-Stroh S, Stadler PF, Zehl M. Design of multistable RNA molecules. RNA. 2001;7:254–265.',
  tpt: 'Metzner P, Schütte C, Vanden-Eijnden E. Transition path theory for Markov jump processes. Multiscale Model Simul. 2009;7:1192–1219.',
  phi: 'Matouschek A, Kellis JT, Serrano L, Fersht AR. Mapping the transition state and pathway of protein folding by protein engineering. Nature. 1989;340:122–126.',
  shapes: 'Giegerich R, Voß B, Rehmsmeier M. Abstract shapes of RNA. Nucleic Acids Res. 2004;32:4843–4851.',
  trust: 'Venna J, Kaski S. Neighborhood preservation in nonlinear projection methods: an experimental study. ICANN 2001, LNCS 2130:485–491.',
  deigan: 'Deigan KE, Li TW, Mathews DH, Weeks KM. Accurate SHAPE-directed RNA structure determination. Proc Natl Acad Sci USA. 2009;106:97–102.',
  bme: 'Bottaro S, Bengtsen T, Lindorff-Larsen K. Integrating molecular simulation and experimental data: a Bayesian/maximum entropy reweighting approach. Methods Mol Biol. 2020;2112:219–240.',
  salt: 'Lorenz R, Hofacker IL, Stadler PF. RNA folding with hard and soft constraints. Algorithms Mol Biol. 2016;11:8.',
  lddt: 'Mariani V, Biasini M, Barbato A, Schwede T. lDDT: a local superposition-free score for comparing protein structures and models using distance difference tests. Bioinformatics. 2013;29:2722–2728.',
};

export function analysisMethods(app, ens) {
  const paras = [], refs = new Set();
  const r = ens._pathway, sc = ens._phiScan;
  const cond = conditionsOf(app, ens);
  if (r) {
    refs.add('vienna').add('turner').add('findpath').add('tpt');
    if (cond.sodium !== 1.021 || cond.ligand) refs.add('salt');
    const minima = r.nodes.filter(n => n.kind === 'minimum').length, states = r.nodes.filter(n => n.kind === 'ensemble').length;
    let t = `Refolding from ${ens.shortId(r.a)} (A) to ${ens.shortId(r.b)} (B) was analysed at the secondary-structure level with ViennaRNA ${r.method?.vienna || ''} (Turner 2004 parameters; ${r.conditions_text}). `;
    if (r.A.dbn === r.B.dbn) t += 'A and B share their nested canonical base pairs, so no secondary-structure barrier separates them. ';
    else if (r.same_basin) t += 'A and B relax into the same local minimum of the secondary-structure landscape, so no barrier separates them. ';
    else {
      t += `Kinetics were computed between the basins of A and B (the local minima they reach by steepest descent; A lies ${fmt(r.A.above_basin ?? 0, 1)} and B ${fmt(r.B.above_basin ?? 0, 1)} kcal/mol above its basin minimum). `;
      t += `Direct refolding paths were computed with the findpath heuristic (${r.method?.path?.replace('findpath (Flamm et al. 2001), ', '') || ''}); routes through ${states} further 2D state${states === 1 ? '' : 's'} of the ensemble and ${minima} local minim${minima === 1 ? 'um' : 'a'} (steepest descent from ${r.method?.samples ?? 0} Boltzmann-sampled structures, seed ${r.method?.seed ?? '–'}) formed a network whose edges carry the findpath saddle energies. Transition path theory on this network, with Arrhenius rates k = k₀·exp(−(S − E)/RT) and k₀ = ${r.k0.toExponential(0)} s⁻¹, gave committors, reactive fluxes and flux-ranked routes; the switching time is the mean first-passage time (MFPT) from the basin of A to that of B on the same network, and the effective barrier is RT·ln(k₀·MFPT). `;
      t += `The barrier of the direct path was ${fmt(r.direct.barrier, 1)} kcal/mol and the effective barrier over all routes ${fmt(r.effective_barrier, 1)} kcal/mol (ΔG(A→B) = ${fmt(r.dG, 1)} kcal/mol). `;
      if (r.trap) t += `A trap was reported because wells below the basin of A held the RNA for ${Math.round((r.trap.time_share ?? 1) * 100)} % of the mean first-passage time (expected time in each node before B is reached, from the fundamental matrix of the network), making the switch more than twice as slow as crossing the best single route. `;
      if (r.bottleneck) t += `The cooperative unit was defined as the base pairs opened between the lowest point before the highest saddle and the saddle (${r.bottleneck.unit.length} bp), and the nucleus as the pairs formed between the saddle and the next minimum (${r.bottleneck.nucleus.length} bp). `;
    }
    if (r.motion) t += `Rigid modules between the two 3D models were identified by merging groups seeded by helices shared by A and B while the joint fit RMSD of their C4′, C1′, P and O3′ atoms stayed below ${r.motion.tau} Å (${r.motion.modules.length} module${r.motion.modules.length === 1 ? '' : 's'}); module motions are given relative to the largest module as rotation angle and screw axis.`;
    paras.push(t.trim());
  }
  if (sc && r) {
    refs.add('phi');
    paras.push(`In-silico φ-values were computed for all ${3 * ens.nRes} single-point mutations as φ = ΔΔG‡/ΔΔG(B−A), with the route kept fixed and its energies re-evaluated for each mutant (findpath search width ${sc.width}); φ is reported where |ΔΔG| ≥ ${sc.phi_min_ddg} kcal/mol.`);
  }
  const L = ens.doc?.ensemble?.landscape || ens._landscape;
  if (L?.quality?.map) {
    refs.add('shapes').add('trust');
    const q = L.quality.map, m = L.quality.mds;
    paras.push(`Ensemble maps place conformations by pairwise C4′ RMSD under a topology constraint: each abstract shape (level 5) of the nested secondary structure occupies its own non-overlapping disk (constrained stress majorisation, oriented by Procrustes to classical MDS). Trustworthiness and continuity of ${q.k}-nearest neighbourhoods were ${fmt(q.trustworthiness, 2)} and ${fmt(q.continuity, 2)}, topology purity ${fmt(q.topology_purity, 2)} (classical MDS: ${fmt(m?.trustworthiness, 2)}, ${fmt(m?.continuity, 2)} and ${fmt(m?.topology_purity, 2)}).`);
  }
  const p = ens.probing;
  if (p?.fit) {
    refs.add('bme');
    const w = p.fit.reweighting;
    paras.push(`${p.kind} reactivities (${p.values.filter(v => v !== null).length} of ${ens.nRes} nucleotides${p.factor ? ', normalised by the 2–8 % rule' : ''}) were compared with each conformation as the area under the ROC curve of reactivity versus unpaired nucleotides (base pairs annotated in 3D)` +
      (w ? `, and the ensemble was reweighted by maximum entropy to fit them (forward model: reactivity linear in the unpaired fraction; θ = ${fmt(w.theta, 2)}, σ = ${w.sigma}; χ²/n ${fmt(w.chi2_prior, 2)} → ${fmt(w.chi2, 2)}; effective number of conformations ${fmt(w.n_eff_conformations, 1)} of ${ens.n}).` : '.'));
    if (cond.useProbing) { refs.add('deigan'); paras.push(`Data-informed state populations used the ${p.kind} reactivities as pseudo-energies (Deigan et al. 2009; m = 1.8, b = −0.6 kcal/mol).`); }
  }
  if (ens._confidence?.info) {
    refs.add('lddt');
    const i = ens._confidence.info;
    paras.push(`Per-nucleotide confidence was predicted with the RCV confidence model (a two-layer neural network on ${i.features?.length ?? 41} structural features per nucleotide, trained on decoys of ${i.training?.targets ?? '–'} RNA structures that share no family with the evaluation benchmark; cross-validated per-nucleotide Pearson r = ${i.cv?.per_nucleotide_r ?? '–'} against lDDT).`);
  }
  if (!r && (cond.sodium !== 1.021 || cond.ligand || cond.temperature !== 37)) paras.push(`Secondary-structure energies were computed at ${conditionsText(cond, ens)}.`);
  return { text: paras.join(' '), references: [...refs].map(k => REFERENCES[k]) };
}
