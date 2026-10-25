// Optional confidence model (rnafold3d/qa.py): predicted local accuracy (lDDT) of every nucleotide of every
// conformation, computed by the local engine in the background when an ensemble opens. Adds the metric
// "Confidence (predicted lDDT)" (mean over the nucleotides; higher = more reliable) and the colour scheme
// "Confidence (per nucleotide)". Switch it off in View ▸ Preferences.
import { registerScheme } from './colors.js';

export const CONF_BANDS = [
  { min: 0.8, color: '#1c5cab', label: 'high (≥ 0.8)' },
  { min: 0.65, color: '#6da7ec', label: 'medium (0.65–0.8)' },
  { min: 0.5, color: '#eda100', label: 'low (0.5–0.65)' },
  { min: -1, color: '#eb6834', label: 'very low (< 0.5)' },
];
export const confColor = v => (v === null || v === undefined ? '#c4c2ba' : CONF_BANDS.find(b => v >= b.min).color);

export class Confidence {
  constructor(app) {
    this.app = app;
    registerScheme({
      id: 'confidence', label: 'Confidence (per nucleotide)',
      residueColors: (ens, ctx) => (ens._confidence ? ens._confidence.conformations[ctx?.conf ?? 0].per_residue.map(confColor) : ens.residues.map(() => '#c4c2ba')),
      legend: ens => ({ title: ens._confidence ? 'Predicted local accuracy (lDDT, confidence model)' : 'Confidence model not computed', items: CONF_BANDS.map(b => ({ color: b.color, label: b.label })) }),
    });
    app.on('ensemble', ens => { if (ens) this.ensure(ens); });
  }

  enabled() { return this.app.prefs.confidence !== false; }

  async ensure(ens) {
    const app = this.app;
    if (!this.enabled() || !app.api.online || ens._confidence || ens._confidenceLoading) return;
    ens._confidenceLoading = true;
    try {
      const body = ens.serverId ? { document: ens.serverId } : { doc: ens.doc };
      const r = await app.api.post('/api/confidence', body);
      if (!r.available) return;
      ens._confidence = r;
      for (let k = 0; k < ens.n; k++) ens.conformations[k].values.confidence = r.conformations[k].global;
      const m = { key: 'confidence', label: 'Confidence (predicted lDDT)', unit: '', group: 'confidence', lower_is_better: false,
        description: `Mean predicted local accuracy of the nucleotides (lDDT, 0–1) from the RCV confidence model${r.info?.cv ? ` (cross-validated per-nucleotide r = ${r.info.cv.per_nucleotide_r})` : ''}. A guide to reliability, not a measurement.` };
      if (ens.metricByKey.confidence) Object.assign(ens.metricByKey.confidence, m); else { ens.metrics.push(m); ens.metricByKey.confidence = m; }
      app.refreshSchemes();
      if (app.ens === ens) { app.emit('values'); if (app.state.colorScheme === 'confidence') app.refreshColors(); }
    } catch { /* optional feature: stay silent */ } finally { ens._confidenceLoading = false; }
  }
}
