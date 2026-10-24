// Scenes: named views (camera, conformation, colours, representation, selection) that can be recalled with a
// smooth camera move — for presentations, figure series and reproducible screenshots. Saved in .rcv files.
import * as THREE from '../vendor/three/three.module.min.js';
import { easeInOut } from './util.js';

export class Scenes {
  constructor(app) {
    this.app = app;
    app.registerSessionPart('scenes', {
      save: ens => (ens?.scenes || []).map(s => ({ ...s })),
      load: (ens, list) => { ens.scenes = (list || []).map(s => ({ ...s })); app.emit('scenes'); },
    });
  }

  list(ens = this.app.ens) { return ens ? (ens.scenes ||= []) : []; }

  save(name) {
    const app = this.app, ens = app.ens;
    if (!ens) throw new Error('no ensemble open');
    const list = this.list();
    name = String(name || `scene${list.length + 1}`).trim();
    const s = {
      name, camera: app.viewer.controls.state, conf: app.state.conf, confId: ens.conformations[app.state.conf].id,
      colorScheme: app.state.colorScheme, representation: app.viewer.settings.representation, style: app.viewer.settings.style,
      background: app.viewer.settings.background, selection: [...app.state.selection], ghosts: app.state.ghosts, putty: app.state.putty,
    };
    const i = list.findIndex(x => x.name === name);
    if (i >= 0) list[i] = s; else list.push(s);
    app.emit('scenes');
    return s;
  }

  remove(name) {
    const list = this.list();
    const i = list.findIndex(x => x.name === name);
    if (i < 0) return false;
    list.splice(i, 1);
    this.app.emit('scenes');
    return true;
  }

  recall(name, { duration = 0.9 } = {}) {
    const app = this.app, ens = app.ens;
    const s = this.list().find(x => x.name === name);
    if (!s) throw new Error(`no scene "${name}" (scenes: ${this.list().map(x => x.name).join(', ') || 'none'})`);
    const k = ens.conformations.findIndex(c => c.id === s.confId);
    const conf = k >= 0 ? k : Math.min(s.conf ?? 0, ens.n - 1);
    if (s.colorScheme && s.colorScheme !== app.state.colorScheme) app.setColorScheme(s.colorScheme);
    if (s.representation && s.representation !== app.viewer.settings.representation) app.setRepresentation(s.representation);
    if (s.style && s.style !== app.viewer.settings.style) app.setStyle(s.style);
    if (s.background && s.background !== app.viewer.settings.background) app.setBackground(s.background);
    if (!!s.ghosts !== app.state.ghosts) app.toggleGhosts();
    if (!!s.putty !== app.state.putty) app.togglePutty();
    app.select(s.selection || [], 'set');
    if (conf !== app.state.conf) app.goTo(conf, { duration });
    this.flyTo(s.camera, duration);
    return s;
  }

  // Animate the camera (rotation slerp, pan and zoom interpolation).
  flyTo(cam, duration = 0.9) {
    const c = this.app.viewer.controls;
    if (!cam) return;
    const q0 = c.quaternion.clone(), q1 = new THREE.Quaternion().fromArray(cam.q);
    const p0 = c.pan.clone(), p1 = new THREE.Vector2().fromArray(cam.pan);
    const z0 = c.zoom, z1 = cam.zoom;
    if (duration <= 0) { c.state = cam; this.app.viewer.needsRender = true; return; }
    const t0 = performance.now();
    c.velocity?.set(0, 0, 0);
    this.app.viewer._addAnimator(() => {
      const t = Math.min(1, (performance.now() - t0) / (duration * 1000));
      const e = easeInOut(t);
      c.quaternion.copy(q0).slerp(q1, e);
      c.pan.lerpVectors(p0, p1, e);
      c.zoom = z0 + (z1 - z0) * e;
      return t < 1;
    });
  }
}
