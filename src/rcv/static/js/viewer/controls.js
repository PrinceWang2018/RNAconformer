// Molecular-viewer mouse/touch controls: the model rotates freely (trackball) under a fixed camera.
//   left drag: rotate · Ctrl/⌘+left drag: roll · right / middle / Shift+left drag: pan · wheel / pinch: zoom
import * as THREE from '../../vendor/three/three.module.min.js';

export class TrackballControls {
  constructor(dom, { onChange, onStart, onEnd } = {}) {
    this.dom = dom;
    this.quaternion = new THREE.Quaternion();
    this.pan = new THREE.Vector2(0, 0);         // view-space offset (Å)
    this.zoom = 1;                              // >1 = closer
    this.onChange = onChange || (() => {});
    this.onStart = onStart || (() => {});
    this.onEnd = onEnd || (() => {});
    this.rotateSpeed = 0.006;
    this.velocity = new THREE.Vector3();        // axis * angle per frame for inertia
    this.inertia = true;
    this.spin = 0;                              // auto-rotate rad/s around view Y
    this.enabled = true;
    this.pointers = new Map();
    this.mode = null;
    this.unitsPerPixel = 0.05;
    this._bind();
  }

  _bind() {
    const d = this.dom;
    d.addEventListener('contextmenu', e => e.preventDefault());
    d.addEventListener('pointerdown', e => this._down(e));
    d.addEventListener('pointermove', e => this._move(e));
    d.addEventListener('pointerup', e => this._up(e));
    d.addEventListener('pointercancel', e => this._up(e));
    d.addEventListener('wheel', e => {
      if (!this.enabled) return;
      e.preventDefault();
      const s = Math.exp(-e.deltaY * (e.deltaMode === 1 ? 0.05 : 0.0015));
      this.zoom = THREE.MathUtils.clamp(this.zoom * s, 0.05, 40);
      this.onChange('zoom');
    }, { passive: false });
  }

  _down(e) {
    if (!this.enabled) return;
    if (this.interceptor && this.pointers.size === 0 && this.interceptor(e)) return;   // e.g. pull tool took the gesture
    this.dom.setPointerCapture(e.pointerId);
    this.pointers.set(e.pointerId, { x: e.clientX, y: e.clientY });
    this.velocity.set(0, 0, 0);
    this.moved = 0;
    if (this.pointers.size === 2) { this.mode = 'pinch'; this._pinch = this._pinchState(); return; }
    if (e.button === 2 || e.button === 1 || (e.button === 0 && e.shiftKey)) this.mode = 'pan';
    else if (e.button === 0 && (e.ctrlKey || e.metaKey)) this.mode = 'roll';
    else this.mode = 'rotate';
    this.last = { x: e.clientX, y: e.clientY, t: performance.now() };
    this.onStart(this.mode);
  }

  _pinchState() {
    const [a, b] = [...this.pointers.values()];
    return { d: Math.hypot(a.x - b.x, a.y - b.y), cx: (a.x + b.x) / 2, cy: (a.y + b.y) / 2 };
  }

  _move(e) {
    if (!this.pointers.has(e.pointerId) || !this.enabled) return;
    this.pointers.set(e.pointerId, { x: e.clientX, y: e.clientY });
    if (this.mode === 'pinch' && this.pointers.size === 2) {
      const s = this._pinchState();
      this.zoom = THREE.MathUtils.clamp(this.zoom * (s.d / Math.max(1, this._pinch.d)), 0.05, 40);
      this.pan.x += (s.cx - this._pinch.cx) * this.unitsPerPixel;
      this.pan.y -= (s.cy - this._pinch.cy) * this.unitsPerPixel;
      this._pinch = s;
      this.onChange('pinch');
      return;
    }
    if (!this.last) return;
    const dx = e.clientX - this.last.x, dy = e.clientY - this.last.y;
    this.moved += Math.abs(dx) + Math.abs(dy);
    const now = performance.now();
    if (this.mode === 'rotate') {
      const angle = Math.hypot(dx, dy) * this.rotateSpeed;
      if (angle > 0) {
        const axis = new THREE.Vector3(dy, dx, 0).normalize();
        this._rotate(axis, angle);
        const dt = Math.max(1, now - this.last.t);
        this.velocity.copy(axis).multiplyScalar(angle * Math.min(1, 16 / dt));
      }
    } else if (this.mode === 'roll') {
      const rect = this.dom.getBoundingClientRect();
      const cx = rect.left + rect.width / 2, cy = rect.top + rect.height / 2;
      const a0 = Math.atan2(this.last.y - cy, this.last.x - cx), a1 = Math.atan2(e.clientY - cy, e.clientX - cx);
      this._rotate(new THREE.Vector3(0, 0, 1), -(a1 - a0));
    } else if (this.mode === 'pan') {
      this.pan.x += dx * this.unitsPerPixel;
      this.pan.y -= dy * this.unitsPerPixel;
    }
    this.last = { x: e.clientX, y: e.clientY, t: now };
    this.onChange(this.mode);
  }

  _up(e) {
    if (!this.pointers.has(e.pointerId)) return;
    this.pointers.delete(e.pointerId);
    try { this.dom.releasePointerCapture(e.pointerId); } catch { /* already released */ }
    if (this.pointers.size === 0) {
      const wasRotate = this.mode === 'rotate';
      if (!wasRotate || performance.now() - (this.last?.t || 0) > 80) this.velocity.set(0, 0, 0);
      this.mode = null;
      this.last = null;
      this.onEnd(wasRotate);
    } else if (this.pointers.size === 1) {
      const p = [...this.pointers.values()][0];
      this.mode = 'rotate';
      this.last = { x: p.x, y: p.y, t: performance.now() };
    }
  }

  _rotate(axis, angle) {
    const q = new THREE.Quaternion().setFromAxisAngle(axis, angle);
    this.quaternion.premultiply(q).normalize();
  }

  rotateBy(axis, angle) { this._rotate(axis.clone().normalize(), angle); this.onChange('api'); }

  // Advance inertia / auto-spin; returns true while still animating.
  update(dt) {
    let active = false;
    if (this.inertia && this.mode === null && this.velocity.lengthSq() > 1e-10) {
      const ang = this.velocity.length();
      this._rotate(this.velocity.clone().normalize(), ang * Math.min(3, dt * 60));
      this.velocity.multiplyScalar(Math.pow(0.92, dt * 60));
      if (this.velocity.lengthSq() < 1e-9) this.velocity.set(0, 0, 0);
      active = true;
    }
    if (this.spin && this.mode === null) {
      this._rotate(new THREE.Vector3(0, 1, 0), this.spin * dt);
      active = true;
    }
    return active;
  }

  reset() { this.quaternion.identity(); this.pan.set(0, 0); this.zoom = 1; this.velocity.set(0, 0, 0); }
  get state() { return { q: this.quaternion.toArray(), pan: this.pan.toArray(), zoom: this.zoom }; }
  set state(s) {
    if (!s) return;
    this.quaternion.fromArray(s.q || [0, 0, 0, 1]);
    this.pan.fromArray(s.pan || [0, 0]);
    this.zoom = s.zoom || 1;
  }
}
