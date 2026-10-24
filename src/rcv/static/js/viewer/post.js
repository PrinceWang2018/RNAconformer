// Post-processing for publication-style rendering: screen-space ambient occlusion from depth,
// depth-based silhouette outlines, optional depth cueing, sRGB output and FXAA.
import * as THREE from '../../vendor/three/three.module.min.js';

const VERT = /* glsl */`
varying vec2 vUv;
void main() { vUv = uv; gl_Position = vec4(position.xy, 0.0, 1.0); }`;

const COMMON = /* glsl */`
uniform sampler2D tDepth;
uniform mat4 invProjection;
uniform vec2 resolution;
vec3 viewPos(vec2 uv) {
  float d = texture2DLodEXT(tDepth, uv, 0.0).x;
  vec4 c = vec4(uv * 2.0 - 1.0, d * 2.0 - 1.0, 1.0);
  vec4 v = invProjection * c;
  return v.xyz / v.w;
}`;

const AO_FRAG = /* glsl */`
${COMMON}
uniform mat4 projection;
uniform float radius;
uniform float intensity;
uniform float bias;
uniform float isOrtho;
varying vec2 vUv;
#define SAMPLES 16
void main() {
  float d = texture2D(tDepth, vUv).x;
  if (d >= 0.99999) { gl_FragColor = vec4(1.0); return; }
  vec2 px = 1.0 / resolution;
  vec3 p = viewPos(vUv);
  vec3 pl = viewPos(vUv - vec2(px.x, 0.0)), pr = viewPos(vUv + vec2(px.x, 0.0));
  vec3 pd = viewPos(vUv - vec2(0.0, px.y)), pu = viewPos(vUv + vec2(0.0, px.y));
  vec3 dx = abs(pr.z - p.z) < abs(p.z - pl.z) ? pr - p : p - pl;
  vec3 dy = abs(pu.z - p.z) < abs(p.z - pd.z) ? pu - p : p - pd;
  vec3 n = normalize(cross(dx, dy));
  float scale = projection[1][1] * 0.5 * resolution.y;
  float rpx = isOrtho > 0.5 ? radius * scale : radius * scale / max(1e-3, -p.z);
  rpx = clamp(rpx, 3.0, 160.0);
  float noise = fract(52.9829189 * fract(dot(gl_FragCoord.xy, vec2(0.06711056, 0.00583715))));
  float sum = 0.0;
  float r2 = radius * radius;
  for (int i = 0; i < SAMPLES; i++) {
    float fi = (float(i) + 0.5) / float(SAMPLES);
    float a = float(i) * 2.39996323 + noise * 6.2831853;
    vec2 off = vec2(cos(a), sin(a)) * sqrt(fi) * rpx * px;
    vec3 q = viewPos(vUv + off);
    vec3 v = q - p;
    float vv = dot(v, v);
    float fall = max(r2 - vv, 0.0) / r2;
    sum += fall * max(dot(v, n) - bias, 0.0) / (vv + 0.05);
  }
  float ao = clamp(1.0 - intensity * 2.0 * sum / float(SAMPLES), 0.0, 1.0);
  gl_FragColor = vec4(ao, ao, ao, 1.0);
}`;

const BLUR_FRAG = /* glsl */`
${COMMON}
uniform sampler2D tAO;
varying vec2 vUv;
void main() {
  vec2 px = 1.0 / resolution;
  float z0 = viewPos(vUv).z;
  float sum = 0.0, wsum = 0.0;
  for (int x = -2; x <= 1; x++) for (int y = -2; y <= 1; y++) {
    vec2 uv = vUv + (vec2(float(x), float(y)) + 0.5) * px;
    float z = viewPos(uv).z;
    float w = 1.0 / (0.02 + abs(z - z0));
    sum += texture2DLodEXT(tAO, uv, 0.0).r * w;
    wsum += w;
  }
  float ao = sum / max(wsum, 1e-5);
  gl_FragColor = vec4(ao, ao, ao, 1.0);
}`;

const COMPOSITE_FRAG = /* glsl */`
${COMMON}
uniform sampler2D tColor;
uniform sampler2D tAO;
uniform float aoStrength;
uniform float outlineWidth;
uniform float outlineStrength;
uniform vec3 outlineColor;
uniform float outlineThreshold;
uniform vec3 bgColor;
uniform float bgAlpha;
uniform float fogStrength;
uniform float fogNear;
uniform float fogFar;
varying vec2 vUv;
float lin2srgb(float c) { return c <= 0.0031308 ? c * 12.92 : 1.055 * pow(c, 1.0 / 2.4) - 0.055; }
float vz(vec2 uv) { float d = texture2DLodEXT(tDepth, uv, 0.0).x; return d >= 0.99999 ? -1e6 : viewPos(uv).z; }
void main() {
  vec4 src = texture2D(tColor, vUv);
  vec3 col = src.rgb;
  float alpha = src.a;
  float d = texture2D(tDepth, vUv).x;
  bool bg = d >= 0.99999;
  float z0 = vz(vUv);
  if (!bg) {
    float ao = texture2D(tAO, vUv).r;
    col *= mix(1.0, ao, aoStrength);
    if (fogStrength > 0.0) {
      float f = clamp((-z0 - fogNear) / max(1e-3, fogFar - fogNear), 0.0, 1.0);
      col = mix(col, bgColor, f * f * fogStrength);
    }
  }
  // background pixels keep the scene colour: it already holds the clear colour plus any
  // translucent overlays (ensemble ghosts) that are drawn without writing depth
  if (outlineStrength > 0.0) {
    vec2 px = outlineWidth / resolution;
    float edge = 0.0;
    for (int i = 0; i < 8; i++) {
      float a = float(i) * 0.785398;
      vec2 o = vec2(cos(a), sin(a)) * px;
      float z = vz(vUv + o);
      float thr = outlineThreshold + 0.012 * abs(z);
      edge += step(thr, z - z0);
    }
    float e = clamp(edge / 2.5, 0.0, 1.0) * outlineStrength;
    col = mix(col, outlineColor, e);
    alpha = max(alpha, e);
  }
  gl_FragColor = vec4(lin2srgb(col.r), lin2srgb(col.g), lin2srgb(col.b), alpha);
}`;

const FXAA_FRAG = /* glsl */`
uniform sampler2D tInput;
uniform vec2 resolution;
varying vec2 vUv;
void main() {
  vec2 inv = 1.0 / resolution;
  vec4 cM = texture2D(tInput, vUv);
  vec4 cNW = texture2D(tInput, vUv + vec2(-1.0, -1.0) * inv);
  vec4 cNE = texture2D(tInput, vUv + vec2(1.0, -1.0) * inv);
  vec4 cSW = texture2D(tInput, vUv + vec2(-1.0, 1.0) * inv);
  vec4 cSE = texture2D(tInput, vUv + vec2(1.0, 1.0) * inv);
  vec3 L = vec3(0.299, 0.587, 0.114);
  float lNW = dot(cNW.rgb, L) + (1.0 - cNW.a), lNE = dot(cNE.rgb, L) + (1.0 - cNE.a);
  float lSW = dot(cSW.rgb, L) + (1.0 - cSW.a), lSE = dot(cSE.rgb, L) + (1.0 - cSE.a);
  float lM = dot(cM.rgb, L) + (1.0 - cM.a);
  float lMin = min(lM, min(min(lNW, lNE), min(lSW, lSE)));
  float lMax = max(lM, max(max(lNW, lNE), max(lSW, lSE)));
  vec2 dir = vec2(-((lNW + lNE) - (lSW + lSE)), ((lNW + lSW) - (lNE + lSE)));
  float reduce = max((lNW + lNE + lSW + lSE) * 0.03125, 1.0 / 128.0);
  float rcp = 1.0 / (min(abs(dir.x), abs(dir.y)) + reduce);
  dir = clamp(dir * rcp, vec2(-8.0), vec2(8.0)) * inv;
  vec4 a = 0.5 * (texture2D(tInput, vUv + dir * (1.0 / 3.0 - 0.5)) + texture2D(tInput, vUv + dir * (2.0 / 3.0 - 0.5)));
  vec4 b = a * 0.5 + 0.25 * (texture2D(tInput, vUv - dir * 0.5) + texture2D(tInput, vUv + dir * 0.5));
  float lB = dot(b.rgb, L) + (1.0 - b.a);
  gl_FragColor = (lB < lMin || lB > lMax) ? a : b;
}`;

const COPY_FRAG = /* glsl */`
uniform sampler2D tInput;
varying vec2 vUv;
void main() { gl_FragColor = texture2D(tInput, vUv); }`;

function rt(w, h, opts = {}) {
  return new THREE.WebGLRenderTarget(w, h, {
    minFilter: THREE.LinearFilter, magFilter: THREE.LinearFilter, depthBuffer: !!opts.depth,
    type: opts.type || THREE.UnsignedByteType, format: THREE.RGBAFormat, colorSpace: THREE.NoColorSpace,
    depthTexture: opts.depth ? new THREE.DepthTexture(w, h, THREE.UnsignedIntType) : undefined,
  });
}

export class PostFX {
  constructor(renderer) {
    this.renderer = renderer;
    this.quadCam = new THREE.OrthographicCamera(-1, 1, 1, -1, 0, 1);
    this.quad = new THREE.Mesh(new THREE.PlaneGeometry(2, 2));
    this.quad.frustumCulled = false;
    this.quadScene = new THREE.Scene();
    this.quadScene.add(this.quad);
    const mk = frag => new THREE.ShaderMaterial({
      vertexShader: VERT, fragmentShader: frag, depthTest: false, depthWrite: false,
      uniforms: {
        tDepth: { value: null }, tColor: { value: null }, tAO: { value: null }, tInput: { value: null },
        invProjection: { value: new THREE.Matrix4() }, projection: { value: new THREE.Matrix4() },
        resolution: { value: new THREE.Vector2(1, 1) }, radius: { value: 5.0 }, intensity: { value: 1.0 },
        bias: { value: 0.15 }, isOrtho: { value: 0 }, aoStrength: { value: 0.8 }, outlineWidth: { value: 1.5 },
        outlineStrength: { value: 1.0 }, outlineColor: { value: new THREE.Color(0x111111) }, outlineThreshold: { value: 1.2 },
        bgColor: { value: new THREE.Color(0xffffff) }, bgAlpha: { value: 1 }, fogStrength: { value: 0 },
        fogNear: { value: 0 }, fogFar: { value: 1 },
      },
    });
    this.matAO = mk(AO_FRAG);
    this.matBlur = mk(BLUR_FRAG);
    this.matComposite = mk(COMPOSITE_FRAG);
    this.matFXAA = mk(FXAA_FRAG);
    this.matCopy = mk(COPY_FRAG);
    this.size = [0, 0];
    this.targets = null;
  }

  setSize(w, h) {
    w = Math.max(1, Math.round(w)); h = Math.max(1, Math.round(h));
    if (this.size[0] === w && this.size[1] === h && this.targets) return;
    this.dispose();
    this.size = [w, h];
    this.targets = {
      scene: rt(w, h, { depth: true, type: THREE.HalfFloatType }),
      ao: rt(w, h), blur: rt(w, h), comp: rt(w, h),
    };
  }

  dispose() {
    if (!this.targets) return;
    for (const t of Object.values(this.targets)) { t.depthTexture?.dispose(); t.dispose(); }
    this.targets = null;
  }

  _pass(mat, target) {
    this.quad.material = mat;
    this.renderer.setRenderTarget(target);
    this.renderer.render(this.quadScene, this.quadCam);
  }

  // Render `scene` with `camera` into `output` (null = canvas) using settings `s`.
  render(scene, camera, output, s) {
    const r = this.renderer;
    const [w, h] = this.size;
    const T = this.targets;
    const bg = new THREE.Color(s.background);          // linear working colour
    r.setRenderTarget(T.scene);
    r.setClearColor(bg, s.transparent ? 0 : 1);
    r.clear(true, true, true);
    r.render(scene, camera);

    const inv = camera.projectionMatrixInverse;
    for (const m of [this.matAO, this.matBlur, this.matComposite]) {
      m.uniforms.tDepth.value = T.scene.depthTexture;
      m.uniforms.invProjection.value.copy(inv);
      m.uniforms.projection.value.copy(camera.projectionMatrix);
      m.uniforms.resolution.value.set(w, h);
    }
    const useAO = s.ao > 0;
    if (useAO) {
      this.matAO.uniforms.radius.value = s.aoRadius;
      this.matAO.uniforms.intensity.value = s.aoIntensity;
      this.matAO.uniforms.isOrtho.value = camera.isOrthographicCamera ? 1 : 0;
      this._pass(this.matAO, T.ao);
      this.matBlur.uniforms.tAO.value = T.ao.texture;
      this._pass(this.matBlur, T.blur);
    }
    const c = this.matComposite.uniforms;
    c.tColor.value = T.scene.texture;
    c.tAO.value = T.blur.texture;
    c.aoStrength.value = useAO ? s.ao : 0;
    c.outlineWidth.value = s.outlineWidth * (s.pixelScale || 1);
    c.outlineStrength.value = s.outline;
    c.outlineColor.value.set(s.outlineColor || '#111111');
    c.outlineThreshold.value = s.outlineThreshold ?? 1.2;
    c.bgColor.value.copy(bg);
    c.bgAlpha.value = s.transparent ? 0 : 1;
    c.fogStrength.value = s.fog || 0;
    c.fogNear.value = s.fogNear || 0;
    c.fogFar.value = s.fogFar || 1;
    if (s.fxaa) {
      this._pass(this.matComposite, T.comp);
      this.matFXAA.uniforms.tInput.value = T.comp.texture;
      this.matFXAA.uniforms.resolution.value.set(w, h);
      this._pass(this.matFXAA, output);
    } else {
      this._pass(this.matComposite, output);
    }
    r.setRenderTarget(null);
  }
}
