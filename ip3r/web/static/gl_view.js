// Minimal WebGL 1 viewer for coarse-grained models (one bead per residue): shaded sphere impostors as point sprites,
// backbone bonds as lines, an orbit camera (drag, touch, wheel / pinch, arrow keys, +/-), picking by the nearest
// projected bead, and an optional slow turn about the four-fold axis. No library.
//
// Coordinates are in the deposit's axis frame (z = four-fold axis, cytosol +z), drawn with +z up the screen.
// const v = createViewer(canvas, { onPick: i => ... });
// v.setModel({ positions: Float32Array(3n), colors: Float32Array(3n), sizes: Float32Array(n), bonds: Uint16Array|Uint32Array(2m) })
// v.updatePositions(p); v.setColors(c); v.setSizes(s); v.spin(on); v.dispose();

const VS_P = `attribute vec3 aPos; attribute vec3 aCol; attribute float aSize;
uniform mat4 uMV; uniform mat4 uP; uniform float uScale; varying vec3 vCol; varying float vDepth;
void main() { vec4 p = uMV * vec4(aPos, 1.0); gl_Position = uP * p;
  gl_PointSize = aSize > 0.0 ? max(1.5, aSize * uScale / -p.z) : 0.0; vCol = aCol; vDepth = -p.z; }`;
const FS_P = `precision mediump float; varying vec3 vCol; varying float vDepth; uniform vec2 uFog;
void main() { vec2 c = gl_PointCoord * 2.0 - 1.0; c.y = -c.y; float r2 = dot(c, c); if (r2 > 1.0) discard;
  vec3 n = vec3(c, sqrt(1.0 - r2)); vec3 L = normalize(vec3(0.4, 0.55, 0.75));
  float diff = max(dot(n, L), 0.0); float spec = pow(max(dot(reflect(-L, n), vec3(0.0, 0.0, 1.0)), 0.0), 28.0);
  vec3 col = vCol * (0.3 + 0.75 * diff) + vec3(0.3) * spec;
  float fog = clamp((vDepth - uFog.x) / (uFog.y - uFog.x), 0.0, 0.7); col = mix(col, vec3(0.02, 0.03, 0.06), fog);
  gl_FragColor = vec4(col, 1.0); }`;
const VS_L = `attribute vec3 aPos; attribute vec3 aCol; uniform mat4 uMV; uniform mat4 uP; varying vec3 vCol;
void main() { gl_Position = uP * uMV * vec4(aPos, 1.0); vCol = aCol; }`;
const FS_L = `precision mediump float; varying vec3 vCol; void main() { gl_FragColor = vec4(vCol * 0.7, 1.0); }`;

function compile(gl, vs, fs) {
  const p = gl.createProgram();
  for (const [type, src] of [[gl.VERTEX_SHADER, vs], [gl.FRAGMENT_SHADER, fs]]) {
    const s = gl.createShader(type); gl.shaderSource(s, src); gl.compileShader(s);
    if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(s));
    gl.attachShader(p, s);
  }
  gl.linkProgram(p);
  if (!gl.getProgramParameter(p, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(p));
  return p;
}

// Column-major 4x4 helpers.
function mul(a, b) {
  const o = new Float32Array(16);
  for (let c = 0; c < 4; c++) for (let r = 0; r < 4; r++) {
    let s = 0; for (let k = 0; k < 4; k++) s += a[k * 4 + r] * b[c * 4 + k]; o[c * 4 + r] = s;
  }
  return o;
}
function perspective(fovy, aspect, near, far) {
  const f = 1 / Math.tan(fovy / 2), o = new Float32Array(16);
  o[0] = f / aspect; o[5] = f; o[10] = (far + near) / (near - far); o[11] = -1; o[14] = 2 * far * near / (near - far);
  return o;
}
function rotX(a) { const c = Math.cos(a), s = Math.sin(a); return new Float32Array([1, 0, 0, 0, 0, c, s, 0, 0, -s, c, 0, 0, 0, 0, 1]); }
function rotY(a) { const c = Math.cos(a), s = Math.sin(a); return new Float32Array([c, 0, -s, 0, 0, 1, 0, 0, s, 0, c, 0, 0, 0, 0, 1]); }
function translate(x, y, z) { const o = new Float32Array([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, x, y, z, 1]); return o; }
// Frame (x, y, z) -> view (x, z, -y): the four-fold axis points up the screen.
const ZUP = new Float32Array([1, 0, 0, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 1]);

export function webglAvailable() {
  try { const c = document.createElement('canvas'); return !!(c.getContext('webgl') || c.getContext('experimental-webgl')); }
  catch { return false; }
}

export function createViewer(canvas, { onPick = null, background = [0.02, 0.03, 0.06] } = {}) {
  const gl = canvas.getContext('webgl', { antialias: true }) || canvas.getContext('experimental-webgl');
  if (!gl) throw new Error('WebGL is not available in this browser');
  const uint32 = !!gl.getExtension('OES_element_index_uint');
  const progP = compile(gl, VS_P, FS_P), progL = compile(gl, VS_L, FS_L);
  const buf = { pos: gl.createBuffer(), col: gl.createBuffer(), size: gl.createBuffer(), idx: gl.createBuffer() };
  const cam = { yaw: 0.6, pitch: 0.12, dist: 400, center: [0, 0, 0], radius: 100, spin: false };
  let n = 0, nIdx = 0, idxType = gl.UNSIGNED_SHORT, positions = new Float32Array(0), raf = 0, alive = true, last = 0;

  function setBuffer(b, data) { gl.bindBuffer(gl.ARRAY_BUFFER, b); gl.bufferData(gl.ARRAY_BUFFER, data, gl.DYNAMIC_DRAW); }
  function fit() {
    let lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
    for (let i = 0; i < n; i++) for (let a = 0; a < 3; a++) { const v = positions[3 * i + a]; if (v < lo[a]) lo[a] = v; if (v > hi[a]) hi[a] = v; }
    cam.center = lo.map((v, a) => (v + hi[a]) / 2);
    let r = 1; for (let i = 0; i < n; i++) r = Math.max(r, Math.hypot(positions[3 * i] - cam.center[0], positions[3 * i + 1] - cam.center[1], positions[3 * i + 2] - cam.center[2]));
    cam.radius = r; cam.dist = r / Math.tan(Math.PI / 9) * 1.02;
  }
  function matrices() {
    const w = canvas.clientWidth || 300, h = canvas.clientHeight || 200;
    const P = perspective(Math.PI * 2 / 9, w / h, Math.max(1, cam.dist - 3 * cam.radius), cam.dist + 3 * cam.radius);
    const fitScale = w < h ? h / w : 1;                         // keep the model inside a narrow (phone) view
    let MV = translate(-cam.center[0], -cam.center[1], -cam.center[2]);
    MV = mul(ZUP, MV); MV = mul(rotY(cam.yaw), MV); MV = mul(rotX(cam.pitch), MV);
    MV = mul(translate(0, 0, -cam.dist * fitScale), MV);
    return { P, MV, w, h };
  }
  function draw() {
    const dpr = Math.min(devicePixelRatio || 1, 2), w = canvas.clientWidth, h = canvas.clientHeight;
    if (canvas.width !== Math.round(w * dpr) || canvas.height !== Math.round(h * dpr)) { canvas.width = Math.round(w * dpr); canvas.height = Math.round(h * dpr); }
    gl.viewport(0, 0, canvas.width, canvas.height);
    gl.clearColor(...background, 1); gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT); gl.enable(gl.DEPTH_TEST);
    if (!n) return;
    const { P, MV } = matrices();
    const bind = (prog, name, b, size) => { const loc = gl.getAttribLocation(prog, name); if (loc < 0) return; gl.bindBuffer(gl.ARRAY_BUFFER, b); gl.enableVertexAttribArray(loc); gl.vertexAttribPointer(loc, size, gl.FLOAT, false, 0, 0); };
    if (nIdx) {
      gl.useProgram(progL); bind(progL, 'aPos', buf.pos, 3); bind(progL, 'aCol', buf.col, 3);
      gl.uniformMatrix4fv(gl.getUniformLocation(progL, 'uMV'), false, MV); gl.uniformMatrix4fv(gl.getUniformLocation(progL, 'uP'), false, P);
      gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, buf.idx); gl.drawElements(gl.LINES, nIdx, idxType, 0);
    }
    gl.useProgram(progP); bind(progP, 'aPos', buf.pos, 3); bind(progP, 'aCol', buf.col, 3); bind(progP, 'aSize', buf.size, 1);
    gl.uniformMatrix4fv(gl.getUniformLocation(progP, 'uMV'), false, MV); gl.uniformMatrix4fv(gl.getUniformLocation(progP, 'uP'), false, P);
    gl.uniform1f(gl.getUniformLocation(progP, 'uScale'), canvas.height / (2 * Math.tan(Math.PI / 9)));
    const d = cam.dist * (canvas.clientWidth < canvas.clientHeight ? canvas.clientHeight / canvas.clientWidth : 1);
    gl.uniform2f(gl.getUniformLocation(progP, 'uFog'), d - cam.radius * 0.2, d + cam.radius * 1.2);
    gl.drawArrays(gl.POINTS, 0, n);
  }
  const reduced = matchMedia('(prefers-reduced-motion: reduce)').matches;
  function loop(t) {
    if (!alive) return;
    if (cam.spin && !reduced && last) cam.yaw += (t - last) * 0.00025;
    last = t; draw(); raf = requestAnimationFrame(loop);
  }
  raf = requestAnimationFrame(loop);

  // Interaction
  let drag = null, pinch = null;
  const down = e => { canvas.setPointerCapture?.(e.pointerId); drag = { x: e.clientX, y: e.clientY, moved: false }; canvas.style.cursor = 'grabbing'; };
  const move = e => {
    if (!drag) return;
    const dx = e.clientX - drag.x, dy = e.clientY - drag.y;
    if (Math.abs(dx) + Math.abs(dy) > 3) drag.moved = true;
    cam.yaw += dx * 0.008; cam.pitch = Math.max(-1.5, Math.min(1.5, cam.pitch + dy * 0.008)); drag.x = e.clientX; drag.y = e.clientY;
  };
  const up = e => { if (drag && !drag.moved && onPick) pick(e); drag = null; canvas.style.cursor = 'grab'; };
  const wheel = e => { e.preventDefault(); cam.dist *= Math.exp(e.deltaY * 0.001); clampDist(); };
  const key = e => {
    const k = { ArrowLeft: () => cam.yaw -= 0.1, ArrowRight: () => cam.yaw += 0.1, ArrowUp: () => cam.pitch -= 0.1, ArrowDown: () => cam.pitch += 0.1,
      '+': () => cam.dist /= 1.15, '=': () => cam.dist /= 1.15, '-': () => cam.dist *= 1.15 }[e.key];
    if (k) { e.preventDefault(); k(); clampDist(); }
  };
  const touch = e => {
    if (e.touches.length === 2) {
      e.preventDefault();
      const d = Math.hypot(e.touches[0].clientX - e.touches[1].clientX, e.touches[0].clientY - e.touches[1].clientY);
      if (pinch) { cam.dist *= pinch / d; clampDist(); } pinch = d;
    } else pinch = null;
  };
  const clampDist = () => { cam.dist = Math.max(cam.radius * 0.6, Math.min(cam.radius * 12, cam.dist)); };
  canvas.addEventListener('pointerdown', down); canvas.addEventListener('pointermove', move);
  canvas.addEventListener('pointerup', up); canvas.addEventListener('pointercancel', () => { drag = null; });
  canvas.addEventListener('wheel', wheel, { passive: false }); canvas.addEventListener('keydown', key);
  canvas.addEventListener('touchmove', touch, { passive: false }); canvas.addEventListener('touchend', () => { pinch = null; });
  canvas.style.touchAction = 'none'; canvas.style.cursor = 'grab';
  if (!canvas.hasAttribute('tabindex')) canvas.tabIndex = 0;

  let sizes = new Float32Array(0);
  function pick(e) {
    const rect = canvas.getBoundingClientRect(), { P, MV, w, h } = matrices(), M = mul(P, MV);
    const px = e.clientX - rect.left, py = e.clientY - rect.top;
    let best = -1, bestD = 14 * 14, bestZ = Infinity;
    for (let i = 0; i < n; i++) {
      if (!(sizes[i] > 0)) continue;
      const x = positions[3 * i], y = positions[3 * i + 1], z = positions[3 * i + 2];
      const cw = M[3] * x + M[7] * y + M[11] * z + M[15];
      const sx = ((M[0] * x + M[4] * y + M[8] * z + M[12]) / cw * 0.5 + 0.5) * w;
      const sy = (0.5 - (M[1] * x + M[5] * y + M[9] * z + M[13]) / cw * 0.5) * h;
      const d2 = (sx - px) ** 2 + (sy - py) ** 2;
      if (d2 < bestD - 4 || (d2 < bestD + 4 && cw < bestZ)) { best = i; bestD = Math.min(bestD, d2); bestZ = cw; }
    }
    onPick(best);
  }

  return {
    setModel({ positions: p, colors, sizes: s, bonds }) {
      n = p.length / 3; positions = Float32Array.from(p); sizes = Float32Array.from(s);
      setBuffer(buf.pos, positions); setBuffer(buf.col, Float32Array.from(colors)); setBuffer(buf.size, sizes);
      const big = n > 65535;
      if (big && !uint32) bonds = null;
      idxType = big ? gl.UNSIGNED_INT : gl.UNSIGNED_SHORT;
      nIdx = bonds ? bonds.length : 0;
      if (nIdx) { gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, buf.idx); gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, big ? Uint32Array.from(bonds) : Uint16Array.from(bonds), gl.STATIC_DRAW); }
      fit();
    },
    updatePositions(p) { positions = Float32Array.from(p); setBuffer(buf.pos, positions); },
    setColors(c) { setBuffer(buf.col, Float32Array.from(c)); },
    setSizes(s) { sizes = Float32Array.from(s); setBuffer(buf.size, sizes); },
    spin(on) { cam.spin = !!on; },
    view(yaw, pitch) { cam.yaw = yaw; cam.pitch = pitch; },
    fit,
    dispose() { alive = false; cancelAnimationFrame(raf); gl.getExtension('WEBGL_lose_context')?.loseContext(); },
  };
}
