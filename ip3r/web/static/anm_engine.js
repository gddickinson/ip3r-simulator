// Anisotropic network model in the browser: a port of ip3r/physics/anm.py (build_hessian, ANM.calc_modes,
// ANM.label_symmetry, ModeSet.collectivity) for a reduced C-alpha network. Pure functions, no DOM; tested by
// anm_engine.test.mjs against tests/fixtures/anm_reduced.json (the repo's Python ANM on the same sites).
//
// Python solves the sparse Hessian by ARPACK shift-invert (sigma = -1e-6) and drops 6 x components zero modes.
// Here: a dense Cholesky factor of H + 1e-6 I, Lanczos on its inverse with full reorthogonalisation, and the
// rigid-body motions of each connected piece projected out (the same modes Python drops), so the k lowest internal
// modes come out directly. Dense storage limits this to ~2,000 sites; the page uses a few hundred.

export const IRREP_CHARACTERS = { A: 1, B: -1, E: 0 };
const SHIFT = 1e-6;

// base64 little-endian int16 -> Float64Array scaled (ip3r/web/codec.py)
export function decodeI16(b64, scale) {
  const bin = typeof atob === 'function' ? atob(b64) : Buffer.from(b64, 'base64').toString('binary');
  const n = bin.length >> 1, out = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    let v = bin.charCodeAt(2 * i) | (bin.charCodeAt(2 * i + 1) << 8);
    if (v & 0x8000) v -= 0x10000;
    out[i] = v * scale;
  }
  return out;
}
export const COORD_SCALE = 0.01;
export const decodeCoords = b64 => decodeI16(b64, COORD_SCALE);

// Residues every subunit resolves, every stride-th; coords (4m x 3, subunit blocks) as reduced_anm.reduced_sites.
export function reducedSites(record, stride) {
  const maps = record.subunits.map(s => {
    const xyz = decodeCoords(s.xyz), m = new Map();
    s.resi.forEach((r, i) => m.set(r, i));
    return { xyz, m };
  });
  let shared = record.subunits[0].resi.filter(r => maps.every(t => t.m.has(r)));
  shared = [...new Set(shared)].sort((a, b) => a - b).filter((_, i) => i % Math.max(1, stride | 0) === 0);
  const coords = new Float64Array(maps.length * shared.length * 3);
  let k = 0;
  for (const t of maps) for (const r of shared) { const i = t.m.get(r); coords.set(t.xyz.subarray(3 * i, 3 * i + 3), 3 * k++); }
  return { coords, residues: shared, nSub: maps.length };
}

// Contact pairs within cutoff (i < j).
export function contacts(coords, cutoff) {
  const n = coords.length / 3, c2 = cutoff * cutoff, pairs = [];
  for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
    const dx = coords[3 * j] - coords[3 * i], dy = coords[3 * j + 1] - coords[3 * i + 1], dz = coords[3 * j + 2] - coords[3 * i + 2];
    const d2 = dx * dx + dy * dy + dz * dz;
    if (d2 <= c2) pairs.push(i, j);
  }
  return pairs;
}

// Connected pieces: a label per site.
export function components(n, pairs) {
  const p = Int32Array.from({ length: n }, (_, i) => i);
  const find = x => { while (p[x] !== x) { p[x] = p[p[x]]; x = p[x]; } return x; };
  for (let k = 0; k < pairs.length; k += 2) { const a = find(pairs[k]), b = find(pairs[k + 1]); if (a !== b) p[a] = b; }
  const lab = new Map(), out = new Int32Array(n);
  for (let i = 0; i < n; i++) { const r = find(i); if (!lab.has(r)) lab.set(r, lab.size); out[i] = lab.get(r); }
  return { label: out, count: lab.size };
}

// Dense 3N x 3N Hessian with inverse-square springs gamma (d0/d)^2 (build_hessian).
export function buildHessian(coords, pairs, gamma, d0) {
  const n = coords.length / 3, N = 3 * n, H = new Float64Array(N * N);
  for (let k = 0; k < pairs.length; k += 2) {
    const i = pairs[k], j = pairs[k + 1];
    const d = [coords[3 * j] - coords[3 * i], coords[3 * j + 1] - coords[3 * i + 1], coords[3 * j + 2] - coords[3 * i + 2]];
    const dist = Math.max(Math.hypot(d[0], d[1], d[2]), 1e-6), s = gamma * (d0 / dist) ** 2;
    for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) {
      const v = -s * d[a] * d[b] / (dist * dist);
      H[(3 * i + a) * N + 3 * j + b] += v; H[(3 * j + a) * N + 3 * i + b] += v;
      H[(3 * i + a) * N + 3 * i + b] -= v; H[(3 * j + a) * N + 3 * j + b] -= v;
    }
  }
  return H;
}

// In-place lower Cholesky of a symmetric positive-definite N x N matrix (row-major).
export function cholesky(A, N) {
  for (let j = 0; j < N; j++) {
    const rj = j * N;
    let s = A[rj + j];
    for (let k = 0; k < j; k++) s -= A[rj + k] * A[rj + k];
    if (!(s > 0)) throw new Error('Hessian + shift is not positive definite');
    const d = Math.sqrt(s);
    A[rj + j] = d;
    for (let i = j + 1; i < N; i++) {
      const ri = i * N;
      let t = A[ri + j];
      for (let k = 0; k < j; k++) t -= A[ri + k] * A[rj + k];
      A[ri + j] = t / d;
    }
  }
  return A;
}

// Solve (L L^T) x = b in place.
function cholSolve(L, N, x) {
  for (let i = 0; i < N; i++) { const ri = i * N; let t = x[i]; for (let k = 0; k < i; k++) t -= L[ri + k] * x[k]; x[i] = t / L[ri + i]; }
  for (let i = N - 1; i >= 0; i--) { let t = x[i]; for (let k = i + 1; k < N; k++) t -= L[k * N + i] * x[k]; x[i] = t / L[i * N + i]; }
  return x;
}

const dot = (a, b) => { let s = 0; for (let i = 0; i < a.length; i++) s += a[i] * b[i]; return s; };
const axpy = (y, a, x) => { for (let i = 0; i < y.length; i++) y[i] += a * x[i]; };

// Orthonormal basis of the rigid-body motions of each connected piece (6 per piece, fewer if degenerate).
export function rigidBasis(coords, label, count) {
  const n = coords.length / 3, basis = [];
  for (let c = 0; c < count; c++) {
    const idx = []; for (let i = 0; i < n; i++) if (label[i] === c) idx.push(i);
    const cen = [0, 0, 0]; for (const i of idx) for (let a = 0; a < 3; a++) cen[a] += coords[3 * i + a] / idx.length;
    for (let a = 0; a < 3; a++) {
      const t = new Float64Array(3 * n), r = new Float64Array(3 * n);
      for (const i of idx) {
        t[3 * i + a] = 1;
        const x = [coords[3 * i] - cen[0], coords[3 * i + 1] - cen[1], coords[3 * i + 2] - cen[2]];
        const e = [0, 0, 0]; e[a] = 1;                      // e x x
        r[3 * i] = e[1] * x[2] - e[2] * x[1]; r[3 * i + 1] = e[2] * x[0] - e[0] * x[2]; r[3 * i + 2] = e[0] * x[1] - e[1] * x[0];
      }
      basis.push(t, r);
    }
  }
  const out = [];                                         // Gram-Schmidt twice, dropping null vectors
  for (const v of basis) {
    for (let pass = 0; pass < 2; pass++) for (const q of out) axpy(v, -dot(v, q), q);
    const nv = Math.sqrt(dot(v, v));
    if (nv > 1e-8) { for (let i = 0; i < v.length; i++) v[i] /= nv; out.push(v); }
  }
  return out;
}

// Cyclic Jacobi eigen-decomposition of a small symmetric matrix (m x m, row-major). Returns {values, vectors (columns)}.
export function jacobiEigen(A0, m) {
  const A = Float64Array.from(A0), V = new Float64Array(m * m);
  for (let i = 0; i < m; i++) V[i * m + i] = 1;
  for (let sweep = 0; sweep < 100; sweep++) {
    let off = 0;
    for (let p = 0; p < m; p++) for (let q = p + 1; q < m; q++) off += A[p * m + q] ** 2;
    if (off < 1e-30) break;
    for (let p = 0; p < m; p++) for (let q = p + 1; q < m; q++) {
      const apq = A[p * m + q];
      if (Math.abs(apq) < 1e-300) continue;
      const theta = (A[q * m + q] - A[p * m + p]) / (2 * apq);
      const t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
      const c = 1 / Math.sqrt(t * t + 1), s = t * c;
      for (let k = 0; k < m; k++) {
        const akp = A[k * m + p], akq = A[k * m + q];
        A[k * m + p] = c * akp - s * akq; A[k * m + q] = s * akp + c * akq;
      }
      for (let k = 0; k < m; k++) {
        const apk = A[p * m + k], aqk = A[q * m + k];
        A[p * m + k] = c * apk - s * aqk; A[q * m + k] = s * apk + c * aqk;
      }
      for (let k = 0; k < m; k++) {
        const vkp = V[k * m + p], vkq = V[k * m + q];
        V[k * m + p] = c * vkp - s * vkq; V[k * m + q] = s * vkp + c * vkq;
      }
    }
  }
  return { values: Array.from({ length: m }, (_, i) => A[i * m + i]), vectors: V };
}

// The k lowest internal modes. Returns {eigenvalues, vectors: [Float64Array 3n], iterations, converged}.
export function lowestModes(H, N, rigid, k, { tol = 1e-10, maxIter = 400 } = {}) {
  const L = Float64Array.from(H);
  for (let i = 0; i < N; i++) L[i * N + i] += SHIFT;
  cholesky(L, N);
  const deflate = v => { for (let pass = 0; pass < 2; pass++) for (const q of rigid) axpy(v, -dot(v, q), q); };
  const Q = [];
  let q = new Float64Array(N);
  for (let i = 0; i < N; i++) q[i] = Math.sin(1.7 * i + 0.3) + 0.5 * Math.cos(0.37 * i);   // fixed start
  deflate(q);
  let nq = Math.sqrt(dot(q, q)); for (let i = 0; i < N; i++) q[i] /= nq;
  const alpha = [], beta = [];
  let result = null;
  const limit = Math.min(maxIter, N - rigid.length);
  for (let j = 0; j < limit; j++) {
    Q.push(q);
    const w = cholSolve(L, N, Float64Array.from(q));
    deflate(w);
    const a = dot(w, q); alpha.push(a);
    for (let pass = 0; pass < 2; pass++) for (const v of Q) axpy(w, -dot(w, v), v);   // full reorthogonalisation
    const b = Math.sqrt(dot(w, w));
    const m = j + 1;
    if ((m >= k + 10 && m % 5 === 0) || m === limit || b < 1e-14) {
      const T = new Float64Array(m * m);
      for (let i = 0; i < m; i++) { T[i * m + i] = alpha[i]; if (i + 1 < m) T[i * m + i + 1] = T[(i + 1) * m + i] = beta[i]; }
      const { values, vectors } = jacobiEigen(T, m);
      const order = values.map((v, i) => i).sort((x, y) => values[y] - values[x]).slice(0, Math.min(k, m));
      const ok = order.every(i => Math.abs(b * vectors[(m - 1) * m + i]) <= tol * Math.abs(values[i]));
      if (ok || m === limit || b < 1e-14) {
        const vecs = order.map(i => { const v = new Float64Array(N); for (let r = 0; r < m; r++) axpy(v, vectors[r * m + i], Q[r]); return v; });
        result = { eigenvalues: order.map(i => 1 / values[i] - SHIFT), vectors: vecs, iterations: m, converged: ok };
        break;
      }
    }
    beta.push(b);
    q = w; for (let i = 0; i < N; i++) q[i] /= b;
  }
  return result;
}

// The C_n generator applied to a per-site field (apply_generator; axis = frame z).
export function applyGenerator(field, nSub = 4, power = 1) {
  const per = field.length / 3 / nSub, ang = 2 * Math.PI * power / nSub, c = Math.cos(ang), s = Math.sin(ang);
  const out = new Float64Array(field.length);
  for (let b = 0; b < nSub; b++) {
    const dst = (b + power) % nSub;
    for (let i = 0; i < per; i++) {
      const si = 3 * (b * per + i), di = 3 * (dst * per + i), x = field[si], y = field[si + 1];
      out[di] = c * x - s * y; out[di + 1] = s * x + c * y; out[di + 2] = field[si + 2];
    }
  }
  return out;
}

export function labelSymmetry(vectors, nSub = 4, tolerance = 0.15) {
  const character = vectors.map(u => dot(u, applyGenerator(u, nSub)) / dot(u, u));
  const symmetry = character.map(ch => {
    for (const [name, chi] of Object.entries(IRREP_CHARACTERS)) if (Math.abs(ch - chi) < tolerance) return name;
    return 'mixed';
  });
  return { character, symmetry };
}

// Brüschweiler's kappa per mode.
export function collectivity(v) {
  const n = v.length / 3, u2 = new Float64Array(n);
  let tot = 0;
  for (let i = 0; i < n; i++) { u2[i] = v[3 * i] ** 2 + v[3 * i + 1] ** 2 + v[3 * i + 2] ** 2; tot += u2[i]; }
  let ent = 0;
  for (let i = 0; i < n; i++) { const p = u2[i] / tot; if (p > 0) ent -= p * Math.log(p); }
  return Math.exp(ent) / n;
}

// Everything the page needs for one reduced network.
export function reducedANM(record, { stride, cutoff, gamma, d0, nModes, tolerance = 0.15 }) {
  const { coords, residues, nSub } = reducedSites(record, stride);
  const n = coords.length / 3, pairs = contacts(coords, cutoff), comp = components(n, pairs);
  const H = buildHessian(coords, pairs, gamma, d0);
  const rigid = rigidBasis(coords, comp.label, comp.count);
  const k = Math.min(nModes, 3 * n - rigid.length);
  const res = lowestModes(H, 3 * n, rigid, k);
  const vectors = res.vectors.map(v => { const s = Math.sqrt(dot(v, v)); return v.map(x => x / s); });
  const sym = labelSymmetry(vectors, nSub, tolerance);
  return { coords, residues, nSub, nSites: n, nPairs: pairs.length / 2, nComponents: comp.count,
    zeroModesDropped: rigid.length, eigenvalues: res.eigenvalues, vectors, ...sym,
    collectivity: vectors.map(collectivity), iterations: res.iterations, converged: res.converged };
}
