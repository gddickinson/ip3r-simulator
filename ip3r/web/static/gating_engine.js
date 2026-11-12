// Steady-state IP3R gating in the browser: ports of ip3r/physics/gating.py (De Young-Keizer / Li-Rinzel),
// gating_mak.py (Mak et al. 1998), park_drive.py (Siekmann 2012 / Cao 2013, stationary by detailed balance) and
// bell.py (peak and half-peak crossings). Constants come from /api/gating (the parameter registry), never from here.
// Tested by gating_engine.test.mjs against tests/fixtures/gating_curves.json. Concentrations in µM.

// --- De Young-Keizer, reduced by Li & Rinzel ---------------------------------------------------------------
export function dykOpen(c, p, g) {
  const m = p / (p + g.d1), n = c / (c + g.d5);
  const q2 = g.d2 * (p + g.d1) / (p + g.d3), h = q2 / (q2 + c);
  return (m * n * h) ** g.subunits;
}

// --- Mak, McBride & Foskett 1998 (Eqs. 1-2) ---------------------------------------------------------------
export function makKinh(p, g) { return g.k_inf / (1 + (g.k_ip3 / p) ** g.h_ip3); }
export function makOpen(c, p, g) {
  const ki = makKinh(p, g);
  const denom = 1 + (g.k_act / c) ** g.h_act + (c / ki) ** g.h_inh;
  return Number.isFinite(denom) ? g.p_max / denom : 0;
}

// --- park/drive ---------------------------------------------------------------------------------------------
const rise = (x, k, h) => x ** h / (x ** h + k ** h);
export function pdIp3Functions(p, g) {
  return {
    V42: g.v42_amp * rise(p, g.v42_k, g.v42_h),
    k42: g.k42_base + g.k42_amp * rise(p, g.k42_k, g.k42_h),
    kn42: g.kn42_base + g.kn42_amp * rise(p, g.kn42_k, g.kn42_h),
    a42: g.a42_amp * rise(p, g.a42_k, g.a42_h),
    V24: g.v24_base + g.v24_amp * (1 - rise(p, g.v24_k, g.v24_h)),
    a24: g.a24_base + g.a24_amp * (1 - rise(p, g.a24_k, g.a24_h)),
  };
}
// Stationary occupancy of C1 C2 C3 C4 O5 O6 (detailed balance on the tree, every state relative to C2).
export function pdStationary(c, p, g) {
  const f = pdIp3Functions(p, g);
  const m24 = rise(c, g.k24, g.n24), h24 = 1 - rise(c, g.kn24, g.nn24);
  const m42 = rise(c, f.k42, g.n42), h42 = 1 - rise(c, f.kn42, g.nn42);
  const q24 = f.a24 + f.V24 * (1 - m24 * h24), q42 = f.a42 + f.V42 * m42 * h42;
  const c4 = q42 > 0 ? q24 / q42 : Infinity;
  let w;
  if (!Number.isFinite(c4)) w = [0, 0, 0, 1, g.q45 / g.q54, 0];
  else w = [g.q21 / g.q12, 1, g.q23 / g.q32, c4, c4 * g.q45 / g.q54, g.q26 / g.q62];
  const s = w.reduce((a, b) => a + b, 0);
  return w.map(x => x / s);
}
export function pdOpen(c, p, g) { const s = pdStationary(c, p, g); return s[4] + s[5]; }
export function pdParkFraction(c, p, g) { const s = pdStationary(c, p, g); return s[3] + s[4]; }
export function pdDriveOpen(g) {
  const w = [g.q21 / g.q12, 1, g.q23 / g.q32, g.q26 / g.q62];
  return w[3] / w.reduce((a, b) => a + b, 0);
}

export const MODEL_FN = { dyk: dykOpen, mak: makOpen, pd: pdOpen };

// --- bell measurement (physics/bell.py) ---------------------------------------------------------------------
const LOG_C = [-4, 4];
function golden(f, lo, hi, tol = 1e-10) {          // minimise f on [lo, hi]
  const r = (Math.sqrt(5) - 1) / 2;
  let a = lo, b = hi, x1 = b - r * (b - a), x2 = a + r * (b - a), f1 = f(x1), f2 = f(x2);
  while (b - a > tol) {
    if (f1 < f2) { b = x2; x2 = x1; f2 = f1; x1 = b - r * (b - a); f1 = f(x1); }
    else { a = x1; x1 = x2; f1 = f2; x2 = a + r * (b - a); f2 = f(x2); }
  }
  return (a + b) / 2;
}
function bisect(f, a, b, tol = 1e-12) {
  let fa = f(a);
  for (let i = 0; i < 200 && b - a > tol; i++) {
    const m = (a + b) / 2, fm = f(m);
    if ((fm < 0) === (fa < 0)) { a = m; fa = fm; } else b = m;
  }
  return (a + b) / 2;
}
// {c_peak, po_peak, c_half_act, c_half_inh (null if the curve never falls to half)}
export function measureBell(f) {
  const g = lc => f(10 ** lc);
  let best = -Infinity, ib = 0;
  const N = 801, step = (LOG_C[1] - LOG_C[0]) / (N - 1);
  for (let i = 0; i < N; i++) { const v = g(LOG_C[0] + i * step); if (v > best) { best = v; ib = i; } }
  const lo = LOG_C[0] + Math.max(ib - 1, 0) * step, hi = LOG_C[0] + Math.min(ib + 1, N - 1) * step;
  const xpk = golden(x => -g(x), lo, hi), pk = g(xpk);
  const half = x => g(x) - pk / 2;
  const xact = bisect(half, LOG_C[0], xpk);
  const xinh = half(LOG_C[1]) < 0 ? bisect(half, xpk, LOG_C[1]) : null;
  return { c_peak: 10 ** xpk, po_peak: pk, c_half_act: 10 ** xact, c_half_inh: xinh === null ? null : 10 ** xinh };
}
export function bellAt(model, p, g) { const fn = MODEL_FN[model]; return measureBell(c => fn(c, p, g)); }

// Log-spaced curve for plotting: {ca, po}
export function curve(model, p, g, { lo = -3, hi = 3, n = 241 } = {}) {
  const fn = MODEL_FN[model], ca = [], po = [];
  for (let i = 0; i < n; i++) { const c = 10 ** (lo + (hi - lo) * i / (n - 1)); ca.push(c); po.push(fn(c, p, g)); }
  return { ca, po };
}
