// Turning a precomputed deposit (/api/structure/PDB) into beads, bonds and colours for gl_view.js.
// The colour rules are the desktop viewer's (ip3r/render/colormaps.py): fixed scales, grey where a value is missing or
// the deposit is not in a human paralog's numbering, never the low end of a scale.
import { decodeCoords } from './anm_engine.js';

// colormaps._RAMP: blue -> pale -> red, and MISSING grey.
const RAMP = [[0.23, 0.30, 0.75], [0.55, 0.69, 0.99], [0.87, 0.87, 0.87], [0.96, 0.60, 0.48], [0.71, 0.02, 0.15]];
export const MISSING = [0.42, 0.43, 0.46];
export const GOLD = [1.0, 0.82, 0.25];
export function ramp(t) {
  if (t === null || t === undefined || !Number.isFinite(t)) return MISSING;
  const x = Math.min(1, Math.max(0, t)) * (RAMP.length - 1), lo = Math.floor(x), hi = Math.min(lo + 1, RAMP.length - 1), f = x - lo;
  return RAMP[lo].map((v, a) => v * (1 - f) + RAMP[hi][a] * f);
}
export const rampCss = t => `rgb(${ramp(t).map(v => Math.round(v * 255)).join(',')})`;
export const rgbCss = c => `rgb(${c.map(v => Math.round(v * 255)).join(',')})`;

export const BEAD = 4.4;          // Å, drawn bead diameter for a C-alpha
const BOND_MAX = 4.3;             // Å, consecutive C-alphas closer than this are joined

// Flatten the four subunits: positions, which subunit and residue each bead is, and backbone bonds.
export function beads(record) {
  const pos = [], sub = [], resi = [], element = [], jsd = [], bfac = [];
  record.subunits.forEach((s, k) => {
    const xyz = decodeCoords(s.xyz);
    s.resi.forEach((r, i) => {
      pos.push(xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2]); sub.push(k); resi.push(r);
      element.push(s.element ? s.element[i] : -1); jsd.push(s.jsd ? s.jsd[i] : null); bfac.push(s.bfactor[i]);
    });
  });
  const positions = Float32Array.from(pos), n = sub.length, bonds = [];
  for (let i = 1; i < n; i++) {
    if (sub[i] !== sub[i - 1] || resi[i] - resi[i - 1] !== 1) continue;
    const d = Math.hypot(pos[3 * i] - pos[3 * i - 3], pos[3 * i + 1] - pos[3 * i - 2], pos[3 * i + 2] - pos[3 * i - 1]);
    if (d < BOND_MAX) bonds.push(i - 1, i);
  }
  const ligands = [];
  for (const l of record.ligands || []) { const x = decodeCoords(l.xyz); for (let i = 0; i < x.length; i++) ligands.push(x[i]); }
  return { positions, sub, resi, element, jsd, bfac, bonds, n, ligands: Float32Array.from(ligands) };
}

// Index of (subunit, residue) -> bead.
export function beadIndex(b) {
  const m = new Map();
  for (let i = 0; i < b.n; i++) m.set(`${b.sub[i]}:${b.resi[i]}`, i);
  return m;
}

function percentile(values, q) {
  const v = values.filter(Number.isFinite).sort((a, b) => a - b);
  if (!v.length) return 0;
  const x = (v.length - 1) * q, lo = Math.floor(x), hi = Math.ceil(x);
  return v[lo] + (v[hi] - v[lo]) * (x - lo);
}

// Colours per bead for a scheme: 'element' | 'subunit' | 'conservation' | 'bfactor' | 'values' (with opts.values, 0..1).
export function colours(b, scheme, palette, opts = {}) {
  const out = new Float32Array(3 * b.n);
  const [lo, hi] = palette.conservation_range;
  const b2 = percentile(b.bfac, 0.02), b98 = percentile(b.bfac, 0.98);
  for (let i = 0; i < b.n; i++) {
    let c;
    if (scheme === 'subunit') c = palette.chains[b.sub[i] % palette.chains.length];
    else if (scheme === 'element') c = b.element[i] >= 0 ? palette.elements[b.element[i]].rgb : MISSING;
    else if (scheme === 'conservation') c = b.jsd[i] === null ? MISSING : ramp((b.jsd[i] - lo) / (hi - lo));
    else if (scheme === 'bfactor') c = ramp((b.bfac[i] - b2) / Math.max(b98 - b2, 1e-6));
    else if (scheme === 'values') c = ramp(opts.values[i]);
    else c = MISSING;
    out.set(c, 3 * i);
  }
  return out;
}

// Assemble the drawable model: protein beads (+ highlighted residues) and IP3 atoms as small gold beads.
export function model(b, cols, { visible = null, highlight = null, showLigands = true } = {}) {
  const nl = showLigands ? b.ligands.length / 3 : 0, n = b.n + nl;
  const positions = new Float32Array(3 * n), colors = new Float32Array(3 * n), sizes = new Float32Array(n);
  positions.set(b.positions); colors.set(cols);
  for (let i = 0; i < b.n; i++) {
    const on = !visible || visible[b.sub[i]];
    sizes[i] = on ? (highlight && highlight.has(i) ? BEAD * 1.7 : BEAD) : 0;
    if (on && highlight && highlight.has(i)) colors.set(GOLD, 3 * i);
  }
  for (let j = 0; j < nl; j++) {
    positions.set(b.ligands.subarray(3 * j, 3 * j + 3), 3 * (b.n + j)); colors.set(GOLD, 3 * (b.n + j)); sizes[b.n + j] = 2.2;
  }
  const bonds = visible ? b.bonds.filter((_, k) => visible[b.sub[b.bonds[k - (k % 2)]]]) : b.bonds;
  return { positions, colors, sizes, bonds };
}

// Residue numbers of the IP3 contacts on each subunit (own-subunit and cross-subunit), as bead indices.
export function contactBeads(record, b) {
  const idx = beadIndex(b), out = new Set();
  record.chains.forEach((ch, k) => {
    const c = record.ip3_contacts[ch];
    if (!c) return;
    for (const lab of [...c.same_subunit]) { const r = +lab.replace(/^\D+/, ''); const i = idx.get(`${k}:${r}`); if (i !== undefined) out.add(i); }
  });
  return out;
}
