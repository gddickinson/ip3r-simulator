// Transition page: 8TKG (resting) -> 8TKF (activated), the restrained morph on the matched C-alpha basis, coloured by
// each residue's displacement on the fixed 0-25 Å scale; the gate along the path; and whether 8TKG's elastic-network
// modes point towards 8TKF. All numbers are the Python code's (python -m ip3r transition 8TKG 8TKF --gate).
import { route, api, app, esc, onLeave } from './core.js';
import { createViewer, webglAvailable } from './gl_view.js';
import { decodeCoords } from './anm_engine.js';
import { ramp, rampCss, BEAD } from './mol.js';
import { plot, onResize } from './plot.js';

const LABELS = { nterm_trefoil: 'β-trefoil (IP3-binding core)', MIR: 'MIR', RIH_N: 'RIH, N', RIH_C: 'RIH, C', RIH_assoc: 'RIH-associated',
  channel: 'pore domain', luminal_loop: 'luminal loop', selectivity_filter: 'selectivity filter', gate: 'gate', linker: 'linkers' };

route('transition', async () => {
  const t = await api('/transition');
  const ov = t.overlap, A = ov.A_subspace, gp = t.gate_path, m = t.morph;
  const frames = m.frames.map(decodeCoords), nRes = t.n_residues, nSub = t.chains.length, n = nRes * nSub;
  const dispCol = new Float32Array(3 * n);
  for (let k = 0; k < nSub; k++) for (let i = 0; i < nRes; i++) dispCol.set(ramp(t.residue_displacement[i] / t.displacement_max), 3 * (k * nRes + i));
  const bonds = [];
  for (let k = 0; k < nSub; k++) for (let i = 1; i < nRes; i++) if (t.residues[i] - t.residues[i - 1] === 1) bonds.push(k * nRes + i - 1, k * nRes + i);
  const means = Object.entries(t.element_means).sort((a, b) => b[1] - a[1]);
  app().innerHTML = `
  <p class="kicker">Transition · ${t.start} → ${t.end}</p>
  <h1 class="glow">From resting to open</h1>
  <p class="lead">The same ${nRes.toLocaleString('en')} residues on each of the four subunits of resting ${t.start} and activated ${t.end},
  the end superposed onto the start on the pore domain, and a restrained interpolation between them. ${esc(m.note)}</p>
  <div class="side">
    <div class="panel ctlpanel">
      <label for="frac">Path fraction <b id="fv">0.00</b></label>
      <input type="range" id="frac" min="0" max="1" step="0.005" value="0">
      <button id="play">Play</button>
      <dl class="kv small" style="margin-top:12px">
        <dt>gate</dt><dd><b class="num" id="gate"></b></dd><dt>filter</dt><dd id="filt"></dd><dt>lining</dt><dd id="lin"></dd></dl>
      <div class="legend"><span><i style="background:${rampCss(0)}"></i>0 Å (still)</span><span><i style="background:${rampCss(0.5)}"></i>${t.displacement_max / 2} Å</span>
        <span><i style="background:${rampCss(1)}"></i>≥ ${t.displacement_max} Å</span></div>
      <p class="small mute">Fit on the pore domain: RMSD ${t.fit_rmsd.toFixed(2)} Å there, ${t.rmsd.toFixed(2)} Å overall.</p>
    </div>
    <div class="viewer"><canvas id="gl" role="img" aria-label="Morph of the IP3 receptor from resting 8TKG to activated 8TKF, coloured by displacement"></canvas>
      <div class="hud">cytosol ↑ · colour = how far each residue moves</div></div>
  </div>
  <div class="two">
    <div class="panel"><h3>Do the resting structure's modes point to the open one?</h3>
      <canvas class="plot" id="ov" role="img" aria-label="Cumulative overlap of the symmetric A modes with the observed displacement"></canvas>
      <div class="legend"><span><i style="background:#ffb547"></i>collective A modes together</span><span><i style="background:#9aa3b5"></i>random direction, same symmetry</span></div>
      <p class="small">The ${A.index.length} collective A modes (all four subunits move alike) among ${t.start}'s lowest ${ov.symmetry.length}
      reach <b class="num">${A.cumulative.at(-1).toFixed(2)}</b> of the observed movement; a random direction with the same symmetry reaches
      <b class="num">${A.null.at(-1).toFixed(2)}</b>. The displacement is ${(100 * ov.irrep_fraction.A).toFixed(1)} % A-symmetric
      (${(100 * ov.rigid_fraction).toFixed(1)} % rigid body removed first). Best single mode #${ov.best + 1}, overlap ${ov.overlap[ov.best].toFixed(2)}.</p></div>
    <div class="panel"><h3>The gate along the path</h3>
      <canvas class="plot" id="gp" role="img" aria-label="Gate and filter radius along the morph"></canvas>
      <div class="legend"><span><i style="background:#ffb547"></i>gate, atoms interpolated</span><span><i style="background:#ff5fa2"></i>gate, side chains rigid</span><span><i style="background:#5fa8ff"></i>filter</span></div>
      <p class="small">Half-way open at <b class="num">${gp.half_open.toFixed(2)}</b> of the path; the gate goes from ${gp.gate[0].toFixed(2)} to
      ${gp.gate.at(-1).toFixed(2)} Å while the filter stays near ${gp.filter[0].toFixed(1)} Å.</p></div>
  </div>
  <div class="panel scroll"><h3>Mean displacement by region (Å)</h3><table><tbody>
    ${means.map(([k, v]) => `<tr><td>${esc(LABELS[k] || k)}</td><td class="num">${v.toFixed(1)}</td><td style="width:55%"><span class="bar" style="width:${Math.min(100, 100 * v / t.displacement_max)}%;background:${rampCss(v / t.displacement_max)}"></span></td></tr>`).join('')}
  </tbody></table></div>
  <details class="panel"><summary>Python report</summary><pre>${esc(ov.report.join('\n'))}</pre></details>`;

  const draw = () => {
    plot(document.getElementById('ov'), { series: [
      { x: A.cumulative.map((_, i) => i + 1), y: A.cumulative, color: '#ffb547', width: 2.4 },
      { x: A.cumulative.map((_, i) => i + 1), y: A.cumulative, color: '#ffb547', points: true },
      { x: A.null.map((_, i) => i + 1), y: A.null, color: '#9aa3b5', dash: [4, 4] }],
      xlabel: 'collective A modes added, lowest first', ylabel: 'cumulative overlap', ylim: [0, 1], xlim: [0.5, A.cumulative.length + 0.5],
      hlines: [{ y: A.ceiling, color: 'rgba(255,255,255,.25)' }] });
    plot(document.getElementById('gp'), { series: [
      { x: gp.fraction, y: gp.gate, color: '#ffb547', width: 2.4 }, { x: gp.fraction, y: gp.rigid_gate, color: '#ff5fa2', dash: [5, 4] },
      { x: gp.fraction, y: gp.filter, color: '#5fa8ff' }],
      xlabel: 'fraction of the path', ylabel: 'radius (Å)', ylim: [0, 7], xlim: [0, 1], vlines: [{ x: gp.half_open, color: 'rgba(255,255,255,.3)' }] });
  };
  draw(); onResize(document.getElementById('ov'), draw);

  const slider = document.getElementById('frac'), playBtn = document.getElementById('play');
  const readout = f => {
    const i = Math.round(f * (gp.fraction.length - 1));
    document.getElementById('fv').textContent = f.toFixed(2);
    document.getElementById('gate').textContent = gp.gate[i].toFixed(2) + ' Å';
    document.getElementById('filt').textContent = gp.filter[i].toFixed(2) + ' Å';
    document.getElementById('lin').textContent = gp.lining[i].join(' ');
  };
  readout(0);
  const canvas = document.getElementById('gl');
  if (!webglAvailable()) { canvas.parentElement.innerHTML = '<p class="err" style="padding:20px">This page needs WebGL.</p>'; return; }
  const v = createViewer(canvas);
  v.setModel({ positions: frames[0], colors: dispCol, sizes: new Float32Array(n).fill(BEAD), bonds });
  const pos = new Float32Array(3 * n);
  const show = f => {
    // interpolate linearly between the two stored frames either side of f (on their own path fractions)
    let j = 0; while (j < m.fraction.length - 2 && m.fraction[j + 1] < f) j++;
    const span = m.fraction[j + 1] - m.fraction[j], u = span > 0 ? Math.min(1, Math.max(0, (f - m.fraction[j]) / span)) : 0;
    const a = frames[j], b = frames[j + 1];
    for (let k = 0; k < pos.length; k++) pos[k] = a[k] + (b[k] - a[k]) * u;
    v.updatePositions(pos); readout(f);
  };
  slider.oninput = () => show(+slider.value);
  let playing = false, dir = 1, raf = 0, lastT = 0;
  const step = ts => {
    if (!playing) return;
    const f = Math.min(1, Math.max(0, +slider.value + dir * (ts - (lastT || ts)) / 4000)); lastT = ts;
    if (f >= 1 || f <= 0) dir = -dir;
    slider.value = f; show(f); raf = requestAnimationFrame(step);
  };
  playBtn.onclick = () => { playing = !playing; playBtn.textContent = playing ? 'Pause' : 'Play'; lastT = 0; if (playing) raf = requestAnimationFrame(step); };
  onLeave(() => { playing = false; cancelAnimationFrame(raf); v.dispose(); });
});
