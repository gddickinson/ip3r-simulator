// Modes page: the elastic-network (ANM) modes of the tetramer with their C4 irreps. The full network (registered
// stride, cutoff and springs) is precomputed by the Python code (python -m ip3r modes PDB); a reduced network is solved
// live in the browser (anm_engine.js, a Web Worker) on any deposit, with stride and cutoff of your choice.
import { route, api, app, esc, onLeave, go } from './core.js';
import { createViewer, webglAvailable } from './gl_view.js';
import { decodeCoords, decodeI16 } from './anm_engine.js';
import { ramp, rampCss } from './mol.js';

const IRREP_COL = { A: 'a', B: 'b', E: 'c', mixed: '' };

function modeRows(m, prefix) {
  const l0 = m.eigenvalues[0];
  return m.eigenvalues.map((lam, i) => `<tr class="click" data-mode="${prefix}${i}" tabindex="0">
    <td>#${i + 1}</td><td><span class="tag ${IRREP_COL[m.symmetry[i]]}">${m.symmetry[i]}</span></td>
    <td class="num">${(lam / l0).toFixed(2)}</td><td class="num">${m.character[i].toFixed(2)}</td>
    <td class="num">${m.collectivity[i].toFixed(2)}${m.collective && !m.collective[i] ? ' <span class="mute">local</span>' : ''}</td></tr>`).join('');
}
const TABLE_HEAD = '<thead><tr><th>Mode</th><th>Irrep</th><th>λ / λ₁</th><th>χ</th><th>κ</th></tr></thead>';

route('modes', async (args) => {
  const [dep, anmP] = await Promise.all([api('/deposits'), api('/anm-parameters')]);
  const pdb = dep.modes.includes((args[0] || '').toUpperCase()) ? args[0].toUpperCase() : dep.modes[0];
  const full = await api('/modes/' + pdb);
  const P = k => anmP[k].value;
  app().innerHTML = `
  <p class="kicker">Modes · anisotropic network model with C4 symmetry</p>
  <h1 class="glow">The receptor's softest motions</h1>
  <p class="lead">Beads on the C-alphas joined by springs (stiffness γ(d₀/d)², within a cutoff) give the protein's large-scale
  normal modes. Because the receptor is a four-fold (C4) tetramer, every mode is one symmetry type: <b>A</b> (all four subunits
  move alike, the only kind that can open a four-fold pore symmetrically), <b>B</b> (neighbours opposite) or <b>E</b> (a degenerate pair).</p>
  <div class="side">
    <div class="panel ctlpanel">
      <label>Full network (precomputed)</label>
      <div class="seg" role="group" aria-label="Deposit">${dep.modes.map(p => `<button data-pdb="${p}" class="${p === pdb ? 'on' : ''}">${p}</button>`).join('')}</div>
      <p class="small mute">${full.n_sites.toLocaleString('en')} sites (stride ${P('anm.stride')} over the residues all four subunits resolve),
      cutoff ${P('anm.cutoff')} Å, d₀ ${P('anm.d0')} Å, ${full.eigenvalues.length} modes. Lowest collective A mode:
      <b>#${full.lowest_collective_A + 1}</b>.</p>
      <label for="amp">Amplitude <b id="ampv">6 Å</b></label><input type="range" id="amp" min="1" max="15" step="0.5" value="6">
      <p class="small" id="showing"></p>
    </div>
    <div class="viewer"><canvas id="gl" role="img" aria-label="The selected normal mode animated on the receptor"></canvas>
      <div class="hud">cytosol ↑ · colour = the mode's amplitude per site (ranged per mode)</div></div>
  </div>
  <div class="two">
    <div class="panel scroll"><h3>${pdb}, full network (Python)</h3><table class="modes">${TABLE_HEAD}<tbody>${modeRows(full, 'f')}</tbody></table>
      <p class="small mute">λ only as a ratio: γ sets an arbitrary stiffness without a B-factor fit. κ = Brüschweiler's collectivity
      (a mode below ${P('anm.min_collectivity')} moves a fragment, not the receptor). ${full.local_artefacts.length ? 'Local artefacts: ' + esc(full.local_artefacts.join('; ')) : ''}</p></div>
    <div class="panel"><h3>Reduced network, solved in your browser</h3>
      <div class="row">
        <div><label for="rdep">Deposit</label><select id="rdep">${dep.deposits.map(d => `<option ${d.pdb_id === pdb ? 'selected' : ''}>${d.pdb_id}</option>`).join('')}</select></div>
        <div><label for="rstride">Stride <b id="rsv">16</b></label><input type="range" id="rstride" min="12" max="32" step="1" value="16"></div>
        <div><label for="rcut">Cutoff <b id="rcv">25 Å</b></label><input type="range" id="rcut" min="12" max="40" step="0.5" value="25"></div>
      </div>
      <button id="solve">Solve</button> <span class="small mute" id="rstatus"></span>
      <div class="scroll"><table class="modes" id="rtable"></table></div>
      <p class="small mute">The same model with every <i>stride</i>-th shared residue (the JavaScript is checked against the repo's Python
      ANM on these sites: eigenvalues to 10⁻⁶, irreps exact). It illustrates, it does not replace the full network: the repo's own check
      finds stride 2 agrees with stride 1 and stride 4 already does not. A coarse network also needs a longer cutoff than the registered
      ${P('anm.cutoff')} Å to stay in one piece; each separate piece adds six zero modes, which are removed, and is reported.
      The defaults here (stride 16, cutoff 25 Å) are the web page's choice, not registered parameters.</p>
    </div>
  </div>`;

  document.querySelectorAll('[data-pdb]').forEach(b => b.onclick = () => go('#/modes/' + b.dataset.pdb));
  const canvas = document.getElementById('gl');
  if (!webglAvailable()) { canvas.parentElement.innerHTML = '<p class="err" style="padding:20px">This page needs WebGL.</p>'; return; }
  const v = createViewer(canvas);
  const ampIn = document.getElementById('amp');
  ampIn.oninput = () => { document.getElementById('ampv').textContent = ampIn.value + ' Å'; };
  let anim = null, raf = 0;
  const bondsFor = (nSub, residues) => {
    const per = residues.length, out = [];
    for (let k = 0; k < nSub; k++) for (let i = 1; i < per; i++) if (residues[i] - residues[i - 1] <= 4) out.push(k * per + i - 1, k * per + i);
    return out;
  };
  const showMode = (sites, vec, residues, nSub, label, beadSize) => {
    const n = sites.length / 3, col = new Float32Array(3 * n), amp = new Float64Array(n);
    let peak = 0; for (let i = 0; i < n; i++) { amp[i] = Math.hypot(vec[3 * i], vec[3 * i + 1], vec[3 * i + 2]); peak = Math.max(peak, amp[i]); }
    for (let i = 0; i < n; i++) col.set(ramp(amp[i] / peak), 3 * i);
    v.setModel({ positions: sites, colors: col, sizes: new Float32Array(n).fill(beadSize), bonds: bondsFor(nSub, residues) });
    anim = { sites, vec: vec.map(x => x / peak), pos: new Float32Array(3 * n) };
    document.getElementById('showing').innerHTML = `Showing <b>${esc(label)}</b>. <span class="legend"><span><i style="background:${rampCss(0)}"></i>still</span><span><i style="background:${rampCss(1)}"></i>moves most</span></span>`;
  };
  const reduced = matchMedia('(prefers-reduced-motion: reduce)').matches;
  const tick = t => {
    if (anim) {
      const s = (reduced ? 1 : Math.sin(t / 450)) * +ampIn.value;
      for (let k = 0; k < anim.pos.length; k++) anim.pos[k] = anim.sites[k] + s * anim.vec[k];
      v.updatePositions(anim.pos);
    }
    raf = requestAnimationFrame(tick);
  };
  raf = requestAnimationFrame(tick);

  const fullSites = decodeCoords(full.sites);
  const pickFull = i => showMode(fullSites, decodeI16(full.vectors[i], full.vector_scale), full.residues, 4,
    `${pdb} full network, mode #${i + 1} (${full.symmetry[i]}, κ ${full.collectivity[i].toFixed(2)})`, 4.4);
  const bindRows = (sel, fn) => document.querySelectorAll(sel).forEach(tr => {
    tr.onclick = () => { document.querySelectorAll('tr.sel').forEach(x => x.classList.remove('sel')); tr.classList.add('sel'); fn(+tr.dataset.mode.slice(1)); };
    tr.onkeydown = e => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); tr.click(); } };
  });
  bindRows('[data-mode^="f"]', pickFull);
  const first = full.lowest_collective_A ?? 0;
  document.querySelector(`[data-mode="f${first}"]`)?.classList.add('sel');
  pickFull(first);

  // Reduced network in a worker.
  const worker = new Worker(new URL('./anm_worker.js', import.meta.url), { type: 'module' });
  let reqId = 0;
  const st = document.getElementById('rstride'), cu = document.getElementById('rcut');
  st.oninput = () => { document.getElementById('rsv').textContent = st.value; };
  cu.oninput = () => { document.getElementById('rcv').textContent = cu.value + ' Å'; };
  const solve = async () => {
    const id = ++reqId, rpdb = document.getElementById('rdep').value, status = document.getElementById('rstatus');
    status.textContent = 'solving…';
    const record = await api('/structure/' + rpdb);
    worker.onmessage = e => {
      if (e.data.id !== reqId) return;
      if (e.data.error) { status.innerHTML = `<span class="err">${esc(e.data.error)}</span>`; return; }
      const r = e.data.result;
      status.textContent = `${r.nSites} sites, ${r.nPairs.toLocaleString('en')} springs, ${r.nComponents} piece${r.nComponents > 1 ? 's' : ''}` +
        ` (${r.zeroModesDropped} zero modes removed), ${r.iterations} Lanczos steps, ${Math.round(r.ms)} ms`;
      const floppy = r.eigenvalues.filter(x => x < 1e-8 * r.eigenvalues.at(-1)).length;
      const tbl = document.getElementById('rtable');
      tbl.innerHTML = TABLE_HEAD + '<tbody>' + modeRows({ eigenvalues: r.eigenvalues.map(x => Math.max(x, 1e-300)), symmetry: r.symmetry, character: r.character, collectivity: r.collectivity }, 'r') + '</tbody>';
      if (floppy) status.textContent += ` · ${floppy} near-zero mode${floppy > 1 ? 's' : ''}: a site held by too few springs`;
      bindRows('#rtable [data-mode]', i => showMode(r.coords, r.vectors[i], r.residues, r.nSub,
        `${rpdb} reduced (stride ${st.value}, cutoff ${cu.value} Å), mode #${i + 1} (${r.symmetry[i]})`, 7));
    };
    worker.postMessage({ id, record, params: { stride: +st.value, cutoff: +cu.value, gamma: P('anm.gamma'), d0: P('anm.d0'), nModes: 12,
      tolerance: P('anm.symmetry_tolerance') } });
  };
  document.getElementById('solve').onclick = solve;
  solve();
  onLeave(() => { cancelAnimationFrame(raf); worker.terminate(); v.dispose(); });
});
