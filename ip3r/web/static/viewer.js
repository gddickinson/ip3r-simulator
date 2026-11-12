// Structure page: one deposit as a C-alpha bead model in WebGL, coloured as the desktop viewer colours it, with the
// measured channel summary. #/viewer/<PDB>?color=element|subunit|conservation|bfactor
import { route, api, app, esc, onLeave, go } from './core.js';
import { createViewer, webglAvailable } from './gl_view.js';
import { beads, colours, model, contactBeads, rgbCss, rampCss, MISSING, GOLD } from './mol.js';

export const SCHEMES = [['element', 'Region'], ['subunit', 'Subunit'], ['conservation', 'Conservation'], ['bfactor', 'B-factor']];

export function legendHtml(scheme, palette, rec) {
  const sw = (c, t) => `<span><i style="background:${c}"></i>${esc(t)}</span>`;
  if (scheme === 'subunit') return rec.chains.map((ch, k) => sw(rgbCss(palette.chains[k]), `subunit ${ch}`)).join('');
  if (scheme === 'element') {
    if (!rec.painted_in) return sw(rgbCss(MISSING), 'not in a human paralog numbering: grey, never guessed');
    const used = new Set(rec.subunits.flatMap(s => s.element || []));
    return palette.elements.filter((_, i) => used.has(i)).map(e => sw(rgbCss(e.rgb), e.label)).join('') + sw(rgbCss(GOLD), 'IP3 contacts and IP3');
  }
  if (scheme === 'conservation') {
    const [lo, hi] = palette.conservation_range;
    return `${sw(rampCss(0), `JSD ≤ ${lo}`)}${sw(rampCss(0.5), ((lo + hi) / 2).toFixed(2))}${sw(rampCss(1), `≥ ${hi} (most conserved)`)}${sw(rgbCss(MISSING), 'not scored')}
      <span class="mute">fixed scale, ip3r_genes S17 "${esc(palette.conservation_layer)}" layer</span>`;
  }
  return `${sw(rampCss(0), 'low')}${sw(rampCss(1), 'high')}<span class="mute">B-factor ranged over this deposit (2nd–98th percentile), as the desktop app</span>`;
}

function summaryHtml(rec) {
  const c = rec.constrictions, num = rec.numbering;
  const con = n => c[n] ? `<b class="num">${c[n].radius.toFixed(2)} Å</b> <span class="mute">${esc(c[n].residues.join(', '))}</span>` : '<span class="mute">not found</span>';
  const contacts = Object.values(rec.ip3_contacts)[0];
  return `<dl class="kv">
    <dt>Deposit</dt><dd><a href="${esc(rec.source.url)}" target="_blank" rel="noopener">PDB ${esc(rec.pdb_id)}</a> · ${esc(rec.paralog)} · ${esc(rec.organism)}</dd>
    <dt>Title</dt><dd>${esc(rec.title)}</dd>
    <dt>State</dt><dd>${esc(rec.state)}${rec.is_control ? ' <span class="tag b">open-state control</span>' : ''} · ${rec.resolution} Å</dd>
    <dt>Numbering</dt><dd>${num ? `${esc(num.paralog)}, ${(100 * num.identity).toFixed(1)} % identical over ${num.n_compared} residues${num.passed ? '' : ' (fails the rule)'}` : 'fits no human paralog'}</dd>
    <dt>Gate radius</dt><dd>${con('gate')}</dd>
    <dt>Filter radius</dt><dd>${con('filter')}</dd>
    <dt>C4 residual</dt><dd>${rec.c4_residual} Å · axis methods agree to ${rec.axis_disagreement_deg}°</dd>
    <dt>IP3 bound</dt><dd>${rec.ligands.length ? `${rec.ligands.length} sites; contacts on each subunit: ${esc(contacts.same_subunit.join(' '))}` : 'none'}</dd>
    <dt>Beads</dt><dd>${rec.subunits.reduce((a, s) => a + s.resi.length, 0).toLocaleString('en')} C-alphas (of ${rec.n_atoms.toLocaleString('en')} atoms)</dd>
  </dl>
  <p class="small mute">Radii are r_min, the distance from the four-fold axis to the nearest heavy-atom centre (ip3r_genes S0's
  definition), measured by <code>measure_channel</code>. Source file SHA-256 <code>${esc(rec.source.sha256.slice(0, 16))}…</code></p>`;
}

route('viewer', async (args, q) => {
  const [dep, palette] = await Promise.all([api('/deposits'), api('/palette')]);
  const ids = dep.deposits.map(d => d.pdb_id);
  const pdb = ids.includes((args[0] || '').toUpperCase()) ? args[0].toUpperCase() : '6DQN';
  const scheme = SCHEMES.some(s => s[0] === q.color) ? q.color : 'element';
  const rec = await api('/structure/' + pdb);
  app().innerHTML = `
  <p class="kicker">Structure · coarse-grained C-alpha model</p>
  <h1 class="glow">PDB ${esc(pdb)} <span class="mute small">${esc(rec.state)}</span></h1>
  <div class="side">
    <div class="panel ctlpanel">
      <label for="dep">Deposit</label>
      <select id="dep">${dep.deposits.map(d => `<option value="${d.pdb_id}" ${d.pdb_id === pdb ? 'selected' : ''}>${d.pdb_id} · ${esc(d.state)}</option>`).join('')}</select>
      <label>Colour by</label>
      <div class="seg" role="group" aria-label="Colour by">${SCHEMES.map(([k, t]) => `<button data-c="${k}" class="${k === scheme ? 'on' : ''}">${t}</button>`).join('')}</div>
      <label>Subunits</label>
      <div class="chk">${rec.chains.map((ch, k) => `<label><input type="checkbox" data-sub="${k}" checked> ${esc(ch)}</label>`).join('')}</div>
      <div class="chk" style="margin-top:8px">
        <label><input type="checkbox" id="contacts" checked> IP3 contacts</label>
        <label><input type="checkbox" id="spin" checked> turn</label>
      </div>
      <button class="alt" id="side">Side view</button><button class="alt" id="top">Top view</button>
      <p class="small mute" id="picked">Click a bead to name it. Drag to rotate, wheel or pinch to zoom, arrow keys too.</p>
    </div>
    <div>
      <div class="viewer"><canvas id="gl" role="img" aria-label="3-D C-alpha model of PDB ${esc(pdb)}, cytosol at the top"></canvas>
        <div class="hud">cytosol ↑ · lumen ↓ · <b>${esc(pdb)}</b></div></div>
      <div class="legend" id="legend"></div>
    </div>
  </div>
  <div class="panel">${summaryHtml(rec)}</div>
  <p class="small mute">Coordinates: wwPDB entry <a href="${esc(rec.source.url)}" target="_blank" rel="noopener">${esc(pdb)}</a> (CC0), one bead per
  C-alpha in the deposit's own axis frame (z along the four-fold axis, cytosol up). The desktop app draws every atom as cartoon,
  tube, spheres or sticks; this page draws the backbone trace only.</p>`;

  document.getElementById('dep').onchange = e => go(`#/viewer/${e.target.value}?color=${scheme}`);
  document.querySelectorAll('[data-c]').forEach(btn => btn.onclick = () => go(`#/viewer/${pdb}?color=${btn.dataset.c}`));
  const canvas = document.getElementById('gl');
  if (!webglAvailable()) { canvas.parentElement.innerHTML = '<p class="err" style="padding:20px">This page needs WebGL, which this browser does not provide.</p>'; return; }
  const b = beads(rec), cols = colours(b, scheme, palette), contacts = contactBeads(rec, b);
  const visible = rec.chains.map(() => true);
  let picked = '';
  const v = createViewer(canvas, { onPick: i => {
    const el = document.getElementById('picked');
    if (i < 0 || i >= b.n) { el.textContent = i >= b.n ? 'IP3 (bound ligand)' : 'Nothing there.'; return; }
    const e = b.element[i] >= 0 ? palette.elements[b.element[i]].label : 'not painted (numbering)';
    picked = `subunit ${rec.chains[b.sub[i]]}, residue ${b.resi[i]} · ${e}${b.jsd[i] !== null ? ` · JSD ${b.jsd[i]}` : ''} · B ${b.bfac[i]}`;
    el.textContent = picked;
  } });
  const redraw = () => v.setModel(model(b, cols, { visible, highlight: document.getElementById('contacts').checked ? contacts : null }));
  redraw(); v.spin(true);
  document.querySelectorAll('[data-sub]').forEach(cb => cb.onchange = () => { visible[+cb.dataset.sub] = cb.checked; redraw(); });
  document.getElementById('contacts').onchange = redraw;
  document.getElementById('spin').onchange = e => v.spin(e.target.checked);
  document.getElementById('side').onclick = () => v.view(0.6, 0.12);
  document.getElementById('top').onclick = () => v.view(0, -Math.PI / 2 + 0.01);
  document.getElementById('legend').innerHTML = legendHtml(scheme, palette, rec);
  onLeave(() => v.dispose());
});
