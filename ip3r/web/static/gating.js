// Gating page: steady-state open probability against cytosolic Ca2+ for the three published IP3R models, computed live
// in the browser (gating_engine.js) from the registered constants (/api/gating), with Python's bell measurements beside
// the JavaScript's. #/gating/<dyk|mak|pd|compare>
import { route, api, app, esc, go } from './core.js';
import { curve, bellAt, pdParkFraction, pdDriveOpen } from './gating_engine.js';
import { plot, legend, onResize } from './plot.js';

const MODELS = [['dyk', 'De Young–Keizer'], ['mak', 'Mak 1998'], ['pd', 'Park/drive'], ['compare', 'All three']];
const COL = { dyk: '#5fe3ff', mak: '#ffb547', pd: '#ff5fa2' };
const SHORT = Object.fromEntries(MODELS);
const fmtC = c => c === null ? '–' : c < 0.1 ? c.toPrecision(3) : c.toFixed(c < 10 ? 3 : 1);

route('gating', async (args) => {
  const G = await api('/gating');
  const model = MODELS.some(m => m[0] === args[0]) ? args[0] : 'dyk';
  const [pLo, pHi] = G.compare_ip3;
  app().innerHTML = `
  <p class="kicker">Gating · steady state at clamped Ca²⁺ and IP3</p>
  <h1 class="glow">The bell-shaped calcium response</h1>
  <p class="lead">A little cytosolic Ca²⁺ helps the receptor open and a lot shuts it, so the open probability rises and falls with Ca²⁺.
  Three published models describe it; IP3 moves the curve differently in each. The curves are computed in your browser from the
  repository's registered constants.</p>
  <div class="seg" role="group" aria-label="Model">${MODELS.map(([k, t]) => `<button data-m="${k}" class="${k === model ? 'on' : ''}">${t}</button>`).join('')}</div>
  <div class="side" style="margin-top:14px">
    <div class="panel ctlpanel">
      ${model === 'compare' ? `<p class="small">Each model at Mak's two IP3 levels, ${pLo} µM (dashed) and ${pHi} µM (solid), each curve scaled
        to its own peak so the shapes compare.</p>` : `
      <label for="ip3">IP3 <b id="ip3v"></b></label><input type="range" id="ip3" min="-2" max="2" step="0.01" value="0">
      <dl class="kv small" id="bell"></dl>`}
      <p class="small mute" id="about"></p>
    </div>
    <div class="panel"><canvas class="plot tall" id="po" role="img" aria-label="Open probability against cytosolic calcium"></canvas><div id="leg"></div></div>
  </div>
  <div class="two">
    <div class="panel scroll"><h3>Browser vs Python: the bell measured</h3><table id="cmp"></table>
      <p class="small mute">Python: <code>physics.bell.measure_bell</code> on the repo's model functions, computed by the server;
      browser: the same search (801-point grid, then refined) on the JavaScript port.</p></div>
    <div class="panel scroll"><h3>Constants</h3><div id="consts"></div></div>
  </div>
  <div class="panel"><h3>Sources</h3><ul class="small">${Object.entries(G.references).map(([k, t]) => `<li><code>${esc(k)}</code> ${esc(t)}</li>`).join('')}</ul></div>`;
  document.querySelectorAll('[data-m]').forEach(b => b.onclick = () => go('#/gating/' + b.dataset.m));

  const keys = model === 'compare' ? ['dyk', 'mak', 'pd'] : [model];
  const about = { dyk: 'Each subunit has an IP3 site and two Ca²⁺ sites (activating, inhibitory); the channel conducts when three of four subunits are active (Li & Rinzel). IP3 moves both flanks.',
    mak: 'A biphasic Hill fit to single IP3R-1 channels; only the inhibitory constant depends on IP3, so IP3 moves the falling flank alone.',
    pd: `Six states in two modes; Ca²⁺ and IP3 act only on the switch between park and drive. A drive-mode channel is open ${(100 * pdDriveOpen(G.models.pd.values)).toFixed(0)} % of the time.`,
    compare: 'De Young–Keizer moves both flanks with IP3; Mak 1998 moves the inhibitory flank only; park/drive moves mostly the inhibitory flank.' };
  document.getElementById('about').textContent = about[model];

  // Python vs JS bell table at the reference IP3 levels.
  const rows = keys.flatMap(m => G.reference_bells[m].map(ref => {
    const js = bellAt(m, ref.ip3, G.models[m].values);
    return `<tr><td>${esc(SHORT[m])}</td><td>${ref.ip3}</td>
      <td class="num">${fmtC(ref.c_peak)} / ${fmtC(js.c_peak)}</td><td class="num">${ref.po_peak.toPrecision(4)} / ${js.po_peak.toPrecision(4)}</td>
      <td class="num">${fmtC(ref.c_half_act)}–${fmtC(ref.c_half_inh)} / ${fmtC(js.c_half_act)}–${fmtC(js.c_half_inh)}</td></tr>`;
  }));
  document.getElementById('cmp').innerHTML = `<thead><tr><th>Model</th><th>IP3 µM</th><th>Peak Ca²⁺ µM (Py / JS)</th><th>Peak P_open</th><th>Half-points µM</th></tr></thead><tbody>${rows.join('')}</tbody>`;
  document.getElementById('consts').innerHTML = keys.map(m => `<h3>${esc(G.models[m].label)}</h3><table><tbody>${Object.entries(G.models[m].values).map(([k, val]) => {
    const meta = G.models[m].meta[k];
    return `<tr><td><code>${esc(meta.key)}</code></td><td class="num">${val}</td><td>${esc(meta.unit)}</td><td class="small" title="${esc(meta.source_note)}">${esc(meta.citation)}</td></tr>`;
  }).join('')}</tbody></table>`).join('');

  const draw = () => {
    const series = [];
    if (model === 'compare') {
      for (const m of keys) for (const [p, dash] of [[pLo, [5, 4]], [pHi, null]]) {
        const c = curve(m, p, G.models[m].values, { lo: -2.5, hi: 3 }), peak = Math.max(...c.po);
        series.push({ x: c.ca, y: c.po.map(v => v / peak), color: COL[m], dash, width: 2, label: dash ? null : SHORT[m] });
      }
      plot(document.getElementById('po'), { series, xlog: true, xlim: [10 ** -2.5, 1e3], ylim: [0, 1.05], xlabel: 'cytosolic Ca²⁺ (µM)', ylabel: 'P_open / own peak' });
    } else {
      const p = 10 ** +document.getElementById('ip3').value, g = G.models[model].values;
      document.getElementById('ip3v').textContent = `${p < 0.1 ? p.toPrecision(2) : p.toFixed(2)} µM`;
      const c = curve(model, p, g, { lo: -3, hi: 3 }), b = bellAt(model, p, g);
      series.push({ x: c.ca, y: c.po, color: COL[model], width: 2.6, label: `P_open at ${p.toPrecision(2)} µM IP3` });
      if (model === 'pd') series.push({ x: c.ca, y: c.ca.map(x => pdParkFraction(x, p, g)), color: '#9aa3b5', dash: [5, 4], label: 'fraction parked' });
      plot(document.getElementById('po'), { series, xlog: true, xlim: [1e-3, 1e3], ylim: [0, model === 'pd' ? 1.02 : Math.max(0.05, b.po_peak * 1.15)],
        xlabel: 'cytosolic Ca²⁺ (µM)', ylabel: 'open probability', vlines: [{ x: b.c_peak, color: 'rgba(255,255,255,.3)' }] });
      document.getElementById('bell').innerHTML = `<dt>peak</dt><dd class="num">${b.po_peak.toPrecision(3)} at ${fmtC(b.c_peak)} µM</dd>
        <dt>half-activation</dt><dd class="num">${fmtC(b.c_half_act)} µM</dd><dt>half-inhibition</dt><dd class="num">${fmtC(b.c_half_inh)} µM</dd>
        <dt>width</dt><dd class="num">${b.c_half_inh ? Math.log10(b.c_half_inh / b.c_half_act).toFixed(2) + ' decades' : '–'}</dd>`;
    }
    document.getElementById('leg').innerHTML = legend(series);
  };
  document.getElementById('ip3')?.addEventListener('input', draw);
  draw(); onResize(document.getElementById('po'), draw);
});
