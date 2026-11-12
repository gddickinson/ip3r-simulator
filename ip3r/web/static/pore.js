// Pore page: the radius of every ITPR3 state's pore along the four-fold axis (ip3r_genes S11's state panel, measured
// by structure.states.state_panel), the 7T3T open-state control dashed, and the gate / filter table.
import { route, api, app, esc } from './core.js';
import { plot, legend, onResize } from './plot.js';

const COLOURS = ['#7aa2ff', '#5fe3ff', '#64d98b', '#c3e15a', '#ffd23f', '#b18cff', '#ff5fa2', '#9aa3b5'];

route('pore', async () => {
  const s = await api('/states');
  const rows = s.rows;
  const series = rows.map((r, i) => ({ x: r.profile.z, y: r.profile.r_min, label: `${r.pdb_id} ${r.state}`,
    color: r.pdb_id === '8TKF' ? '#ff5f5f' : r.control ? '#c8cbd4' : COLOURS[i % COLOURS.length],
    width: r.pdb_id === '8TKF' ? 3 : 1.6, dash: r.control ? [6, 4] : null }));
  const open = rows.filter(r => r.gate && r.gate.radius > 4);
  app().innerHTML = `
  <p class="kicker">Pore · the ITPR3 state panel</p>
  <h1 class="glow">Only the activated structures are open</h1>
  <p class="lead">Each deposit of human ITPR3 that ip3r_genes chose for a gating state is measured the same way: the four-fold
  axis by superposing each subunit onto its neighbour, the pore radius along it, the narrowest luminal point (the selectivity
  filter) and the narrowest cytosolic point (the gate). Profiles are aligned on each deposit's filter; the lumen is to the left.</p>
  <div class="panel">
    <canvas class="plot tall" id="prof" role="img" aria-label="Pore radius against distance from the filter for every ITPR3 state"></canvas>
    ${legend(series)}
  </div>
  <div class="two">
    <div class="panel"><h3>Radius at the gate and the filter</h3>
      <canvas class="plot" id="bars" role="img" aria-label="Gate and filter radius per deposit"></canvas>
      <div class="legend"><span><i style="background:#ffb547"></i>gate</span><span><i style="background:#5fa8ff"></i>filter</span></div></div>
    <div class="panel scroll"><table>
      <thead><tr><th>PDB</th><th>State</th><th>Res. (Å)</th><th>IP3</th><th>Gate r (Å)</th><th>Filter r (Å)</th><th>Gate lining</th></tr></thead>
      <tbody>${rows.map(r => `<tr class="click" data-pdb="${r.pdb_id}"><td><a href="#/viewer/${r.pdb_id}">${r.pdb_id}</a></td><td>${esc(r.state)}${r.control ? ' <span class="tag b">control</span>' : ''}</td>
        <td>${r.resolution}</td><td>${r.ip3_bound ? 'yes' : 'no'}</td><td class="num">${r.gate ? r.gate.radius.toFixed(2) : '–'}</td>
        <td class="num">${r.filter ? r.filter.radius.toFixed(2) : '–'}</td><td class="small">${r.gate ? esc(r.gate.residues.join(', ')) : ''}</td></tr>`).join('')}</tbody>
    </table></div>
  </div>
  <div class="callout a">Open (gate radius above 4 Å): ${open.map(r => `<b>${r.pdb_id}</b> (${esc(r.state)}, ${r.gate.radius.toFixed(2)} Å)`).join(' and ')}.
  Every other state's gate is between ${Math.min(...rows.filter(r => !open.includes(r)).map(r => r.gate.radius)).toFixed(2)} and
  ${Math.max(...rows.filter(r => !open.includes(r)).map(r => r.gate.radius)).toFixed(2)} Å, while the filter barely changes.</div>
  <p class="small mute">${esc(s.definition)}. Computed by <code>python -m ip3r states</code> (and <code>measure_channel</code> for 7T3T)
  and served precomputed; the README's state table and Figure 2 come from the same code.</p>`;
  const draw = () => {
    plot(document.getElementById('prof'), { series, xlabel: 'distance from the filter along the axis (Å) · lumen ← → cytosol',
      ylabel: 'pore radius r_min (Å)', ylim: [0, 14] });
    const xs = rows.map((_, i) => i);
    plot(document.getElementById('bars'), { series: [
      { x: xs, y: rows.map(r => r.gate?.radius), color: '#ffb547', points: true, r: 5 },
      { x: xs, y: rows.map(r => r.filter?.radius), color: '#5fa8ff', points: true, r: 5 }],
      xticks: rows.map((r, i) => ({ v: i, label: r.pdb_id })), xlabel: 'deposit, narrowest gate first', ylabel: 'radius (Å)', ylim: [0, 7], xlim: [-0.5, rows.length - 0.5] });
  };
  draw();
  onResize(document.getElementById('prof'), draw);
});
