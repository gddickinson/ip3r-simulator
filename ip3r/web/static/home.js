// Home page: hero with the IP3-bound receptor turning in WebGL, counts from the API, entry cards.
import { route, api, app, esc, onLeave } from './core.js';
import { createViewer, webglAvailable } from './gl_view.js';
import { beads, colours, model, contactBeads } from './mol.js';

const DOORS = [
  ['#/viewer/6DQN', 'Structure', 'Eight cryo-EM deposits of human ITPR3 as C-alpha models: colour by functional region, subunit, conservation or B-factor; the IP3 contacts in gold.'],
  ['#/pore', 'Pore', 'The pore radius along the four-fold axis for every ITPR3 gating state: only activated 8TKF and the 7T3T control are open at the gate.'],
  ['#/transition', 'Transition', 'Resting 8TKG morphing into activated 8TKF, coloured by how far each residue moves, with the gate opening along the path.'],
  ['#/modes', 'Modes', 'Elastic-network modes labelled by C4 symmetry: the full network from the Python code, and a reduced one solved live in your browser.'],
  ['#/gating', 'Gating', 'Three published models of how Ca²⁺ and IP3 open the channel: the bell-shaped Ca²⁺ response, computed in your browser.'],
  ['#/about', 'About', 'What the simulator is, what this web version shows and leaves out, the PDB entries and papers it uses, and how to run it.'],
];

route('home', async () => {
  const [info, dep, palette, rec] = await Promise.all([api('/info'), api('/deposits'), api('/palette'), api('/structure/6DQN')]);
  const open = dep.deposits.filter(d => d.constrictions.gate && d.constrictions.gate.radius > 4).map(d => d.pdb_id);
  app().innerHTML = `
  <section class="hero two">
    <div>
      <p class="kicker">Ca²⁺ release channel · cryo-EM · elastic network · gating</p>
      <h1 class="glow">IP3 Receptor Structural Simulator</h1>
      <p class="lead">${esc(info.description)}</p>
      <p class="small mute">The desktop app (PyQt6 + OpenGL) loads every deposit atom by atom and runs the permeation, puff and spark models;
      this web version shows a polished subset: the structures as C-alpha models, the pore of each gating state, the resting-to-open
      morph, the elastic-network modes and the three gating models.</p>
      <a class="btn" href="#/viewer/6DQN">Open the structure</a><a class="btn alt" href="#/transition">Watch it open</a>
    </div>
    <figure class="fig"><div class="viewer hero-view"><canvas id="hero" role="img" aria-label="Human ITPR3 with IP3 bound (PDB 6DQN) turning about its four-fold axis"></canvas></div>
      <figcaption class="small mute">PDB 6DQN, human ITPR3 with IP3 bound, coloured by functional region; the IP3 contacts and IP3 in gold.
      Cytosol at the top, the membrane pore at the bottom.</figcaption></figure>
  </section>
  <div class="tiles">
    <div class="tile"><b>${info.deposits_shown}</b><span>deposits in the web viewer</span></div>
    <div class="tile"><b>${info.registry.ip3r + info.registry.ryr}</b><span>deposits in the desktop app (${info.registry.ip3r} IP3R, ${info.registry.ryr} RyR1)</span></div>
    <div class="tile"><b>${open.length}</b><span>open at the gate (${open.join(', ')})</span></div>
    <div class="tile"><b>${info.parameters}</b><span>registered, cited parameters</span></div>
    <div class="tile"><b>${info.references}</b><span>sources behind them</span></div>
  </div>
  <div class="grid wide">${DOORS.map(([href, t, d]) => `<a class="card door" href="${href}"><b>${t}</b><span>${d}</span></a>`).join('')}</div>`;
  const canvas = document.getElementById('hero');
  if (!webglAvailable()) { canvas.parentElement.innerHTML = '<p class="mute" style="padding:20px">The 3-D view needs WebGL.</p>'; return; }
  const b = beads(rec), v = createViewer(canvas);
  v.setModel(model(b, colours(b, 'element', palette), { highlight: contactBeads(rec, b) }));
  v.spin(true);
  onLeave(() => v.dispose());
});
