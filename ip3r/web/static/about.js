// About page: what the simulator is, what the web version computes where, the PDB entries and papers, how to run it.
import { route, api, app, esc, isStatic } from './core.js';

// The PDB entries shown and the papers that describe them (README reference list, items 7-9).
const PAPERS = {
  paknejad2018: 'Paknejad N, Hite RK (2018). Structural basis for the regulation of inositol trisphosphate receptors by Ca²⁺ and IP3. Nat Struct Mol Biol 25:660–668. doi:10.1038/s41594-018-0089-6',
  schmitz2022: 'Schmitz EA, Takahashi H, Karakas E (2022). Structural basis for activation and gating of IP3 receptors. Nat Commun 13:1408. doi:10.1038/s41467-022-29073-2',
  paknejad2023: 'Paknejad N, Sapuru V, Hite RK (2023). Structural titration reveals Ca²⁺-dependent conformational landscape of the IP3 receptor. Nat Commun 14:6897. doi:10.1038/s41467-023-42707-3',
};
const PAPER_OF = { '6DQJ': 'paknejad2018', '6DQN': 'paknejad2018', '7T3P': 'schmitz2022', '7T3T': 'schmitz2022',
  '8TKF': 'paknejad2023', '8TKG': 'paknejad2023', '8TKH': 'paknejad2023', '8TLA': 'paknejad2023' };
const METHODS = [
  'Berman HM, et al. (2000). The Protein Data Bank. Nucleic Acids Res 28:235–242. doi:10.1093/nar/28.1.235',
  'Kabsch W (1976). A solution for the best rotation to relate two sets of vectors. Acta Cryst A32:922–923. doi:10.1107/S0567739476001873',
  'Atilgan AR, et al. (2001). Anisotropy of fluctuation dynamics of proteins with an elastic network model. Biophys J 80:505–515. doi:10.1016/S0006-3495(01)76033-X',
  'Tama F, Sanejouand YH (2001). Conformational change of proteins arising from normal mode calculations. Protein Eng 14:1–6. doi:10.1093/protein/14.1.1',
  'De Young GW, Keizer J (1992). PNAS 89:9895–9899. doi:10.1073/pnas.89.20.9895 · Li YX, Rinzel J (1994). J Theor Biol 166:461–473. doi:10.1006/jtbi.1994.1041',
  'Mak DOD, McBride S, Foskett JK (1998). PNAS 95:15821–15825. doi:10.1073/pnas.95.26.15821',
  'Siekmann I, et al. (2012). Biophys J 103:658–668. doi:10.1016/j.bpj.2012.07.016 · Cao P, et al. (2013). Biophys J 105:1133–1142. doi:10.1016/j.bpj.2013.07.034',
];
const doi = t => esc(t).replace(/doi:(\S+?)(?=\s|$|·)/g, (m, d) => `<a href="https://doi.org/${d}" target="_blank" rel="noopener">doi:${d}</a>`);

route('about', async () => {
  const [info, dep] = await Promise.all([api('/info'), api('/deposits')]);
  app().innerHTML = `
  <p class="kicker">About</p>
  <h1 class="glow">IP3R Structural Simulator</h1>
  <p class="lead">${esc(info.description)} It re-derives that project's published results from their original inputs and shows them on
  the receptor's structure. Every number a calculation uses is a registered, cited parameter (${info.parameters} of them).</p>
  <h2>What this web version shows</h2>
  <table><thead><tr><th>Page</th><th>Computed where</th></tr></thead><tbody>
    <tr><td><a href="#/viewer">Structure</a></td><td>C-alpha traces, annotation and channel measurements precomputed by <code>measure_channel</code>; drawn by a small WebGL viewer in your browser.</td></tr>
    <tr><td><a href="#/pore">Pore</a></td><td>Precomputed by <code>python -m ip3r states</code>.</td></tr>
    <tr><td><a href="#/transition">Transition</a></td><td>Precomputed by <code>python -m ip3r transition 8TKG 8TKF --gate</code> (morph frames kept every other one, interpolated between them for display).</td></tr>
    <tr><td><a href="#/modes">Modes</a></td><td>Full network precomputed by <code>python -m ip3r modes</code>; the reduced network solved live in JavaScript, tested against the Python ANM.</td></tr>
    <tr><td><a href="#/gating">Gating</a></td><td>Live in JavaScript from the registered constants, tested against the Python models.</td></tr>
  </tbody></table>
  <p>Left out of the web version: the atomic representations (cartoon, sticks, surfaces), the RyR1 control structures and the 17 further
  IP3R deposits, the findings checks against ip3r_genes (they read that project's tables), the permeation and selectivity models in 1-D and
  3-D, Ca²⁺ oscillations, stochastic puffs and sparks, AlphaFold fills, superposition, sessions and movies. All of them are in the desktop
  app. ${isStatic() ? 'This static site computes nothing in Python; the local web app serves the same pages from the Python code.' : ''}</p>
  <h2>Structures</h2>
  <p>Atomic coordinates from the wwPDB (CC0), reduced here to one bead per C-alpha and stored at 0.01 Å; each file's SHA-256 is recorded
  with it. Residue numbers are canonical human UniProt numbers (ITPR3, Q14573); annotation is painted only where the deposit is in that numbering.</p>
  <div class="scroll"><table><thead><tr><th>PDB</th><th>State</th><th>Resolution</th><th>Described in</th></tr></thead><tbody>
    ${dep.deposits.map(d => `<tr><td><a href="https://www.rcsb.org/structure/${d.pdb_id}" target="_blank" rel="noopener">${d.pdb_id}</a></td>
      <td>${esc(d.state)}</td><td>${d.resolution} Å</td><td class="small">${doi(PAPERS[PAPER_OF[d.pdb_id]] || '')}</td></tr>`).join('')}
  </tbody></table></div>
  <h2>Methods and models</h2>
  <ul class="small">${METHODS.map(m => `<li>${doi(m)}</li>`).join('')}</ul>
  <h2>Run it yourself</h2>
  <pre>git clone https://github.com/gddickinson/ip3r-simulator
cd ip3r-simulator
make env &amp;&amp; make fetch          # conda env ip3r_sim; ~86 MB of structures into ref/
python -m ip3r                  # the desktop app (loads 6DQN)

pip install -r requirements-web.txt
python run_web.py               # this web app at http://127.0.0.1:8194
python scripts/make_web_data.py # rebuild its precomputed results (needs ref/)</pre>
  <h2>Licence and credits</h2>
  <p>Code: MIT licence. By George Dickinson; part of the <a href="https://github.com/gddickinson/ip3r_genes" target="_blank" rel="noopener">ip3r_genes</a>
  project. Cite the repository through its CITATION.cff. The web interface was added with the rest of George's repositories; its numbers
  are the Python code's.</p>`;
});
