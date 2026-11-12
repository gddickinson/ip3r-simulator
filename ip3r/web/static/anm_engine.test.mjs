// node ip3r/web/static/anm_engine.test.mjs — the JS elastic network against the repo's Python ANM on the same sites
// (tests/fixtures/anm_reduced.json, written by scripts/make_anm_fixture.py from ip3r/web/precomputed/).
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import assert from 'node:assert/strict';
import { reducedANM, applyGenerator, jacobiEigen, cholesky } from './anm_engine.js';

const here = dirname(fileURLToPath(import.meta.url));
const root = join(here, '..', '..', '..');
const fixture = JSON.parse(readFileSync(join(root, 'tests', 'fixtures', 'anm_reduced.json'), 'utf8'));
let checks = 0;
const close = (a, b, rtol, atol, what) => { assert.ok(Math.abs(a - b) <= atol + rtol * Math.abs(b), `${what}: ${a} vs ${b}`); checks++; };

// Small pieces by hand.
{
  const { values } = jacobiEigen([2, 1, 1, 2], 2);
  close(Math.min(...values), 1, 1e-12, 1e-12, 'jacobi'); close(Math.max(...values), 3, 1e-12, 1e-12, 'jacobi');
  const A = Float64Array.from([4, 2, 2, 3]); cholesky(A, 2);
  close(A[0], 2, 1e-12, 0, 'chol'); close(A[2], 1, 1e-12, 0, 'chol'); close(A[3], Math.SQRT2, 1e-12, 0, 'chol');
  // Four times the generator is the identity; once moves block 0 onto block 1 turned by 90 degrees.
  const f = Float64Array.from({ length: 24 }, (_, i) => Math.sin(i + 1));
  let g = f; for (let p = 0; p < 4; p++) g = applyGenerator(g, 4);
  f.forEach((v, i) => close(g[i], v, 1e-12, 1e-12, 'S^4'));
  const once = applyGenerator(f, 4);
  close(once[6], -f[1], 1e-12, 1e-12, 'S x'); close(once[7], f[0], 1e-12, 1e-12, 'S y');
}

for (const c of fixture.cases) {
  const record = JSON.parse(readFileSync(join(root, 'ip3r', 'web', 'precomputed', `structure_${c.pdb_id}.json`), 'utf8'));
  const t0 = Date.now();
  const r = reducedANM(record, { stride: c.stride, cutoff: c.cutoff, gamma: c.gamma, d0: c.d0, nModes: c.eigenvalues.length });
  const ms = Date.now() - t0;
  assert.equal(r.nSites, c.n_sites); assert.deepEqual(r.residues, c.residues); assert.equal(r.nComponents, c.n_components);
  assert.ok(r.converged, 'Lanczos converged');
  c.eigenvalues.forEach((v, i) => close(r.eigenvalues[i], v, 1e-6, 1e-12, `${c.pdb_id} lambda ${i + 1}`));
  c.symmetry.forEach((s, i) => { assert.equal(r.symmetry[i], s, `${c.pdb_id} mode ${i + 1} irrep`); checks++; });
  c.character.forEach((v, i) => close(r.character[i], v, 0, 2e-3, `${c.pdb_id} chi ${i + 1}`));
  c.collectivity.forEach((v, i) => close(r.collectivity[i], v, 0, 2e-3, `${c.pdb_id} kappa ${i + 1}`));
  // Each Python mode lies in the span of the JS modes with (nearly) its eigenvalue: E pairs are only defined as a pair.
  c.vectors.forEach((pv, i) => {
    const nrm = Math.hypot(...pv);
    let proj = 0;
    r.vectors.forEach((jv, j) => {
      if (Math.abs(r.eigenvalues[j] - c.eigenvalues[i]) > 1e-3 * c.eigenvalues[i]) return;
      let d = 0; for (let k = 0; k < pv.length; k++) d += pv[k] * jv[k];
      proj += (d / nrm) ** 2;
    });
    close(Math.sqrt(proj), 1, 0, 1e-5, `${c.pdb_id} mode ${i + 1} subspace`);
  });
  console.log(`${c.pdb_id} stride ${c.stride} cutoff ${c.cutoff}: ${r.nSites} sites, ${r.iterations} Lanczos steps, ${ms} ms`);
}
console.log(`anm_engine: ${checks} checks passed`);
