// node ip3r/web/static/gating_engine.test.mjs — the JS gating models against the repo's Python
// (tests/fixtures/gating_curves.json, written by scripts/make_anm_fixture.py with the registered constants).
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import assert from 'node:assert/strict';
import { MODEL_FN, bellAt, pdParkFraction } from './gating_engine.js';

const here = dirname(fileURLToPath(import.meta.url));
const fx = JSON.parse(readFileSync(join(here, '..', '..', '..', 'tests', 'fixtures', 'gating_curves.json'), 'utf8'));
let checks = 0;
const close = (a, b, rtol, atol, what) => { assert.ok(Math.abs(a - b) <= atol + rtol * Math.abs(b), `${what}: ${a} vs ${b}`); checks++; };

for (const [model, byIp3] of Object.entries(fx.curves)) {
  const g = fx.constants[model];
  for (const [p, po] of Object.entries(byIp3)) {
    fx.ca.forEach((c, i) => close(MODEL_FN[model](c, +p, g), po[i], 1e-9, 1e-15, `${model} P_open(c=${c}, p=${p})`));
  }
  for (const ref of fx.bells[model]) {
    const b = bellAt(model, ref.ip3, g);
    for (const k of ['c_peak', 'c_half_act']) close(b[k], ref[k], 1e-5, 0, `${model} ${k} at ${ref.ip3}`);
    close(b.po_peak, ref.po_peak, 1e-8, 0, `${model} po_peak at ${ref.ip3}`);
    if (ref.c_half_inh === null) assert.equal(b.c_half_inh, null); else close(b.c_half_inh, ref.c_half_inh, 1e-5, 0, `${model} c_half_inh`);
  }
}
fx.ca.forEach((c, i) => close(pdParkFraction(c, 1.0, fx.constants.pd), fx.park_fraction_ip3_1[i], 1e-9, 1e-15, 'park fraction'));
console.log(`gating_engine: ${checks} checks passed`);
