// Static-site data layer: answers api() calls from the JSON files exported by scripts/build_static.py.
// Every GET the pages use (/info, /deposits, /palette, /structure/<PDB>, /modes/<PDB>, /transition, /states, /gating,
// /anm-parameters) was produced by the Python server at build time. Nothing needs Python at request time.
const files = new Map();
async function file(rel) {
  if (!files.has(rel)) files.set(rel, fetch('api/' + rel + '.json').then(r => { if (!r.ok) throw new Error(`${r.status} ${rel}`); return r.json(); }));
  return files.get(rel);
}

export async function staticApi(path, body) {
  if (body) throw new Error('This needs the local app (Python server).');
  return file(path.split('?')[0].replace(/^\/+/, ''));
}
