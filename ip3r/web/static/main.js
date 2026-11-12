// Entry point: registers the lazily loaded pages, sets the footer mode text and starts the router.
import { start, route, lazy, app, isStatic, $ } from './core.js';
import { startBackground } from './bg.js';

const PAGES = { 'home.js': ['home'], 'viewer.js': ['viewer'], 'pore.js': ['pore'], 'transition.js': ['transition'],
  'modes.js': ['modes'], 'gating.js': ['gating'], 'about.js': ['about'] };
lazy(Object.fromEntries(Object.entries(PAGES).flatMap(([file, names]) => names.map(n => [n, () => import('./' + file)]))));

route('notfound', () => { app().innerHTML = '<h1>Not found</h1><p><a href="#/">Back to the start</a></p>'; });

$('#mode').textContent = isStatic()
  ? 'static site: viewer, reduced elastic network and gating models run in your browser; structure measurements precomputed by the Python code'
  : 'local app: Python server on this machine';
startBackground($('#bg'));
start();
