// Web Worker: runs the reduced elastic network (anm_engine.reducedANM) off the page's main thread.
import { reducedANM } from './anm_engine.js';

self.onmessage = e => {
  const { id, record, params } = e.data;
  try {
    const t0 = performance.now();
    const r = reducedANM(record, params);
    r.ms = performance.now() - t0;
    self.postMessage({ id, result: r });
  } catch (err) {
    self.postMessage({ id, error: String(err && err.message || err) });
  }
};
