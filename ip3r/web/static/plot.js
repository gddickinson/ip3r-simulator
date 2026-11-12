// Small canvas line plot (no library): series of (x, y), optional log x, axis ticks, legend, marker lines.
// plot(canvas, { series: [{ x, y, color, width, dash, label, points }], xlabel, ylabel, xlog, xlim, ylim, vlines, hlines })

const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim() || '#94a3c6';

function niceTicks(lo, hi, count = 5) {
  const span = hi - lo || 1, step0 = span / count, mag = 10 ** Math.floor(Math.log10(step0));
  const step = [1, 2, 5, 10].map(m => m * mag).find(s => span / s <= count) || 10 * mag;
  const out = []; for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9 * span; v += step) out.push(+v.toPrecision(10));
  return out;
}
const fmt = v => Math.abs(v) >= 1e4 || (Math.abs(v) < 1e-2 && v !== 0) ? v.toExponential(0) : String(+v.toPrecision(4));

export function plot(canvas, o) {
  const dpr = Math.min(devicePixelRatio || 1, 2), W = canvas.clientWidth || 600, H = canvas.clientHeight || 260;
  canvas.width = Math.round(W * dpr); canvas.height = Math.round(H * dpr);
  const ctx = canvas.getContext('2d'); ctx.setTransform(dpr, 0, 0, dpr, 0, 0); ctx.clearRect(0, 0, W, H);
  const m = { l: 56, r: 12, t: 12, b: 40 }, pw = W - m.l - m.r, ph = H - m.t - m.b;
  const tx = v => o.xlog ? Math.log10(v) : v;
  const xs = o.series.flatMap(s => s.x.map(tx)).filter(Number.isFinite), ys = o.series.flatMap(s => s.y).filter(Number.isFinite);
  const [x0, x1] = o.xlim ? o.xlim.map(tx) : [Math.min(...xs), Math.max(...xs)];
  const [y0, y1] = o.ylim || [Math.min(0, ...ys), Math.max(...ys) * 1.05 || 1];
  const X = v => m.l + (tx(v) - x0) / (x1 - x0 || 1) * pw, Y = v => m.t + ph - (v - y0) / (y1 - y0 || 1) * ph;
  ctx.font = '11px ui-monospace, Menlo, monospace'; ctx.fillStyle = css('--mute'); ctx.strokeStyle = 'rgba(148,163,198,.18)'; ctx.lineWidth = 1;
  const xt = o.xlog ? niceTicks(x0, x1, 6).filter(Number.isInteger).map(e => 10 ** e) : niceTicks(x0, x1, Math.max(3, Math.floor(pw / 80)));
  ctx.textAlign = 'center';
  const xticks = o.xticks || xt.map(v => ({ v, label: fmt(v) }));
  for (const { v, label } of xticks) { const x = X(v); ctx.beginPath(); ctx.moveTo(x, m.t); ctx.lineTo(x, m.t + ph); ctx.stroke(); ctx.fillText(label, x, m.t + ph + 14); }
  ctx.textAlign = 'right';
  for (const v of niceTicks(y0, y1, 5)) { const y = Y(v); ctx.beginPath(); ctx.moveTo(m.l, y); ctx.lineTo(m.l + pw, y); ctx.stroke(); ctx.fillText(fmt(v), m.l - 6, y + 4); }
  ctx.textAlign = 'center'; ctx.fillText(o.xlabel || '', m.l + pw / 2, H - 6);
  ctx.save(); ctx.translate(12, m.t + ph / 2); ctx.rotate(-Math.PI / 2); ctx.fillText(o.ylabel || '', 0, 0); ctx.restore();
  ctx.save(); ctx.beginPath(); ctx.rect(m.l, m.t, pw, ph); ctx.clip();
  for (const l of o.vlines || []) { ctx.strokeStyle = l.color; ctx.setLineDash(l.dash || [4, 4]); ctx.beginPath(); ctx.moveTo(X(l.x), m.t); ctx.lineTo(X(l.x), m.t + ph); ctx.stroke(); }
  for (const l of o.hlines || []) { ctx.strokeStyle = l.color; ctx.setLineDash(l.dash || [4, 4]); ctx.beginPath(); ctx.moveTo(m.l, Y(l.y)); ctx.lineTo(m.l + pw, Y(l.y)); ctx.stroke(); }
  for (const s of o.series) {
    ctx.strokeStyle = s.color; ctx.fillStyle = s.color; ctx.lineWidth = s.width || 1.6; ctx.setLineDash(s.dash || []);
    if (s.points) { s.x.forEach((x, i) => { if (Number.isFinite(s.y[i])) { ctx.beginPath(); ctx.arc(X(x), Y(s.y[i]), s.r || 3, 0, 7); ctx.fill(); } }); continue; }
    ctx.beginPath(); let pen = false;
    s.x.forEach((x, i) => { const y = s.y[i]; if (!Number.isFinite(y) || !Number.isFinite(tx(x))) { pen = false; return; } pen ? ctx.lineTo(X(x), Y(y)) : ctx.moveTo(X(x), Y(y)); pen = true; });
    ctx.stroke();
  }
  ctx.restore(); ctx.setLineDash([]);
  ctx.strokeStyle = 'rgba(148,163,198,.5)'; ctx.strokeRect(m.l, m.t, pw, ph);
  return { X, Y };
}

// HTML legend for series with labels.
export function legend(series) {
  return `<div class="legend">${series.filter(s => s.label).map(s =>
    `<span><i style="background:${s.color}${s.dash ? ';opacity:.6' : ''}"></i>${s.label}</span>`).join('')}</div>`;
}

// Redraw when the canvas is resized; stops (and never redraws) once the page has replaced it.
export function onResize(canvas, draw) {
  const ro = new ResizeObserver(() => { if (canvas.isConnected) draw(); else ro.disconnect(); });
  ro.observe(canvas);
}
