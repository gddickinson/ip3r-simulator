// Subtle animated background: slow drifting Ca2+ "puffs" (glowing dots). Off under prefers-reduced-motion.
export function startBackground(canvas) {
  if (!canvas || matchMedia('(prefers-reduced-motion: reduce)').matches) return;
  const ctx = canvas.getContext('2d');
  const dots = Array.from({ length: 26 }, (_, i) => ({ x: Math.random(), y: Math.random(), r: 30 + Math.random() * 90,
    p: Math.random() * 6.28, s: 0.0002 + Math.random() * 0.0004, h: i % 3 }));
  const hues = ['255,95,162', '255,181,71', '95,227,255'];
  const size = () => { canvas.width = innerWidth; canvas.height = innerHeight; };
  size(); addEventListener('resize', size);
  const tick = t => {
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    for (const d of dots) {
      const a = 0.05 * (0.5 + 0.5 * Math.sin(d.p + t * d.s * 6));
      const x = (d.x + Math.sin(t * d.s + d.p) * 0.03) * canvas.width, y = (d.y + Math.cos(t * d.s + d.p) * 0.03) * canvas.height;
      const g = ctx.createRadialGradient(x, y, 0, x, y, d.r);
      g.addColorStop(0, `rgba(${hues[d.h]},${a})`); g.addColorStop(1, `rgba(${hues[d.h]},0)`);
      ctx.fillStyle = g; ctx.fillRect(x - d.r, y - d.r, 2 * d.r, 2 * d.r);
    }
    requestAnimationFrame(tick);
  };
  requestAnimationFrame(tick);
}
