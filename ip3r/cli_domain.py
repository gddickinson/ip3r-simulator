"""The puff microdomain (Round 7.2), ``register(sub)``:

    python -m ip3r microdomain --ip3 0.1 -n 10       # puffs read from F/F0
    python -m ip3r microdomain --clamp none          # the store free
    python -m ip3r microdomain --scan ah42           # IPI shape vs h42 recovery
    python -m ip3r microdomain --scan n [--scale S]  # amplitude vs cluster size
    python -m ip3r microdomain --scan store          # the sustained-open question
"""

from __future__ import annotations

__all__ = ["register"]

_CLAMP_WORDS = {"store": "store clamped", "none": "store free",
                "bath": "store and cytosol clamped"}


def _stats_line(s) -> str:
    ipi = f"{s.ipis.mean():.2f}" if len(s.ipis) else "-"
    return (f"puffs {s.n_puffs:4d} ({s.rate:.2f}/s)  blip dF/F0 {s.blip_mean:.2f}  "
            f"amp {s.amp_mean:.2f} ({s.amp_in_blips:.1f} blips, {s.ca_mean:.2f} uM, "
            f"{s.open_mean:.1f} open)  IPI {ipi} s  CV {s.cv:.2f}  "
            f"lam {s.lam:.2f} xi {s.xi:.2f}  LR {s.refractory_lr:.1f}")


def _microdomain(args) -> int:
    from .parameters import PARAMETERS as _P
    from .physics import puff_domain_scans as scans
    from .physics.puff_stats import puff_stats
    from .physics.puffs_domain import DomainPuffParams, simulate_cluster_domain

    if args.scan == "ah42":
        print(f"IPI vs h42 recovery a_h42 (N = {args.n}, IP3 {args.ip3} uM); "
              "Cao 2013 Fig. 4: lam 0.1-0.5, xi 0.5-2.2 /s, CV 0.65-0.95")
        for a, s in scans.ah42_scan(args.ip3, args.n, duration=args.duration):
            print(f"a_h42 {a:6.3f}  " + _stats_line(s))
        return 0
    if args.scan == "n":
        scale = args.scale if args.scale is not None else _P.value("domain.blip_scale")
        print(f"Amplitude vs N (IP3 {args.ip3} uM, release x{scale:g}); Cao 2013 "
              "Figs. 8/S9: dF/F0 bends near N = 12, Ca2+ does not")
        for n, s in scans.n_scan(args.ip3, duration=args.duration, scale=scale):
            print(f"N {n:3d}  " + _stats_line(s))
        return 0
    if args.scan == "store":
        print(f"Open fraction at raised release (IP3 {args.ip3} uM, N 20); "
              "conditions: " + "; ".join(f"{k} = {v}" for k, v in
                                          scans.CONDITION_NOTES.items()))
        for r in scans.store_scan(args.ip3, duration=args.duration):
            print(f"x{r.scale:5.2f}  {r.coupling:.2f} uM/open  free store "
                  f"{r.store_start:.0f} -> min {r.store_min:.0f} uM")
            for c, o in r.runs.items():
                cyt = r.cytosol_max.get(c)
                extra = f"  cytosol max {cyt:.2f} uM" if cyt is not None else ""
                print(f"    {c:24s} open {100 * o.open_fraction:5.1f} %  "
                      f"active {100 * o.active:5.1f} %  Fano {o.fano:.2f}{extra}")
        return 0
    dp = DomainPuffParams(n_channels=args.n, clamp=args.clamp)
    if args.scale is not None:
        dp.domain = dp.domain.scaled(args.scale)
    duration = args.duration or _P.value("domain.ipi_duration")
    tr = simulate_cluster_domain(args.ip3, duration, args.seed, dp)
    print(f"park/drive cluster in a microdomain: N {args.n}, IP3 {args.ip3} uM, "
          f"{_CLAMP_WORDS[args.clamp]}, {duration:g} s")
    print(f"rest: c {tr.rest.c:.3f} uM, store {tr.rest.cs:.0f} uM; open fraction "
          f"{tr.n_open.mean() / args.n:.3f}; microdomain max {tr.ca.max():.2f} uM")
    print(_stats_line(puff_stats(tr)))
    return 0


def register(sub) -> None:
    p = sub.add_parser("microdomain", help="park/drive puffs in Cao's "
                       "microdomain with fluo-4 (IPIs, amplitudes, the store)")
    p.add_argument("--ip3", type=float, default=0.1)
    p.add_argument("-n", type=int, default=10, help="receptors in the cluster")
    p.add_argument("--clamp", choices=("store", "none", "bath"), default="store")
    p.add_argument("--scale", type=float, default=None,
                   help="multiply the per-receptor release k_ipr")
    p.add_argument("--scan", choices=("ah42", "n", "store"), default=None)
    p.add_argument("--duration", type=float, default=None,
                   help="s per run (default: domain.ipi_duration / store_duration)")
    p.add_argument("--seed", type=int, default=0)
    p.set_defaults(fn=_microdomain)
