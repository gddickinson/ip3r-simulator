"""Command line: ``python -m ip3r <command>``.

Everything the GUI computes is available headless, which is what makes it
testable and scriptable:

    python -m ip3r                  # the GUI (same as `gui`)
    python -m ip3r fetch            # download every registry structure
    python -m ip3r info 6DQN        # measure a deposit (axis, pore, IP3 sites)
    python -m ip3r checks           # re-derive the ip3r_genes findings
    python -m ip3r states           # pore of every ITPR3 gating state (--paralog RYR1)
    python -m ip3r unitary          # K+ conductance of each state (--paralog RYR1)
    python -m ip3r selectivity 8TKF # P_Cl:P_K, P_Ca:P_K, i_Ca vs Vais 2010
    python -m ip3r mutants | ryr-gating | sparks | spark-termination  # RyR1
    python -m ip3r modes 6DQN       # elastic-network modes with C4 irreps
    python -m ip3r transition 8TKG 8TKF   # displacement, morph, mode overlap
    python -m ip3r gating           # the bell curve at several IP3 levels (--model mak|pd)
    python -m ip3r oscillate --ip3 0.5 [--window]
    python -m ip3r puffs --ip3 0.2 [--model park-drive] [--scan]
    python -m ip3r microdomain [--scan ah42|n|store]   # puffs read from fluo-4
    python -m ip3r graft 8TKG [--mode full] [--calibrate] [--long]   # AlphaFold fills (cli_graft)
    python -m ip3r params           # every registered number and its source
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

__all__ = ["main"]

_ICON = {"confirmed": "✓", "discrepancy": "✗", "not_run": "·", "error": "!"}


def _info(args) -> int:
    from .io import loader
    from .structure.channel import measure_channel
    loader.ALLOW_FETCH = args.fetch
    s = measure_channel(loader.load(args.pdb))
    num = s.numbering
    print(f"{s.name}: subunit rotation {s.superposition_angle:.2f}°, C4 residual "
          f"{s.c4_residual:.3f} Å, axes differ by {s.axis_disagreement_deg:.3f}°")
    print(f"  numbering: {num.paralog + f' ({num.identity:.1%})' if num else 'no human paralog'}"
          + (f"; mismatch segments {list(num.mismatch_segments)}" if num and num.mismatch_segments else ""))
    print(f"  pore-domain span z {s.span[0]:.1f}..{s.span[1]:.1f} Å ({s.span_source})")
    for c in s.constrictions.values():
        print(f"  {c.name:6s} r = {c.radius:.2f} Å at z = {c.z:+.1f} Å: {', '.join(c.residues)}")
    for k, v in s.ip3_contacts.items():
        print(f"  IP3 on {k}: {', '.join(f'{a}{b}' for a, b in v['same_subunit'])}"
              + (f"; cross-subunit {v['other_subunit']}" if v["other_subunit"] else ""))
    return 0


def _checks(args) -> int:
    from .analysis.checks import run_checks
    from .io import loader
    loader.ALLOW_FETCH = args.fetch
    results = run_checks(ids=set(args.id) if args.id else None, paper=args.paper,
                         progress=lambda i, n, c: print(f"  [{i + 1}/{n}] {c.id}",
                                                        file=sys.stderr, end="\r"))
    print(" " * 60, file=sys.stderr, end="\r")
    counts: dict[str, int] = {}
    for r in results:
        counts[r.status] = counts.get(r.status, 0) + 1
        print(f"{_ICON[r.status]} {r.status:11s} {r.check.id:32s} [{r.check.kind}]")
        if args.verbose or r.status != "confirmed":
            print(f"      published: {r.outcome.published}")
            print(f"      found:     {r.outcome.found or r.outcome.detail}")
    print("\n" + ", ".join(f"{v} {k}" for k, v in sorted(counts.items())))
    if args.json:
        Path(args.json).write_text(json.dumps([{
            "id": r.check.id, "paper": r.check.paper, "kind": r.check.kind,
            "claim": r.check.claim, "status": r.status,
            "published": r.outcome.published, "found": r.outcome.found,
            "detail": r.outcome.detail} for r in results], indent=1))
    if args.figures:
        _figures(results, Path(args.figures))
    if counts.get("error"):
        return 2
    return 1 if args.strict and counts.get("discrepancy") else 0


def _figures(results, out: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from .analysis import exhibits
    out.mkdir(parents=True, exist_ok=True)
    for r in results:
        fig, ax = plt.subplots(figsize=(5, 3.2))
        if exhibits.draw(ax, r.check.id, r.outcome):
            fig.tight_layout()
            fig.savefig(out / f"{r.check.id}.png", dpi=150)
        plt.close(fig)


def _modes(args) -> int:
    from .io import loader
    from .physics.anm import ANM, tetramer_sites
    from .structure.symmetry import tetramer_frame
    loader.ALLOW_FETCH = args.fetch
    st = loader.load(args.pdb)
    fr = tetramer_frame(st)
    from .physics.network_checks import describe_local, local_modes
    coords, res = tetramer_sites(st, fr)
    anm = ANM(coords, axis=fr.axis)
    ms = anm.label_symmetry(anm.calc_modes(args.n))
    print(f"{st.name}: {len(coords)} sites, {len(res)} per subunit")
    kappa = ms.collectivity()
    for i in range(ms.n_modes):
        print(f"  mode {i + 1:3d}  λ = {ms.eigenvalues[i]:.4e}  {ms.symmetry[i]:5s} "
              f"χ = {ms.character[i]:+.3f}  κ {kappa[i]:.2f}")
    local = local_modes(ms, res, tetramer_sites(st, fr, stride=1)[1])
    for line in describe_local(local):
        print(f"local network artefact: {line}")
    return 0


def _transition(args) -> int:
    from .io import loader
    from .physics.transition_modes import transition_overlap
    from .structure.morph import morph
    from .structure.transition import prepare_transition
    loader.ALLOW_FETCH = args.fetch
    tr = prepare_transition(loader.load(args.start), loader.load(args.end), args.fit)
    print(f"fit on {tr.fit}: {tr.meta['n_fit_sites']} sites, RMSD {tr.fit_rmsd:.2f} Å; "
          f"subunits {''.join(tr.chains_start)} -> {''.join(tr.chains_end)}"
          + ("" if tr.meta["correspondence_determined"] else
             " (all cyclic correspondences fit alike: C4-symmetric pair)"))
    print("mean residue displacement by element (Å):")
    for name, d in sorted(tr.element_means().items(), key=lambda kv: -kv[1]):
        print(f"  {name:20s} {d:6.2f}")
    m = morph(tr.start, tr.end, args.method)
    print(f"morph: {m.summary()}")
    if args.gate:
        _gate_table(loader.load(args.start), loader.load(args.end), tr, m)
    ov = transition_overlap(tr, args.reference, stride=args.stride, n_modes=args.n)
    print("\n".join(ov.report()))
    kappa = ov.modes.collectivity()
    for i in range(ov.modes.n_modes):
        print(f"  mode {i + 1:3d} {ov.modes.symmetry[i]:5s} κ {kappa[i]:.2f}  overlap {ov.overlap[i]:.3f}"
              f"  cumulative {ov.cumulative[i]:.3f}  (null {ov.null_cumulative[i]:.3f})")
    if args.cutoff_scan:
        _cutoff_table(tr, args)
    if args.stride_check:
        from .physics.network_checks import stride_agreement
        for s, v in stride_agreement(tr, reference=args.reference, n_modes=args.n).items():
            print(f"stride {s}: RMSIP with stride 1 over {ov.modes.n_modes} modes {v:.3f}")
    return 0


def _cutoff_table(tr, args) -> None:
    from .physics.network_checks import cutoff_scan
    print("cutoff  local  lowest A  overlap  collective A  cumulative  (null)")
    for r in cutoff_scan(tr, stride=args.stride, reference=args.reference, n_modes=args.n):
        a = "   -" if r.lowest_a is None else f"#{r.lowest_a + 1:3d}"
        print(f"{r.cutoff:6.1f}  {r.n_local:5d}  {a:>8s}  {r.lowest_a_overlap:7.3f}"
              f"  {r.collective_a:12.3f}  {r.cumulative:10.3f}  ({r.null:.3f})")


def _gate_table(st_start, st_end, tr, m) -> None:
    from .structure.morph_pore import atom_path, gate_path
    g = gate_path(atom_path(st_start, st_end, tr), tr, m)
    print(f"gate along the morph ({g.meta['n_atoms']:,} heavy atoms matched, "
          f"{g.meta['unmatched']} unmatched): half-way at fraction {g.half_open():.2f}, "
          f"overshoot {g.overshoot():.2f} Å")
    print(f"  {'t':>4s} {'gate':>6s} {'rigid':>6s} {'filter':>6s} {'chord':>6s}  lining")
    for i in range(len(g.fraction)):
        print(f"  {g.fraction[i]:4.2f} {g.gate[i]:6.2f} {g.rigid_gate[i]:6.2f} "
              f"{g.filter[i]:6.2f} {g.lining_error[i]:6.2f}  {' '.join(g.lining[i])}")


def _states(args) -> int:
    from .io import loader
    from .structure.states import state_panel
    loader.ALLOW_FETCH = args.fetch
    print(f"{'PDB':5s} {'state':24s} {'res':>5s} {'IP3':>4s} {'gate r':>7s} {'filter r':>8s}  gate lining")
    for r in state_panel(args.paralog):
        g = r.summary.constrictions.get("gate")
        print(f"{r.pdb_id:5s} {r.state:24s} {r.resolution:5.2f} {'yes' if r.ip3_bound else 'no':>4s} "
              f"{r.radius('gate'):7.2f} {r.radius('filter'):8.2f}  {', '.join(g.residues) if g else ''}")
    return 0


def _unitary(args) -> int:
    from .io import loader
    from .physics.unitary import published, unitary_panel
    loader.ALLOW_FETCH = args.fetch
    rows = unitary_panel(args.paralog, sweep=True)
    from .parameters import PARAMETERS as _P
    bath = rows[0].bath if rows and rows[0].bath else _P.value("permeation.bath_concentration")
    print(f"{args.paralog}: symmetric {bath * 1000:.0f} mM KCl")
    print("K+ conductance in symmetric KCl (series = closed form; neutral / "
          "charged = drift-diffusion without / with the lining side chains; "
          "paired = charged less salt-bridged groups)")
    for u in rows:
        print(u.row())
        if u.neutral.is_conducting:
            lo, hi = u.sweep["neutral"]
            clo, chi = u.sweep["charged"]
            plo, phi = u.sweep["paired"]
            print(f"      sweep (diffusivity x ion radius): neutral {lo:.0f}-{hi:.0f}, "
                  f"charged {clo:.0f}-{chi:.0f}, paired {plo:.0f}-{phi:.0f} pS; Debye "
                  f"{u.charged.meta['debye_length_A']:.1f} A; in-pore peak "
                  f"{u.charged.meta.get('peak_in_pore_M', float('nan')):.1f} M"
                  + (" (above the packing ceiling)"
                     if u.charged.meta.get("exceeds_packing_limit") else ""))
            print("      lining charges: " + ", ".join(
                f"{lab} x{n} at z {z:+.0f}" for lab, n, z in u.charge.residues()))
            if u.paired_charge.bridged:
                print("      salt-bridged, dropped when paired: " + ", ".join(
                    b.label() for b in u.paired_charge.bridged))
    print("measured: " + "; ".join(f"{k} {v:.0f} pS"
                                   for k, v in published(args.paralog).items()))
    return 0


def _gating(args) -> int:
    from .physics import gating, gating_mak, park_drive
    model = {"mak": gating_mak, "pd": park_drive}.get(args.model, gating)
    for p in args.ip3:
        b = model.bell_at(p)
        print(f"IP3 {p:6.3f} µM: peak P_open {b.po_peak:.4f} at Ca2+ "
              f"{b.c_peak:.3f} µM; half-activation {b.c_half_act:.3f}, "
              f"half-inhibition {b.c_half_inh:.3f} µM")
    for name, (a, i) in gating_mak.compare_flanks().items():
        print(f"{name}: IP3 low -> high moves half-activation {a:.3f}x, "
              f"half-inhibition {i:.2f}x")
    return 0


def _oscillate(args) -> int:
    from .physics.calcium import oscillation_metrics, oscillation_window, simulate
    if args.window:
        lo, hi, _ = oscillation_window()
        print(f"Ca2+ oscillates for IP3 in [{lo:.2f}, {hi:.2f}] µM (grid 0.01)")
        return 0
    m = oscillation_metrics(simulate(args.ip3, t_end=args.t_end))
    print(json.dumps(m, indent=1))
    return 0


def _puffs(args) -> int:
    from .physics import puff_compare as pc
    if args.scan:
        rows = pc.coupling_scan(args.ip3, args.duration, args.seed)
        print(f"IP3 {args.ip3} µM, {args.duration:g} s, seed {args.seed}; "
              f"large = peak ≥ {rows[pc.MODELS[0]][0]['large_at']} channels")
        print(f"{'model':11s} {'µM/open':>8s} {'Fano':>5s} {'open':>6s} "
              f"{'blips':>6s} {'multi':>6s} {'large':>6s} {'/s':>5s}")
        for model, scan in rows.items():
            for r in scan:
                print(f"{model:11s} {r['coupling']:8.3f} {r['fano']:5.2f} "
                      f"{r['open_fraction']:6.3f} {r['blips']:6d} {r['multi']:6d} "
                      f"{r['large']:6d} {r['large_per_s']:5.2f}")
        return 0
    pp = pc.params_for(args.model, args.coupling)
    print(json.dumps(pc.coupling_effect(args.ip3, args.duration, args.seed, pp,
                                        model=args.model), indent=1))
    return 0


def _params(args) -> int:
    from .parameters import PARAMETERS
    for r in PARAMETERS.provenance_rows():
        print(f"{r['key']:28s} {r['value']:>10g} {r['unit']:10s} {r['kind']:10s} {r['citation']}")
    return 0


def _fetch(args) -> int:
    from .io.fetch import fetch_all
    from .io.predictions import fetch_predictions
    status = fetch_all(args.pdb or None)
    bad = {k: v for k, v in status.items() if v != "ok"}
    print(f"{len(status) - len(bad)} of {len(status)} structures present")
    if not args.pdb:
        models = fetch_predictions()
        failed = [k for k, v in models.items() if v.startswith("failed")]
        print(f"{len(models) - len(failed)} of {len(models)} AlphaFold models present")
        bad = bad or failed
    return 1 if bad else 0


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if not argv or argv[0] == "gui" or argv[0].startswith("--"):
        from .ui.app import main as gui
        return gui(argv[1:] if argv and argv[0] == "gui" else argv)
    ap = argparse.ArgumentParser(prog="python -m ip3r")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("fetch")
    p.add_argument("pdb", nargs="*")
    p.set_defaults(fn=_fetch)
    p = sub.add_parser("info")
    p.add_argument("pdb")
    p.add_argument("--fetch", action="store_true")
    p.set_defaults(fn=_info)
    p = sub.add_parser("checks")
    p.add_argument("--paper")
    p.add_argument("--id", action="append")
    p.add_argument("--fetch", action="store_true", help="download missing structures")
    p.add_argument("--json")
    p.add_argument("--figures", help="write exhibit PNGs to this directory")
    p.add_argument("--strict", action="store_true", help="exit 1 on any discrepancy")
    p.add_argument("-v", "--verbose", action="store_true")
    p.set_defaults(fn=_checks)
    p = sub.add_parser("states")
    p.add_argument("--paralog", default="ITPR3")
    p.add_argument("--fetch", action="store_true")
    p.set_defaults(fn=_states)
    p = sub.add_parser("unitary", help="K+ conductance of every state from "
                       "its pore profile")
    p.add_argument("--paralog", default="ITPR3")
    p.add_argument("--fetch", action="store_true")
    p.set_defaults(fn=_unitary)
    from .cli_perm import register as _register_perm
    _register_perm(sub)
    from .cli_ryr import register as _register_ryr
    _register_ryr(sub)
    from .cli_domain import register as _register_domain
    _register_domain(sub)
    from .cli_shortfall import register as _register_shortfall
    _register_shortfall(sub)
    from .cli_graft import register as _register_graft
    _register_graft(sub)
    p = sub.add_parser("modes")
    p.add_argument("pdb")
    p.add_argument("-n", type=int, default=None)
    p.add_argument("--fetch", action="store_true")
    p.set_defaults(fn=_modes)
    p = sub.add_parser("transition", help="two states of one paralog: "
                       "displacement, morph and ANM overlap")
    p.add_argument("start", nargs="?", default="8TKG")
    p.add_argument("end", nargs="?", default="8TKF")
    p.add_argument("--fit", choices=("pore", "global"), default="pore")
    p.add_argument("--method", choices=("restrained", "linear"), default="restrained")
    p.add_argument("--reference", choices=("start", "end"), default="start",
                   help="whose elastic network is solved")
    p.add_argument("--stride", type=int, default=None)
    p.add_argument("-n", type=int, default=None)
    p.add_argument("--gate", action="store_true",
                   help="also measure the gate on every frame (atoms interpolated)")
    p.add_argument("--cutoff-scan", action="store_true",
                   help="re-solve the network over the registered cutoff grid")
    p.add_argument("--stride-check", action="store_true",
                   help="RMSIP of strides 2-4 against stride 1 (slow: ~10 s)")
    p.add_argument("--fetch", action="store_true")
    p.set_defaults(fn=_transition)
    p = sub.add_parser("gating")
    p.add_argument("--ip3", type=float, nargs="+", default=[0.1, 0.3, 1.0, 10.0])
    p.add_argument("--model", choices=["dyk", "mak", "pd"], default="dyk")
    p.set_defaults(fn=_gating)
    p = sub.add_parser("oscillate")
    p.add_argument("--ip3", type=float, default=0.5)
    p.add_argument("--t-end", type=float, default=200.0)
    p.add_argument("--window", action="store_true")
    p.set_defaults(fn=_oscillate)
    p = sub.add_parser("puffs")
    p.add_argument("--ip3", type=float, default=0.2)
    p.add_argument("--coupling", type=float, default=None)
    p.add_argument("--model", choices=["dyk", "park-drive"], default="dyk")
    p.add_argument("--scan", action="store_true",
                   help="both receptors over the registered coupling scan")
    p.add_argument("--duration", type=float, default=20.0)
    p.add_argument("--seed", type=int, default=0)
    p.set_defaults(fn=_puffs)
    p = sub.add_parser("params")
    p.set_defaults(fn=_params)
    args = ap.parse_args(argv)
    return args.fn(args)
