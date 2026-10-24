"""Command-line interface: ``ceatofoam <command> ...``.

Commands:
  design  INPUT            print the CEA data, mixture and design point
  case    INPUT CASE_DIR   write the OpenFOAM case (and ceatofoam.json)
  post    CASE_DIR         thrust and mass flow of a finished case
  sweep   ...              CEA chamber temperature and Isp over O/F
  interactive CASE_DIR     ask for the inputs on the console, then 'case'
"""

import argparse
import json
import sys

from . import __version__
from .units import to_si


def _print_design(d):
    t = d["target"]
    mix = d["mixture"]
    c = d["contour"]
    print(f"Chamber: pc = {d['pc']/1e5:.4g} bar, Tc = {d['tc']:.1f} K")
    print(f"Ambient: p = {d['p_amb']:.1f} Pa, T = {d['t_amb']:.2f} K")
    for name, pt in d["cea"].items():
        print(
            f"CEA {name:11s}: Isp = {pt.isp:7.2f} s, c* = {pt.cstar:7.1f}"
            f" m/s, gamma_s(c) = {pt.gammas_c:.4f}, M(c) = {pt.mw_c:.3f}"
        )
    print(
        f"Frozen mixture: M = {mix.mw:.4f} kg/kmol, gamma(Tc) ="
        f" {mix.gamma(d['tc']):.4f}, dropped species {mix.dropped}"
    )
    print(f"Design ({t.model}):")
    for key, val in t.as_dict().items():
        if key != "model":
            print(f"  {key:12s} {val:.6g}")
    print("Contour:")
    print(f"  cylinder length {c.cylinder_length:.6g} m, throat at"
          f" x = {c.x_throat:.6g} m, exit at x = {c.x_exit:.6g} m")
    print(f"  contraction ratio {c.contraction_ratio:.4g}, area ratio"
          f" {c.area_ratio:.4g}")


def cmd_design(args):
    from .build import load, make_design

    _print_design(make_design(load(args.input)))


def cmd_case(args):
    from .build import load, make_case

    rec = make_case(load(args.input), args.case_dir)
    print(
        f"Wrote {args.case_dir}: {rec['mesh']['cells']} cells, design"
        f" thrust {rec['design']['thrust']:.6g} N, mdot"
        f" {rec['design']['mdot']:.6g} kg/s. Run ./Allrun in the case."
    )


def cmd_post(args):
    from .post import write

    print(json.dumps(write(args.case_dir, args.window), indent=2))


def cmd_sweep(args):
    import numpy as np

    from .atmosphere import ambient
    from .cea import sweep_mixture_ratio

    pc = to_si(args.pressure, "pressure")
    _, p_amb = ambient(to_si(args.altitude, "length"))
    ratios = np.arange(args.of_min, args.of_max + 0.5*args.of_step,
                       args.of_step)
    pts = sweep_mixture_ratio(args.fuel, args.oxidizer, pc, pc/p_amb,
                              ratios, args.analysis)
    rows = [(of, p.tc, p.isp) for of, p in zip(ratios, pts) if p]
    print("O/F      Tc [K]     Isp [s]")
    for of, tc, isp in rows:
        print(f"{of:5.2f}  {tc:9.1f}  {isp:9.2f}")
    if rows:
        of_t = max(rows, key=lambda r: r[1])
        of_i = max(rows, key=lambda r: r[2])
        print(f"Peak chamber temperature {of_t[1]:.1f} K at O/F = {of_t[0]:.2f}")
        print(f"Peak Isp {of_i[2]:.2f} s at O/F = {of_i[0]:.2f}")
    if args.plot:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        of, tc, isp = (np.array(v) for v in zip(*rows))
        fig, ax = plt.subplots(1, 2, figsize=(9, 3.5))
        ax[0].plot(of, tc, "k-")
        ax[0].set_xlabel("O/F")
        ax[0].set_ylabel("Chamber temperature [K]")
        ax[1].plot(of, isp, "k-")
        ax[1].set_xlabel("O/F")
        ax[1].set_ylabel("Isp [s]")
        fig.tight_layout()
        fig.savefig(args.plot, dpi=200)
        print(f"Wrote {args.plot}")


def _ask(prompt, default=None):
    text = input(f"{prompt}" + (f" [{default}]" if default else "") + ": ")
    return text.strip() or default


def cmd_interactive(args):
    """The console dialogue of ceaToFoam 1.x, writing a case at the end."""
    from .build import make_case

    cfg = {
        "propellants": {
            "fuel": _ask("Fuel (RP-1, H2(L), CH4(L), CH6N2(L), ...)", "RP-1"),
            "oxidizer": _ask("Oxidizer (O2(L), N2O4(L), O2)", "O2(L)"),
        },
        "chamber": {
            "pressure": _ask("Chamber pressure (e.g. 500 psi)", "500 psi"),
            "mixture_ratio": float(_ask("O/F ratio", "2.4")),
        },
        "design": {
            "thrust": _ask("Design thrust (e.g. 500 lbf)", "500 lbf"),
            "altitude": _ask("Design altitude (e.g. 0 ft)", "0 ft"),
            "sizing": _ask("Sizing model (frozen, equilibrium,"
                           " constant-gamma)", "frozen"),
        },
        "geometry": {
            "characteristic_length": _ask("L* (e.g. 40 in)", "40 in"),
            "convergent_half_angle": _ask("Convergent half-angle [deg]",
                                          "30"),
            "divergent_half_angle": _ask("Divergent half-angle [deg]", "15"),
        },
    }
    rec = make_case(cfg, args.case_dir)
    print(f"Wrote {args.case_dir} ({rec['mesh']['cells']} cells)")


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="ceatofoam",
        description="NASA CEA chamber thermochemistry to OpenFOAM"
        " rocket-nozzle cases.",
    )
    parser.add_argument("--version", action="version",
                        version=f"ceatofoam {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("design", help="print the design point")
    p.add_argument("input")
    p.set_defaults(func=cmd_design)

    p = sub.add_parser("case", help="write an OpenFOAM case")
    p.add_argument("input")
    p.add_argument("case_dir")
    p.set_defaults(func=cmd_case)

    p = sub.add_parser("post", help="evaluate a finished case")
    p.add_argument("case_dir")
    p.add_argument("--window", type=float, default=0.1,
                   help="fraction of iterations averaged (default 0.1)")
    p.set_defaults(func=cmd_post)

    p = sub.add_parser("sweep", help="CEA over a range of O/F")
    p.add_argument("--fuel", required=True)
    p.add_argument("--oxidizer", required=True)
    p.add_argument("--pressure", required=True, help="e.g. '500 psi'")
    p.add_argument("--altitude", default="0", help="e.g. '0 ft'")
    p.add_argument("--of-min", type=float, default=1.0)
    p.add_argument("--of-max", type=float, default=5.0)
    p.add_argument("--of-step", type=float, default=0.05)
    p.add_argument("--analysis", default="equilibrium",
                   choices=["equilibrium", "frozen"])
    p.add_argument("--plot", help="write a figure to this file")
    p.set_defaults(func=cmd_sweep)

    p = sub.add_parser("interactive", help="console dialogue, then a case")
    p.add_argument("case_dir")
    p.set_defaults(func=cmd_interactive)

    args = parser.parse_args(argv)
    args.func(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
