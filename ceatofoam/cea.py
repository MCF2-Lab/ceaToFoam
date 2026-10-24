"""NASA CEA rocket problems through the CEA_Wrap package.

CEA (Gordon and McBride, NASA RP-1311, 1994/1996) is run twice for one
operating point: once with shifting equilibrium and once frozen at the
chamber composition (nfz = 1). The chamber state of both runs is identical;
the frozen run gives the one-dimensional reference that is consistent with
the frozen-mixture CFD model.

On Windows CEA_Wrap ships its executable. On Linux and macOS the bundled
Fortran source must be compiled once and its location exported as
CEA_EXE_LOCATION (see README).
"""

from dataclasses import dataclass, field

from .units import G0

#: Default propellant inlet temperatures [K]. For the cryogenic and
#: fixed-temperature condensed species these are the reference temperatures
#: of the species in the CEA thermodynamic library (thermo.inp, McBride
#: et al. 1996 for O2(L), H2(L), CH4(L), N2O4(L), C3H8(L); Mehta et al.
#: AIAA 95-2962 for RP-1; TRC 12/93 for CH6N2(L)); CEA assigns their
#: enthalpy at that temperature. The others are storable liquids at 298.15 K
#: inside the range of their thermo.inp fit.
PROPELLANT_TEMPERATURES = {
    "RP-1": 298.15,
    "H2(L)": 20.270,
    "CH4(L)": 111.643,
    "CH6N2(L)": 298.15,
    "N2H4(L)": 298.15,
    "C2H5OH(L)": 298.15,
    "C3H8(L)": 231.076,
    "O2(L)": 90.170,
    "N2O4(L)": 298.15,
    "O2": 298.15,
}


@dataclass
class CeaPoint:
    """Chamber state and one-dimensional performance of one CEA run.

    Pressures in Pa, temperatures in K, velocities in m/s, Isp in s,
    molecular weight in kg/kmol, viscosity in Pa s.
    """

    analysis: str
    pc: float
    of: float
    pressure_ratio: float
    tc: float
    mw_c: float
    gammas_c: float
    gamma_c: float
    cp_c: float
    visc_c: float
    composition: dict = field(default_factory=dict)
    #: CEA M,(1/n) at chamber, throat and exit (equal to the molecular
    #: weight without condensed species).
    m_stations: tuple = ()
    t_e: float = float("nan")
    p_e: float = float("nan")
    mach_e: float = float("nan")
    isp: float = float("nan")
    ivac: float = float("nan")
    cstar: float = float("nan")
    area_ratio: float = float("nan")

    @property
    def exit_velocity(self):
        """Exit velocity [m/s]; Isp of CEA is the ideal (matched) value."""
        return self.isp*G0


def _materials(fuel, oxidizer, t_fuel=None, t_oxidizer=None):
    import CEA_Wrap

    def temperature(name, value):
        if value is not None:
            return float(value)
        if name not in PROPELLANT_TEMPERATURES:
            raise ValueError(
                f"no default temperature for {name!r}; give it explicitly"
            )
        return PROPELLANT_TEMPERATURES[name]

    return [
        CEA_Wrap.Fuel(fuel, temperature(fuel, t_fuel)),
        CEA_Wrap.Oxidizer(oxidizer, temperature(oxidizer, t_oxidizer)),
    ]


def run_rocket(
    fuel, oxidizer, pc, of, pressure_ratio, analysis="equilibrium",
    t_fuel=None, t_oxidizer=None,
):
    """Run one CEA rocket problem and return a :class:`CeaPoint`.

    ``pc`` is the chamber pressure [Pa], ``pressure_ratio`` is pc/pe (the
    CEA ``pip`` option) and ``analysis`` is ``"equilibrium"`` or ``"frozen"``
    (frozen at the chamber, nfz = 1).
    """
    import CEA_Wrap

    if analysis not in ("equilibrium", "frozen"):
        raise ValueError(f"unknown CEA analysis {analysis!r}")
    options = {"nfz": 1} if analysis == "frozen" else {}
    problem = CEA_Wrap.RocketProblem(
        pressure=pc/1.0e5,
        pressure_units="bar",
        o_f=of,
        pip=pressure_ratio,
        analysis_type=analysis,
        materials=_materials(fuel, oxidizer, t_fuel, t_oxidizer),
        **options,
    )
    out = problem.run()
    if out.c_condensed:
        raise ValueError(
            "CEA predicts condensed products in the chamber; the single-phase"
            " frozen mixture model does not apply"
        )
    return CeaPoint(
        analysis=analysis,
        pc=pc,
        of=of,
        pressure_ratio=pressure_ratio,
        tc=out.c_t,
        mw_c=out.c_mw,
        gammas_c=out.c_gammas,
        gamma_c=out.c_gamma,
        cp_c=out.c_cp*1.0e3,
        visc_c=out.c_visc,
        composition={k: v for k, v in out.prod_c.items() if v > 0.0},
        m_stations=(out.c_m, out.t_m, out.m),
        t_e=out.t,
        p_e=out.p*1.0e5,
        mach_e=out.mach,
        isp=out.isp,
        ivac=out.ivac,
        cstar=out.cstar,
        area_ratio=out.ae if out.ae > 0.0 else float("nan"),
    )


def sweep_mixture_ratio(
    fuel, oxidizer, pc, pressure_ratio, ratios, analysis="equilibrium",
    t_fuel=None, t_oxidizer=None,
):
    """Run CEA over the mixture ratios ``ratios``; return a list of points.

    A ratio at which CEA fails or predicts condensed chamber products is
    returned as None so that the caller can plot gaps.
    """
    points = []
    for of in ratios:
        try:
            points.append(
                run_rocket(
                    fuel, oxidizer, pc, float(of), pressure_ratio, analysis,
                    t_fuel, t_oxidizer,
                )
            )
        except (RuntimeError, ValueError):
            points.append(None)
    return points
