"""Die Simulation selbst: Minimierung, Äquilibrierung, Produktion (OpenMM).

OpenMM wird erst in den Funktionen importiert, die es brauchen. So laufen die
Tests und alle anderen Module ohne OpenMM, und `doctor` kann berichten, dass es
fehlt, statt mit einem ImportError abzubrechen.

Die Planung (wie viele Schritte, welche Plattform, welche Atome festhalten) ist
von der Ausführung getrennt. Der planende Teil ist reine Rechnerei und wird
getestet; der ausführende braucht eine Grafikkarte.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from .config import ProtocolConfig
from .gates import GateReport, GateResult, gate, run_gates
from .run import Run

#: Reihenfolge, in der eine Plattform gewählt wird, wenn "auto" eingestellt ist
PLATTFORM_VORZUG = ("CUDA", "HIP", "OpenCL", "CPU")

#: Cutoff für die nichtbindenden Wechselwirkungen
CUTOFF_NM = 0.9
#: Masse, auf die Wasserstoffe bei HMR gesetzt werden
HMR_MASSE_AMU = 1.5
#: Reibung des Thermostaten
REIBUNG_PRO_PS = 1.0
#: Wie oft Energie und Temperatur ins Logbuch geschrieben werden
PROTOKOLL_PS = 10.0

TRAJEKTORIE = "produktion.dcd"
ZUSTAND = "produktion.xml"
PRUEFPUNKT = "produktion.chk"
MESSWERTE = "produktion.csv"


class SimulateError(RuntimeError):
    pass


# ---------- Planung (ohne OpenMM) ----------

@dataclass
class Stufe:
    name: str
    schritte: int
    temperatur_k: float
    druck_bar: float | None
    restraints_kj_nm2: float
    zeitschritt_fs: float

    @property
    def ps(self) -> float:
        return round(self.schritte * self.zeitschritt_fs / 1000.0, 3)


@dataclass
class Plan:
    minimierung_schritte: int
    aequilibrierung: list[Stufe]
    produktion: Stufe
    speichern_alle: int           # Schritte zwischen zwei Bildern der Trajektorie
    protokoll_alle: int           # Schritte zwischen zwei Zeilen Messwerte
    zufallsstartwert: int
    hmr: bool
    stufen: list[Stufe] = field(default_factory=list)

    def __post_init__(self):
        self.stufen = [*self.aequilibrierung, self.produktion]

    @property
    def gesamt_ns(self) -> float:
        return round(sum(s.ps for s in self.stufen) / 1000.0, 4)


def _schritte(ps: float, zeitschritt_fs: float) -> int:
    return max(1, int(round(ps * 1000.0 / zeitschritt_fs)))


def plan_erstellen(protokoll: ProtocolConfig) -> Plan:
    """Rechnet Zeiten in Schritte um.

    Die Äquilibrierung läuft bewusst mit 2 fs, auch wenn die Produktion 4 fs
    benutzt: Solange noch Restraints wirken und sich die Box setzt, ist der
    kleine Zeitschritt die sicherere Wahl.
    """
    p = protokoll.produktion
    aequi_dt = min(2.0, p.zeitschritt_fs)
    aequilibrierung = [
        Stufe(name=stufe.name, schritte=_schritte(stufe.ps, aequi_dt),
              temperatur_k=stufe.temperatur_k, druck_bar=stufe.druck_bar,
              restraints_kj_nm2=stufe.restraints_kj_nm2, zeitschritt_fs=aequi_dt)
        for stufe in protokoll.aequilibrierung
    ]
    produktion = Stufe(name="produktion", schritte=_schritte(p.ns * 1000.0, p.zeitschritt_fs),
                       temperatur_k=p.temperatur_k, druck_bar=p.druck_bar,
                       restraints_kj_nm2=0.0, zeitschritt_fs=p.zeitschritt_fs)
    return Plan(
        minimierung_schritte=protokoll.minimierung.schritte,
        aequilibrierung=aequilibrierung,
        produktion=produktion,
        speichern_alle=_schritte(p.speichern_ps, p.zeitschritt_fs),
        protokoll_alle=_schritte(PROTOKOLL_PS, p.zeitschritt_fs),
        zufallsstartwert=protokoll.zufallsstartwert,
        hmr=p.hmr,
    )


def waehle_plattform(gewuenscht: str, verfuegbar: list[str]) -> str:
    """Welche Plattform benutzt wird. 'auto' nimmt die schnellste vorhandene."""
    if gewuenscht != "auto":
        if gewuenscht not in verfuegbar:
            raise SimulateError(
                f"Die Plattform {gewuenscht} gibt es in dieser OpenMM-Installation nicht. "
                f"Vorhanden: {', '.join(verfuegbar) or 'keine'}. Bei einer AMD-Karte läuft "
                "OpenMM unter nativem Windows über OpenCL, unter WSL2 meist gar nicht."
            )
        return gewuenscht
    for name in PLATTFORM_VORZUG:
        if name in verfuegbar:
            return name
    raise SimulateError("OpenMM meldet keine einzige Rechenplattform.")


# ---------- Prüfpunkte auf die Simulation ----------

@dataclass
class SimContext:
    messwerte: list[dict]
    soll_temperatur_k: float
    stufe: str = "produktion"
    toleranz_k: float = 15.0


def _spalte(messwerte: list[dict], name: str) -> list[float]:
    werte = []
    for zeile in messwerte:
        for schluessel, wert in zeile.items():
            if name.lower() in schluessel.lower():
                try:
                    werte.append(float(wert))
                except (TypeError, ValueError):
                    pass
                break
    return werte


@gate("simulation.laeuft", "Simulation hat Messwerte geliefert")
def gate_messwerte(ctx: SimContext) -> GateResult:
    if len(ctx.messwerte) < 2:
        return GateResult("simulation.laeuft", "Simulation hat Messwerte geliefert", "fehler",
                          f"Nur {len(ctx.messwerte)} Zeile(n) in {MESSWERTE}. Die Simulation "
                          "ist abgebrochen, bevor sie etwas gerechnet hat.")
    return GateResult("simulation.laeuft", "Simulation hat Messwerte geliefert", "ok",
                      f"{len(ctx.messwerte)} Messpunkte aufgezeichnet.")


@gate("simulation.energie", "Energie bleibt endlich")
def gate_energie(ctx: SimContext) -> GateResult:
    import math
    energien = _spalte(ctx.messwerte, "Total Energy")
    schlecht = [e for e in energien if math.isnan(e) or math.isinf(e)]
    if schlecht or not energien:
        return GateResult("simulation.energie", "Energie bleibt endlich", "fehler",
                          "Die Gesamtenergie ist nicht endlich (NaN oder unendlich). Das System "
                          "ist explodiert: meist zu großer Zeitschritt, schlechte Startstruktur "
                          "oder fehlerhafte Parameter.")
    return GateResult("simulation.energie", "Energie bleibt endlich", "ok",
                      f"Gesamtenergie zwischen {min(energien):.0f} und {max(energien):.0f} kJ/mol.")


@gate("simulation.temperatur", "Temperatur stimmt")
def gate_temperatur(ctx: SimContext) -> GateResult:
    temperaturen = _spalte(ctx.messwerte, "Temperature")
    if not temperaturen:
        return GateResult("simulation.temperatur", "Temperatur stimmt", "warnung",
                          "In den Messwerten steht keine Temperatur.")
    # die zweite Hälfte zählt; am Anfang heizt das System noch auf
    haelfte = temperaturen[len(temperaturen) // 2:]
    mittel = sum(haelfte) / len(haelfte)
    if abs(mittel - ctx.soll_temperatur_k) > ctx.toleranz_k:
        return GateResult(
            "simulation.temperatur", "Temperatur stimmt", "fehler",
            f"Die mittlere Temperatur liegt bei {mittel:.1f} K statt {ctx.soll_temperatur_k:.0f} K. "
            "Prüfe Thermostat, Zeitschritt und ob die Äquilibrierung lang genug war.",
            {"mittel_k": round(mittel, 1)},
        )
    return GateResult("simulation.temperatur", "Temperatur stimmt", "ok",
                      f"Mittlere Temperatur {mittel:.1f} K.", {"mittel_k": round(mittel, 1)})


@gate("simulation.dichte", "Dichte plausibel")
def gate_dichte(ctx: SimContext) -> GateResult:
    dichten = _spalte(ctx.messwerte, "Density")
    if not dichten:
        return GateResult("simulation.dichte", "Dichte plausibel", "warnung",
                          "In den Messwerten steht keine Dichte (NVT-Lauf?).")
    haelfte = dichten[len(dichten) // 2:]
    mittel = sum(haelfte) / len(haelfte)
    if not 0.9 < mittel < 1.15:
        return GateResult("simulation.dichte", "Dichte plausibel", "fehler",
                          f"Die Dichte liegt bei {mittel:.3f} g/ml. Wässrige Systeme liegen bei "
                          "etwa 1,0; das deutet auf ein Problem mit Barostat oder Box hin.",
                          {"mittel_g_ml": round(mittel, 3)})
    return GateResult("simulation.dichte", "Dichte plausibel", "ok",
                      f"Mittlere Dichte {mittel:.3f} g/ml.", {"mittel_g_ml": round(mittel, 3)})


SIM_GATES = [gate_messwerte, gate_energie, gate_temperatur, gate_dichte]


def lies_messwerte(pfad: Path) -> list[dict]:
    """Die CSV-Datei, die OpenMM mitschreibt."""
    import csv
    if not pfad.is_file():
        return []
    with pfad.open(encoding="utf-8", errors="replace") as datei:
        text = datei.read()
    zeilen = [z.lstrip("#").strip() for z in text.splitlines() if z.strip()]
    if len(zeilen) < 2:
        return []
    return list(csv.DictReader(zeilen))


# ---------- Ausführung (braucht OpenMM) ----------

def _plattformen() -> list[str]:
    import openmm
    return [openmm.Platform.getPlatform(i).getName()
            for i in range(openmm.Platform.getNumPlatforms())]


def _restraint_kraft(prmtop, positionen, staerke_kj_nm2: float):
    """Hält die Schweratome des Soluts an ihrem Platz, Wasser und Ionen nicht."""
    import openmm
    from openmm import unit
    from .prepare import WASSER_RESTE

    kraft = openmm.CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2")
    kraft.addGlobalParameter("k", staerke_kj_nm2 * unit.kilojoule_per_mole / unit.nanometer**2)
    for name in ("x0", "y0", "z0"):
        kraft.addPerParticleParameter(name)
    ionen_und_wasser = WASSER_RESTE | {"K+", "Na+", "Cl-", "MG", "CA", "ZN", "Li+", "IOD"}
    n = 0
    for atom in prmtop.topology.atoms():
        if atom.residue.name.strip().upper() in ionen_und_wasser:
            continue
        if atom.element is not None and atom.element.symbol == "H":
            continue
        kraft.addParticle(atom.index, positionen[atom.index].value_in_unit(unit.nanometer))
        n += 1
    return kraft, n


def simulieren(lauf: Run, fortsetzen: bool = True) -> GateReport:
    """Minimierung, Äquilibrierung und Produktion. Braucht OpenMM und eine Grafikkarte."""
    try:
        import openmm
        from openmm import app, unit
    except ImportError:
        raise SimulateError(
            "OpenMM ist in dieser Umgebung nicht installiert. Die Simulation läuft unter "
            "nativem Windows (environment-windows.yml), nicht in WSL2."
        ) from None

    bericht_vorher = lauf.lies_bericht("vorbereiten")
    if not (bericht_vorher and bericht_vorher.get("bestanden")):
        raise SimulateError(
            "Für diesen Lauf gibt es kein bestandenes 'vorbereiten'. Erst "
            "'python -m oligomd vorbereiten', dann simulieren."
        )
    prmtop_pfad = lauf.ordner / bericht_vorher["prmtop"]
    inpcrd_pfad = lauf.ordner / bericht_vorher["inpcrd"]

    protokoll = lauf.konfiguration.protokoll
    plan = plan_erstellen(protokoll)
    arbeit = lauf.pfad("simulation")
    arbeit.mkdir(parents=True, exist_ok=True)

    prmtop = app.AmberPrmtopFile(str(prmtop_pfad))
    inpcrd = app.AmberInpcrdFile(str(inpcrd_pfad))
    system = prmtop.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=CUTOFF_NM * unit.nanometer,
        constraints=app.HBonds,
        rigidWater=True,
        hydrogenMass=(HMR_MASSE_AMU * unit.amu) if plan.hmr else None,
    )
    plattform_name = waehle_plattform(protokoll.plattform, _plattformen())
    plattform = openmm.Platform.getPlatformByName(plattform_name)
    lauf.log(f"Simulation auf Plattform {plattform_name}, "
             f"{system.getNumParticles()} Teilchen, {plan.gesamt_ns} ns geplant")

    integrator = openmm.LangevinMiddleIntegrator(
        plan.aequilibrierung[0].temperatur_k * unit.kelvin,
        REIBUNG_PRO_PS / unit.picosecond,
        plan.aequilibrierung[0].zeitschritt_fs * unit.femtosecond)
    integrator.setRandomNumberSeed(plan.zufallsstartwert)
    simulation = app.Simulation(prmtop.topology, system, integrator, plattform)
    simulation.context.setPositions(inpcrd.positions)
    if inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

    # Minimierung
    lauf.log(f"Minimierung über {plan.minimierung_schritte} Schritte")
    simulation.minimizeEnergy(maxIterations=plan.minimierung_schritte)
    minimiert = simulation.context.getState(getPositions=True).getPositions()
    with (arbeit / "minimiert.pdb").open("w", encoding="utf-8") as datei:
        app.PDBFile.writeFile(prmtop.topology, minimiert, datei)

    # Äquilibrierung: je Stufe ein eigener Kontext, weil Restraints und Barostat wechseln
    positionen, boxvektoren = minimiert, simulation.context.getState().getPeriodicBoxVectors()
    for stufe in plan.aequilibrierung:
        positionen, boxvektoren = _stufe_laufen(
            lauf, arbeit, prmtop, stufe, positionen, boxvektoren, plan)

    # Produktion
    bericht = _produktion(lauf, arbeit, prmtop, plan, positionen, boxvektoren, fortsetzen)
    lauf.schreibe_bericht("simulieren", {
        "bestanden": bericht.ok,
        "plattform": plattform_name,
        "teilchen": system.getNumParticles(),
        "ns": plan.produktion.ps / 1000.0,
        "trajektorie": str((arbeit / TRAJEKTORIE).relative_to(lauf.ordner)),
        **bericht.als_dict(),
    })
    lauf.log("Simulieren " + ("bestanden" if bericht.ok else "NICHT bestanden"))
    return bericht


def _stufe_laufen(lauf: Run, arbeit: Path, prmtop, stufe: Stufe, positionen, boxvektoren, plan: Plan):
    import openmm
    from openmm import app, unit

    system = prmtop.createSystem(
        nonbondedMethod=app.PME, nonbondedCutoff=CUTOFF_NM * unit.nanometer,
        constraints=app.HBonds, rigidWater=True,
        hydrogenMass=(HMR_MASSE_AMU * unit.amu) if plan.hmr else None)
    if stufe.restraints_kj_nm2 > 0:
        kraft, n = _restraint_kraft(prmtop, positionen, stufe.restraints_kj_nm2)
        system.addForce(kraft)
        lauf.log(f"Stufe {stufe.name}: {n} Atome festgehalten mit "
                 f"{stufe.restraints_kj_nm2} kJ/mol/nm²")
    if stufe.druck_bar is not None:
        system.addForce(openmm.MonteCarloBarostat(
            stufe.druck_bar * unit.bar, stufe.temperatur_k * unit.kelvin))

    integrator = openmm.LangevinMiddleIntegrator(
        stufe.temperatur_k * unit.kelvin, REIBUNG_PRO_PS / unit.picosecond,
        stufe.zeitschritt_fs * unit.femtosecond)
    integrator.setRandomNumberSeed(plan.zufallsstartwert)
    simulation = app.Simulation(prmtop.topology, system, integrator)
    simulation.context.setPositions(positionen)
    simulation.context.setPeriodicBoxVectors(*boxvektoren)
    simulation.context.setVelocitiesToTemperature(stufe.temperatur_k * unit.kelvin,
                                                  plan.zufallsstartwert)
    simulation.reporters.append(app.StateDataReporter(
        str(arbeit / f"{stufe.name}.csv"), plan.protokoll_alle, step=True, time=True,
        potentialEnergy=True, totalEnergy=True, temperature=True, density=True, speed=True))
    lauf.log(f"Stufe {stufe.name}: {stufe.schritte} Schritte à {stufe.zeitschritt_fs} fs "
             f"= {stufe.ps} ps")
    simulation.step(stufe.schritte)
    zustand = simulation.context.getState(getPositions=True)
    return zustand.getPositions(), zustand.getPeriodicBoxVectors()


def _produktion(lauf: Run, arbeit: Path, prmtop, plan: Plan, positionen, boxvektoren,
                fortsetzen: bool) -> GateReport:
    import openmm
    from openmm import app, unit

    stufe = plan.produktion
    system = prmtop.createSystem(
        nonbondedMethod=app.PME, nonbondedCutoff=CUTOFF_NM * unit.nanometer,
        constraints=app.HBonds, rigidWater=True,
        hydrogenMass=(HMR_MASSE_AMU * unit.amu) if plan.hmr else None)
    system.addForce(openmm.MonteCarloBarostat(stufe.druck_bar * unit.bar,
                                              stufe.temperatur_k * unit.kelvin))
    integrator = openmm.LangevinMiddleIntegrator(
        stufe.temperatur_k * unit.kelvin, REIBUNG_PRO_PS / unit.picosecond,
        stufe.zeitschritt_fs * unit.femtosecond)
    integrator.setRandomNumberSeed(plan.zufallsstartwert)
    simulation = app.Simulation(prmtop.topology, system, integrator)

    pruefpunkt = arbeit / PRUEFPUNKT
    offen = "a"
    if fortsetzen and pruefpunkt.is_file():
        simulation.loadCheckpoint(str(pruefpunkt))
        getan = simulation.currentStep
        lauf.log(f"Produktion wird fortgesetzt ab Schritt {getan}")
    else:
        simulation.context.setPositions(positionen)
        simulation.context.setPeriodicBoxVectors(*boxvektoren)
        simulation.context.setVelocitiesToTemperature(stufe.temperatur_k * unit.kelvin,
                                                      plan.zufallsstartwert)
        getan, offen = 0, "w"

    rest = max(0, stufe.schritte - getan)
    simulation.reporters.append(app.DCDReporter(str(arbeit / TRAJEKTORIE), plan.speichern_alle,
                                                append=(offen == "a")))
    simulation.reporters.append(app.StateDataReporter(
        str(arbeit / MESSWERTE), plan.protokoll_alle, step=True, time=True,
        potentialEnergy=True, totalEnergy=True, temperature=True, density=True,
        progress=True, remainingTime=True, speed=True,
        totalSteps=stufe.schritte, append=(offen == "a")))
    simulation.reporters.append(app.CheckpointReporter(str(pruefpunkt), plan.speichern_alle * 10))

    lauf.log(f"Produktion: {rest} Schritte à {stufe.zeitschritt_fs} fs = {stufe.ps} ps")
    simulation.step(rest)
    simulation.saveState(str(arbeit / ZUSTAND))

    return run_gates(SIM_GATES, SimContext(
        messwerte=lies_messwerte(arbeit / MESSWERTE),
        soll_temperatur_k=stufe.temperatur_k))
