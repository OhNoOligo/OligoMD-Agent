"""Struktur → fertiges System für die Simulation (tleap).

Ablauf:
  1. Struktur lesen, auf die gewählten Ketten eindampfen, Prüfpunkte anwenden
  2. als PDB schreiben, ohne Wasserstoffe - tleap setzt sie mit den Namen des
     Kraftfelds selbst
  3. tleap, erster Durchgang: solvatisieren und zählen, wie viele Wasser es wurden
  4. daraus die Ionenzahlen ausrechnen
  5. tleap, zweiter Durchgang: solvatisieren, Ionen setzen, prmtop und inpcrd schreiben
  6. die Ausgabe von tleap lesen und die Prüfpunkte darauf anwenden

Zwei Durchgänge, weil tleap nicht rechnen kann: Die Zahl der Salzionen hängt
davon ab, wie viele Wassermoleküle in der Box gelandet sind.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from .config import SystemConfig
from .gates import GateReport, GateSettings, StructureContext, check_structure, run_gates
from .run import Run
from .structure import Structure, inventory, net_charge, read_structure, write_pdb
from .tools import ToolError, run_tool

#: Mol Wasser pro Liter. Daraus wird aus der Salzkonzentration die Ionenzahl.
WASSER_PRO_LITER = 55.56

#: leaprc-Dateien je Wassermodell und der zugehörige Boxname in tleap
WASSERMODELLE = {
    "OPC": ("leaprc.water.opc", "OPCBOX"),
    "TIP3P": ("leaprc.water.tip3p", "TIP3PBOX"),
    "TIP4P-Ew": ("leaprc.water.tip4pew", "TIP4PEWBOX"),
    "SPC/E": ("leaprc.water.spce", "SPCBOX"),
}

#: Schreibweisen, die Nutzer eingeben, auf die Restnamen von Amber abgebildet
IONEN_NAMEN = {
    "MG": "MG", "MG2+": "MG", "MG2": "MG", "MAGNESIUM": "MG",
    "CA": "CA", "CA2+": "CA", "ZN": "ZN", "ZN2+": "ZN", "MN": "MN",
    "K": "K+", "K+": "K+", "NA": "Na+", "NA+": "Na+", "LI": "Li+",
    "CL": "Cl-", "CL-": "Cl-", "BR": "Br-", "I": "IOD",
}

#: Ladung der Ionen unter ihrem Amber-Restnamen
AMBER_ION_LADUNG = {"K+": 1, "Na+": 1, "Li+": 1, "Rb+": 1, "Cs+": 1,
                    "Cl-": -1, "Br-": -1, "IOD": -1, "F-": -1,
                    "MG": 2, "CA": 2, "ZN": 2, "MN": 2, "FE2": 2}

#: Restnamen, die tleap für Wasser vergibt
WASSER_RESTE = {"WAT", "HOH", "SOL"}

PRMTOP = "system.prmtop"
INPCRD = "system.inpcrd"
SOLVATISIERT = "system_solvatisiert.pdb"
EINGABE_PDB = "eingabe.pdb"


class PrepareError(RuntimeError):
    pass


def amber_ion(name: str) -> str:
    """'Mg2+' → 'MG'. Unbekanntes wird unverändert durchgereicht; tleap meckert dann."""
    return IONEN_NAMEN.get(name.strip().upper(), name.strip())


def ion_ladung(amber_name: str) -> int:
    """Ladung eines Ions unter seinem Amber-Restnamen; 0, wenn unbekannt."""
    return AMBER_ION_LADUNG.get(amber_name, 0)


def extra_ionen(system: SystemConfig, n_wasser: int) -> dict[str, int]:
    """Zusätzliche Ionen als {Amber-Name: Anzahl}, Konzentration in Zahl umgerechnet."""
    out: dict[str, int] = {}
    for extra in system.ionen.zusaetzlich:
        name = amber_ion(extra.ion)
        anzahl = (extra.anzahl if extra.anzahl is not None
                  else max(1, int(round(extra.mol_l * n_wasser / WASSER_PRO_LITER))))
        out[name] = out.get(name, 0) + anzahl
    return out


def gegenionen(ladung_solut: int, extra: dict[str, int]) -> int:
    """Wie viele Gegenionen die Box braucht - die Zusatzionen zählen mit.

    Zehn Mg2+ bringen +20 mit. Wer erst neutralisiert und sie danach hineinwirft,
    hat am Ende eine geladene Box, und unter periodischen Randbedingungen wird
    die Elektrostatik davon unsinnig.
    """
    gesamt = ladung_solut + sum(ion_ladung(name) * n for name, n in extra.items())
    return abs(gesamt)


def salz_ionenpaare(n_wasser: int, salz_mol_l: float) -> int:
    """Zahl der Salzpaare für eine Konzentration, über das Verhältnis zum Wasser."""
    return int(round(salz_mol_l * n_wasser / WASSER_PRO_LITER))


# ---------- tleap-Skript ----------

@dataclass
class LeapPlan:
    """Alles, was in ein tleap-Skript muss. Rein datenhaltig, damit testbar."""
    system: SystemConfig
    pdb: str = EINGABE_PDB
    neutralisieren: int = 0          # Zahl der Gegenionen; 0 = tleap rechnet selbst
    salzpaare: int = 0
    nur_solvatisieren: bool = False  # erster Durchgang: zählen, nicht speichern
    extra_ionen: dict[str, int] = field(default_factory=dict)
    ausgleich_ion: str | None = None  # None = das Gegenion aus der Konfiguration
    zusatz_lib: list[str] = field(default_factory=list)
    zusatz_frcmod: list[str] = field(default_factory=list)


def leap_script(plan: LeapPlan) -> str:
    """Baut das tleap-Skript. Jede Zeile ist absichtlich sichtbar und nachlesbar."""
    s = plan.system
    if s.kraftfeld.wasser not in WASSERMODELLE:
        raise PrepareError(
            f"Unbekanntes Wassermodell {s.kraftfeld.wasser}. Möglich: "
            + ", ".join(WASSERMODELLE)
        )
    wasser_leaprc, boxname = WASSERMODELLE[s.kraftfeld.wasser]

    zeilen = [
        "# Von OligoMD erzeugt. Nicht von Hand ändern - der Lauf-Ordner ist die Wahrheit.",
        f"source leaprc.DNA.{s.kraftfeld.dna}",
        f"source leaprc.RNA.{s.kraftfeld.rna}",
        f"source {wasser_leaprc}",
    ]
    for lib in plan.zusatz_lib:
        zeilen.append(f"loadoff {lib}")
    for frcmod in plan.zusatz_frcmod:
        zeilen.append(f"loadamberparams {frcmod}")

    zeilen += [
        "",
        f"mol = loadpdb {plan.pdb}",
        "check mol",
        "",
        "# Box: Puffer in Ångström, die Konfiguration rechnet in Nanometern",
    ]
    puffer_a = round(s.box.puffer_nm * 10.0, 2)
    befehl = "solvateoct" if s.box.typ == "oktaeder" else "solvatebox"
    zeilen.append(f"{befehl} mol {boxname} {puffer_a}")

    if plan.nur_solvatisieren:
        # Erster Durchgang: nur zählen, wie viele Wasser hineinpassen
        zeilen += ["", f"savepdb mol {SOLVATISIERT}", "quit", ""]
        return "\n".join(zeilen)

    gegenion = amber_ion(s.ionen.neutralisieren_mit)
    # Reihenfolge ist wichtig: erst die Zusatzionen, dann ausgleichen. Sonst
    # trägt die Box am Ende die Ladung der Zusatzionen.
    for name, anzahl in plan.extra_ionen.items():
        ladung = ion_ladung(name)
        zeilen += [
            "",
            f"# frei in der Lösung, ohne Vorgabe der Bindestelle: {anzahl} x {name} "
            f"(bringt {ladung * anzahl:+d} Ladung mit)",
            f"addionsrand mol {name} {anzahl}",
        ]
    ausgleich = amber_ion(plan.ausgleich_ion or s.ionen.neutralisieren_mit)
    zeilen += [
        "",
        f"# Ladung ausgleichen: {plan.neutralisieren} x {ausgleich}",
        f"addions mol {ausgleich} {plan.neutralisieren}",
    ]
    if plan.salzpaare > 0:
        zeilen += [
            f"# Salz: {plan.salzpaare} Paare für {s.ionen.salz_mol_l} mol/l",
            f"addionsrand mol {gegenion} {plan.salzpaare} Cl- {plan.salzpaare}",
        ]
    zeilen += [
        "",
        "check mol",
        f"saveamberparm mol {PRMTOP} {INPCRD}",
        f"savepdb mol {SOLVATISIERT}",
        "quit",
        "",
    ]
    return "\n".join(zeilen)


# ---------- tleap-Ausgabe lesen ----------

_FATAL = re.compile(r"(FATAL|Fatal Error|Could not open|No such file)", re.I)
_FEHLENDE_PARAMETER = re.compile(
    r"(Could not find (?:bond|angle|torsion) parameter|"
    r"^\s*\*\*.*(?:missing|not found)|"
    r"Parameter file was not saved|"
    r"Unknown residue)", re.I | re.M)
# tleap schreibt die Ladung in zwei Formen, je nach Befehl:
#   "Total unperturbed charge:  -50.000000"   (saveamberparm)
#   "WARNING: The unperturbed charge of the unit: -50.000000 is not zero."  (check)
_LADUNG = re.compile(
    r"(?:Total unperturbed charge|unperturbed charge of the unit)\s*[:(]?\s*(-?[\d.]+)", re.I)
_ERZEUGT = re.compile(r"Added\s+(\d+)\s+residues", re.I)


@dataclass
class LeapErgebnis:
    ausgabe: str
    rueckgabe: int
    restladung: float | None = None
    fatale: list[str] = field(default_factory=list)
    fehlende_parameter: list[str] = field(default_factory=list)

    @property
    def erfolg(self) -> bool:
        return self.rueckgabe == 0 and not self.fatale and not self.fehlende_parameter


def lies_leap_ausgabe(ausgabe: str, rueckgabe: int) -> LeapErgebnis:
    """tleap meldet Fehler im Text, nicht im Rückgabewert. Also Text lesen."""
    zeilen = ausgabe.splitlines()
    fatale = [z.strip() for z in zeilen if _FATAL.search(z)]
    fehlend = [z.strip() for z in zeilen if _FEHLENDE_PARAMETER.search(z)]
    ladung = None
    treffer = _LADUNG.findall(ausgabe)
    if treffer:
        try:
            ladung = float(treffer[-1])
        except ValueError:
            ladung = None
    return LeapErgebnis(ausgabe, rueckgabe, ladung, fatale, fehlend)


def zaehle_wasser(pdb: Path) -> int:
    """Wassermoleküle in einer von tleap geschriebenen PDB-Datei.

    Gezählt werden die Sauerstoffe, nicht die Restnummern: Letztere laufen im
    PDB-Format bei 9999 um, und eine solvatisierte Box hat leicht mehr Wasser.
    """
    if not pdb.is_file():
        return 0
    anzahl = 0
    for zeile in pdb.read_text(encoding="utf-8", errors="replace").splitlines():
        if not zeile.startswith(("ATOM  ", "HETATM")):
            continue
        if zeile[17:20].strip().upper() not in WASSER_RESTE:
            continue
        name = zeile[12:16].strip().upper()
        element = (zeile[76:78].strip().upper() if len(zeile) >= 78 else "")
        if name in ("O", "OW", "OH2") or (element == "O" and not name.startswith("H")):
            anzahl += 1
    return anzahl


# ---------- Prüfpunkte auf das gebaute System ----------

@dataclass
class SystemContext:
    ergebnis: LeapErgebnis
    n_wasser: int
    erwartete_ladung: int
    prmtop: Path
    inpcrd: Path
    ionen_bestellt: dict[str, int]
    einstellungen: GateSettings = field(default_factory=GateSettings)


def _gate(gate_id: str, titel: str):
    from .gates import gate as _dekorator
    return _dekorator(gate_id, titel)


@_gate("system.tleap", "tleap ohne Abbruch")
def gate_leap_lief(ctx: SystemContext):
    from .gates import GateResult
    if ctx.ergebnis.fatale:
        return GateResult("system.tleap", "tleap ohne Abbruch", "fehler",
                          "tleap hat abgebrochen: " + ctx.ergebnis.fatale[0],
                          {"meldungen": ctx.ergebnis.fatale[:10]})
    if ctx.ergebnis.rueckgabe != 0:
        return GateResult("system.tleap", "tleap ohne Abbruch", "fehler",
                          f"tleap endete mit Rückgabewert {ctx.ergebnis.rueckgabe}. "
                          "Die ganze Ausgabe steht in tleap.log im Lauf-Ordner.")
    return GateResult("system.tleap", "tleap ohne Abbruch", "ok", "tleap ist sauber durchgelaufen.")


@_gate("system.parameter", "Alle Parameter vorhanden")
def gate_parameter(ctx: SystemContext):
    from .gates import GateResult
    if ctx.ergebnis.fehlende_parameter:
        erste = ctx.ergebnis.fehlende_parameter[:5]
        return GateResult(
            "system.parameter", "Alle Parameter vorhanden", "fehler",
            f"{len(ctx.ergebnis.fehlende_parameter)} Meldung(en) über fehlende Parameter oder "
            f"unbekannte Reste: {erste[0]}. Ohne vollständige Parameter wird nicht simuliert.",
            {"meldungen": erste},
        )
    return GateResult("system.parameter", "Alle Parameter vorhanden", "ok",
                      "tleap hat für jeden Rest Parameter gefunden.")


@_gate("system.dateien", "prmtop und inpcrd geschrieben")
def gate_dateien(ctx: SystemContext):
    from .gates import GateResult
    fehlend = [p.name for p in (ctx.prmtop, ctx.inpcrd) if not p.is_file()]
    if fehlend:
        return GateResult("system.dateien", "prmtop und inpcrd geschrieben", "fehler",
                          f"tleap hat {', '.join(fehlend)} nicht geschrieben.")
    groessen = {p.name: p.stat().st_size for p in (ctx.prmtop, ctx.inpcrd)}
    return GateResult("system.dateien", "prmtop und inpcrd geschrieben", "ok",
                      f"Topologie und Koordinaten liegen vor ({groessen[PRMTOP] // 1024} KB).",
                      {"groessen": groessen})


@_gate("system.neutral", "System ist neutral")
def gate_neutral(ctx: SystemContext):
    from .gates import GateResult
    ladung = ctx.ergebnis.restladung
    if ladung is None:
        return GateResult("system.neutral", "System ist neutral", "warnung",
                          "tleap hat keine Restladung gemeldet; prüfe tleap.log von Hand.")
    if abs(ladung) > 0.01:
        return GateResult(
            "system.neutral", "System ist neutral", "fehler",
            f"Das System trägt noch {ladung:+.2f} Elementarladungen. Erwartet waren "
            f"{ctx.erwartete_ladung:+d} vor dem Ausgleich. Mit geladener Box wird die "
            "Elektrostatik unter periodischen Randbedingungen unsinnig.",
            {"restladung": ladung},
        )
    ionen = ", ".join(f"{n}x {name}" for name, n in sorted(ctx.ionen_bestellt.items()))
    return GateResult("system.neutral", "System ist neutral", "ok",
                      f"Restladung {ladung:+.2f}. Das Solut trägt {ctx.erwartete_ladung:+d}, "
                      f"ausgeglichen mit: {ionen}.",
                      {"restladung": ladung, "ionen": ctx.ionen_bestellt})


@_gate("system.wasser", "Wasserbox gefüllt")
def gate_wasser(ctx: SystemContext):
    from .gates import GateResult
    if ctx.n_wasser < 1000:
        return GateResult("system.wasser", "Wasserbox gefüllt", "fehler",
                          f"Nur {ctx.n_wasser} Wassermoleküle. Das ist zu wenig; prüfe "
                          "box.puffer_nm und ob die Solvatisierung überhaupt lief.",
                          {"wasser": ctx.n_wasser})
    return GateResult("system.wasser", "Wasserbox gefüllt", "ok",
                      f"{ctx.n_wasser} Wassermoleküle in der Box.", {"wasser": ctx.n_wasser})


SYSTEM_GATES = [gate_leap_lief, gate_parameter, gate_dateien, gate_neutral, gate_wasser]


# ---------- der Schritt selbst ----------

def lade_struktur(lauf: Run) -> Structure:
    s = lauf.konfiguration.system.struktur
    struktur = read_structure(lauf.struktur_datei(), s.modell)
    return struktur.subset(ketten=s.ketten, ohne_restnamen=set(s.entfernen))


def vorbereiten(lauf: Run) -> GateReport:
    """Führt den ganzen Schritt aus und schreibt den Bericht in den Lauf-Ordner."""
    system = lauf.konfiguration.system
    arbeit = lauf.pfad("vorbereitet")
    arbeit.mkdir(parents=True, exist_ok=True)

    # 1 und 2: Struktur lesen und prüfen
    struktur = lade_struktur(lauf)
    info = inventory(struktur)
    ladung = net_charge(struktur)
    lauf.log(f"Struktur gelesen: {info['atome']} Atome, Ladung {ladung:+d}")

    bericht = check_structure(StructureContext(
        structure=struktur,
        erwartete_ketten=system.struktur.ketten,
        erlaubte_sonderreste=set(system.struktur.behalten_hetatm),
        modifizierter_rest=system.modifikation.name if system.modifikation else None,
    ))
    if not bericht.ok:
        lauf.schreibe_bericht("vorbereiten", {"bestanden": False, "struktur": info,
                                              **bericht.als_dict()})
        return bericht

    write_pdb(struktur, arbeit / EINGABE_PDB, ohne_wasserstoff=True,
              kommentar=f"{lauf.struktur_datei().name}, Modell {system.struktur.modell}, "
                        f"Ketten {', '.join(system.struktur.ketten)}")

    # 3: erster tleap-Durchgang, nur zum Zählen
    plan = LeapPlan(system=system, nur_solvatisieren=True)
    (arbeit / "leap_zaehlen.in").write_text(leap_script(plan), encoding="utf-8")
    try:
        lauf_1 = run_tool("tleap", ["-f", "leap_zaehlen.in"], arbeit)
    except ToolError as exc:
        raise PrepareError(str(exc)) from None
    n_wasser = zaehle_wasser(arbeit / SOLVATISIERT)
    lauf.log(f"Solvatisierung: {n_wasser} Wassermoleküle")
    if n_wasser == 0:
        ergebnis = lies_leap_ausgabe(lauf_1.ausgabe, lauf_1.rueckgabe)
        bericht = run_gates(SYSTEM_GATES, SystemContext(
            ergebnis, 0, ladung, arbeit / PRMTOP, arbeit / INPCRD, {}))
        lauf.schreibe_bericht("vorbereiten", {"bestanden": False, **bericht.als_dict()})
        return bericht

    # 4: Ionen ausrechnen. Erst die Zusatzionen, denn sie verschieben die Ladung.
    extra = extra_ionen(system, n_wasser)
    ausgleich_anzahl = gegenionen(ladung, extra)
    gesamt = ladung + sum(ion_ladung(n) * a for n, a in extra.items())
    ausgleich_ion = system.ionen.neutralisieren_mit if gesamt < 0 else "Cl-"
    paare = salz_ionenpaare(n_wasser, system.ionen.salz_mol_l)

    bestellt = dict(extra)
    name_ausgleich = amber_ion(ausgleich_ion)
    bestellt[name_ausgleich] = bestellt.get(name_ausgleich, 0) + ausgleich_anzahl
    if paare:
        kation = amber_ion(system.ionen.neutralisieren_mit)
        bestellt[kation] = bestellt.get(kation, 0) + paare
        bestellt["Cl-"] = bestellt.get("Cl-", 0) + paare
    lauf.log(f"Ionen: {extra} zusätzlich, {ausgleich_anzahl} x {name_ausgleich} zum Ausgleich "
             f"(Solut {ladung:+d}, mit Zusatz {gesamt:+d}), {paare} Salzpaare")

    # 5: zweiter Durchgang, diesmal vollständig
    plan = LeapPlan(system=system, neutralisieren=ausgleich_anzahl, salzpaare=paare,
                    extra_ionen=extra, ausgleich_ion=ausgleich_ion,
                    zusatz_lib=[str(p) for p in system.kraftfeld.zusatz_lib],
                    zusatz_frcmod=[str(p) for p in system.kraftfeld.zusatz_frcmod])
    (arbeit / "leap.in").write_text(leap_script(plan), encoding="utf-8")
    try:
        lauf_2 = run_tool("tleap", ["-f", "leap.in"], arbeit)
    except ToolError as exc:
        raise PrepareError(str(exc)) from None

    # 6: prüfen
    ergebnis = lies_leap_ausgabe(lauf_2.ausgabe, lauf_2.rueckgabe)
    n_wasser = zaehle_wasser(arbeit / SOLVATISIERT) or n_wasser
    bericht = run_gates(SYSTEM_GATES, SystemContext(
        ergebnis=ergebnis, n_wasser=n_wasser, erwartete_ladung=ladung,
        prmtop=arbeit / PRMTOP, inpcrd=arbeit / INPCRD, ionen_bestellt=bestellt))

    lauf.schreibe_bericht("vorbereiten", {
        "bestanden": bericht.ok,
        "struktur": info,
        "ladung_solut": ladung,
        "wasser": n_wasser,
        "ionen": bestellt,
        "prmtop": str((arbeit / PRMTOP).relative_to(lauf.ordner)),
        "inpcrd": str((arbeit / INPCRD).relative_to(lauf.ordner)),
        **bericht.als_dict(),
    })
    lauf.log("Vorbereiten " + ("bestanden" if bericht.ok else "NICHT bestanden"))
    return bericht
