"""Prüfpunkte: die Tests, die zwischen den Schritten bestehen müssen.

Die Prüfpunkte stehen hier im Code, nicht in der Konfiguration. Einstellbar sind
nur Schwellen, und auch die nur in dokumentierten Grenzen (GateSettings).
Ein Prüfpunkt mit Status "fehler" hält die Pipeline an; weder Nutzer noch Agent
können ihn abschalten. Das ist Absicht: Ein übergangener Prüfpunkt erzeugt eine
Simulation, die durchläuft und falsch ist.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Literal, Optional

from pydantic import BaseModel, Field

from .structure import (
    Structure, clashes, inventory, missing_atoms, net_charge, unknown_residues,
)

Status = Literal["ok", "warnung", "fehler"]


class GateSettings(BaseModel):
    """Die einzigen Stellschrauben. Grenzen bewusst eng."""
    kollision_toleranz_a: float = Field(0.4, ge=0.2, le=0.8,
                                        description="Å unter der VdW-Summe, ab der es eine Kollision ist")
    kollisionen_erlaubt: int = Field(0, ge=0, le=10,
                                     description="geduldete Kollisionen, normalerweise 0")
    rmsd_grenze_a: float = Field(8.0, ge=2.0, le=20.0,
                                 description="RMSD, ab dem die Struktur als davongelaufen gilt")
    torsion_abweichung_kj: float = Field(4.0, ge=1.0, le=20.0,
                                         description="erlaubte Abweichung Kraftfeld gegen Quantenchemie")


@dataclass(frozen=True)
class GateResult:
    id: str
    titel: str
    status: Status
    meldung: str
    details: dict = field(default_factory=dict)

    @property
    def blockiert(self) -> bool:
        return self.status == "fehler"


@dataclass
class GateReport:
    ergebnisse: list[GateResult]

    @property
    def ok(self) -> bool:
        return not any(r.blockiert for r in self.ergebnisse)

    @property
    def fehler(self) -> list[GateResult]:
        return [r for r in self.ergebnisse if r.status == "fehler"]

    @property
    def warnungen(self) -> list[GateResult]:
        return [r for r in self.ergebnisse if r.status == "warnung"]

    def als_text(self) -> str:
        symbol = {"ok": "ok     ", "warnung": "Warnung", "fehler": "FEHLER "}
        zeilen = [f"[{symbol[r.status]}] {r.titel}: {r.meldung}" for r in self.ergebnisse]
        return "\n".join(zeilen)

    def als_dict(self) -> dict:
        return {
            "bestanden": self.ok,
            "pruefungen": [
                {"id": r.id, "titel": r.titel, "status": r.status,
                 "meldung": r.meldung, "details": r.details}
                for r in self.ergebnisse
            ],
        }


class GateError(RuntimeError):
    """Wird geworfen, wenn die Pipeline trotz blockierender Prüfpunkte weiterlaufen soll."""

    def __init__(self, bericht: GateReport):
        self.bericht = bericht
        namen = ", ".join(r.titel for r in bericht.fehler)
        super().__init__(f"Prüfpunkte nicht bestanden: {namen}\n\n{bericht.als_text()}")


# ---------- Rahmen ----------

@dataclass
class StructureContext:
    """Alles, was die Struktur-Prüfpunkte brauchen."""
    structure: Structure
    erwartete_ketten: list[str]
    erlaubte_sonderreste: set[str] = field(default_factory=set)
    modifizierter_rest: Optional[str] = None
    einstellungen: GateSettings = field(default_factory=GateSettings)


Gate = Callable[[StructureContext], GateResult]


def run_gates(gates: list[Gate], ctx) -> GateReport:
    ergebnisse = []
    for gate in gates:
        try:
            ergebnisse.append(gate(ctx))
        except Exception as exc:  # ein kaputter Prüfpunkt ist ein Fehler, kein Freifahrtschein
            ergebnisse.append(GateResult(
                id=getattr(gate, "gate_id", gate.__name__),
                titel=getattr(gate, "gate_titel", gate.__name__),
                status="fehler",
                meldung=f"Der Prüfpunkt selbst ist abgestürzt: {type(exc).__name__}: {exc}",
            ))
    return GateReport(ergebnisse)


def gate(gate_id: str, titel: str):
    def decorate(func: Gate) -> Gate:
        func.gate_id = gate_id  # type: ignore[attr-defined]
        func.gate_titel = titel  # type: ignore[attr-defined]
        return func
    return decorate


# ---------- Struktur-Prüfpunkte ----------

@gate("struktur.ketten", "Ketten vorhanden")
def gate_chains(ctx: StructureContext) -> GateResult:
    vorhanden = set(ctx.structure.chains())
    fehlend = [k for k in ctx.erwartete_ketten if k not in vorhanden]
    if fehlend:
        return GateResult(
            "struktur.ketten", "Ketten vorhanden", "fehler",
            f"Die Struktur hat keine Kette {', '.join(fehlend)}. "
            f"Vorhanden sind: {', '.join(sorted(vorhanden))}.",
            {"fehlend": fehlend, "vorhanden": sorted(vorhanden)},
        )
    return GateResult("struktur.ketten", "Ketten vorhanden", "ok",
                      f"{len(ctx.erwartete_ketten)} Ketten gefunden: "
                      f"{', '.join(ctx.erwartete_ketten)}.")


@gate("struktur.fehlende_atome", "Keine fehlenden Schweratome")
def gate_missing_atoms(ctx: StructureContext) -> GateResult:
    fehlend = missing_atoms(ctx.structure)
    if fehlend:
        beispiel = "; ".join(f"{e['rest']}: {', '.join(e['fehlt'])}" for e in fehlend[:5])
        rest = f" und {len(fehlend) - 5} weitere" if len(fehlend) > 5 else ""
        return GateResult(
            "struktur.fehlende_atome", "Keine fehlenden Schweratome", "fehler",
            f"{len(fehlend)} Rest(e) haben fehlende Schweratome ({beispiel}{rest}). "
            "Ergänze sie mit PDBFixer, bevor du weitermachst.",
            {"reste": fehlend},
        )
    return GateResult("struktur.fehlende_atome", "Keine fehlenden Schweratome", "ok",
                      "Alle Standardnukleotide sind vollständig.")


@gate("struktur.unbekannte_reste", "Alle Reste bekannt")
def gate_unknown_residues(ctx: StructureContext) -> GateResult:
    unbekannt = unknown_residues(ctx.structure)
    erlaubt = {r.strip().upper() for r in ctx.erlaubte_sonderreste}
    if ctx.modifizierter_rest:
        erlaubt.add(ctx.modifizierter_rest.strip().upper())
    offen = [r for r in unbekannt if r.split()[0].upper() not in erlaubt]
    if offen:
        return GateResult(
            "struktur.unbekannte_reste", "Alle Reste bekannt", "fehler",
            f"Für diese Reste gibt es kein Standard-Kraftfeld: {', '.join(offen)}. "
            "Entweder entfernst du sie (struktur.entfernen), erklärst sie für behaltenswert "
            "(struktur.behalten_hetatm) oder du brauchst eigene Parameter.",
            {"reste": offen},
        )
    meldung = "Nur Standardreste." if not unbekannt else \
        f"Sonderreste, aber angekündigt: {', '.join(unbekannt)}."
    return GateResult("struktur.unbekannte_reste", "Alle Reste bekannt", "ok", meldung)


@gate("struktur.ladung", "Gesamtladung bestimmbar")
def gate_charge(ctx: StructureContext) -> GateResult:
    ladung = net_charge(ctx.structure)
    unbekannt = unknown_residues(ctx.structure)
    if unbekannt:
        return GateResult(
            "struktur.ladung", "Gesamtladung bestimmbar", "warnung",
            f"Aus den Standardresten ergibt sich {ladung:+d}. Die Ladung von "
            f"{', '.join(unbekannt)} ist dabei nicht enthalten und muss von Hand dazu.",
            {"ladung_standard": ladung, "offen": unbekannt},
        )
    return GateResult("struktur.ladung", "Gesamtladung bestimmbar", "ok",
                      f"Gesamtladung {ladung:+d}. So viele Gegenionen braucht die Box.",
                      {"ladung": ladung})


@gate("struktur.kollisionen", "Keine Atomüberlappungen")
def gate_clashes(ctx: StructureContext) -> GateResult:
    treffer = clashes(ctx.structure, ctx.einstellungen.kollision_toleranz_a)
    erlaubt = ctx.einstellungen.kollisionen_erlaubt
    if len(treffer) > erlaubt:
        schlimmste = treffer[0]
        return GateResult(
            "struktur.kollisionen", "Keine Atomüberlappungen", "fehler",
            f"{len(treffer)} Atompaar(e) überlappen, erlaubt sind {erlaubt}. "
            f"Am engsten: {schlimmste['a']} und {schlimmste['b']} mit "
            f"{schlimmste['abstand']} Å statt mindestens {schlimmste['grenze']} Å.",
            {"anzahl": len(treffer), "paare": treffer[:20]},
        )
    return GateResult("struktur.kollisionen", "Keine Atomüberlappungen", "ok",
                      "Kein Atompaar kommt sich zu nah.")


@gate("struktur.wasserstoffe", "Wasserstoffe vorhanden")
def gate_hydrogens(ctx: StructureContext) -> GateResult:
    info = inventory(ctx.structure)
    if info["wasserstoffe"] == 0:
        return GateResult(
            "struktur.wasserstoffe", "Wasserstoffe vorhanden", "warnung",
            "Die Struktur enthält keine Wasserstoffe. PDBFixer muss sie ergänzen, "
            "bevor das Kraftfeld zugewiesen wird.",
            {"wasserstoffe": 0},
        )
    return GateResult("struktur.wasserstoffe", "Wasserstoffe vorhanden", "ok",
                      f"{info['wasserstoffe']} Wasserstoffe in der Datei.",
                      {"wasserstoffe": info["wasserstoffe"]})


STRUCTURE_GATES: list[Gate] = [
    gate_chains,
    gate_missing_atoms,
    gate_unknown_residues,
    gate_charge,
    gate_clashes,
    gate_hydrogens,
]


def check_structure(ctx: StructureContext) -> GateReport:
    return run_gates(STRUCTURE_GATES, ctx)
