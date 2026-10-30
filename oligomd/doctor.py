"""Umgebung prüfen: `python -m oligomd doctor`.

Das ist der Abnahmetest für Phase 0. Er sagt für jede Abhängigkeit, ob sie da
ist, welche Version läuft und wofür sie gebraucht wird. Nichts davon wird
installiert, nur berichtet.
"""

from __future__ import annotations

import importlib
import platform
import sys
from dataclasses import dataclass

from .tools import TOOLS, find_tool

#: Python-Pakete: Name, ab welcher Phase nötig, wofür
PACKAGES = [
    ("numpy", 0, "Zahlen"),
    ("rdkit", 0, "Strukturen lesen, Atomradien"),
    ("yaml", 0, "Konfiguration"),
    ("pydantic", 0, "Konfiguration prüfen"),
    ("openmm", 1, "die Simulation selbst"),
    ("openmmforcefields", 1, "AMBER-Kraftfelder in OpenMM"),
    ("pdbfixer", 1, "fehlende Atome und Wasserstoffe ergänzen"),
    ("MDAnalysis", 3, "Auswertung der Trajektorien"),
    ("claude_agent_sdk", 4, "Agenten-Schicht"),
]


@dataclass
class Finding:
    name: str
    status: str  # "ok" | "fehlt"
    info: str
    phase: int
    zweck: str


def _version(modul) -> str:
    for attribut in ("__version__", "version", "VERSION"):
        wert = getattr(modul, attribut, None)
        if isinstance(wert, str):
            return wert
    return "?"


def check_packages() -> list[Finding]:
    befunde = []
    for name, phase, zweck in PACKAGES:
        try:
            modul = importlib.import_module(name)
        except Exception:
            befunde.append(Finding(name, "fehlt", "nicht installiert", phase, zweck))
        else:
            befunde.append(Finding(name, "ok", _version(modul), phase, zweck))
    return befunde


def check_tools() -> list[Finding]:
    befunde = []
    for name, zweck in TOOLS.items():
        pfad = find_tool(name)
        phase = 2 if name != "tleap" else 1
        befunde.append(Finding(
            name, "ok" if pfad else "fehlt", str(pfad) if pfad else "nicht im PATH", phase, zweck,
        ))
    return befunde


def check_openmm_platforms() -> list[Finding]:
    """Welche Rechenplattformen OpenMM sieht. Ohne GPU wird es nichts."""
    try:
        import openmm
    except Exception:
        return [Finding("OpenMM-Plattformen", "fehlt",
                        "OpenMM ist nicht installiert", 1, "GPU-Beschleunigung")]
    namen = [openmm.Platform.getPlatform(i).getName()
             for i in range(openmm.Platform.getNumPlatforms())]
    schnell = [n for n in namen if n in ("CUDA", "HIP", "OpenCL")]
    return [Finding(
        "OpenMM-Plattformen",
        "ok" if schnell else "fehlt",
        ", ".join(namen) if namen else "keine",
        1,
        "CUDA wäre ideal; nur CPU ist für dieses System zu langsam",
    )]


def report(max_phase: int = 4) -> tuple[list[Finding], bool]:
    befunde = check_packages() + check_tools() + check_openmm_platforms()
    relevant = [b for b in befunde if b.phase <= max_phase]
    bereit = all(b.status == "ok" for b in relevant)
    return befunde, bereit


def als_text(befunde: list[Finding], bis_phase: int) -> str:
    zeilen = [
        f"Python   {sys.version.split()[0]}  ({platform.system()} {platform.machine()})",
        "",
        f"{'Status':8} {'Name':22} {'Phase':6} Info",
        "-" * 78,
    ]
    for b in sorted(befunde, key=lambda x: (x.phase, x.name.lower())):
        marke = "ok" if b.status == "ok" else ("FEHLT" if b.phase <= bis_phase else "offen")
        zeilen.append(f"{marke:8} {b.name:22} {b.phase:<6} {b.info}")
    zeilen.append("")

    fehlend = [b for b in befunde if b.status != "ok" and b.phase <= bis_phase]
    if fehlend:
        zeilen.append(f"Für Phase {bis_phase} fehlen noch: "
                      + ", ".join(b.name for b in fehlend))
        zeilen.append("")
        for b in fehlend:
            zeilen.append(f"  {b.name}: {b.zweck}")
        zeilen.append("")
        zeilen.append("Installation: conda env create -f environment.yml && conda activate oligomd")
        zeilen.append("modXNA getrennt, siehe README.")
    else:
        zeilen.append(f"Alles da, was Phase {bis_phase} braucht.")
    return "\n".join(zeilen)
