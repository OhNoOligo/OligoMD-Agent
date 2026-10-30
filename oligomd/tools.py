"""Externe Programme finden und aufrufen.

Immer mit vollem Pfad und Argumentliste, nie über die Shell. Welche Programme
vorhanden sind, beantwortet `find_tool`; ob die Umgebung insgesamt taugt,
beantwortet `oligomd.doctor`.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

#: Umgebungsvariable, mit der ein Pfad erzwungen werden kann, z. B. OLIGOMD_TLEAP
ENV_PREFIX = "OLIGOMD_"

#: Programme, die die Pipeline benutzt, mit der Phase, ab der sie nötig sind
TOOLS = {
    "tleap": "Kraftfeld zuweisen und Topologie bauen (AmberTools)",
    "antechamber": "Ladungen und Atomtypen für neue Reste (AmberTools)",
    "parmchk2": "fehlende Bindungsparameter ergänzen (AmberTools)",
    "pdb4amber": "PDB für Amber aufräumen (AmberTools)",
    "cpptraj": "Trajektorien verarbeiten (AmberTools)",
    "modxna": "modifizierte Nukleotide parametrisieren (modXNA)",
}


class ToolError(RuntimeError):
    pass


def find_tool(name: str) -> Path | None:
    """Pfad zum Programm: erst OLIGOMD_<NAME>, dann PATH."""
    override = os.environ.get(f"{ENV_PREFIX}{name.upper()}")
    if override:
        pfad = Path(override)
        return pfad if pfad.is_file() else None
    gefunden = shutil.which(name)
    return Path(gefunden) if gefunden else None


def require_tool(name: str) -> Path:
    pfad = find_tool(name)
    if pfad is None:
        zweck = TOOLS.get(name, "")
        hinweis = f" ({zweck})" if zweck else ""
        raise ToolError(
            f"{name}{hinweis} wurde nicht gefunden. Aktiviere die conda-Umgebung "
            f"(conda activate oligomd) oder setze {ENV_PREFIX}{name.upper()} auf den vollen Pfad."
        )
    return pfad


@dataclass
class ToolRun:
    befehl: list[str]
    rueckgabe: int
    ausgabe: str
    arbeitsordner: Path

    @property
    def erfolg(self) -> bool:
        return self.rueckgabe == 0


def run_tool(name: str, args: list[str], cwd: Path, eingabe: str | None = None,
             timeout: float = 3600) -> ToolRun:
    """Startet ein Programm und gibt stdout und stderr zusammen zurück."""
    pfad = require_tool(name)
    cwd = Path(cwd)
    cwd.mkdir(parents=True, exist_ok=True)
    try:
        proc = subprocess.run(
            [str(pfad), *args], cwd=cwd, input=eingabe, capture_output=True,
            text=True, timeout=timeout, shell=False, encoding="utf-8", errors="replace",
        )
    except subprocess.TimeoutExpired:
        raise ToolError(f"{name} hat nach {int(timeout)} s nicht geantwortet.") from None
    except OSError as exc:
        raise ToolError(f"{name} ließ sich nicht starten: {exc}") from None
    ausgabe = (proc.stdout or "") + (proc.stderr or "")
    (cwd / f"{name}.log").write_text(ausgabe, encoding="utf-8")
    return ToolRun([str(pfad), *args], proc.returncode, ausgabe, cwd)
