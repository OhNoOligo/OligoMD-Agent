"""Lauf-Ordner: eingefrorene Konfiguration, Logbuch, Berichte.

Der Lauf-Ordner ist die Wahrheit. Was hier liegt, beschreibt den Lauf
vollständig und wird nie nachträglich überschrieben. Ändert sich später die
Konfiguration im Repository, bleibt der Lauf davon unberührt.
"""

from __future__ import annotations

import json
import re
import shutil
from dataclasses import dataclass
from datetime import date, datetime, timezone
from pathlib import Path

from .config import RunConfig, load_run

RUNS = Path("runs")
KONFIG = "konfiguration"
LOG = "logbuch.txt"


class RunError(RuntimeError):
    pass


def slug(text: str, max_len: int = 40) -> str:
    text = text.strip().lower()
    for alt, neu in (("ä", "ae"), ("ö", "oe"), ("ü", "ue"), ("ß", "ss")):
        text = text.replace(alt, neu)
    text = re.sub(r"[^a-z0-9]+", "-", text).strip("-")
    return (text[:max_len].rstrip("-")) or "lauf"


def _jetzt() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


@dataclass
class Run:
    """Ein Lauf auf der Platte."""
    ordner: Path
    konfiguration: RunConfig

    # ---------- anlegen und öffnen ----------

    @classmethod
    def neu(cls, name: str, system_datei: Path, protokoll_datei: Path,
            basis: Path = RUNS) -> "Run":
        """Legt runs/<datum>_<name>/ an und friert die Konfiguration dort ein."""
        konfiguration = load_run(Path(system_datei), Path(protokoll_datei))
        basis = Path(basis)
        basis.mkdir(parents=True, exist_ok=True)
        grundname = f"{date.today():%Y-%m-%d}_{slug(name)}"
        ordner, n = basis / grundname, 2
        while ordner.exists():
            ordner = basis / f"{grundname}-{n}"
            n += 1
        (ordner / KONFIG).mkdir(parents=True)
        shutil.copyfile(system_datei, ordner / KONFIG / "system.yaml")
        shutil.copyfile(protokoll_datei, ordner / KONFIG / "protocol.yaml")

        lauf = cls(ordner, load_run(ordner / KONFIG / "system.yaml",
                                    ordner / KONFIG / "protocol.yaml"))
        lauf.log(f"Lauf '{name}' angelegt aus {system_datei} und {protokoll_datei}")
        return lauf

    @classmethod
    def oeffnen(cls, ordner: Path) -> "Run":
        ordner = Path(ordner)
        system = ordner / KONFIG / "system.yaml"
        if not system.is_file():
            raise RunError(
                f"{ordner} sieht nicht wie ein Lauf aus: {KONFIG}/system.yaml fehlt."
            )
        return cls(ordner, load_run(system, ordner / KONFIG / "protocol.yaml"))

    @classmethod
    def letzter(cls, basis: Path = RUNS) -> "Run":
        ordner = sorted((p for p in Path(basis).iterdir()
                         if p.is_dir() and (p / KONFIG / "system.yaml").is_file()),
                        key=lambda p: p.name)
        if not ordner:
            raise RunError(f"Unter {basis} liegt noch kein Lauf. Lege einen mit "
                           "'python -m oligomd vorbereiten --name <name>' an.")
        return cls.oeffnen(ordner[-1])

    # ---------- Pfade ----------

    @property
    def name(self) -> str:
        return self.ordner.name

    def pfad(self, *teile: str) -> Path:
        """Pfad im Lauf-Ordner; Zwischenordner werden angelegt."""
        ziel = self.ordner.joinpath(*teile)
        ziel.parent.mkdir(parents=True, exist_ok=True)
        return ziel

    def struktur_datei(self) -> Path:
        """Die Eingabestruktur, relativ zum Repository aufgelöst."""
        pdb = self.konfiguration.system.struktur.pdb
        if pdb.is_absolute():
            return pdb
        # Lauf-Ordner ist runs/<name>/, das Repository liegt zwei Ebenen höher
        kandidaten = [Path.cwd() / pdb, self.ordner.parent.parent / pdb]
        for kandidat in kandidaten:
            if kandidat.is_file():
                return kandidat
        raise RunError(
            f"Die Strukturdatei {pdb} wurde nicht gefunden. Gesucht wurde unter: "
            + ", ".join(str(k) for k in kandidaten)
        )

    # ---------- Logbuch und Berichte ----------

    def log(self, text: str) -> None:
        with (self.ordner / LOG).open("a", encoding="utf-8") as datei:
            datei.write(f"{_jetzt()}  {text}\n")

    def schreibe_bericht(self, name: str, daten: dict) -> Path:
        """Maschinenlesbarer Bericht eines Schritts, z. B. 'vorbereiten'."""
        ziel = self.pfad("berichte", f"{name}.json")
        ziel.write_text(
            json.dumps({"schritt": name, "zeitpunkt": _jetzt(), **daten},
                       indent=2, ensure_ascii=False),
            encoding="utf-8",
        )
        return ziel

    def lies_bericht(self, name: str) -> dict | None:
        pfad = self.ordner / "berichte" / f"{name}.json"
        if not pfad.is_file():
            return None
        try:
            return json.loads(pfad.read_text(encoding="utf-8", errors="replace"))
        except ValueError:
            return None

    def schritt_fertig(self, name: str) -> bool:
        bericht = self.lies_bericht(name)
        return bool(bericht and bericht.get("bestanden"))
