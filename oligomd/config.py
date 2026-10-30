"""Konfiguration laden und prüfen.

Zwei Dateien beschreiben einen Lauf vollständig:
  system.yaml    Was gerechnet wird: Struktur, Ketten, Kraftfeld, Box, Ionen.
  protocol.yaml  Wie gerechnet wird: Minimierung, Äquilibrierung, Produktion.

Beide werden beim Start in den Lauf-Ordner kopiert. Ab dann gilt die Kopie:
Ändert jemand später die Datei im Repository, bleibt der Lauf davon unberührt.
"""

from __future__ import annotations

from pathlib import Path
from typing import Literal, Optional

import yaml
from pydantic import BaseModel, Field, ValidationError, field_validator, model_validator


class ConfigError(ValueError):
    """Verständliche Meldung für die Oberfläche und die Kommandozeile."""


# ---------- system.yaml ----------

class StructureConfig(BaseModel):
    pdb: Path = Field(description="PDB-Datei, relativ zum Projektordner oder absolut")
    modell: int = Field(1, ge=1, description="Bei NMR-Strukturen: welches Modell (1-basiert)")
    ketten: list[str] = Field(min_length=1, description="Ketten, die übernommen werden")
    entfernen: list[str] = Field(default_factory=list, description="Restnamen, die wegfallen, z. B. HOH")
    behalten_hetatm: list[str] = Field(default_factory=list, description="HETATM-Reste, die bleiben, z. B. MG")

    @field_validator("ketten")
    @classmethod
    def _unique(cls, value: list[str]) -> list[str]:
        if len(set(value)) != len(value):
            raise ValueError("Jede Kette darf nur einmal vorkommen.")
        return value


class ForceFieldConfig(BaseModel):
    dna: str = Field("OL21", description="DNA-Kraftfeld, z. B. OL21 oder bsc1")
    rna: str = Field("OL3", description="RNA-Kraftfeld, z. B. OL3")
    wasser: Literal["OPC", "TIP3P", "TIP4P-Ew", "SPC/E"] = "OPC"
    ionen: str = Field("Joung-Cheatham", description="Ionensatz, muss zum Wassermodell passen")
    zusatz_frcmod: list[Path] = Field(default_factory=list, description="eigene frcmod-Dateien")
    zusatz_lib: list[Path] = Field(default_factory=list, description="eigene Bibliotheken (.lib/.off)")


class BoxConfig(BaseModel):
    typ: Literal["oktaeder", "wuerfel"] = "oktaeder"
    puffer_nm: float = Field(1.2, ge=0.8, le=3.0, description="Abstand Solut zur Boxwand")


class ExtraIon(BaseModel):
    """Ein Ion, das zusätzlich in die Lösung kommt, z. B. katalytisches Mg2+.

    Entweder eine feste Zahl oder eine Konzentration. Die Konzentration ist der
    ehrlichere Weg, weil sie nicht von der Boxgröße abhängt; die Zahl ergibt
    sich dann aus der Wassermenge.
    """
    ion: str
    anzahl: Optional[int] = Field(None, ge=1)
    mol_l: Optional[float] = Field(None, gt=0, le=1.0)

    @model_validator(mode="after")
    def _genau_eins(self) -> "ExtraIon":
        if (self.anzahl is None) == (self.mol_l is None):
            raise ValueError(
                f"Ion {self.ion}: Gib entweder 'anzahl' oder 'mol_l' an, nicht beides und nicht keins."
            )
        return self


class IonConfig(BaseModel):
    neutralisieren_mit: str = Field("K+", description="Gegenion für die Ladung des Rückgrats")
    salz_mol_l: float = Field(0.15, ge=0.0, le=2.0)
    zusaetzlich: list[ExtraIon] = Field(default_factory=list,
                                        description="z. B. katalytisches Mg2+, bewusst gesetzt")


class ResidueRef(BaseModel):
    kette: str
    nummer: int


class ModificationConfig(BaseModel):
    rest: ResidueRef = Field(description="Welcher Rest ersetzt wird, z. B. Kette A, Nummer 23")
    definition: Path = Field(description="YAML unter config/residues/ mit der Beschreibung")
    name: str = Field(description="Dreibuchstaben-Restname, muss zur Bibliothek passen",
                      min_length=1, max_length=4)


class SystemConfig(BaseModel):
    struktur: StructureConfig
    kraftfeld: ForceFieldConfig = Field(default_factory=ForceFieldConfig)
    box: BoxConfig = Field(default_factory=BoxConfig)
    ionen: IonConfig = Field(default_factory=IonConfig)
    modifikation: Optional[ModificationConfig] = None

    @model_validator(mode="after")
    def _modification_in_chain(self) -> "SystemConfig":
        if self.modifikation and self.modifikation.rest.kette not in self.struktur.ketten:
            raise ValueError(
                f"Die Modifikation sitzt in Kette {self.modifikation.rest.kette}, "
                f"die aber nicht in struktur.ketten steht ({', '.join(self.struktur.ketten)})."
            )
        return self


# ---------- protocol.yaml ----------

class MinimizationConfig(BaseModel):
    schritte: int = Field(5000, ge=100)


class EquilibrationStage(BaseModel):
    name: str
    ps: float = Field(gt=0)
    temperatur_k: float = Field(300.0, gt=0)
    druck_bar: Optional[float] = Field(None, description="leer = NVT, sonst NPT")
    restraints_kj_nm2: float = Field(0.0, ge=0, description="Haltekraft auf das Solut")


class ProductionConfig(BaseModel):
    ns: float = Field(gt=0)
    zeitschritt_fs: float = Field(4.0, gt=0, le=5.0)
    hmr: bool = Field(True, description="Wasserstoff-Massen-Umverteilung; nötig ab 4 fs")
    temperatur_k: float = Field(300.0, gt=0)
    druck_bar: float = 1.0
    speichern_ps: float = Field(10.0, gt=0)

    @model_validator(mode="after")
    def _hmr_needed(self) -> "ProductionConfig":
        if self.zeitschritt_fs > 2.5 and not self.hmr:
            raise ValueError(
                f"Zeitschritt {self.zeitschritt_fs} fs ohne HMR ist nicht stabil. "
                "Setze hmr: true oder höchstens 2 fs."
            )
        return self


class ProtocolConfig(BaseModel):
    minimierung: MinimizationConfig = Field(default_factory=MinimizationConfig)
    aequilibrierung: list[EquilibrationStage] = Field(min_length=1)
    produktion: ProductionConfig
    zufallsstartwert: int = Field(42, description="damit zwei Läufe gleich ausgehen")
    plattform: Literal["auto", "CUDA", "HIP", "OpenCL", "CPU"] = "auto"

    @model_validator(mode="after")
    def _ends_unrestrained(self) -> "ProtocolConfig":
        if self.aequilibrierung[-1].restraints_kj_nm2 > 0:
            raise ValueError(
                "Die letzte Stufe der Äquilibrierung muss ohne Restraints laufen, "
                "sonst startet die Produktion aus einer festgehaltenen Struktur."
            )
        return self


# ---------- Laden ----------

def _load_yaml(path: Path) -> dict:
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        raise ConfigError(f"{path} ließ sich nicht lesen: {exc}") from None
    try:
        data = yaml.safe_load(text)
    except yaml.YAMLError as exc:
        raise ConfigError(f"{path} ist kein gültiges YAML: {exc}") from None
    if not isinstance(data, dict):
        raise ConfigError(f"{path} muss eine Zuordnung von Schlüsseln zu Werten enthalten.")
    return data


def _pretty(path: Path, exc: ValidationError) -> str:
    zeilen = []
    for error in exc.errors():
        wo = ".".join(str(p) for p in error["loc"]) or "(Wurzel)"
        zeilen.append(f"  {wo}: {error['msg']}")
    return f"{path} ist nicht gültig:\n" + "\n".join(zeilen)


def load_system(path: Path) -> SystemConfig:
    try:
        return SystemConfig.model_validate(_load_yaml(path))
    except ValidationError as exc:
        raise ConfigError(_pretty(path, exc)) from None


def load_protocol(path: Path) -> ProtocolConfig:
    try:
        return ProtocolConfig.model_validate(_load_yaml(path))
    except ValidationError as exc:
        raise ConfigError(_pretty(path, exc)) from None


class RunConfig(BaseModel):
    """Beides zusammen, plus die Pfade, aus denen es kam."""
    system: SystemConfig
    protokoll: ProtocolConfig
    system_datei: Path
    protokoll_datei: Path

    @property
    def gesamtdauer_ns(self) -> float:
        aequi = sum(s.ps for s in self.protokoll.aequilibrierung) / 1000.0
        return round(aequi + self.protokoll.produktion.ns, 3)


def load_run(system_datei: Path, protokoll_datei: Path) -> RunConfig:
    return RunConfig(
        system=load_system(system_datei),
        protokoll=load_protocol(protokoll_datei),
        system_datei=system_datei,
        protokoll_datei=protokoll_datei,
    )
