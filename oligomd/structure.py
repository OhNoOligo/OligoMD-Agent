"""PDB lesen und prüfen, ohne MD-Programme.

Nur rdkit und numpy, damit dieses Modul und seine Tests ohne OpenMM und
AmberTools laufen. Es beantwortet die Fragen, die vor jeder Simulation geklärt
sein müssen: Welche Reste sind da, fehlen Atome, stimmt die Ladung, überlappt
irgendwo etwas.
"""

from __future__ import annotations

import gzip
import xml.etree.ElementTree as ET
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

#: Namensraum der PDBML-Dateien von der RCSB (…/7PDU.xml.gz)
PDBML_NS = "{http://pdbml.pdb.org/schema/pdbx-v50.xsd}"
XSI_NIL = "{http://www.w3.org/2001/XMLSchema-instance}nil"


class StructureError(ValueError):
    pass


# ---------- Erwartete Atome ----------

# Rückgrat eines Nukleotids (PDB-Namen, Version 3). P/OP1/OP2 fehlen am 5'-Ende.
BACKBONE_5P = ("P", "OP1", "OP2")
BACKBONE = ("O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'")
RIBOSE_EXTRA = ("O2'",)

BASES = {
    "A": ("N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4"),
    "G": ("N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4"),
    "C": ("N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"),
    "T": ("N1", "C2", "O2", "N3", "C4", "O4", "C5", "C7", "C6"),
    "U": ("N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"),
}

DNA_RESIDUES = {"DA": "A", "DG": "G", "DC": "C", "DT": "T"}
RNA_RESIDUES = {"A": "A", "G": "G", "C": "C", "U": "U"}

# Formale Ladung häufiger Ionen (PDB-Restname)
ION_CHARGES = {"MG": 2, "CA": 2, "ZN": 2, "MN": 2, "NA": 1, "K": 1, "LI": 1,
               "CL": -1, "BR": -1, "IOD": -1}

WATER = {"HOH", "WAT", "H2O", "TIP3", "SOL"}

#: Elemente, die Wasserstoffbrücken bilden. Zwei davon dürfen sich nahe kommen:
#: eine Watson-Crick-Paarung liegt bei 2,7 bis 3,0 Å zwischen den Schweratomen.
HBOND_ELEMENTS = {"N", "O", "F", "S"}


def expected_atoms(resname: str, is_first_in_chain: bool) -> set[str] | None:
    """Erwartete Schweratome eines Standardnukleotids; None, wenn unbekannt."""
    resname = resname.strip().upper()
    if resname in DNA_RESIDUES:
        base, ribose = DNA_RESIDUES[resname], ()
    elif resname in RNA_RESIDUES:
        base, ribose = RNA_RESIDUES[resname], RIBOSE_EXTRA
    else:
        return None
    atoms = set(BACKBONE) | set(ribose) | set(BASES[base])
    if not is_first_in_chain:
        atoms |= set(BACKBONE_5P)
    return atoms


# ---------- Datenmodell ----------

@dataclass(frozen=True)
class Atom:
    index: int
    name: str
    element: str
    chain: str
    resnum: int
    resname: str
    icode: str
    xyz: tuple[float, float, float]

    @property
    def is_hydrogen(self) -> bool:
        return self.element == "H"


@dataclass
class Residue:
    chain: str
    number: int
    name: str
    icode: str = ""
    atoms: list[Atom] = field(default_factory=list)

    @property
    def key(self) -> str:
        return f"{self.chain}{self.number}{self.icode}".strip()

    @property
    def label(self) -> str:
        return f"{self.name} {self.chain}{self.number}{self.icode}".strip()

    @property
    def kind(self) -> str:
        name = self.name.strip().upper()
        if name in DNA_RESIDUES:
            return "dna"
        if name in RNA_RESIDUES:
            return "rna"
        if name in WATER:
            return "wasser"
        if name in ION_CHARGES:
            return "ion"
        return "sonstiges"

    def atom(self, name: str) -> Atom | None:
        return next((a for a in self.atoms if a.name == name), None)


@dataclass
class Structure:
    atoms: list[Atom]
    quelle: Path | None = None
    modell: int = 1

    @property
    def coords(self) -> np.ndarray:
        return np.array([a.xyz for a in self.atoms], dtype=float)

    @property
    def has_hydrogens(self) -> bool:
        return any(a.is_hydrogen for a in self.atoms)

    def residues(self) -> list[Residue]:
        out: dict[tuple, Residue] = {}
        for atom in self.atoms:
            key = (atom.chain, atom.resnum, atom.icode)
            residue = out.get(key)
            if residue is None:
                residue = out[key] = Residue(atom.chain, atom.resnum, atom.resname, atom.icode)
            residue.atoms.append(atom)
        return list(out.values())

    def chains(self) -> list[str]:
        return list(dict.fromkeys(a.chain for a in self.atoms))

    def residue(self, chain: str, number: int) -> Residue | None:
        return next((r for r in self.residues() if r.chain == chain and r.number == number), None)

    def subset(self, ketten: list[str] | None = None,
               ohne_restnamen: set[str] | None = None) -> "Structure":
        ohne = {n.strip().upper() for n in (ohne_restnamen or set())}
        atoms = [a for a in self.atoms
                 if (ketten is None or a.chain in ketten)
                 and a.resname.strip().upper() not in ohne]
        return Structure(atoms, self.quelle, self.modell)


# ---------- Lesen ----------

def _element_from(symbol: str, name: str) -> str:
    symbol = (symbol or "").strip()
    if symbol:
        return symbol.capitalize()
    # Fallback: erstes Zeichen des Atomnamens, Ziffern weg
    stripped = name.strip().lstrip("0123456789")
    return stripped[:1].upper() if stripped else "X"


def _suffixes(path: Path) -> list[str]:
    return [s.lower() for s in Path(path).suffixes]


def _read_text(path: Path) -> str:
    """Datei lesen, .gz wird unterwegs ausgepackt."""
    path = Path(path)
    try:
        if _suffixes(path)[-1:] == [".gz"]:
            with gzip.open(path, "rt", encoding="utf-8", errors="replace") as datei:
                return datei.read()
        return path.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        raise StructureError(f"{path} ließ sich nicht lesen: {exc}") from None


def is_pdbml(path: Path) -> bool:
    """PDBML (XML) statt klassischem PDB? Erkennt auch .xml.gz."""
    return ".xml" in _suffixes(path)


def models_in(path: Path) -> int:
    """Zahl der Modelle in der Datei; 1, wenn es keine Modelle gibt."""
    path = Path(path)
    if is_pdbml(path):
        nummern = {n for n, _ in _pdbml_atoms(path)}
        return max(len(nummern), 1)
    text = _read_text(path)
    count = sum(1 for line in text.splitlines() if line.startswith("MODEL "))
    return max(count, 1)


def read_structure(path: Path, modell: int = 1) -> Structure:
    """PDB oder PDBML lesen, auch als .gz. Der Dateiname entscheidet."""
    path = Path(path)
    if is_pdbml(path):
        return read_pdbml(path, modell)
    if ".cif" in _suffixes(path) or ".mmcif" in _suffixes(path):
        raise StructureError(
            f"{path.name} ist mmCIF. Das liest OligoMD noch nicht. Lade die Struktur "
            "von der RCSB als PDBML (.xml.gz) oder als PDB-Datei herunter."
        )
    return read_pdb(path, modell)


# ---------- PDBML (XML von der RCSB) ----------

def _text_of(element, tag: str) -> str | None:
    kind = element.find(f"{PDBML_NS}{tag}")
    if kind is None or kind.get(XSI_NIL) == "true":
        return None
    wert = (kind.text or "").strip()
    return wert or None


def _pdbml_atoms(path: Path):
    """(Modellnummer, atom_site-Element) für jedes Atom, speicherschonend."""
    path = Path(path)
    opener = gzip.open if _suffixes(path)[-1:] == [".gz"] else open
    try:
        with opener(path, "rb") as datei:
            for ereignis, element in ET.iterparse(datei, events=("end",)):
                if element.tag != f"{PDBML_NS}atom_site":
                    continue
                nummer = _text_of(element, "pdbx_PDB_model_num") or "1"
                yield int(nummer), element
                element.clear()
    except ET.ParseError as exc:
        raise StructureError(f"{path.name} ist kein gültiges XML: {exc}") from None
    except OSError as exc:
        raise StructureError(f"{path} ließ sich nicht lesen: {exc}") from None


def read_pdbml(path: Path, modell: int = 1) -> Structure:
    """Ein Modell aus einer PDBML-Datei (…/7PDU.xml.gz von der RCSB).

    Es gelten die auth_*-Felder, also die Nummerierung aus der Veröffentlichung,
    nicht die internen label_*-Felder.
    """
    path = Path(path)
    atoms: list[Atom] = []
    gesehen: set[int] = set()
    index = 0
    for nummer, element in _pdbml_atoms(path):
        gesehen.add(nummer)
        if nummer != modell:
            continue
        altloc = _text_of(element, "label_alt_id")
        if altloc not in (None, "A"):
            continue  # nur die erste Alternativlage
        try:
            xyz = (float(_text_of(element, "Cartn_x")),
                   float(_text_of(element, "Cartn_y")),
                   float(_text_of(element, "Cartn_z")))
        except (TypeError, ValueError):
            raise StructureError(
                f"{path.name}: Atom {_text_of(element, 'id')} hat keine lesbaren Koordinaten."
            ) from None
        name = _text_of(element, "auth_atom_id") or _text_of(element, "label_atom_id") or "?"
        resnum = _text_of(element, "auth_seq_id") or _text_of(element, "label_seq_id")
        atoms.append(Atom(
            index=index,
            name=name,
            element=(_text_of(element, "type_symbol") or name[:1]).capitalize(),
            chain=_text_of(element, "auth_asym_id") or _text_of(element, "label_asym_id") or " ",
            resnum=int(resnum) if resnum else 0,
            resname=(_text_of(element, "auth_comp_id")
                     or _text_of(element, "label_comp_id") or "UNK"),
            icode=_text_of(element, "pdbx_PDB_ins_code") or "",
            xyz=xyz,
        ))
        index += 1

    if not gesehen:
        raise StructureError(f"{path.name} enthält keine atom_site-Einträge.")
    if modell not in gesehen:
        raise StructureError(
            f"{path.name} enthält {len(gesehen)} Modell(e) "
            f"({min(gesehen)} bis {max(gesehen)}), angefragt war Nummer {modell}."
        )
    if not atoms:
        raise StructureError(f"{path.name} enthält keine Atome für Modell {modell}.")
    return Structure(atoms, path, modell)


# ---------- Klassisches PDB ----------

def read_pdb(path: Path, modell: int = 1) -> Structure:
    """Liest genau ein Modell einer PDB-Datei. ATOM und HETATM, mit Wasserstoffen."""
    path = Path(path)
    text = _read_text(path)

    gesamt = models_in(path)
    if modell < 1 or modell > gesamt:
        raise StructureError(
            f"{path} enthält {gesamt} Modell(e), angefragt war Nummer {modell}."
        )

    atoms: list[Atom] = []
    aktuelles = 1
    index = 0
    for line in text.splitlines():
        if line.startswith("MODEL "):
            try:
                aktuelles = int(line[10:14])
            except ValueError:
                aktuelles += 1
            continue
        if line.startswith("ENDMDL"):
            continue
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if gesamt > 1 and aktuelles != modell:
            continue
        if len(line) < 54:
            raise StructureError(f"Zu kurze Koordinatenzeile in {path.name}: {line!r}")
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            raise StructureError(f"Unlesbare Koordinaten in {path.name}: {line!r}") from None
        altloc = line[16]
        if altloc not in (" ", "A"):
            continue  # nur die erste Alternativlage
        name = line[12:16].strip()
        atoms.append(Atom(
            index=index,
            name=name,
            element=_element_from(line[76:78] if len(line) >= 78 else "", name),
            chain=line[21].strip() or " ",
            resnum=int(line[22:26]),
            resname=line[17:20].strip(),
            icode=line[26].strip(),
            xyz=xyz,
        ))
        index += 1

    if not atoms:
        raise StructureError(f"{path} enthält keine Atome für Modell {modell}.")
    return Structure(atoms, path, modell)


# ---------- Auswerten ----------

def inventory(structure: Structure) -> dict:
    """Übersicht: Ketten, Reste je Art, Ionen, Wasser."""
    nach_art: dict[str, int] = defaultdict(int)
    ionen: dict[str, int] = defaultdict(int)
    ketten: dict[str, int] = defaultdict(int)
    unbekannt: list[str] = []
    for residue in structure.residues():
        nach_art[residue.kind] += 1
        ketten[residue.chain] += 1
        if residue.kind == "ion":
            ionen[residue.name.strip().upper()] += 1
        elif residue.kind == "sonstiges":
            unbekannt.append(residue.label)
    return {
        "atome": len(structure.atoms),
        "wasserstoffe": sum(1 for a in structure.atoms if a.is_hydrogen),
        "reste": dict(nach_art),
        "ketten": dict(ketten),
        "ionen": dict(ionen),
        "unbekannte_reste": unbekannt,
    }


def net_charge(structure: Structure) -> int:
    """Gesamtladung: ein Phosphat = −1, dazu die Ionen.

    Für Standard-Nukleinsäuren bei pH 7 verlässlich. Unbekannte Reste zählen 0;
    deren Ladung muss die Prüfung an anderer Stelle erfragen.
    """
    ladung = 0
    for residue in structure.residues():
        if residue.kind == "ion":
            ladung += ION_CHARGES.get(residue.name.strip().upper(), 0)
        elif residue.kind in ("dna", "rna") and residue.atom("P") is not None:
            ladung -= 1
    return ladung


def first_residues(structure: Structure) -> set[tuple[str, int, str]]:
    """Der jeweils erste Rest jeder Kette (dort fehlt das 5'-Phosphat planmäßig)."""
    erste: dict[str, Residue] = {}
    for residue in structure.residues():
        if residue.kind not in ("dna", "rna"):
            continue
        vorhanden = erste.get(residue.chain)
        if vorhanden is None or residue.number < vorhanden.number:
            erste[residue.chain] = residue
    return {(r.chain, r.number, r.icode) for r in erste.values()}


def missing_atoms(structure: Structure) -> list[dict]:
    """Fehlende Schweratome in Standardnukleotiden, Rest für Rest."""
    erste = first_residues(structure)
    fehlend = []
    for residue in structure.residues():
        erwartet = expected_atoms(residue.name, (residue.chain, residue.number, residue.icode) in erste)
        if erwartet is None:
            continue
        vorhanden = {a.name for a in residue.atoms if not a.is_hydrogen}
        luecke = sorted(erwartet - vorhanden)
        if luecke:
            fehlend.append({"rest": residue.label, "fehlt": luecke})
    return fehlend


def unknown_residues(structure: Structure) -> list[str]:
    """Reste, für die es kein Standard-Kraftfeld gibt (Modifikationen, Liganden)."""
    return sorted({r.label for r in structure.residues() if r.kind == "sonstiges"})


# ---------- Kollisionen ----------

def _vdw_radii(structure: Structure) -> np.ndarray:
    table = Chem.GetPeriodicTable()
    radien = []
    for atom in structure.atoms:
        try:
            radien.append(table.GetRvdw(atom.element))
        except RuntimeError:
            radien.append(1.7)
    return np.array(radien)


def _neighbor_pairs(coords: np.ndarray, cutoff: float):
    """Paare näher als cutoff, über ein Raumgitter statt aller Paare."""
    if len(coords) < 2:
        return
    zellen: dict[tuple[int, int, int], list[int]] = defaultdict(list)
    keys = np.floor(coords / cutoff).astype(int)
    for i, key in enumerate(map(tuple, keys)):
        zellen[key].append(i)
    nachbarn = [(dx, dy, dz) for dx in (-1, 0, 1) for dy in (-1, 0, 1) for dz in (-1, 0, 1)]
    for key, indices in zellen.items():
        for offset in nachbarn:
            other = zellen.get((key[0] + offset[0], key[1] + offset[1], key[2] + offset[2]))
            if not other:
                continue
            for i in indices:
                for j in other:
                    if j > i:
                        yield i, j


def clashes(structure: Structure, toleranz: float = 0.4, nur_schwer: bool = True,
            max_treffer: int = 50, hbond_min: float = 2.2) -> list[dict]:
    """Atompaare, die sich zu nahe kommen.

    Zwei Regeln, weil eine nicht reicht:
    - Können beide Atome eine Wasserstoffbrücke bilden (N, O, F, S), gilt erst
      ein Abstand unter `hbond_min` als Kollision. Sonst würde jede
      Watson-Crick-Paarung als Fehler gemeldet.
    - Sonst: die Summe der VdW-Radien minus `toleranz`.

    Atome desselben oder benachbarter Reste werden übersprungen: Dort sind kurze
    Abstände normal (kovalente Bindungen).
    """
    atoms = [a for a in structure.atoms if not (nur_schwer and a.is_hydrogen)]
    if len(atoms) < 2:
        return []
    teil = Structure(atoms, structure.quelle, structure.modell)
    coords = teil.coords
    radien = _vdw_radii(teil)
    cutoff = float(2 * radien.max() - toleranz)

    treffer = []
    for i, j in _neighbor_pairs(coords, cutoff):
        a, b = atoms[i], atoms[j]
        if a.chain == b.chain and abs(a.resnum - b.resnum) <= 1:
            continue
        wasserstoffbruecke = a.element in HBOND_ELEMENTS and b.element in HBOND_ELEMENTS
        grenze = hbond_min if wasserstoffbruecke else radien[i] + radien[j] - toleranz
        d = float(np.linalg.norm(coords[i] - coords[j]))
        if d < grenze:
            treffer.append({
                "a": f"{a.resname} {a.chain}{a.resnum} {a.name}",
                "b": f"{b.resname} {b.chain}{b.resnum} {b.name}",
                "abstand": round(d, 2),
                "grenze": round(float(grenze), 2),
                "ueberlappung": round(float(grenze) - d, 2),
                "art": "wasserstoffbruecke" if wasserstoffbruecke else "kontakt",
            })
    treffer.sort(key=lambda t: -t["ueberlappung"])
    return treffer[:max_treffer]
