"""Gegen echte Koordinaten prüfen, nicht nur gegen gebaute Testfälle.

tests/data/7pdu_paar.xml.gz enthält ein Watson-Crick-Paar aus der NMR-Struktur
7PDU (Modell 1, DG A28 und C B39), unverändert von der RCSB. Daten aus der
Protein Data Bank sind gemeinfrei (CC0).
"""

from pathlib import Path

from oligomd.gates import StructureContext, check_structure
from oligomd.structure import clashes, inventory, missing_atoms, net_charge, read_structure

PAAR = Path(__file__).parent / "data" / "7pdu_paar.xml.gz"


def test_liest_echte_pdbml_datei():
    s = read_structure(PAAR)
    info = inventory(s)
    assert info["reste"] == {"dna": 1, "rna": 1}
    assert [r.label for r in s.residues()] == ["DG A28", "C B39"]
    assert info["wasserstoffe"] > 0
    # Das DG trägt sein O6 - genau das Atom, das später modifiziert wird
    assert s.residue("A", 28).atom("O6") is not None
    # Das RNA-Cytidin hat die 2'-OH-Gruppe, die DNA nicht
    assert s.residue("B", 39).atom("O2'") is not None
    assert s.residue("A", 28).atom("O2'") is None


def test_watson_crick_paar_ist_keine_kollision():
    """Ohne die Wasserstoffbrücken-Regel würde jede Basenpaarung als Fehler gelten."""
    s = read_structure(PAAR)
    assert clashes(s) == []
    # Zum Beleg: ohne die Regel meldet dieselbe Struktur Dutzende "Kollisionen"
    assert len(clashes(s, hbond_min=99.0, max_treffer=500)) > 10


def test_echte_reste_sind_vollstaendig_und_geladen():
    s = read_structure(PAAR)
    assert missing_atoms(s) == []
    assert net_charge(s) == -2  # beide Reste tragen ein Phosphat


def test_pruefpunkte_lassen_echte_struktur_durch():
    bericht = check_structure(StructureContext(
        structure=read_structure(PAAR), erwartete_ketten=["A", "B"]))
    assert bericht.ok, bericht.als_text()
    assert not bericht.warnungen, bericht.als_text()
