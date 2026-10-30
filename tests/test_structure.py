"""Struktur lesen und bewerten, ohne OpenMM und ohne AmberTools."""

import pytest

from oligomd.structure import (
    StructureError, clashes, expected_atoms, inventory, missing_atoms, models_in,
    net_charge, read_pdb, unknown_residues,
)

from .conftest import build_pdb, pdb_line


def test_reads_atoms_chains_and_residues(duplex_pdb):
    s = read_pdb(duplex_pdb)
    info = inventory(s)
    assert info["reste"]["dna"] == 5
    assert info["ketten"] == {"A": 3, "B": 2}
    assert s.chains() == ["A", "B"]
    assert s.residue("A", 2).name == "DA"
    assert s.residue("A", 9) is None


def test_first_residue_has_no_phosphate(duplex_pdb):
    s = read_pdb(duplex_pdb)
    assert s.residue("A", 1).atom("P") is None      # 5'-Ende, planmäßig
    assert s.residue("A", 2).atom("P") is not None
    assert missing_atoms(s) == []                   # und das gilt nicht als Lücke


def test_net_charge_counts_phosphates_and_ions(tmp_path):
    text = build_pdb([("DG", "A", 1), ("DA", "A", 2), ("DC", "A", 3)],
                     extra=[pdb_line(900, "MG", "MG", "A", 101, (80.0, 80.0, 80.0), "MG")])
    pfad = tmp_path / "s.pdb"
    pfad.write_text(text, encoding="utf-8")
    s = read_pdb(pfad)
    # zwei Phosphate (Rest 2 und 3) = -2, ein Mg = +2
    assert net_charge(s) == 0
    assert inventory(s)["ionen"] == {"MG": 1}


def test_missing_heavy_atom_is_found(tmp_path):
    text = build_pdb([("DG", "A", 1), ("DG", "A", 2)],
                     weglassen={("A", 2): ("O6", "N7")})
    pfad = tmp_path / "luecke.pdb"
    pfad.write_text(text, encoding="utf-8")
    fehlend = missing_atoms(read_pdb(pfad))
    assert len(fehlend) == 1
    assert fehlend[0]["rest"].startswith("DG A2")
    assert fehlend[0]["fehlt"] == ["N7", "O6"]


def test_unknown_residue_is_reported(tmp_path):
    text = build_pdb([("DG", "A", 1)],
                     extra=[pdb_line(900, "C1", "G6X", "A", 2, (60.0, 60.0, 60.0), "C")])
    pfad = tmp_path / "mod.pdb"
    pfad.write_text(text, encoding="utf-8")
    s = read_pdb(pfad)
    assert unknown_residues(s) == ["G6X A2"]
    assert missing_atoms(s) == []  # unbekannte Reste werden nicht auf Vollständigkeit geprüft


def test_models(tmp_path):
    text = build_pdb([("DG", "A", 1)], modelle=3)
    pfad = tmp_path / "nmr.pdb"
    pfad.write_text(text, encoding="utf-8")
    assert models_in(pfad) == 3
    assert read_pdb(pfad, 2).modell == 2
    with pytest.raises(StructureError, match="3 Modell"):
        read_pdb(pfad, 7)


def test_clean_structure_has_no_clashes(duplex_pdb):
    assert clashes(read_pdb(duplex_pdb)) == []


def test_clash_is_found(duplex_pdb, tmp_path):
    # Ein fremdes Atom 0,5 Å neben ein Atom von Kette A legen
    ziel = read_pdb(duplex_pdb).residue("A", 1).atom("C1'")
    x, y, z = ziel.xyz
    text = duplex_pdb.read_text(encoding="utf-8").replace(
        "END\n", pdb_line(900, "C1", "LIG", "C", 1, (x + 0.5, y, z), "C") + "\nEND\n")
    pfad = tmp_path / "clash.pdb"
    pfad.write_text(text, encoding="utf-8")

    treffer = clashes(read_pdb(pfad))
    assert treffer, "Die Überlappung hätte auffallen müssen"
    assert all(t["abstand"] < t["grenze"] for t in treffer)
    beteiligt = {t["a"] for t in treffer} | {t["b"] for t in treffer}
    assert any(e.startswith("LIG C1") for e in beteiligt)
    # Ohne die fremde Kette ist wieder alles sauber
    assert clashes(read_pdb(pfad).subset(ketten=["A", "B"])) == []


def test_subset_filters_chains_and_residue_names(duplex_pdb):
    s = read_pdb(duplex_pdb)
    nur_a = s.subset(ketten=["A"])
    assert nur_a.chains() == ["A"] and len(nur_a.residues()) == 3
    ohne_dg = s.subset(ohne_restnamen={"DG"})
    assert all(r.name != "DG" for r in ohne_dg.residues())


def test_expected_atoms_table():
    mit_p = expected_atoms("DG", is_first_in_chain=False)
    ohne_p = expected_atoms("DG", is_first_in_chain=True)
    assert "P" in mit_p and "P" not in ohne_p
    assert "O6" in mit_p and "O2'" not in mit_p      # DNA: keine 2'-OH
    assert "O2'" in expected_atoms("G", is_first_in_chain=False)   # RNA: schon
    assert expected_atoms("G6X", is_first_in_chain=False) is None  # unbekannt
