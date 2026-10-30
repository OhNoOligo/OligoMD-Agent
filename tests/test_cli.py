"""Kommandozeile: python -m oligomd …"""

import yaml

from oligomd.__main__ import main

from .conftest import build_pdb


def test_info(duplex_pdb, capsys):
    assert main(["info", str(duplex_pdb)]) == 0
    ausgabe = capsys.readouterr().out
    assert "Atome" in ausgabe and "Ladung" in ausgabe


def test_info_json_und_unbekanntes_modell(duplex_pdb, capsys):
    assert main(["info", str(duplex_pdb), "--json"]) == 0
    import json
    daten = json.loads(capsys.readouterr().out)
    assert daten["reste"]["dna"] == 5 and daten["modelle_in_datei"] == 1

    assert main(["info", str(duplex_pdb), "--modell", "4"]) == 1
    assert "FEHLER" in capsys.readouterr().err


def test_doctor_phase_0_laeuft_ohne_openmm(capsys):
    assert main(["doctor", "--phase", "0"]) == 0
    assert "Alles da, was Phase 0 braucht" in capsys.readouterr().out


def test_pruefe_bestanden_und_durchgefallen(tmp_path, capsys):
    (tmp_path / "config").mkdir()
    (tmp_path / "eingaben").mkdir()
    pdb = tmp_path / "eingaben" / "test.pdb"
    pdb.write_text(build_pdb([("DG", "A", 1), ("DA", "A", 2)]), encoding="utf-8")

    system = tmp_path / "config" / "system.yaml"
    system.write_text(yaml.safe_dump({
        "struktur": {"pdb": "eingaben/test.pdb", "ketten": ["A"]},
    }), encoding="utf-8")
    protokoll = tmp_path / "config" / "protocol.yaml"
    protokoll.write_text(yaml.safe_dump({
        "aequilibrierung": [{"name": "npt", "ps": 100, "restraints_kj_nm2": 0}],
        "produktion": {"ns": 1},
    }), encoding="utf-8")

    assert main(["pruefe", "--system", str(system), "--protokoll", str(protokoll)]) == 0
    assert "Bestanden." in capsys.readouterr().out

    # Kette, die es nicht gibt: Rückgabewert 2, kein Absturz
    system.write_text(yaml.safe_dump({
        "struktur": {"pdb": "eingaben/test.pdb", "ketten": ["A", "Q"]},
    }), encoding="utf-8")
    assert main(["pruefe", "--system", str(system), "--protokoll", str(protokoll)]) == 2
    assert "NICHT bestanden" in capsys.readouterr().out


def test_pruefe_meldet_kaputte_konfiguration(tmp_path, capsys):
    system = tmp_path / "system.yaml"
    system.write_text("struktur: {ketten: []}", encoding="utf-8")
    protokoll = tmp_path / "protocol.yaml"
    protokoll.write_text(yaml.safe_dump({
        "aequilibrierung": [{"name": "npt", "ps": 1, "restraints_kj_nm2": 0}],
        "produktion": {"ns": 1},
    }), encoding="utf-8")
    assert main(["pruefe", "--system", str(system), "--protokoll", str(protokoll)]) == 1
    assert "struktur" in capsys.readouterr().err
