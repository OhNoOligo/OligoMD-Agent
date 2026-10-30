"""Prüfpunkte: Sie müssen durchlassen, was sauber ist, und anhalten, was es nicht ist."""

import pytest

from oligomd.gates import (
    GateReport, GateResult, GateSettings, StructureContext, check_structure, gate, run_gates,
)
from oligomd.structure import read_pdb

from .conftest import build_pdb, pdb_line


def ctx_fuer(pfad, ketten=("A", "B"), **kwargs) -> StructureContext:
    return StructureContext(structure=read_pdb(pfad), erwartete_ketten=list(ketten), **kwargs)


def test_sauberer_duplex_besteht(duplex_pdb):
    bericht = check_structure(ctx_fuer(duplex_pdb))
    assert bericht.ok, bericht.als_text()
    assert not bericht.fehler
    # Ohne Wasserstoffe gibt es eine Warnung, aber die blockiert nicht
    assert [w.id for w in bericht.warnungen] == ["struktur.wasserstoffe"]


def test_fehlende_kette_blockiert(duplex_pdb):
    bericht = check_structure(ctx_fuer(duplex_pdb, ketten=("A", "B", "Z")))
    assert not bericht.ok
    assert "Z" in bericht.fehler[0].meldung


def test_fehlendes_atom_blockiert_und_nennt_den_rest(tmp_path):
    text = build_pdb([("DG", "A", 1), ("DG", "A", 2)], weglassen={("A", 2): ("O6",)})
    pfad = tmp_path / "luecke.pdb"
    pfad.write_text(text, encoding="utf-8")
    bericht = check_structure(ctx_fuer(pfad, ketten=("A",)))
    assert not bericht.ok
    fehler = next(f for f in bericht.fehler if f.id == "struktur.fehlende_atome")
    assert "DG A2" in fehler.meldung and "O6" in fehler.meldung


def test_unbekannter_rest_blockiert_bis_er_angekuendigt_ist(tmp_path):
    text = build_pdb([("DG", "A", 1)],
                     extra=[pdb_line(900, "C1", "G6X", "A", 2, (500.0, 500.0, 500.0), "C")])
    pfad = tmp_path / "mod.pdb"
    pfad.write_text(text, encoding="utf-8")

    bericht = check_structure(ctx_fuer(pfad, ketten=("A",)))
    assert not bericht.ok
    assert "G6X" in next(f for f in bericht.fehler if f.id == "struktur.unbekannte_reste").meldung

    # Als Modifikation angekündigt: geht durch, aber die Ladung bleibt offen
    bericht = check_structure(ctx_fuer(pfad, ketten=("A",), modifizierter_rest="G6X"))
    assert bericht.ok
    ladung = next(r for r in bericht.ergebnisse if r.id == "struktur.ladung")
    assert ladung.status == "warnung" and "G6X" in ladung.meldung


def test_ladung_wird_berichtet(duplex_pdb):
    bericht = check_structure(ctx_fuer(duplex_pdb))
    ladung = next(r for r in bericht.ergebnisse if r.id == "struktur.ladung")
    assert ladung.status == "ok"
    assert ladung.details["ladung"] == -3   # 5 Reste, 2 davon am 5'-Ende ohne Phosphat


def test_kollision_blockiert_und_toleranz_wirkt(duplex_pdb, tmp_path):
    ziel = read_pdb(duplex_pdb).residue("A", 1).atom("C1'")
    x, y, z = ziel.xyz
    text = duplex_pdb.read_text(encoding="utf-8").replace(
        "END\n", pdb_line(900, "C1", "LIG", "C", 1, (x + 2.6, y, z), "C") + "\nEND\n")
    pfad = tmp_path / "clash.pdb"
    pfad.write_text(text, encoding="utf-8")

    streng = check_structure(ctx_fuer(pfad, ketten=("A", "B", "C"),
                                      erlaubte_sonderreste={"LIG"}))
    assert not streng.ok
    assert "Å" in next(f for f in streng.fehler if f.id == "struktur.kollisionen").meldung

    # Grosszuegigere Toleranz laesst denselben Abstand durch
    locker = check_structure(ctx_fuer(
        pfad, ketten=("A", "B", "C"), erlaubte_sonderreste={"LIG"},
        einstellungen=GateSettings(kollision_toleranz_a=0.8)))
    assert locker.ok, locker.als_text()


def test_schwellen_haben_grenzen():
    """Niemand soll die Kollisionsprüfung durch eine absurde Toleranz aushebeln."""
    with pytest.raises(ValueError):
        GateSettings(kollision_toleranz_a=5.0)
    with pytest.raises(ValueError):
        GateSettings(kollisionen_erlaubt=10_000)


def test_abgestuerzter_pruefpunkt_gilt_als_fehler():
    @gate("test.kaputt", "Kaputter Prüfpunkt")
    def kaputt(ctx):
        raise RuntimeError("uups")

    bericht = run_gates([kaputt], None)
    assert not bericht.ok
    assert "abgestürzt" in bericht.fehler[0].meldung


def test_bericht_als_dict_und_text():
    bericht = GateReport([
        GateResult("a", "Erster", "ok", "alles gut"),
        GateResult("b", "Zweiter", "fehler", "kaputt"),
    ])
    assert not bericht.ok
    daten = bericht.als_dict()
    assert daten["bestanden"] is False
    assert [p["status"] for p in daten["pruefungen"]] == ["ok", "fehler"]
    assert "FEHLER" in bericht.als_text()
