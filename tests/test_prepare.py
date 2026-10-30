"""System bauen: tleap-Skript, Ausgabe lesen, Prüfpunkte. Ohne echtes tleap."""

import os
import stat
import sys
from pathlib import Path

import pytest
import yaml

from oligomd.config import SystemConfig
from oligomd.prepare import (
    LeapPlan, PrepareError, amber_ion, extra_ionen, gegenionen, ion_ladung, leap_script,
    lies_leap_ausgabe, salz_ionenpaare, vorbereiten, zaehle_wasser,
)
from oligomd.run import Run

from .conftest import build_pdb

FAKE = Path(__file__).parent / "fake_tleap.py"


# ---------- reine Funktionen ----------

def system_config(**kwargs) -> SystemConfig:
    daten = {
        "struktur": {"pdb": "eingaben/x.pdb", "ketten": ["A"]},
        "kraftfeld": {"dna": "OL21", "rna": "OL3", "wasser": "OPC"},
        "box": {"typ": "oktaeder", "puffer_nm": 1.2},
        "ionen": {"neutralisieren_mit": "K+", "salz_mol_l": 0.15},
    }
    daten.update(kwargs)
    return SystemConfig.model_validate(daten)


def test_leap_skript_erster_durchgang_speichert_nichts():
    skript = leap_script(LeapPlan(system=system_config(), nur_solvatisieren=True))
    assert "source leaprc.DNA.OL21" in skript
    assert "source leaprc.RNA.OL3" in skript
    assert "source leaprc.water.opc" in skript
    assert "solvateoct mol OPCBOX 12.0" in skript
    assert "saveamberparm" not in skript   # erster Durchgang zählt nur
    assert "addions" not in skript


def test_leap_skript_zweiter_durchgang_setzt_ionen():
    plan = LeapPlan(system=system_config(), neutralisieren=50, salzpaare=24,
                    extra_ionen={"MG": 10})
    skript = leap_script(plan)
    assert "addionsrand mol MG 10" in skript
    assert "addions mol K+ 50" in skript
    assert "addionsrand mol K+ 24 Cl- 24" in skript
    assert "saveamberparm mol system.prmtop system.inpcrd" in skript
    # Erst die Zusatzionen, dann ausgleichen - sonst bleibt die Box geladen
    assert skript.index("addionsrand mol MG") < skript.index("addions mol K+ 50")


def test_zusatzionen_verschieben_die_ladungsbilanz():
    """Zehn Mg2+ bringen +20 mit. Das muss der Ausgleich berücksichtigen."""
    assert gegenionen(-50, {}) == 50
    assert gegenionen(-50, {"MG": 10}) == 30        # -50 + 20
    assert gegenionen(-50, {"MG": 30}) == 10        # -50 + 60, jetzt positiv
    assert gegenionen(-2, {"MG": 2}) == 2           # -2 + 4 = +2
    assert ion_ladung("MG") == 2 and ion_ladung("Cl-") == -1 and ion_ladung("XYZ") == 0


def test_zusatzionen_aus_konzentration():
    system = system_config(ionen={"neutralisieren_mit": "K+", "salz_mol_l": 0.15,
                                  "zusaetzlich": [{"ion": "Mg2+", "mol_l": 0.05}]})
    # 0,05 mol/l bei 18000 Wasser: 0,05 * 18000 / 55,56
    assert extra_ionen(system, 18000) == {"MG": 16}
    assert extra_ionen(system, 100) == {"MG": 1}    # nie auf null abrunden


def test_leap_skript_wuerfel_und_zusatzdateien():
    plan = LeapPlan(system=system_config(box={"typ": "wuerfel", "puffer_nm": 1.0}),
                    zusatz_lib=["g6x.lib"], zusatz_frcmod=["g6x.frcmod"])
    skript = leap_script(plan)
    assert "solvatebox mol OPCBOX 10.0" in skript
    assert "loadoff g6x.lib" in skript
    assert "loadamberparams g6x.frcmod" in skript
    # Zusatzdateien müssen vor dem Laden der Struktur kommen
    assert skript.index("loadoff") < skript.index("loadpdb")


def test_unbekanntes_wassermodell_wird_abgelehnt():
    system = system_config()
    system.kraftfeld.wasser = "Quatsch"
    with pytest.raises(PrepareError, match="Wassermodell"):
        leap_script(LeapPlan(system=system))


@pytest.mark.parametrize("eingabe,erwartet", [
    ("Mg2+", "MG"), ("MG", "MG"), ("K+", "K+"), ("Na", "Na+"), ("Cl-", "Cl-"), ("XYZ", "XYZ"),
])
def test_ionennamen(eingabe, erwartet):
    assert amber_ion(eingabe) == erwartet


def test_salzmenge_haengt_an_der_wassermenge():
    assert salz_ionenpaare(9000, 0.15) == 24
    assert salz_ionenpaare(9000, 0.0) == 0
    assert salz_ionenpaare(0, 0.15) == 0


def test_leap_ausgabe_lesen():
    ok = lies_leap_ausgabe("Loading ...\nTotal unperturbed charge:  0.000000\n\tQuit", 0)
    assert ok.erfolg and ok.restladung == 0.0

    fatal = lies_leap_ausgabe("FATAL:  Atom does not have a type.", 1)
    assert not fatal.erfolg and fatal.fatale

    fehlend = lies_leap_ausgabe("Unknown residue: G6X   number: 22", 0)
    assert not fehlend.erfolg and fehlend.fehlende_parameter

    # Die zweite Schreibweise von tleap muss auch erkannt werden
    geladen = lies_leap_ausgabe(
        "WARNING: The unperturbed charge of the unit: -50.000000 is not zero.", 0)
    assert geladen.restladung == -50.0


def test_wasser_zaehlen(tmp_path):
    pdb = tmp_path / "s.pdb"
    pdb.write_text(
        "ATOM      1  O   WAT W   1       0.000   0.000   0.000  1.00  0.00           O\n"
        "ATOM      2  H1  WAT W   1       1.000   0.000   0.000  1.00  0.00           H\n"
        "ATOM      3  O   WAT W   2       5.000   0.000   0.000  1.00  0.00           O\n"
        "ATOM      4  P   DG  A   1      10.000   0.000   0.000  1.00  0.00           P\n",
        encoding="utf-8")
    assert zaehle_wasser(pdb) == 2
    assert zaehle_wasser(tmp_path / "gibtsnicht.pdb") == 0


def test_wasser_zaehlen_ueber_9999_reste(tmp_path):
    """Restnummern laufen im PDB-Format bei 9999 um; gezählt werden Sauerstoffe."""
    pdb = tmp_path / "gross.pdb"
    zeilen = []
    for i in range(12000):
        nr = (i % 9999) + 1
        for name, element in (("O", "O"), ("H1", "H"), ("H2", "H")):
            zeilen.append(f"ATOM  {i % 100000:5d}  {name:<3s} WAT W{nr:4d}    "
                          f"{0.0:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00          {element:>2s}")
    pdb.write_text("\n".join(zeilen) + "\n", encoding="utf-8")
    assert zaehle_wasser(pdb) == 12000


# ---------- ganzer Schritt mit Ersatz-tleap ----------

@pytest.fixture
def projekt(tmp_path, monkeypatch):
    """Ein Miniprojekt mit Struktur, Konfiguration und Ersatz-tleap."""
    if sys.platform == "win32":
        pytest.skip("Ersatzskript ohne .exe nur unter Linux/macOS")
    monkeypatch.chdir(tmp_path)
    (tmp_path / "eingaben").mkdir()
    (tmp_path / "config").mkdir()
    (tmp_path / "eingaben" / "test.pdb").write_text(
        build_pdb([("DG", "A", 1), ("DA", "A", 2), ("DC", "A", 3)]), encoding="utf-8")
    (tmp_path / "config" / "system.yaml").write_text(yaml.safe_dump({
        "struktur": {"pdb": "eingaben/test.pdb", "ketten": ["A"]},
        "ionen": {"neutralisieren_mit": "K+", "salz_mol_l": 0.15,
                  "zusaetzlich": [{"ion": "MG", "anzahl": 2}]},
    }), encoding="utf-8")
    (tmp_path / "config" / "protocol.yaml").write_text(yaml.safe_dump({
        "aequilibrierung": [{"name": "npt", "ps": 100, "restraints_kj_nm2": 0}],
        "produktion": {"ns": 1},
    }), encoding="utf-8")

    tleap = tmp_path / "tleap"
    tleap.write_text(f'#!/bin/sh\nexec "{sys.executable}" "{FAKE}" "$@"\n', encoding="utf-8")
    tleap.chmod(tleap.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setenv("OLIGOMD_TLEAP", str(tleap))
    monkeypatch.setenv("FAKE_TLEAP_MODUS", "ok")
    return tmp_path


def neuer_lauf(name="test") -> Run:
    return Run.neu(name, Path("config/system.yaml"), Path("config/protocol.yaml"))


def test_vorbereiten_laeuft_durch(projekt):
    lauf = neuer_lauf()
    bericht = vorbereiten(lauf)
    assert bericht.ok, bericht.als_text()

    arbeit = lauf.ordner / "vorbereitet"
    assert (arbeit / "system.prmtop").is_file()
    assert (arbeit / "system.inpcrd").is_file()
    assert (arbeit / "eingabe.pdb").is_file()
    # Die Eingabe für tleap enthält keine Wasserstoffe; die setzt tleap selbst
    assert " H" not in (arbeit / "eingabe.pdb").read_text(encoding="utf-8")[-2000:]

    bericht_json = lauf.lies_bericht("vorbereiten")
    assert bericht_json["bestanden"] is True
    assert bericht_json["wasser"] == 9000
    assert bericht_json["ladung_solut"] == -2
    # Solut -2, dazu 2 Mg2+ (+4) macht +2: der Ausgleich sind Chloride, nicht Kalium.
    # 24 Salzpaare kommen oben drauf.
    assert bericht_json["ionen"] == {"MG": 2, "Cl-": 26, "K+": 24}
    assert lauf.schritt_fertig("vorbereiten")


def test_konfiguration_wird_eingefroren(projekt):
    lauf = neuer_lauf()
    eingefroren = lauf.ordner / "konfiguration" / "system.yaml"
    assert eingefroren.is_file()
    # Konfiguration im Projekt ändern: der Lauf bleibt, wie er war
    (projekt / "config" / "system.yaml").write_text(yaml.safe_dump({
        "struktur": {"pdb": "eingaben/anders.pdb", "ketten": ["Z"]}}), encoding="utf-8")
    wieder = Run.oeffnen(lauf.ordner)
    assert wieder.konfiguration.system.struktur.ketten == ["A"]


def test_zwei_laeufe_ueberschreiben_sich_nicht(projekt):
    a, b = neuer_lauf("gleich"), neuer_lauf("gleich")
    assert a.ordner != b.ordner
    assert Run.letzter().ordner == b.ordner


def test_fehlende_parameter_halten_an(projekt, monkeypatch):
    monkeypatch.setenv("FAKE_TLEAP_MODUS", "fehlende_parameter")
    bericht = vorbereiten(neuer_lauf())
    assert not bericht.ok
    fehler = next(f for f in bericht.fehler if f.id == "system.parameter")
    assert "G6X" in str(fehler.details)


def test_abbruch_von_tleap_haelt_an(projekt, monkeypatch):
    monkeypatch.setenv("FAKE_TLEAP_MODUS", "fatal")
    bericht = vorbereiten(neuer_lauf())
    assert not bericht.ok
    assert any(f.id in ("system.tleap", "system.wasser") for f in bericht.fehler)


def test_geladenes_system_haelt_an(projekt, monkeypatch):
    monkeypatch.setenv("FAKE_TLEAP_MODUS", "geladen")
    bericht = vorbereiten(neuer_lauf())
    assert not bericht.ok
    fehler = next(f for f in bericht.fehler if f.id == "system.neutral")
    assert "-50" in fehler.meldung


def test_zu_wenig_wasser_haelt_an(projekt, monkeypatch):
    monkeypatch.setenv("FAKE_TLEAP_WASSER", "10")
    bericht = vorbereiten(neuer_lauf())
    assert not bericht.ok
    assert any(f.id == "system.wasser" for f in bericht.fehler)


def test_kaputte_struktur_stoppt_vor_tleap(projekt):
    (projekt / "eingaben" / "test.pdb").write_text(
        build_pdb([("DG", "A", 1), ("DG", "A", 2)], weglassen={("A", 2): ("O6",)}),
        encoding="utf-8")
    lauf = neuer_lauf()
    bericht = vorbereiten(lauf)
    assert not bericht.ok
    assert any(f.id == "struktur.fehlende_atome" for f in bericht.fehler)
    # tleap wurde gar nicht erst gestartet
    assert not (lauf.ordner / "vorbereitet" / "leap.in").exists()


def test_fehlendes_tleap_meldet_sich_verstaendlich(projekt, monkeypatch):
    monkeypatch.setenv("OLIGOMD_TLEAP", str(projekt / "gibtsnicht"))
    with pytest.raises(PrepareError, match="tleap"):
        vorbereiten(neuer_lauf())
