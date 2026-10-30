"""Konfiguration laden und ablehnen, was nicht stimmt."""

from pathlib import Path

import pytest
import yaml

from oligomd.config import ConfigError, load_protocol, load_run, load_system

REPO = Path(__file__).resolve().parents[1]


def schreibe(tmp_path, name, daten) -> Path:
    pfad = tmp_path / name
    pfad.write_text(yaml.safe_dump(daten, allow_unicode=True), encoding="utf-8")
    return pfad


@pytest.fixture
def system_daten():
    return {
        "struktur": {"pdb": "eingaben/7pdu.pdb", "modell": 1, "ketten": ["A", "B"],
                     "entfernen": ["HOH"]},
        "kraftfeld": {"dna": "OL21", "rna": "OL3", "wasser": "OPC"},
        "box": {"typ": "oktaeder", "puffer_nm": 1.2},
        "ionen": {"neutralisieren_mit": "K+", "salz_mol_l": 0.15},
    }


@pytest.fixture
def protokoll_daten():
    return {
        "minimierung": {"schritte": 5000},
        "aequilibrierung": [
            {"name": "nvt", "ps": 100, "temperatur_k": 300, "restraints_kj_nm2": 1000},
            {"name": "npt", "ps": 500, "temperatur_k": 300, "druck_bar": 1.0,
             "restraints_kj_nm2": 0},
        ],
        "produktion": {"ns": 100, "zeitschritt_fs": 4, "hmr": True, "speichern_ps": 10},
        "zufallsstartwert": 42,
    }


def test_mitgelieferte_konfiguration_ist_gueltig():
    """Die Beispiele im Repository müssen immer laden."""
    lauf = load_run(REPO / "config" / "system.yaml", REPO / "config" / "protocol.yaml")
    assert lauf.system.struktur.ketten == ["A", "B"]
    assert lauf.system.modifikation is None
    assert lauf.gesamtdauer_ns == pytest.approx(101.6)


def test_defaults_greifen(tmp_path, system_daten):
    del system_daten["kraftfeld"], system_daten["box"], system_daten["ionen"]
    system = load_system(schreibe(tmp_path, "s.yaml", system_daten))
    assert system.kraftfeld.wasser == "OPC"
    assert system.box.puffer_nm == 1.2
    assert system.ionen.salz_mol_l == 0.15


def test_grosser_zeitschritt_ohne_hmr_wird_abgelehnt(tmp_path, protokoll_daten):
    protokoll_daten["produktion"]["hmr"] = False
    with pytest.raises(ConfigError, match="HMR"):
        load_protocol(schreibe(tmp_path, "p.yaml", protokoll_daten))


def test_aequilibrierung_muss_frei_enden(tmp_path, protokoll_daten):
    protokoll_daten["aequilibrierung"][-1]["restraints_kj_nm2"] = 500
    with pytest.raises(ConfigError, match="ohne Restraints"):
        load_protocol(schreibe(tmp_path, "p.yaml", protokoll_daten))


def test_modifikation_muss_in_einer_gewaehlten_kette_sitzen(tmp_path, system_daten):
    system_daten["modifikation"] = {
        "rest": {"kette": "Z", "nummer": 23},
        "definition": "config/residues/x.yaml", "name": "G6X",
    }
    with pytest.raises(ConfigError, match="Kette Z"):
        load_system(schreibe(tmp_path, "s.yaml", system_daten))

    system_daten["modifikation"]["rest"]["kette"] = "A"
    system = load_system(schreibe(tmp_path, "s2.yaml", system_daten))
    assert system.modifikation.name == "G6X"


def test_fehlermeldung_nennt_feld_und_grund(tmp_path, system_daten):
    system_daten["box"]["puffer_nm"] = 0.1
    with pytest.raises(ConfigError) as exc:
        load_system(schreibe(tmp_path, "s.yaml", system_daten))
    assert "box.puffer_nm" in str(exc.value)


def test_kaputtes_yaml_und_fehlende_datei(tmp_path):
    (tmp_path / "kaputt.yaml").write_text("das: [ist\nkein: yaml", encoding="utf-8")
    with pytest.raises(ConfigError, match="gültiges YAML"):
        load_system(tmp_path / "kaputt.yaml")
    with pytest.raises(ConfigError, match="ließ sich nicht lesen"):
        load_system(tmp_path / "gibtsnicht.yaml")
