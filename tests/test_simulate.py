"""Planung und Prüfpunkte der Simulation. Läuft ohne OpenMM und ohne Grafikkarte."""

import pytest
import yaml

from oligomd.config import load_protocol
from oligomd.simulate import (
    SimContext, SimulateError, lies_messwerte, plan_erstellen, waehle_plattform, SIM_GATES,
)
from oligomd.gates import run_gates


def protokoll(tmp_path, **aender):
    daten = {
        "minimierung": {"schritte": 5000},
        "aequilibrierung": [
            {"name": "nvt", "ps": 100, "temperatur_k": 300, "restraints_kj_nm2": 1000},
            {"name": "npt", "ps": 500, "temperatur_k": 300, "druck_bar": 1.0,
             "restraints_kj_nm2": 0},
        ],
        "produktion": {"ns": 100, "zeitschritt_fs": 4, "hmr": True, "speichern_ps": 10},
        "zufallsstartwert": 42,
    }
    daten.update(aender)
    pfad = tmp_path / "p.yaml"
    pfad.write_text(yaml.safe_dump(daten), encoding="utf-8")
    return load_protocol(pfad)


def test_zeiten_werden_zu_schritten(tmp_path):
    plan = plan_erstellen(protokoll(tmp_path))
    # Produktion: 100 ns bei 4 fs
    assert plan.produktion.schritte == 25_000_000
    assert plan.produktion.ps == 100_000.0
    # alle 10 ps ein Bild
    assert plan.speichern_alle == 2500
    assert plan.gesamt_ns == pytest.approx(100.6)


def test_aequilibrierung_laeuft_mit_kleinerem_zeitschritt(tmp_path):
    """Solange Restraints wirken und sich die Box setzt, sind 2 fs die sichere Wahl."""
    plan = plan_erstellen(protokoll(tmp_path))
    assert all(s.zeitschritt_fs == 2.0 for s in plan.aequilibrierung)
    assert plan.produktion.zeitschritt_fs == 4.0
    assert plan.aequilibrierung[0].schritte == 50_000      # 100 ps bei 2 fs
    assert plan.aequilibrierung[0].restraints_kj_nm2 == 1000


def test_kleiner_zeitschritt_wird_nicht_vergroessert(tmp_path):
    p = protokoll(tmp_path, produktion={"ns": 1, "zeitschritt_fs": 1, "hmr": False,
                                        "speichern_ps": 1})
    plan = plan_erstellen(p)
    assert plan.produktion.zeitschritt_fs == 1.0
    assert all(s.zeitschritt_fs == 1.0 for s in plan.aequilibrierung)
    assert plan.hmr is False


@pytest.mark.parametrize("gewuenscht,verfuegbar,erwartet", [
    ("auto", ["Reference", "CPU", "OpenCL"], "OpenCL"),
    ("auto", ["Reference", "CPU", "OpenCL", "CUDA"], "CUDA"),
    ("auto", ["Reference", "CPU", "HIP", "OpenCL"], "HIP"),
    ("auto", ["Reference", "CPU"], "CPU"),
    ("OpenCL", ["CPU", "OpenCL"], "OpenCL"),
])
def test_plattformwahl(gewuenscht, verfuegbar, erwartet):
    assert waehle_plattform(gewuenscht, verfuegbar) == erwartet


def test_fehlende_plattform_erklaert_sich():
    with pytest.raises(SimulateError, match="AMD"):
        waehle_plattform("CUDA", ["CPU", "OpenCL"])
    with pytest.raises(SimulateError, match="keine einzige"):
        waehle_plattform("auto", [])


# ---------- Prüfpunkte auf die Messwerte ----------

KOPF = ('#"Progress (%)","Step","Time (ps)","Potential Energy (kJ/mole)",'
        '"Total Energy (kJ/mole)","Temperature (K)","Density (g/mL)"')


def messwerte_datei(tmp_path, zeilen):
    pfad = tmp_path / "produktion.csv"
    pfad.write_text(KOPF + "\n" + "\n".join(zeilen) + "\n", encoding="utf-8")
    return pfad


def test_messwerte_lesen(tmp_path):
    pfad = messwerte_datei(tmp_path, ['1,100,0.4,-100.0,-90.0,299.5,1.001',
                                      '2,200,0.8,-101.0,-91.0,300.5,0.999'])
    werte = lies_messwerte(pfad)
    assert len(werte) == 2
    assert werte[0]["Temperature (K)"] == "299.5"
    assert lies_messwerte(tmp_path / "gibtsnicht.csv") == []


def gute_werte(n=10, temperatur=300.0, dichte=1.0, energie=-1000.0):
    return [{"Step": str(i * 100), "Total Energy (kJ/mole)": str(energie),
             "Temperature (K)": str(temperatur), "Density (g/mL)": str(dichte)}
            for i in range(n)]


def test_saubere_simulation_besteht():
    bericht = run_gates(SIM_GATES, SimContext(gute_werte(), 300.0))
    assert bericht.ok, bericht.als_text()


def test_abgebrochene_simulation_faellt_durch():
    bericht = run_gates(SIM_GATES, SimContext([], 300.0))
    assert not bericht.ok
    assert any(f.id == "simulation.laeuft" for f in bericht.fehler)


def test_explodiertes_system_faellt_durch():
    werte = gute_werte()
    werte[5]["Total Energy (kJ/mole)"] = "nan"
    bericht = run_gates(SIM_GATES, SimContext(werte, 300.0))
    assert not bericht.ok
    fehler = next(f for f in bericht.fehler if f.id == "simulation.energie")
    assert "explodiert" in fehler.meldung


def test_falsche_temperatur_faellt_durch():
    bericht = run_gates(SIM_GATES, SimContext(gute_werte(temperatur=350.0), 300.0))
    assert not bericht.ok
    assert any(f.id == "simulation.temperatur" for f in bericht.fehler)


def test_aufheizphase_zaehlt_nicht():
    """Am Anfang ist das System kalt. Gewertet wird die zweite Hälfte."""
    werte = gute_werte(n=10, temperatur=300.0)
    for zeile in werte[:5]:
        zeile["Temperature (K)"] = "120.0"
    bericht = run_gates(SIM_GATES, SimContext(werte, 300.0))
    assert bericht.ok, bericht.als_text()


def test_falsche_dichte_faellt_durch():
    bericht = run_gates(SIM_GATES, SimContext(gute_werte(dichte=0.4), 300.0))
    assert not bericht.ok
    fehler = next(f for f in bericht.fehler if f.id == "simulation.dichte")
    assert "0.400" in fehler.meldung


def test_fehlende_dichte_warnt_nur():
    werte = [{"Step": "1", "Total Energy (kJ/mole)": "-1000", "Temperature (K)": "300"}
             for _ in range(4)]
    bericht = run_gates(SIM_GATES, SimContext(werte, 300.0))
    assert bericht.ok
    assert [w.id for w in bericht.warnungen] == ["simulation.dichte"]
