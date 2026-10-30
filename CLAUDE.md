# OligoMD-Agent – Projektregeln

Pipeline für Molekulardynamik an Nukleinsäure-Komplexen mit modifizierten Bausteinen,
mit einer dünnen Agenten-Schicht darüber. Zielsystem: DNAzym-RNA-Komplexe.
Diese Regeln gelten in allen Phasen.

## Die oberste Regel

**Der Agent schlägt vor, deterministischer Code entscheidet.**

Alles Wissenschaftliche steht in Konfigurationsdateien und normalem Python-Code.
Der Agent wählt Werkzeuge aus, liest Fehler und schreibt Berichte. Er erfindet
keine Parameter, überspringt keine Prüfpunkte und ändert nichts in `runs/`.

## Feste Regeln

- **Prüfpunkte (`oligomd/gates.py`) stehen im Code, nicht in der Konfiguration.**
  Die Konfiguration darf nur Schwellen innerhalb dokumentierter Grenzen setzen.
  Ein Prüfpunkt mit Status `fehler` hält die Pipeline an. Niemand schaltet ihn ab.
- **Der Lauf-Ordner ist die Wahrheit, nicht die Konfiguration im Repository.**
  Jeder Lauf liegt unter `runs/<JJJJ-MM-TT>_<name>/` mit eingefrorener Konfiguration,
  Logs, Prüfpunktbericht und Ergebnissen. Ein Lauf wird nie nachträglich überschrieben.
- **Reproduzierbarkeit vor Bequemlichkeit.** Jeder Zufallsstartwert steht in der
  Konfiguration. Zwei Läufe mit gleicher Konfiguration müssen gleich ausgehen.
- **Externe Programme immer mit vollem Pfad und Argumentliste aufrufen**,
  `subprocess.run` ohne `shell=True`. Pfade zu tleap, modXNA usw. kommen aus
  `oligomd/tools.py`, nie fest verdrahtet.
- **Ausgabedateien immer mit `encoding="utf-8", errors="replace"` lesen.**
- Keine stillen Ersetzungen: Fehlt ein Parameter, ein Atom oder ein Ion, bricht
  die Pipeline ab und sagt, was fehlt. Nie durch etwas Ähnliches ersetzen.

## Architektur

```
config/       system.yaml, protocol.yaml, residues/*.yaml  (Eingaben)
oligomd/      der deterministische Kern
  config.py     Konfiguration laden und prüfen (pydantic)
  structure.py  PDB lesen, Reste, Ladung, Kollisionen (rdkit + numpy)
  gates.py      Prüfpunkte
  tools.py      externe Programme finden und aufrufen
  doctor.py     Umgebung prüfen
  prepare.py    Struktur → solvatisiertes System      (Phase 1)
  params.py     modXNA / antechamber / tleap          (Phase 2)
  qm.py         Anbindung an ORCA-Ergebnisse          (Phase 2)
  simulate.py   OpenMM                                (Phase 1)
  analyze.py    MDAnalysis                            (Phase 3)
  report.py     Bericht je Lauf
agent/        dünne @tool-Hüllen um die Module oben    (Phase 4)
runs/         je Lauf ein Ordner, nicht im Repository
```

- **Schwere Abhängigkeiten nur dort importieren, wo sie gebraucht werden.**
  `openmm`, `MDAnalysis` und `pdbfixer` niemals auf Modulebene in `config.py`,
  `structure.py` oder `gates.py`: Diese Module und ihre Tests müssen ohne sie laufen.
- **Tests laufen ohne OpenMM, ohne AmberTools und ohne GPU** (`pytest`), genauso wie
  die Tests von ORCA Studio ohne ORCA laufen. Externe Programme werden in Tests
  durch kleine Skripte ersetzt.

## Umgebung

- Alles rund um MD läuft in **WSL2 (Ubuntu)**, conda-Umgebung aus `environment.yml`.
- **ORCA Studio bleibt auf nativem Windows.** Quantenchemische Ergebnisse werden
  als Dateien übergeben, es gibt keine direkte Kopplung der beiden Projekte.
- `python -m oligomd doctor` prüft die Umgebung und ist der Abnahmetest für Phase 0.

## Benutzertexte

Deutsch, Du-Form. Fehlermeldungen sagen, was zu tun ist, nicht nur was kaputt ist.
Jede Prüfpunkt-Meldung nennt den konkreten Rest, das Atom oder den Zahlenwert.
