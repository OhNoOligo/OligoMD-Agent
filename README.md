# OligoMD-Agent

Pipeline für Molekulardynamik an Nukleinsäure-Komplexen mit modifizierten
Bausteinen. Gebaut für DNAzym-RNA-Komplexe: Struktur vorbereiten, Kraftfeld
zuweisen, simulieren, auswerten – reproduzierbar und mit harten Prüfpunkten
zwischen den Schritten. Darüber liegt später eine dünne Agenten-Schicht.

**Stand: Phase 0 (Gerüst und Prüfpunkte).** Konfiguration, Strukturprüfung und
Umgebungstest laufen. Simulation, Parametrisierung und Agent kommen in den
nächsten Phasen.

## Grundidee

> Der Agent schlägt vor, deterministischer Code entscheidet.

Alles Wissenschaftliche steht in Konfigurationsdateien und normalem Python-Code.
Zwischen den Schritten stehen **Prüfpunkte**, die bestehen müssen: fehlende
Atome, falsche Ladung, überlappende Atome, fehlende Kraftfeldparameter. Fällt
einer durch, hält die Pipeline an. Weder du noch der Agent können das abschalten.

Der Grund: Ein übergangener Prüfpunkt erzeugt keine Fehlermeldung, sondern eine
Simulation, die sauber durchläuft und falsch ist.

## Installation

Die MD-Werkzeuge laufen nicht nativ unter Windows. Deshalb **WSL2**:

```bash
# in WSL2 (Ubuntu)
conda env create -f environment.yml
conda activate oligomd
pip install -e .
```

**modXNA** (für modifizierte Nukleotide, ab Phase 2) wird getrennt installiert.
Es ist ein BASH-Werkzeug und liegt nicht auf conda-forge:
<https://modxna.chpc.utah.edu/> · [Code](https://github.com/modxna/modxna) ·
[Tutorial-Papier](https://livecomsjournal.org/index.php/livecoms/article/view/v6i1e4545)

Prüfen, ob alles da ist:

```bash
python -m oligomd doctor --phase 1
```

Die Ausgabe sagt für jede Abhängigkeit, ob sie fehlt und wofür sie gebraucht
wird. `--phase` gibt an, wie weit du kommen willst: Phase 0 braucht nur Python,
Phase 1 zusätzlich OpenMM und tleap, Phase 2 die AmberTools und modXNA.

ORCA bleibt davon unberührt: Die quantenchemischen Rechnungen laufen weiter in
ORCA Studio unter Windows, die Ergebnisse werden als Dateien übergeben.

## Benutzung

```bash
# Was steckt in einer PDB-Datei?
python -m oligomd info eingaben/7pdu.pdb

# Struktur und Konfiguration gegen die Prüfpunkte halten
python -m oligomd pruefe --system config/system.yaml --protokoll config/protocol.yaml
```

`pruefe` gibt 0 zurück, wenn alles besteht, 2 bei einem blockierenden Prüfpunkt
und 1, wenn die Konfiguration selbst kaputt ist. Mit `--json` kommt alles
maschinenlesbar heraus; so liest es später der Agent.

## Konfiguration

Zwei Dateien beschreiben einen Lauf vollständig:

| Datei | Inhalt |
|---|---|
| `config/system.yaml` | Was: Struktur, Ketten, Kraftfeld, Box, Ionen, Modifikation |
| `config/protocol.yaml` | Wie: Minimierung, Äquilibrierung, Produktion, Zufallsstartwert |

Beide werden beim Start in den Lauf-Ordner kopiert. Ab dann gilt die Kopie:
Änderst du später die Datei im Repository, bleibt der alte Lauf unberührt.

Die Konfiguration wehrt sich gegen typische Fehler, zum Beispiel gegen einen
Zeitschritt von 4 fs ohne Wasserstoff-Massen-Umverteilung oder gegen eine
Äquilibrierung, die mit festgehaltenem Solut endet.

## Prüfpunkte

Stand Phase 0, alle in `oligomd/gates.py`:

| ID | Prüft | Bei Verstoß |
|---|---|---|
| `struktur.ketten` | Alle konfigurierten Ketten sind da | Fehler |
| `struktur.fehlende_atome` | Keine Lücken in Standardnukleotiden | Fehler |
| `struktur.unbekannte_reste` | Jeder Rest hat ein Kraftfeld oder ist angekündigt | Fehler |
| `struktur.ladung` | Gesamtladung bestimmbar | Warnung, wenn Sonderreste offen sind |
| `struktur.kollisionen` | Kein Atompaar unter der VdW-Grenze | Fehler |
| `struktur.wasserstoffe` | Wasserstoffe vorhanden | Warnung |

Einstellbar sind nur Schwellen, und auch die nur in engen Grenzen
(`GateSettings`). Die Prüfpunkte selbst stehen im Code.

## Phasen

| Phase | Inhalt | Stand |
|---|---|---|
| 0 | Gerüst, Konfiguration, Strukturprüfung, `doctor` | **fertig** |
| 1 | Pipeline für den unmodifizierten Komplex bis zur Produktion | offen |
| 2 | Parameter für die Modifikation (ORCA + modXNA + AmberTools) | offen |
| 3 | Auswertung (MDAnalysis) | offen |
| 4 | Agenten-Schicht (Claude Agent SDK) | offen |
| 5 | Reihenuntersuchung mehrerer Modifikationen | offen |

Nach Phase 2 gibt es verwertbare Ergebnisse, auch ohne Agent.

## Entwicklung

```bash
pytest          # läuft ohne OpenMM, ohne AmberTools, ohne GPU
ruff check .
```

Die Tests ersetzen externe Programme durch kleine Skripte und bauen ihre
PDB-Dateien selbst (`tests/conftest.py`). Wer einen Prüfpunkt hinzufügt,
schreibt einen Test dazu, der ihn auslöst.

Projektregeln: `CLAUDE.md`.
