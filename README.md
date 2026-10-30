# OligoMD-Agent

Pipeline für Molekulardynamik an Nukleinsäure-Komplexen mit modifizierten
Bausteinen. Gebaut für DNAzym-RNA-Komplexe: Struktur vorbereiten, Kraftfeld
zuweisen, simulieren, auswerten – reproduzierbar und mit harten Prüfpunkten
zwischen den Schritten. Darüber liegt später eine dünne Agenten-Schicht.

**Stand: Phase 1 (Pipeline für den unmodifizierten Komplex).** Konfiguration,
Strukturprüfung, Systembau mit tleap und die Simulation mit OpenMM stehen,
jeweils mit Prüfpunkten dazwischen. Parametrisierung der Modifikation und die
Agenten-Schicht kommen als Nächstes.

## Grundidee

> Der Agent schlägt vor, deterministischer Code entscheidet.

Alles Wissenschaftliche steht in Konfigurationsdateien und normalem Python-Code.
Zwischen den Schritten stehen **Prüfpunkte**, die bestehen müssen: fehlende
Atome, falsche Ladung, überlappende Atome, fehlende Kraftfeldparameter. Fällt
einer durch, hält die Pipeline an. Weder du noch der Agent können das abschalten.

Der Grund: Ein übergangener Prüfpunkt erzeugt keine Fehlermeldung, sondern eine
Simulation, die sauber durchläuft und falsch ist.

## Installation: zwei Umgebungen

Die Arbeit teilt sich auf zwei Hälften auf, weil kein System beides kann:

| Wo | Wofür | Umgebung |
|---|---|---|
| **WSL2 (Ubuntu)** | Struktur vorbereiten, Parameter bauen (AmberTools, modXNA) | `environment.yml` |
| **Natives Windows** | die Simulation auf der Grafikkarte (OpenMM) | `environment-windows.yml` |

AmberTools läuft nicht unter Windows, braucht aber keine Grafikkarte. OpenMM
läuft unter Windows und kommt dort an eine **AMD-Karte** heran, unter WSL2 in
der Regel nicht. Beide Hälften sehen denselben Lauf-Ordner, unter WSL2 als
`/mnt/c/...` und unter Windows als `C:\...`.

```bash
# in WSL2 (Ubuntu)
conda env create -f environment.yml
conda activate oligomd
pip install -e .
```

```powershell
# in Windows (Miniforge-Eingabeaufforderung)
conda env create -f environment-windows.yml
conda activate oligomd-win
pip install -e .
```

### AMD-Karte

OpenMM kennt vier Rechenplattformen. Für AMD ist **OpenCL** der verlässliche
Weg; den Treiber bringt die normale Radeon-Software mit. **HIP** ist schneller,
verlangt aber ein passendes HIP-SDK und ist laut
[OpenMM-Diskussionen](https://github.com/openmm/openmm/discussions/4933)
wählerisch. Fang mit OpenCL an (`plattform: OpenCL` in `protocol.yaml`) und
probiere HIP erst, wenn alles andere läuft.

`python -m oligomd doctor` zeigt dir, welche Plattformen OpenMM bei dir sieht.

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
# Was steckt in der Strukturdatei? (PDB oder PDBML, auch gezippt)
python -m oligomd info eingaben/7pdu.xml.gz

# Struktur und Konfiguration gegen die Prüfpunkte halten, ohne etwas anzulegen
python -m oligomd pruefe

# in WSL2: Lauf anlegen und das System bauen (tleap) -> prmtop und inpcrd
python -m oligomd vorbereiten --name 7pdu-referenz

# in Windows: rechnen (OpenMM auf der Grafikkarte)
python -m oligomd simulieren

# Welche Läufe gibt es und wie weit sind sie?
python -m oligomd laeufe
```

`simulieren` setzt an einem Prüfpunkt fort, wenn es einen gibt; mit `--neu`
fängt es von vorn an. Ein langer Lauf übersteht so einen Abbruch.

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

## Der Lauf-Ordner

Jeder Lauf liegt unter `runs/<datum>_<name>/` und beschreibt sich selbst:

```text
runs/2026-10-30_7pdu-referenz/
├─ konfiguration/      system.yaml und protocol.yaml, eingefroren
├─ vorbereitet/        eingabe.pdb, leap.in, tleap.log, system.prmtop, system.inpcrd
├─ simulation/         minimiert.pdb, produktion.dcd, produktion.csv, produktion.chk
├─ berichte/           vorbereiten.json, simulieren.json
└─ logbuch.txt
```

Die Konfiguration wird beim Anlegen hineinkopiert. Änderst du danach
`config/system.yaml`, bleibt der Lauf unberührt. Ein Lauf wird nie
überschrieben; gleicher Name ergibt einen zweiten Ordner.

## Prüfpunkte

Alle in `oligomd/gates.py`, `oligomd/prepare.py` und `oligomd/simulate.py`:

| ID | Prüft | Bei Verstoß |
|---|---|---|
| `struktur.ketten` | Alle konfigurierten Ketten sind da | Fehler |
| `struktur.fehlende_atome` | Keine Lücken in Standardnukleotiden | Fehler |
| `struktur.unbekannte_reste` | Jeder Rest hat ein Kraftfeld oder ist angekündigt | Fehler |
| `struktur.ladung` | Gesamtladung bestimmbar | Warnung, wenn Sonderreste offen sind |
| `struktur.kollisionen` | Kein Atompaar überlappt deutlich | Fehler ab 0,5 Å, darunter Warnung |
| `struktur.wasserstoffe` | Wasserstoffe vorhanden | Warnung |
| `system.tleap` | tleap ohne Abbruch | Fehler |
| `system.parameter` | Kein Rest ohne Kraftfeldparameter | Fehler |
| `system.dateien` | prmtop und inpcrd geschrieben | Fehler |
| `system.neutral` | Box trägt keine Restladung | Fehler |
| `system.wasser` | Wasserbox gefüllt | Fehler |
| `simulation.laeuft` | Messwerte vorhanden | Fehler |
| `simulation.energie` | Energie endlich, kein NaN | Fehler |
| `simulation.temperatur` | Mittlere Temperatur im Sollbereich | Fehler |
| `simulation.dichte` | Dichte zwischen 0,9 und 1,15 g/ml | Fehler |

Einstellbar sind nur Schwellen, und auch die nur in engen Grenzen
(`GateSettings`). Die Prüfpunkte selbst stehen im Code.

Die Kollisionsprüfung kennt zwei Besonderheiten echter Strukturen: Atome, die
beide eine **Wasserstoffbrücke** bilden können (N, O, F, S), dürfen sich bis
2,2 Å nähern, sonst gälte jede Watson-Crick-Paarung als Fehler. Und **leichte
Kontakte** bis 0,5 Å Überlappung sind nur eine Warnung, weil die Minimierung
sie wegräumt. Die sieben Modelle von 7PDU liegen alle unter 0,35 Å.

## Strukturdateien

Gelesen werden **PDB** und **PDBML** (`.xml`, `.xml.gz`), beides auch gezippt.
Von der [RCSB](https://www.rcsb.org/) lädst du am einfachsten die
PDBML-Variante herunter. mmCIF kann das Programm noch nicht.

Es gelten die `auth_*`-Felder, also die Nummerierung aus der Veröffentlichung.
Bei 7PDU heißt das: Kette A sind die DNA-Reste 1 bis 33, Kette B die RNA-Reste
34 bis 52.

## Phasen

| Phase | Inhalt | Stand |
|---|---|---|
| 0 | Gerüst, Konfiguration, Strukturprüfung, `doctor` | **fertig** |
| 1 | Pipeline für den unmodifizierten Komplex bis zur Produktion | **fertig** |
| 2 | Parameter für die Modifikation (ORCA + modXNA + AmberTools) | offen |
| 3 | Auswertung (MDAnalysis) | offen |
| 4 | Agenten-Schicht (Claude Agent SDK) | offen |
| 5 | Reihenuntersuchung mehrerer Modifikationen | offen |

Nach Phase 2 gibt es verwertbare Ergebnisse, auch ohne Agent.

## Ionen und Ladung

Die Ladung einer Box muss null sein, sonst wird die Elektrostatik unter
periodischen Randbedingungen unsinnig. Dabei zählen **Zusatzionen mit**: 16 Mg²⁺
bringen +32 mit, also braucht ein Solut mit −50 nur noch 18 Gegenionen statt 50.
OligoMD rechnet das aus, schreibt die Zahlen ins tleap-Skript und prüft
hinterher nach, was tleap gemeldet hat.

Zusatzionen lassen sich als feste Zahl (`anzahl`) oder als Konzentration
(`mol_l`) angeben. Die Konzentration ist der ehrlichere Weg, weil sie nicht von
der Boxgröße abhängt.

## Entwicklung

```bash
pytest          # läuft ohne OpenMM, ohne AmberTools, ohne GPU
ruff check .
```

Die Tests ersetzen externe Programme durch kleine Skripte und bauen ihre
PDB-Dateien selbst (`tests/conftest.py`). Wer einen Prüfpunkt hinzufügt,
schreibt einen Test dazu, der ihn auslöst.

Projektregeln: `CLAUDE.md`.
