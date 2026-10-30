"""Kommandozeile: python -m oligomd <befehl>

Befehle in Phase 0:
  doctor       Umgebung prüfen
  info         Was steckt in einer Strukturdatei
  pruefe       Struktur und Konfiguration gegen die Prüfpunkte halten
  vorbereiten  Lauf anlegen, System bauen (tleap): prmtop und inpcrd
  simulieren   Minimierung, Äquilibrierung, Produktion (OpenMM)
  laeufe       Welche Läufe es gibt und wie weit sie sind
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .config import ConfigError, load_run
from .doctor import als_text, report
from .gates import GateSettings, StructureContext, check_structure
from .structure import StructureError, inventory, models_in, net_charge, read_structure


def cmd_vorbereiten(args) -> int:
    from .prepare import PrepareError, vorbereiten
    from .run import Run, RunError
    try:
        lauf = Run.neu(args.name, Path(args.system), Path(args.protokoll))
    except (ConfigError, RunError) as exc:
        print(f"FEHLER: {exc}", file=sys.stderr)
        return 1
    print(f"Lauf: {lauf.ordner}")
    try:
        bericht = vorbereiten(lauf)
    except (PrepareError, RunError, StructureError) as exc:
        print(f"FEHLER: {exc}", file=sys.stderr)
        return 1
    if args.json:
        print(json.dumps(bericht.als_dict(), indent=2, ensure_ascii=False))
    else:
        print(bericht.als_text())
        print()
        print("Bestanden. Weiter mit: python -m oligomd simulieren" if bericht.ok
              else "NICHT bestanden. Die Pipeline hält hier an.")
    return 0 if bericht.ok else 2


def cmd_simulieren(args) -> int:
    from .run import Run, RunError
    from .simulate import SimulateError, simulieren
    try:
        lauf = Run.oeffnen(Path(args.lauf)) if args.lauf else Run.letzter()
    except (ConfigError, RunError) as exc:
        print(f"FEHLER: {exc}", file=sys.stderr)
        return 1
    print(f"Lauf: {lauf.ordner}")
    try:
        bericht = simulieren(lauf, fortsetzen=not args.neu)
    except (SimulateError, RunError) as exc:
        print(f"FEHLER: {exc}", file=sys.stderr)
        return 1
    if args.json:
        print(json.dumps(bericht.als_dict(), indent=2, ensure_ascii=False))
    else:
        print(bericht.als_text())
        print("\nBestanden." if bericht.ok else "\nNICHT bestanden.")
    return 0 if bericht.ok else 2


def cmd_laeufe(args) -> int:
    from .run import RUNS, Run
    basis = Path(args.ordner)
    if not basis.is_dir():
        print(f"Unter {basis} liegt noch kein Lauf.")
        return 0
    zeilen = []
    for pfad in sorted(p for p in basis.iterdir() if p.is_dir()):
        try:
            lauf = Run.oeffnen(pfad)
        except Exception:
            continue
        schritte = [f"{name}: {'fertig' if lauf.schritt_fertig(name) else 'offen'}"
                    for name in ("vorbereiten", "simulieren")]
        zeilen.append(f"{pfad.name:40s} {'  '.join(schritte)}")
    print("\n".join(zeilen) if zeilen else f"Unter {basis} liegt noch kein Lauf.")
    return 0


def cmd_doctor(args) -> int:
    befunde, _ = report(args.phase)
    print(als_text(befunde, args.phase))
    fehlend = [b for b in befunde if b.status != "ok" and b.phase <= args.phase]
    return 1 if fehlend else 0


def cmd_info(args) -> int:
    pfad = Path(args.pdb)
    try:
        struktur = read_structure(pfad, args.modell)
    except StructureError as exc:
        print(f"FEHLER: {exc}", file=sys.stderr)
        return 1
    daten = inventory(struktur)
    daten["modelle_in_datei"] = models_in(pfad)
    daten["gelesenes_modell"] = args.modell
    daten["gesamtladung_standardreste"] = net_charge(struktur)
    if args.json:
        print(json.dumps(daten, indent=2, ensure_ascii=False))
        return 0
    print(f"{pfad.name}, Modell {args.modell} von {daten['modelle_in_datei']}")
    print(f"  Atome            {daten['atome']} (davon {daten['wasserstoffe']} Wasserstoffe)")
    print(f"  Ketten           " + ", ".join(f"{k}: {n} Reste" for k, n in daten["ketten"].items()))
    print(f"  Reste nach Art   " + ", ".join(f"{k}: {n}" for k, n in daten["reste"].items()))
    if daten["ionen"]:
        print(f"  Ionen            " + ", ".join(f"{k}: {n}" for k, n in daten["ionen"].items()))
    if daten["unbekannte_reste"]:
        print(f"  Sonderreste      " + ", ".join(daten["unbekannte_reste"]))
    print(f"  Ladung (Standard) {daten['gesamtladung_standardreste']:+d}")
    return 0


def cmd_pruefe(args) -> int:
    try:
        lauf = load_run(Path(args.system), Path(args.protokoll))
    except ConfigError as exc:
        print(f"FEHLER: {exc}", file=sys.stderr)
        return 1

    pdb = lauf.system.struktur.pdb
    if not pdb.is_absolute():
        pdb = Path(args.system).resolve().parent.parent / pdb
    try:
        struktur = read_structure(pdb, lauf.system.struktur.modell)
    except StructureError as exc:
        print(f"FEHLER: {exc}", file=sys.stderr)
        return 1

    struktur = struktur.subset(
        ketten=lauf.system.struktur.ketten,
        ohne_restnamen=set(lauf.system.struktur.entfernen),
    )
    ctx = StructureContext(
        structure=struktur,
        erwartete_ketten=lauf.system.struktur.ketten,
        erlaubte_sonderreste=set(lauf.system.struktur.behalten_hetatm),
        modifizierter_rest=lauf.system.modifikation.name if lauf.system.modifikation else None,
        einstellungen=GateSettings(),
    )
    bericht = check_structure(ctx)
    if args.json:
        print(json.dumps(bericht.als_dict(), indent=2, ensure_ascii=False))
    else:
        print(bericht.als_text())
        print()
        print(f"Gesamtdauer laut Protokoll: {lauf.gesamtdauer_ns} ns")
        print("Bestanden." if bericht.ok else
              "NICHT bestanden. Die Pipeline würde hier anhalten.")
    return 0 if bericht.ok else 2


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="oligomd", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    unter = parser.add_subparsers(dest="befehl", required=True)

    p = unter.add_parser("doctor", help="Umgebung prüfen")
    p.add_argument("--phase", type=int, default=1, choices=[0, 1, 2, 3, 4],
                   help="bis zu welcher Phase alles da sein muss (Vorgabe 1)")
    p.set_defaults(func=cmd_doctor)

    p = unter.add_parser("info", help="Was steckt in einer PDB-Datei")
    p.add_argument("pdb")
    p.add_argument("--modell", type=int, default=1)
    p.add_argument("--json", action="store_true")
    p.set_defaults(func=cmd_info)

    p = unter.add_parser("vorbereiten", help="Lauf anlegen und System mit tleap bauen")
    p.add_argument("--name", default="lauf", help="Name des Lauf-Ordners")
    p.add_argument("--system", default="config/system.yaml")
    p.add_argument("--protokoll", default="config/protocol.yaml")
    p.add_argument("--json", action="store_true")
    p.set_defaults(func=cmd_vorbereiten)

    p = unter.add_parser("simulieren", help="Minimierung, Äquilibrierung, Produktion")
    p.add_argument("--lauf", default=None, help="Lauf-Ordner; ohne Angabe der letzte")
    p.add_argument("--neu", action="store_true",
                   help="von vorn beginnen statt an einem Prüfpunkt fortsetzen")
    p.add_argument("--json", action="store_true")
    p.set_defaults(func=cmd_simulieren)

    p = unter.add_parser("laeufe", help="Welche Läufe es gibt und wie weit sie sind")
    p.add_argument("--ordner", default="runs")
    p.set_defaults(func=cmd_laeufe)

    p = unter.add_parser("pruefe", help="Struktur gegen die Prüfpunkte halten")
    p.add_argument("--system", default="config/system.yaml")
    p.add_argument("--protokoll", default="config/protocol.yaml")
    p.add_argument("--json", action="store_true")
    p.set_defaults(func=cmd_pruefe)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
