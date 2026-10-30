"""Ersatz für tleap in den Tests: schreibt dieselben Dateien, meldet dasselbe.

Über die Umgebungsvariable FAKE_TLEAP_MODUS steuerbar:
  ok                  normaler Durchlauf
  fehlende_parameter  meldet einen unbekannten Rest
  fatal               bricht ab
  kein_wasser         solvatisiert nicht
  geladen             lässt das System geladen
"""
import os
import sys

modus = os.environ.get("FAKE_TLEAP_MODUS", "ok")
skript = open(sys.argv[sys.argv.index("-f") + 1], encoding="utf-8").read()
n_wasser = int(os.environ.get("FAKE_TLEAP_WASSER", "9000"))

print("-I: Adding /fake/dat/leap/prep to search path.")
for zeile in skript.splitlines():
    if zeile.startswith("source"):
        print(f"----- Source: {zeile.split()[1]}")
print("Loading PDB file: ./eingabe.pdb")
print("  total atoms in file: 1077")

if modus == "fatal":
    print("FATAL:  Atom .R<DG 23>.A<XX 7> does not have a type.")
    sys.exit(1)
if modus == "fehlende_parameter":
    print("Unknown residue: G6X   number: 22   type: Nonterminal")
    print("Could not find bond parameter for: CT - XX")

solvatisiert = "solvateoct" in skript or "solvatebox" in skript
if solvatisiert and modus != "kein_wasser":
    print(f"Added {n_wasser} residues.")
    with open("system_solvatisiert.pdb", "w", encoding="utf-8") as datei:
        datei.write("REMARK fake\n")
        nr = 1
        for i in range(n_wasser):
            for name, element in (("O", "O"), ("H1", "H"), ("H2", "H")):
                datei.write(
                    f"ATOM  {nr % 100000:5d}  {name:<3s} WAT W{(i % 9999) + 1:4d}    "
                    f"{i * 0.01:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00          {element:>2s}\n")
                nr += 1
        datei.write("END\n")

if "saveamberparm" in skript:
    ladung = -50.0 if modus == "geladen" else 0.0
    print(f"Total unperturbed charge:  {ladung:.6f}")
    print(f"Total perturbed charge:    {ladung:.6f}")
    if modus != "fehlende_parameter":
        for name, inhalt in (("system.prmtop", "%VERSION  fake\n%FLAG TITLE\n"),
                             ("system.inpcrd", "fake\n  1077\n")):
            with open(name, "w", encoding="utf-8") as datei:
                datei.write(inhalt * 200)
print("\tQuit")
