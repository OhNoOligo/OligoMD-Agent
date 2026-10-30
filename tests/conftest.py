"""Hilfsmittel für die Tests: kleine PDB-Dateien ohne externe Programme bauen."""

from __future__ import annotations

import pytest

from oligomd.structure import BACKBONE, BACKBONE_5P, BASES, DNA_RESIDUES


def pdb_line(serial: int, name: str, resname: str, chain: str, resnum: int,
             xyz: tuple[float, float, float], element: str = "") -> str:
    """Eine ATOM-Zeile im PDB-Format, spaltengenau."""
    element = element or name.lstrip("0123456789")[:1]
    x, y, z = xyz
    return (
        f"ATOM  {serial:5d} {name:^4s} {resname:>3s} {chain:1s}{resnum:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          {element:>2s}"
    )


def build_residue(resname: str, chain: str, resnum: int, origin: tuple[float, float, float],
                  serial_start: int, erstes: bool = False,
                  weglassen: tuple[str, ...] = ()) -> tuple[list[str], int]:
    """Vollständiges Standardnukleotid; die Atome liegen auf einem weiten Raster."""
    base = DNA_RESIDUES[resname]
    namen = list(BACKBONE) + list(BASES[base])
    if not erstes:
        namen = list(BACKBONE_5P) + namen
    namen = [n for n in namen if n not in weglassen]

    zeilen, serial = [], serial_start
    for i, name in enumerate(namen):
        # 6 Å Raster: weit genug auseinander, damit nichts als Kollision zählt und
        # Testatome eindeutig nur einem Nachbarn nahe kommen
        versatz = (origin[0] + (i % 4) * 6.0, origin[1] + ((i // 4) % 4) * 6.0,
                   origin[2] + (i // 16) * 6.0)
        zeilen.append(pdb_line(serial, name, resname, chain, resnum, versatz))
        serial += 1
    return zeilen, serial


def build_pdb(reste: list[tuple[str, str, int]], modelle: int = 1,
              weglassen: dict[tuple[str, int], tuple[str, ...]] | None = None,
              extra: list[str] | None = None) -> str:
    """PDB-Text aus (Restname, Kette, Nummer). Reste liegen weit auseinander."""
    weglassen = weglassen or {}
    erste_je_kette: dict[str, int] = {}
    for resname, chain, resnum in reste:
        if chain not in erste_je_kette or resnum < erste_je_kette[chain]:
            erste_je_kette[chain] = resnum

    bloecke = []
    for modell in range(1, modelle + 1):
        zeilen = [f"MODEL     {modell:4d}"] if modelle > 1 else []
        serial = 1
        for k, (resname, chain, resnum) in enumerate(reste):
            origin = (k * 60.0, (ord(chain) - ord("A")) * 300.0, 0.0)
            neu, serial = build_residue(
                resname, chain, resnum, origin, serial,
                erstes=(resnum == erste_je_kette[chain]),
                weglassen=weglassen.get((chain, resnum), ()),
            )
            zeilen += neu
        zeilen += extra or []
        if modelle > 1:
            zeilen.append("ENDMDL")
        bloecke.append("\n".join(zeilen))
    return "\n".join(bloecke) + "\nEND\n"


@pytest.fixture
def duplex_pdb(tmp_path):
    """Sauberer kleiner Doppelstrang: Kette A mit 3, Kette B mit 2 Resten."""
    text = build_pdb([
        ("DG", "A", 1), ("DA", "A", 2), ("DC", "A", 3),
        ("DT", "B", 1), ("DG", "B", 2),
    ])
    pfad = tmp_path / "duplex.pdb"
    pfad.write_text(text, encoding="utf-8")
    return pfad
