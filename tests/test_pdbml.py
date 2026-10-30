"""PDBML lesen: so kommen Strukturen von der RCSB (…/7PDU.xml.gz)."""

import gzip

import pytest

from oligomd.structure import (
    StructureError, inventory, is_pdbml, missing_atoms, models_in, net_charge, read_structure,
)

KOPF = ('<?xml version="1.0" encoding="UTF-8" ?>\n'
        '<PDBx:datablock datablockName="TEST"\n'
        '   xmlns:PDBx="http://pdbml.pdb.org/schema/pdbx-v50.xsd"\n'
        '   xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">\n'
        '   <PDBx:atom_siteCategory>\n')
FUSS = "   </PDBx:atom_siteCategory>\n</PDBx:datablock>\n"


def atom_xml(nr, name, resname, chain, resnum, xyz, element, modell=1, altloc=None):
    x, y, z = xyz
    alt = (f"<PDBx:label_alt_id>{altloc}</PDBx:label_alt_id>" if altloc
           else '<PDBx:label_alt_id xsi:nil="true" />')
    return f"""      <PDBx:atom_site id="{nr}">
         <PDBx:Cartn_x>{x}</PDBx:Cartn_x>
         <PDBx:Cartn_y>{y}</PDBx:Cartn_y>
         <PDBx:Cartn_z>{z}</PDBx:Cartn_z>
         <PDBx:auth_asym_id>{chain}</PDBx:auth_asym_id>
         <PDBx:auth_atom_id>{name}</PDBx:auth_atom_id>
         <PDBx:auth_comp_id>{resname}</PDBx:auth_comp_id>
         <PDBx:auth_seq_id>{resnum}</PDBx:auth_seq_id>
         <PDBx:group_PDB>ATOM</PDBx:group_PDB>
         {alt}
         <PDBx:pdbx_PDB_model_num>{modell}</PDBx:pdbx_PDB_model_num>
         <PDBx:type_symbol>{element}</PDBx:type_symbol>
      </PDBx:atom_site>
"""


def schreibe_xml(pfad, atome, gezippt=False):
    text = KOPF + "".join(atome) + FUSS
    if gezippt:
        with gzip.open(pfad, "wt", encoding="utf-8") as datei:
            datei.write(text)
    else:
        pfad.write_text(text, encoding="utf-8")
    return pfad


@pytest.fixture
def mini_xml(tmp_path):
    atome = []
    nr = 1
    for modell in (1, 2):
        for i, (name, element) in enumerate([("P", "P"), ("O5'", "O"), ("C1'", "C"), ("O6", "O")]):
            atome.append(atom_xml(nr, name, "DG", "A", 1, (i * 5.0, modell * 1.0, 0.0),
                                  element, modell))
            nr += 1
    return schreibe_xml(tmp_path / "mini.xml", atome)


def test_erkennt_pdbml_am_namen(tmp_path):
    assert is_pdbml(tmp_path / "7pdu.xml")
    assert is_pdbml(tmp_path / "7pdu.xml.gz")
    assert not is_pdbml(tmp_path / "7pdu.pdb")


def test_liest_modelle_und_auth_nummerierung(mini_xml):
    assert models_in(mini_xml) == 2
    s = read_structure(mini_xml, 2)
    assert s.modell == 2 and len(s.atoms) == 4
    assert s.residue("A", 1).name == "DG"
    assert s.residue("A", 1).atom("O6") is not None
    # Modell 2 liegt 1 Å höher als Modell 1
    assert s.residue("A", 1).atom("P").xyz[1] == pytest.approx(2.0)


def test_gzip_wird_ausgepackt(tmp_path):
    atome = [atom_xml(1, "P", "DG", "A", 1, (0.0, 0.0, 0.0), "P")]
    pfad = schreibe_xml(tmp_path / "mini.xml.gz", atome, gezippt=True)
    assert read_structure(pfad).atoms[0].name == "P"


def test_nur_erste_alternativlage(tmp_path):
    atome = [
        atom_xml(1, "P", "DG", "A", 1, (0.0, 0.0, 0.0), "P", altloc="A"),
        atom_xml(2, "P", "DG", "A", 1, (0.3, 0.0, 0.0), "P", altloc="B"),
    ]
    s = read_structure(schreibe_xml(tmp_path / "alt.xml", atome))
    assert len(s.atoms) == 1


def test_unbekanntes_modell_meldet_den_bereich(mini_xml):
    with pytest.raises(StructureError, match="2 Modell"):
        read_structure(mini_xml, 9)


def test_kaputtes_xml_und_mmcif(tmp_path):
    (tmp_path / "kaputt.xml").write_text("<PDBx:datablock>", encoding="utf-8")
    with pytest.raises(StructureError, match="gültiges XML"):
        read_structure(tmp_path / "kaputt.xml")
    (tmp_path / "s.cif").write_text("data_test", encoding="utf-8")
    with pytest.raises(StructureError, match="mmCIF"):
        read_structure(tmp_path / "s.cif")


def test_dna_und_rna_werden_unterschieden(tmp_path):
    atome = [
        atom_xml(1, "P", "DG", "A", 1, (0.0, 0.0, 0.0), "P"),
        atom_xml(2, "P", "G", "B", 2, (20.0, 0.0, 0.0), "P"),
        atom_xml(3, "MG", "MG", "C", 3, (40.0, 0.0, 0.0), "MG"),
    ]
    s = read_structure(schreibe_xml(tmp_path / "mix.xml", atome))
    info = inventory(s)
    assert info["reste"] == {"dna": 1, "rna": 1, "ion": 1}
    assert info["ionen"] == {"MG": 1}
    assert net_charge(s) == 0          # zwei Phosphate −2, ein Mg +2
    assert len(missing_atoms(s)) == 2  # beide Nukleotide sind nur Fragmente
