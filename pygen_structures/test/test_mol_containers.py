from pygen_structures.mol_containers import Atom, Structure, Molecule

PDB_LINE = (
    "ATOM      1  N   ALA     1      -0.992  -0.750  -0.465"
    "  1.00  0.00      PROT N1+\n"
)

PSF_LINE = (
    "         1 PROT     1        ALA      N        NH3     "
    "-0.300000       14.0070           0\n"
)

def test_atom_to_pdb():
    atom = Atom(
        atom_serial = 1,
        atom_name = "N",
        residue_name = "ALA",
        residue_number = 1,
        x = -0.992,
        y = -0.750,
        z = -0.465,
        segment_id = "PROT",
        formal_charge = 1
    )
    assert atom.to_pdb_line() == PDB_LINE

def test_atom_to_psf():
    atom = Atom(
        atom_serial = 1,
        atom_name = "N",
        residue_name = "ALA",
        residue_number = 1,
        segment_id = "PROT",
        atom_type = "NH3",
        partial_charge = -0.3,
        mass = 14.07
    )

def test_atom_from_pdb():
    atom = Atom(
        atom_serial = 1,
        atom_name = "N",
        residue_name = "ALA",
        residue_number = 1,
        x = -0.992,
        y = -0.750,
        z = -0.465,
        segment_id = "PROT",
        formal_charge = 1
    )

    assert atom.from_pdb_line(PDB_LINE).__dict__ == atom.__dict__

