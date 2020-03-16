import os
import warnings
from pygen_structures.convenience_functions import (
    load_charmm_dir,
    pdb_to_mol
)
from pygen_structures import __main__ as cmd_interface

def test_arg_parsing():
    argv = ["HEY", "-o", "HEY_out", "--histidine", "HSP"]
    args = cmd_interface.parse_args(argv)

    assert(args.sequence == "HEY")
    assert(args.segid == "PROT")
    assert(args.patches == None)
    assert(args.toppar == None)
    assert(args.verify == True)
    assert(args.output == "HEY_out")
    assert(args.histidine == "HSP")
    assert(args.use_charmm_names == False)

    argv = [
        "-u", "HSE-TRP-LYS", "-o", "HWK", "--patches", "CT2", "LAST",
        "-v", "--segid", "HWK"
    ]
    args = cmd_interface.parse_args(argv)
    assert(args.sequence == "HSE-TRP-LYS")
    assert(args.segid == "HWK")
    assert(args.patches == ["CT2", "LAST"])
    assert(args.toppar == None)
    assert(args.verify == False)
    assert(args.output == "HWK")
    assert(args.histidine == "HSE")
    assert(args.use_charmm_names == True)

def test_molecule_creation():
    argv = [
        "-u", "AGLC-BFRU-AGAL", "-o", "RAFF",
        "--patches", "RAFF", "0", "1", "2",
        "--segid", "RAFF", "--name", "Raffinose"
    ]
    cmd_interface.main(argv)

    assert(os.path.exists("RAFF.psf"))
    assert(os.path.exists("RAFF.pdb"))

    rtf, prm = load_charmm_dir()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        molecule = pdb_to_mol("RAFF.pdb", rtf, patches={"RAFF": (0, 1, 2)})

    assert(molecule.name == "Raffinose")
    assert(molecule.segment == "RAFF")
    assert(molecule.check_parameters(prm))

    ref_atoms = {
        (0, 'C1'),
        (0, 'H1'),
        (0, 'O1'),
        (0, 'C5'),
        (0, 'H5'),
        (0, 'O5'),
        (0, 'C2'),
        (0, 'H2'),
        (0, 'O2'),
        (0, 'HO2'),
        (0, 'C3'),
        (0, 'H3'),
        (0, 'O3'),
        (0, 'HO3'),
        (0, 'C4'),
        (0, 'H4'),
        (0, 'O4'),
        (0, 'HO4'),
        (0, 'C6'),
        (0, 'H61'),
        (0, 'H62'),
        (1, 'O5'),
        (1, 'C2'),
        (1, 'C5'),
        (1, 'H5'),
        (1, 'C6'),
        (1, 'H61'),
        (1, 'H62'),
        (1, 'O6'),
        (1, 'HO6'),
        (1, 'C1'),
        (1, 'H11'),
        (1, 'H12'),
        (1, 'O1'),
        (1, 'HO1'),
        (1, 'C3'),
        (1, 'H3'),
        (1, 'O3'),
        (1, 'HO3'),
        (1, 'C4'),
        (1, 'H4'),
        (1, 'O4'),
        (1, 'HO4'),
        (2, 'C1'),
        (2, 'H1'),
        (2, 'O1'),
        (2, 'C5'),
        (2, 'H5'),
        (2, 'O5'),
        (2, 'C2'),
        (2, 'H2'),
        (2, 'O2'),
        (2, 'HO2'),
        (2, 'C3'),
        (2, 'H3'),
        (2, 'O3'),
        (2, 'HO3'),
        (2, 'C4'),
        (2, 'H4'),
        (2, 'O4'),
        (2, 'HO4'),
        (2, 'C6'),
        (2, 'H61'),
        (2, 'H62'),
        (2, 'O6'),
        (2, 'HO6')
    }

    ref_bonds = {
        ((0, 'C1'), (0, 'O1')),
        ((0, 'C1'), (0, 'H1')),
        ((0, 'C1'), (0, 'O5')),
        ((0, 'C1'), (0, 'C2')),
        ((0, 'C2'), (0, 'H2')),
        ((0, 'C2'), (0, 'O2')),
        ((0, 'O2'), (0, 'HO2')),
        ((0, 'C2'), (0, 'C3')),
        ((0, 'C3'), (0, 'H3')),
        ((0, 'C3'), (0, 'O3')),
        ((0, 'O3'), (0, 'HO3')),
        ((0, 'C3'), (0, 'C4')),
        ((0, 'C4'), (0, 'H4')),
        ((0, 'C4'), (0, 'O4')),
        ((0, 'O4'), (0, 'HO4')),
        ((0, 'C4'), (0, 'C5')),
        ((0, 'C5'), (0, 'H5')),
        ((0, 'C5'), (0, 'C6')),
        ((0, 'C6'), (0, 'H61')),
        ((0, 'C6'), (0, 'H62')),
        ((0, 'C5'), (0, 'O5')),
        ((0, 'O1'), (1, 'C2')),
        ((1, 'O5'), (1, 'C2')),
        ((1, 'C2'), (1, 'C1')),
        ((1, 'C2'), (1, 'C3')),
        ((1, 'C3'), (1, 'H3')),
        ((1, 'C3'), (1, 'O3')),
        ((1, 'O3'), (1, 'HO3')),
        ((1, 'C3'), (1, 'C4')),
        ((1, 'C4'), (1, 'H4')),
        ((1, 'C4'), (1, 'O4')),
        ((1, 'O4'), (1, 'HO4')),
        ((1, 'C4'), (1, 'C5')),
        ((1, 'C5'), (1, 'H5')),
        ((1, 'C5'), (1, 'C6')),
        ((1, 'C5'), (1, 'O5')),
        ((1, 'C6'), (1, 'H61')),
        ((1, 'C6'), (1, 'H62')),
        ((1, 'C6'), (1, 'O6')),
        ((1, 'O6'), (1, 'HO6')),
        ((1, 'C1'), (1, 'H11')),
        ((1, 'C1'), (1, 'H12')),
        ((1, 'C1'), (1, 'O1')),
        ((1, 'O1'), (1, 'HO1')),
        ((2, 'C1'), (2, 'O1')),
        ((2, 'C1'), (2, 'H1')),
        ((2, 'C1'), (2, 'O5')),
        ((2, 'C1'), (2, 'C2')),
        ((2, 'C2'), (2, 'H2')),
        ((2, 'C2'), (2, 'O2')),
        ((2, 'O2'), (2, 'HO2')),
        ((2, 'C2'), (2, 'C3')),
        ((2, 'C3'), (2, 'H3')),
        ((2, 'C3'), (2, 'O3')),
        ((2, 'O3'), (2, 'HO3')),
        ((2, 'C3'), (2, 'C4')),
        ((2, 'C4'), (2, 'H4')),
        ((2, 'C4'), (2, 'O4')),
        ((2, 'O4'), (2, 'HO4')),
        ((2, 'C4'), (2, 'C5')),
        ((2, 'C5'), (2, 'H5')),
        ((2, 'C5'), (2, 'C6')),
        ((2, 'C6'), (2, 'H61')),
        ((2, 'C6'), (2, 'H62')),
        ((2, 'C6'), (2, 'O6')),
        ((2, 'O6'), (2, 'HO6')),
        ((2, 'C5'), (2, 'O5')),
        ((2, 'O1'), (0, 'C6')),
    }

    atoms = set()
    for atom in molecule.atoms:
        atoms.add((atom.residue_number - 1, atom.atom_name))
    assert(atoms == ref_atoms)

    bonds = set()
    for residue in molecule.residues:
        for bond in residue.bonds:
            if bond in ref_bonds:
                bonds.add(bond)
            else:
                bonds.add((bond[1], bond[0]))
    assert(bonds == ref_bonds)
    os.remove('RAFF.psf')
    os.remove('RAFF.pdb')