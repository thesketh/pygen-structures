import os
import pickle
from pygen_structures import(
    load_charmm_dir,
    code_to_mol,
    sequence_to_mol,
    Molecule,
)
from pygen_structures.convenience_functions import pdb_to_mol

RTF, PRM = load_charmm_dir()


def test_load_charmm_dir():
    pickle_dir = os.path.split(__file__)[0]
    pickle_path = os.path.join(pickle_dir, 'charmm_dir.pkl')

    with open(pickle_path, 'rb') as pickle_file:
        rtf, prm = pickle.load(pickle_file)
    assert(rtf.residues.keys() == RTF.residues.keys())
    assert(rtf.patches.keys() == RTF.patches.keys())
    assert(prm.bonds == PRM.bonds)
    assert(prm.angles == PRM.angles)
    assert(prm.dihedrals == PRM.dihedrals)
    assert(prm.impropers == PRM.impropers)
    assert(prm.cross_maps == PRM.cross_maps)


def test_code_to_mol():
    code = 'YEET'
    molecule = code_to_mol(code, RTF)

    sequence = ['TYR', 'GLU', 'GLU', 'THR']
    residues = []
    for index, residue_name in enumerate(sequence):
        residue = RTF.residues[residue_name].to_residue(index)
        residues.append(residue)
    ref_molecule = Molecule('YEET', residues, RTF)
    ref_molecule.finalize()

    for residue, ref_residue in zip(molecule.residues, ref_molecule.residues):
        assert(residue.name == ref_residue.name)
        assert(set(residue.atoms) == set(ref_residue.atoms))
        assert(set(residue.bonds) == set(ref_residue.bonds))
        assert(set(residue.impropers) == set(ref_residue.impropers))
        assert(set(residue.cross_maps) == set(ref_residue.cross_maps))


def test_code_to_mol_terminal_patches():
    code = 'NNEU-YEET-CT2'
    molecule = code_to_mol(code, RTF)

    sequence = ['TYR', 'GLU', 'GLU', 'THR']
    patches = {'NNEU': ['FIRST'], 'CT2': ['LAST']}

    residues = []
    for index, residue_name in enumerate(sequence):
        residue = RTF.residues[residue_name].to_residue(index)
        residues.append(residue)
    ref_molecule = Molecule('YEET', residues, RTF, patches)
    ref_molecule.finalize()

    for residue, ref_residue in zip(molecule.residues, ref_molecule.residues):
        assert(residue.name == ref_residue.name)
        assert(set(residue.atoms) == set(ref_residue.atoms))
        assert(set(residue.bonds) == set(ref_residue.bonds))
        assert(set(residue.impropers) == set(ref_residue.impropers))
        assert(set(residue.cross_maps) == set(ref_residue.cross_maps))


def test_code_to_mol_patches():
    code = 'CAAC'
    patches = {'DISU': [0, -1]}
    molecule = code_to_mol(code, RTF, patches)

    sequence = ['CYS', 'ALA', 'ALA', 'CYS']
    residues = []
    for index, residue_name in enumerate(sequence):
        residue = RTF.residues[residue_name].to_residue(index)
        residues.append(residue)
    ref_molecule = Molecule('CAAC', residues, RTF, patches)
    ref_molecule.finalize()

    for residue, ref_residue in zip(molecule.residues, ref_molecule.residues):
        assert(residue.name == ref_residue.name)
        assert(set(residue.atoms) == set(ref_residue.atoms))
        assert(set(residue.bonds) == set(ref_residue.bonds))
        assert(set(residue.impropers) == set(ref_residue.impropers))
        assert(set(residue.cross_maps) == set(ref_residue.cross_maps))


def test_sequence_to_mol():
    sequence = ['TYR', 'GLU', 'GLU', 'THR']
    molecule = sequence_to_mol(sequence, RTF)

    residues = []
    for index, residue_name in enumerate(sequence):
        residue = RTF.residues[residue_name].to_residue(index)
        residues.append(residue)
    ref_molecule = Molecule('YEET', residues, RTF)
    ref_molecule.finalize()

    for residue, ref_residue in zip(molecule.residues, ref_molecule.residues):
        assert(residue.name == ref_residue.name)
        assert(set(residue.atoms) == set(ref_residue.atoms))
        assert(set(residue.bonds) == set(ref_residue.bonds))
        assert(set(residue.impropers) == set(ref_residue.impropers))
        assert(set(residue.cross_maps) == set(ref_residue.cross_maps))


def test_pdb_to_mol():
    code = 'AF'

    ref_molecule = code_to_mol(code, RTF)
    ref_molecule.to_pdb_file('AF.pdb')

    molecule = pdb_to_mol('AF.pdb', RTF)
    os.remove('AF.pdb')

    for residue, ref_residue in zip(molecule.residues, ref_molecule.residues):
        assert(residue.name == ref_residue.name)
        assert(set(residue.atoms) == set(ref_residue.atoms))
        assert(set(residue.bonds) == set(ref_residue.bonds))
        assert(set(residue.impropers) == set(ref_residue.impropers))
        assert(set(residue.cross_maps) == set(ref_residue.cross_maps))
