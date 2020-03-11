"""
Contains the Structure class. This is used to represent an arbitrary
3-dimensional chemical species, which may be a complete molecule or
a fragment.

3D coordinates can be generated, and atoms can be fixed.
"""

from typing import List, Tuple, Dict
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D
import numpy as np

from pygen_structures.mol_containers.atom import Atom
from pygen_structures._functions_const import (
    bonds_to_adjacency,
    BOND_ORDER_TO_TYPE,
    get_expected_valence
)


class Structure:
    """
    An arbitrary chemical structure. This may be a complete molecule,
    but could also be an isolated residue fragment.
    """
    def __init__(
            self,
            atoms: List[Atom],
            bonds: List[Tuple[int, int]],
            improper_ics: (Dict[Tuple[int, ...], float], None) = None,
            fixed_atoms: (List[int], None) = None,
            set_charges: bool = True,
            use_etkdg: bool = True,
    ):
        self.atoms = atoms
        self.adjacency_matrix = bonds_to_adjacency(bonds)
        self.bonds = bonds
        if not fixed_atoms:
            fixed_atoms = dict()
        self.fixed_atoms = fixed_atoms
        if not improper_ics:
            improper_ics = dict()
        self.improper_ics = improper_ics
        self.mol: (None, Chem.RWMol) = None
        self.use_etkdg = use_etkdg

        self.set_formal_charges()
        self.to_mol()
        self.reorder_atoms()
        self.generate_coordinates()
        if set_charges is False:
            for rd_atom in self.mol.GetAtoms():
                rd_atom.SetFormalCharge(0)

    def set_formal_charges(self) -> None:
        """
        Set Atom formal charges from valence and bond order.
        """
        for atom, connections in zip(self.atoms, self.adjacency_matrix):
            total_valence = np.sum(connections).item()
            expected_valence = get_expected_valence(atom.element_symbol)
            formal_charge = total_valence - expected_valence
            if formal_charge:
                atom.formal_charge = formal_charge

    def to_mol(self) -> Chem.RWMol:
        """
        Create an RDKit mol from the atoms and adjacency matrix.
        All stereocenters are set to clockwise initially.
        """
        if self.mol:
            return self.mol

        # noinspection PyArgumentList
        mol = Chem.RWMol()
        for atom in self.atoms:
            residue_index = atom.residue_number - 1
            rd_atom = Chem.Atom(atom.element_symbol)
            rd_atom.SetProp("atom_name", atom.atom_name)
            rd_atom.SetProp("atom_type", atom.atom_type)
            rd_atom.SetDoubleProp("partial_charge", atom.partial_charge)
            rd_atom.SetProp("residue_name", atom.residue_name)
            rd_atom.SetIntProp("residue_index", residue_index)
            atom_id = "{}:{}".format(residue_index, atom.atom_name)
            rd_atom.SetProp("atom_id", atom_id)
            rd_atom.SetFormalCharge(atom.formal_charge)
            mol.AddAtom(rd_atom)

        for bond in np.argwhere(self.adjacency_matrix):
            idx_x, idx_y = bond.tolist()
            if idx_y < idx_x:
                continue
            bond_order = self.adjacency_matrix[idx_x, idx_y]
            bond_type = BOND_ORDER_TO_TYPE[bond_order]
            mol.AddBond(idx_x, idx_y, bond_type)

        Chem.SanitizeMol(mol)
        chiral_centers = Chem.FindMolChiralCenters(mol, True, True)
        chiral_centers = [x[0] for x in chiral_centers]
        for atom_index in chiral_centers:
            rd_atom = mol.GetAtomWithIdx(atom_index)
            rd_atom.SetChiralTag(Chem.CHI_TETRAHEDRAL_CW)

        self.mol = mol
        return mol

    def generate_coordinates(self) -> None:
        """
        Generate a 3D conformation for the Structure using RDKit.
        This may take a few attempts to set stereochemistry correctly.
        """
        mol = self.to_mol()

        fixed_coords = {}
        for atom_index in self.fixed_atoms:
            atom = self.atoms[atom_index]
            coordinates = (atom.x, atom.y, atom.z)
            fixed_coords[atom_index] = Point3D(*coordinates)
        if self.use_etkdg:
            embed_parameters = AllChem.ETKDGv2()
        else:
            # noinspection PyArgumentList
            embed_parameters = AllChem.EmbedParameters()
            embed_parameters.useBasicKnowledge = True
            embed_parameters.ignoreSmoothingFailures = True
            embed_parameters.verbose = True
        if fixed_coords:
            embed_parameters.coordMap = fixed_coords
        AllChem.EmbedMolecule(mol, embed_parameters)

        centers_to_invert, previous_centers = set(), set()
        chiral_centers = Chem.FindMolChiralCenters(mol, True, True)
        chiral_centers = [x[0] for x in chiral_centers]

        n_cycles = 0
        while n_cycles < 100:
            n_cycles += 1
            # noinspection PyArgumentList
            conformer = mol.GetConformer()
            for (i, j, k, l), improper in self.improper_ics.items():
                if k not in chiral_centers:
                    continue
                measured_impr = AllChem.GetDihedralDeg(conformer, i, j, k, l)
                if (
                        (improper > 0 and measured_impr < 0) or
                        (improper < 0 and measured_impr > 0)
                ):
                    centers_to_invert.add(k)
            if not centers_to_invert:
                break
            # In case the chirality of one center is dependent upon another
            # This doesn't happen very often.
            if centers_to_invert == previous_centers:
                # Randomly pop one center and randomly add neighboring
                # chiral centers to the mix.
                center = centers_to_invert.pop()
                rd_atom = mol.GetAtomWithIdx(center)
                for neighbor in rd_atom.GetNeighbors():
                    neighbor_index = neighbor.GetIdx()
                    if neighbor_index in chiral_centers:
                        if np.random.random() >= 0.5:
                            centers_to_invert.add(neighbor_index)
                if not centers_to_invert:
                    centers_to_invert.add(center)
            for atom_index in centers_to_invert:
                rd_atom = mol.GetAtomWithIdx(atom_index)
                rd_atom.InvertChirality()
            previous_centers, centers_to_invert = centers_to_invert, set()
            AllChem.EmbedMolecule(mol, embed_parameters)
        else:
            raise ValueError("Unable to set stereochemistry correctly.")

        # noinspection PyArgumentList
        conformer = mol.GetConformer()
        for atom_index, atom in enumerate(self.atoms):
            position = conformer.GetAtomPosition(atom_index)
            atom.x = position.x
            atom.y = position.y
            atom.z = position.z

    def reorder_atoms(self) -> None:
        """
        Reorder the atoms and adjacency matrix to put the hydrogen
        atoms after their parent heavy atom. This is probably an easier
        layout for a person to follow, should they look at the file.
        """
        new_atom_order: List[int] = []
        for rd_atom in self.mol.GetAtoms():
            if rd_atom.GetSymbol() == "H":
                continue
            index = rd_atom.GetIdx()
            new_atom_order.append(index)
            neighbors = rd_atom.GetNeighbors()
            hydrogens = [n for n in neighbors if n.GetSymbol() == "H"]
            hydrogen_indices = [n.GetIdx() for n in hydrogens]
            new_atom_order.extend(hydrogen_indices)
        self.mol = Chem.RenumberAtoms(self.mol, new_atom_order)

        self.atoms[:] = [self.atoms[x] for x in new_atom_order]
        for new_index, atom in enumerate(self.atoms):
            atom.atom_serial = new_index + 1

        for row in self.adjacency_matrix:
            row[:] = [row[x] for x in new_atom_order]
        for column in self.adjacency_matrix.T:
            column[:] = [column[x] for x in new_atom_order]

        for old_improper_indices in list(self.improper_ics.keys()):
            (i, j, k, l) = old_improper_indices
            i, j = new_atom_order[i], new_atom_order[j]
            k, l = new_atom_order[k], new_atom_order[l]
            improper = self.improper_ics.pop(old_improper_indices)
            self.improper_ics[(i, j, k, l)] = improper
