"""
Module containing the Molecule class.

This class is used to create CHARMM structures from a list of RESI
names and a dict of PRES patch names to a list of indices of the atoms
they are to be applied to (or "FIRST"/"LAST").
"""

from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
from typing import Tuple, List, Dict

from charmm_containers import *
from atom import Atom
from _functions_const import *


class Molecule:
    def __init__(
        self,
        name: str = "Unset",
        residues: (list, None) = None,
        topology: (CHARMMResidueTopologyFile, None) = None,
        patches: (Dict[str, List[int]], None) = None,
        segment: (str, None) = None,
    ):
        self.name = name
        if not segment:
            self.segment = name[:4].upper()
        self.residues = residues if residues else list()
        self.patches = patches if patches else dict()
        self.topology = topology if topology else CHARMMResidueTopologyFile()

        # Set during finalization
        self.atoms: (List[Atom], None) = None
        self.adjacency_matrix: (np.ndarray, None) = None
        self.impropers: (List[Tuple[int, ...]], None) = None
        self.cross_maps: (List[Tuple[int, ...]], None) = None
        self.mol: (Chem.RWMol, None) = None
        self.topology_files = set()

        self._id_to_index: (Dict[str, int], None) = None
        self._finalized: bool = False


    def finalize(self) -> None:
        """
        Finalize the molecule, getting it to a state where it can
        be represented directly as a PDB file or PSF file.
        """
        if self._finalized:
            raise ValueError("Molecule already finalized.")
        self._apply_patches()
        self._build_atoms()
        self._remove_dangling_bonds()
        self._build_adjacency_matrix()
        self._set_formal_charges()
        self._reorder_atoms()
        self._build_adjacency_matrix()  # Reorder the connection table
        self._fill_impropers_cross_maps()
        self._generate_coordinates()
        self._finalized = True

    def _apply_patches(self) -> None:
        """
        Apply the patches to the residues. For proteins, for instance,
        C and N termini usually need to be deprotonated and protonated.
        """
        first_patch, last_patch = False, False
        patches = {}
        for patch_name, locations in self.patches.items():
            residues = []
            for location in locations:
                if location == "FIRST":
                    location = 0
                    first_patch = True
                elif location == "LAST":
                    location = -1
                    last_patch = True
                residue = self.residues[location]
                residues.append(residue)
            patches[patch_name] = residues
        if not first_patch:
            first_residue = self.residues[0]
            first_patch = first_residue.first
            patches[first_patch] = [first_residue]
        if not last_patch:
            last_residue = self.residues[-1]
            last_patch = last_residue.last
            patches[last_patch] = [last_residue]

        for patch_name, residues in patches.items():
            if patch_name == "NONE":
                continue
            patch = self.topology.patches[patch_name]
            self.topology_files.add(patch.rtf_file_name)
            patch.apply(*residues)

    def _build_atoms(self) -> None:
        """
        Build self.atoms, the list of Atoms, from the data present in
        each residue.
        """
        self.atoms, self._id_to_index = list(), dict()
        atom_index = 0
        for residue in self.residues:
            self.topology_files.add(residue.rtf_file_name)
            for atom_data in residue.atoms:
                atom_name, atom_type, charge = atom_data
                atom_id = "{}:{}".format(residue.index, atom_name)
                self._id_to_index[atom_id] = atom_index
                atom_name = atom_id.split(":")[1]
                mass = self.topology.masses[atom_type]
                atom = Atom(
                    atom_serial=atom_index + 1,
                    atom_name=atom_name,
                    residue_name=residue.name,
                    residue_number=residue.index + 1,
                    segment_id=self.segment,
                    atom_type=atom_type,
                    partial_charge=charge,
                    mass=mass,
                )
                self.atoms.append(atom)
                atom_index += 1

    def _remove_dangling_bonds(self) -> None:
        """
        Remove all bonds, impropers, cross-maps, and ICs involving
        atoms which don't exist.
        """
        first_residue = self.residues[0]
        last_residue = self.residues[-1]
        final_index = last_residue.index
        for residue in self.residues:
            bonds, impropers, cross_maps, ics = [], [], [], []
            for bond in residue.bonds:
                for atom_id in bond:
                    if atom_id not in self._id_to_index:
                        break
                else:
                    bonds.append(bond)
            for improper in residue.impropers:
                for atom_id in improper:
                    if atom_id not in self._id_to_index:
                        break
                else:
                    impropers.append(improper)
            for cross_map in residue.cross_maps:
                for atom_id in cross_map:
                    if atom_id not in self._id_to_index:
                        break
                else:
                    cross_maps.append(cross_map)
            for ic in residue.ics:
                for atom_id in ic[:4]:
                    if atom_id.replace("*", "") not in self._id_to_index:
                        break
                else:
                    ics.append(ic)
            residue.bonds = bonds
            residue.impropers = impropers
            residue.cross_maps = cross_maps
            residue.ics = ics

    def _build_adjacency_matrix(self) -> None:
        """
        Build an adjacency matrix from the bonds in each residue.

        This is an array of integers which indicate bond order.

        E.g. self.adjacency_matrix[x, y] == 2 would indicate a double
        bond between atoms with index x and y.
        """
        n_atoms = len(self.atoms)
        self.adjacency_matrix = np.zeros((n_atoms, n_atoms), dtype=int)
        bonds = []
        for residue in self.residues:
            for bond in residue.bonds:
                bond = [self._id_to_index[atom_id] for atom_id in bond]
                bonds.append(bond)

        for idx_x, idx_y in bonds:
            self.adjacency_matrix[idx_x, idx_y] += 1
            self.adjacency_matrix[idx_y, idx_x] += 1

    def _set_formal_charges(self) -> None:
        """
        Set Atom formal charges from valence and bond order.
        """
        for atom, connections in zip(self.atoms, self.adjacency_matrix):
            total_valence = np.sum(connections).item()
            expected_valence = get_expected_valence(atom.element_symbol)
            formal_charge = total_valence - expected_valence
            if formal_charge:
                atom.formal_charge = formal_charge

    def _reorder_atoms(self) -> None:
        """
        Reorder the atoms so that hydrogen atoms come after their
        parent heavy atom.
        """
        mol = self._to_rwmol()

        new_atom_order = []
        for rd_atom in mol.GetAtoms():
            if rd_atom.GetSymbol() == "H":
                continue
            index = rd_atom.GetIdx()
            new_atom_order.append(index)
            neighbors = rd_atom.GetNeighbors()
            hydrogens = [n for n in neighbors if n.GetSymbol() == "H"]
            hydrogen_indices = [n.GetIdx() for n in hydrogens]
            new_atom_order.extend(hydrogen_indices)
        self.mol = Chem.RenumberAtoms(mol, new_atom_order)
        self.atoms = [self.atoms[x] for x in new_atom_order]
        for new_index, atom in enumerate(self.atoms):
            atom.atom_serial = new_index + 1
            residue_index = atom.residue_number - 1
            atom_name = atom.atom_name
            atom_id = "{}:{}".format(residue_index, atom_name)
            self._id_to_index[atom_id] = new_index
        self._build_adjacency_matrix()

    def _fill_impropers_cross_maps(self) -> None:
        """
        Fill self.impropers and self.cross_maps from residue data.
        """
        impropers, cross_maps = [], []
        for residue in self.residues:
            for improper in residue.impropers:
                impropers.append([self._id_to_index[x] for x in improper])
            for cross_map in residue.cross_maps:
                cross_maps.append([self._id_to_index[x] for x in cross_map])
        self.impropers, self.cross_maps = impropers, cross_maps

    def _generate_coordinates(self) -> None:
        """
        Generate a 3D conformation for the Molecule using RDKit.
        This may take a few attempts to set stereochemistry correctly.
        """
        mol = self._to_rwmol()
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv2())

        centers_to_invert, previous_centers = set(), set()
        n_cycles = 0

        chiral_centers = Chem.FindMolChiralCenters(mol, True, True)
        chiral_centers = [x[0] for x in chiral_centers]
        Chem.AssignStereochemistry(mol)

        while True:
            if n_cycles == 100:
                raise ValueError("Unable to set stereochemistry correctly.")
            n_cycles += 1
            conformer = mol.GetConformer()
            for residue in self.residues:
                residue_index = residue.index
                for ic in residue.ics:
                    if "*" not in ic[2]:
                        continue
                    ic_atom_ids = ic[:4].copy()
                    ic_atom_ids[2] = ic_atom_ids[2].replace("*", "")
                    i, j, k, l = [self._id_to_index[x] for x in ic_atom_ids]

                    if k not in chiral_centers:
                        continue
                    improper = ic[6]
                    measured_improper = AllChem.GetDihedralDeg(conformer, i, j, k, l)
                    if (
                        (improper > 0 and measured_improper < 0) or
                        (improper < 0 and measured_improper > 0)
                    ):
                        centers_to_invert.add(k)
            if not centers_to_invert:
                break
            # In case the chirality of one center is dependent upon another
            # This doesn't happen very often.
            if centers_to_invert == previous_centers:
                # Randomly pop one center and possibly add neighboring
                # chiral centers to the mix.
                center = centers_to_invert.pop()
                rd_atom = mol.GetAtomWithIdx(center)
                for neighbor in rd_atom.GetNeighbors():
                    neighbor_index = neighbor.GetIdx()
                    if neighbor_index in chiral_centers:
                        if np.random.random >= 0.5:
                            centers_to_invert.add(neighbor_index)
                if not centers_to_invert:
                    centers_to_invert.add(center)
            for atom_index in centers_to_invert:
                rd_atom = mol.GetAtomWithIdx(atom_index)
                rd_atom.InvertChirality()
            previous_centers, centers_to_invert = centers_to_invert, set()
            AllChem.EmbedMolecule(mol, AllChem.ETKDGv2())

        conformation = mol.GetConformer()
        for atom_index, atom in enumerate(self.atoms):
            position = conformation.GetAtomPosition(atom_index)
            atom.x = position.x
            atom.y = position.y
            atom.z = position.z

    def _to_rwmol(self) -> Chem.RWMol:
        """
        Create an RDKit mol from the atoms and adjacency matrix.
        All stereocenters are set to clockwise initially.
        """
        if not self.atoms or self.adjacency_matrix is None:
            raise ValueError("Mol must be finalized first.")
        if self.mol:
            return self.mol

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

        mol.SetProp("_name", self.name)
        self.mol = mol
        return mol

    def get_conect_records(self) -> List[str]:
        """
        Build the PDB CONECT records from the adjacency matrix.
        """
        conect_records = []
        for x_serial, row in enumerate(self.adjacency_matrix, 1):
            connections = np.flatnonzero(row) + 1
            fmt = "{: 5d}" if x_serial < 100000 else "{: 5X}"
            first_atom = fmt.format(x_serial)
            connection_list = ["CONECT", first_atom]
            for connection in connections:
                fmt = "{: 5d}" if connection < 100000 else "{: 5X}"
                connection_list.append(fmt.format(connection))
            connection_list.append("\n")
            conect_records.append("".join(connection_list))
        return conect_records

    def to_pdb_block(self) -> str:
        """
        Build the PDB records and return the PDB as a string.
        """
        if not self._finalized:
            raise ValueError("Molecule not finalized")

        header_records = []
        if self.name:
            header_records.append("COMPND    {}\n".format(self.name))
        header_records.append("AUTHOR    Generated procedurally.\n")

        atom_records = [atom.to_pdb_line() for atom in self.atoms]
        conect_records = self.get_conect_records()
        master_fmt = "MASTER    " + "    0" * 8 + "{:>5}    0{:>5}    0\n"
        master = master_fmt.format(len(atom_records), len(conect_records))
        header_block = "".join(header_records)
        atom_block = "".join(atom_records)
        conect_block = "".join(conect_records)
        return header_block + atom_block + conect_block + master + "END\n"

    def to_pdb_file(self, pdb_path: str):
        """
        Write the molecule to a PDB file.
        """
        with open(pdb_path, "w") as pdb_file:
            pdb_file.write(self.to_pdb_block())

    def to_psf_block(self) -> str:
        """
        Build the PSF records and return the PSF as a string.
        """
        if not self._finalized:
            raise ValueError("Molecule not finalized")

        def format_psf_block(
            data: List[Tuple[int, ...]],
            n_per_line: int
        ) -> str:
            try:
                n_per_item = len(data[0])
            except:
                return ""
            format_string = " ".join([" {:7d}" for _ in range(n_per_item)])
            items_per_line = n_per_line // n_per_item
            lines = []
            for index, item in enumerate(data):
                if not index % n_per_line:
                    if index != 0:
                        lines[-1] = "".join(lines[-1])
                    lines.append([])
                item = [x + 1 for x in item]
                lines[-1].append(format_string.format(*item))
            else:
                if isinstance(lines[-1], list):
                    lines[-1] = "".join(lines[-1])
            return "\n".join(lines) + "\n"

        header = (
            "\n".join(
                [
                    "PSF CMAP",
                    "",
                    "         2 !NTITLE",
                    "* Generated procedurally.",
                    "* Molecule: {}".format(self.name),
                ]
            )
            + "\n"
        )
        bonds, angles, dihedrals = adjacency_to_dof(self.adjacency_matrix)
        atom_head = " {:7d} !NATOM".format(len(self.atoms))
        atom_block = "".join([atom.to_psf_line() for atom in self.atoms])
        bond_head = " {:7d} !NBOND: bonds".format(len(bonds))
        bond_block = format_psf_block(bonds, n_per_line=4)
        angle_head = " {:7d} !NTHETA: angles".format(len(angles))
        angle_block = format_psf_block(angles, n_per_line=3)
        dihedral_head = " {:7d} !NPHI: dihedrals".format(len(dihedrals))
        dihedral_block = format_psf_block(dihedrals, n_per_line=2)
        improper_head = " {:7d} !NIMPHI: impropers".format(len(self.impropers))
        improper_block = format_psf_block(self.impropers, n_per_line=2)
        cmap_head = " {:7d} !NCRTERM: cross-terms".format(len(self.cross_maps))
        cmap_block = format_psf_block(self.cross_maps, n_per_line=1)

        psf_blocks = [
            header,
            atom_head,
            atom_block,
            bond_head,
            bond_block,
            angle_head,
            angle_block,
            dihedral_head,
            dihedral_block,
            improper_head,
            improper_block,
            cmap_head,
            cmap_block,
        ]
        return "\n".join(psf_blocks) + "\n"

    def to_psf_file(self, psf_path):
        """
        Write the molecule to a PSF file.
        """
        with open(psf_path, "w") as psf_file:
            psf_file.write(self.to_psf_block())

    def check_parameters(self, parameter_set: CHARMMParameterFile) -> bool:
        """
        Check if parameters needed by the molecule are present in a
        parameter set. Returns True if all parameters are present,
        otherwise false.
        """
        return parameter_set.check_parameters(self)

    @classmethod
    def from_sequence(
        cls,
        sequence: List[str],
        topology: CHARMMResidueTopologyFile,
        patches:  (Dict[str, List[int]], None) = None,
    ):
        """
        Create the Molecule instance from a list of residues, a
        CHARMM topology file and a list containing patches to be
        applied.
        """
        name = "".join(x[0] for x in sequence)
        instance = cls(name, patches=patches, topology=topology)

        for residue_index, residue_name in enumerate(sequence):
            residue_definition = topology.residues[residue_name]
            residue = residue_definition.to_residue(residue_index)
            instance.residues.append(residue)
        instance.finalize()
        return instance
