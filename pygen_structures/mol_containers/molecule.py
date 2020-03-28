"""
Module containing the Molecule class.

This class is used to create CHARMM structures from a list of RESI
names, a `CHARMMResidueTopologyFile` and a dict mapping PRES names to
indices of the residues they are to be applied to (or "FIRST"/"LAST").

Manipulation of 3D coordinates is performed in the `Structure` class
so that 3D functionality can be applied to residue objects.
"""
from pygen_structures.charmm_containers import *
from pygen_structures._functions_const import *
from pygen_structures.version import version
from pygen_structures.mol_containers.atom import Atom
from pygen_structures.mol_containers.structure import Structure

import numpy as np
import warnings
from typing import Tuple, List, Dict


class Molecule:
    """
    A complete CHARMM molecule. This can be written to PSF and/or PDB

    :param name: The name of the molecule, used in the PSF and PDB.
    :param residues: a list of ``CHARMMResidue``
    :param topology: the ``CHARMMResidueTopologyFile`` to provide the\
    patches.
    :param patches: a dict of patch names to a list of the residue\
    indices the patch is to be applied to (or "FIRST"/"LAST") magic\
    strings.
    :param segment: the segment ID to be used when the ``Atom`` records\
    are created.
    :param fixed_atoms: a mapping of ``atom_id`` to the fixed coords
    :param use_etkdg: In the structure, when generating coords,\
    use empirical distance generation (if ``True``) or a traditional\
    distance generation approach.

    The following attributes are set during finalization:

    :param atoms: a list of ``Atom``
    :param impropers: a list of tuples with ``Atom`` indices for the PSF
    :param cross_maps: a list of tuples with ``Atom`` indices for the PSF
    :param topology_files: a set of paths to the .rtf/.str files

    These are private attributes:

    :param _structure: the underlying ``Structure`` object
    :param _id_to_index: a mapping of ``atom_id`` to index.
    :param _finalized: ``True`` if ``Molecule.finalize()`` has been called.

    """
    def __init__(
        self,
        name: str = "Unset",
        residues: (list, None) = None,
        topology: (CHARMMResidueTopologyFile, None) = None,
        patches: (Dict[str, List[int]], None) = None,
        segment: (str, None) = None,
        fixed_atoms: (
                Dict[
                    Tuple[int, str],
                    Tuple[float, float, float]
                ], None) = None,
        use_etkdg: bool = True
    ):
        self.name = name
        self.residues = residues if residues else list()
        self.patches = patches if patches else dict()
        self.topology = topology if topology else CHARMMResidueTopologyFile()

        self.segment = segment if segment else "U"
        if len(self.segment) > 4:
            err = (
                "Segment {} is longer than 4 characters and will be truncated."
            )
            warnings.warn(err.format(self.segment))

        self.fixed_atoms = fixed_atoms if fixed_atoms else {}
        self.use_etkdg = use_etkdg

        # Set during finalization
        self.atoms: (List[Atom], None) = None
        self.impropers: (List[Tuple[int, ...]], None) = None
        self.cross_maps: (List[Tuple[int, ...]], None) = None
        self.topology_files = set()

        self._structure: (Structure, None) = None
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
        self._structure = self.to_structure()
        self._fill_impropers_cross_maps()
        self._finalized = True

    def _apply_patches(self) -> None:
        """
        Apply the patches to the residues. For proteins, for instance,
        C and N termini usually need to be deprotonated and protonated
        respectively.
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
                atom_id = (residue.index, atom_name)
                self._id_to_index[atom_id] = atom_index
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
                if atom_id in self.fixed_atoms.keys():
                    atom.x, atom.y, atom.z = self.fixed_atoms[atom_id]
                self.atoms.append(atom)
                atom_index += 1

    def _remove_dangling_bonds(self) -> None:
        """
        Remove all bonds, impropers, cross-maps, and ICs involving
        atoms which don't exist.
        """
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
                for res_index, atom_name in ic[:4]:
                    if atom_name.replace("*", "") not in self._id_to_index:
                        break
                else:
                    ics.append(ic)
            residue.bonds = bonds
            residue.impropers = impropers
            residue.cross_maps = cross_maps
            residue.ics = ics

    def to_structure(self) -> Structure:
        """
        Create a ``Structure`` using information from ``atoms`` and
        ``bonds``

        :return: the generated `Structure`

        """
        if self._structure is not None:
            return self._structure
        atoms, use_etkdg = self.atoms, self.use_etkdg
        fixed_atoms = []
        for atom_id in self.fixed_atoms.keys():
            if atom_id in self._id_to_index:
                fixed_atoms.append(self._id_to_index[atom_id])

        bonds, improper_ics = [], {}
        for residue in self.residues:
            for (atom_id_x, atom_id_y) in residue.bonds:
                bond = (
                    self._id_to_index[atom_id_x],
                    self._id_to_index[atom_id_y]
                )
                bonds.append(bond)
            for ic in residue.ics:
                res_index, atom_name = ic[2]
                if "*" not in atom_name[2]:
                    continue
                atom_name = atom_name.replace("*", "")
                k = self._id_to_index[(res_index, atom_name)]
                i, j, _, l = [self._id_to_index[a_id] for a_id in ic[:4]]
                improper = ic[6]
                improper_ics[(i, j, k, l)] = improper

        structure = Structure(
            atoms, bonds, improper_ics, fixed_atoms, True, self.use_etkdg
        )
        self.atoms = structure.atoms
        for index, atom in enumerate(self.atoms):
            atom_id = (atom.residue_number - 1, atom.atom_name)
            self._id_to_index[atom_id] = index
        return structure

    def to_mol(self):
        if not self._finalized:
            raise ValueError("Molecule must be finalized first")
        return self._structure.mol

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

    def get_adjacency_matrix(self) -> np.ndarray:
        return self._structure.adjacency_matrix

    def get_conect_records(self) -> List[str]:
        """
        Build the PDB CONECT records from the adjacency matrix.
        """
        conect_records = []
        adjacency_matrix = self._structure.adjacency_matrix
        for x_serial, row in enumerate(adjacency_matrix, 1):
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

        :return: PDB block as string

        """
        if not self._finalized:
            raise ValueError("Molecule not finalized")

        header_records = []
        if self.name:
            header_records.append("COMPND    {}\n".format(self.name))
        header_records.append("AUTHOR    pygen-structures v{}\n".format(version))
        header_records.append("REMARK  42\n")
        header_records.append("REMARK  42 TOPOLOGY FILES USED\n")
        for top_file in self.topology_files:
            header_records.append("REMARK  42     {}\n".format(top_file))
        atom_records = [atom.to_pdb_line() for atom in self.atoms]
        conect_records = self.get_conect_records()
        master_fmt = "MASTER    " + "    0" * 8 + "{:>5}    0{:>5}    0\n"
        master = master_fmt.format(len(atom_records), len(conect_records))
        header_block = "".join(header_records)
        atom_block = "".join(atom_records)
        conect_block = "".join(conect_records)
        return header_block + atom_block + conect_block + master + "END\n"

    def to_pdb_file(self, pdb_path: str) -> None:
        """
        Write the molecule to a PDB file.

        :param pdb_path: path to PDB file.

        """
        with open(pdb_path, "w") as pdb_file:
            pdb_file.write(self.to_pdb_block())

    def to_psf_block(self) -> str:
        """
        Build the PSF records and return the PSF as a string.

        :return: PSF block as string

        """
        if not self._finalized:
            raise ValueError("Molecule not finalized")

        def format_psf_block(
            data: List[Tuple[int, ...]],
            n_per_line: int
        ) -> str:
            try:
                n_per_item = len(data[0])
            except (ValueError, IndexError, AttributeError):
                return ""
            format_string = "".join(["{:10d}" for _ in range(n_per_item)])
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

        doc_title = ["PSF EXT CMAP XPLOR", ""]
        header_records = [
            "* Generated procedurally by pygen-structures v{}".format(
                version
            ),
            "* Molecule: {}".format(self.name),
            "* Topology files used:",
        ]
        for top_file in self.topology_files:
            header_records.append("* - {}".format(top_file))
        doc_title.append("{:>10} !NTITLE".format(len(header_records)))
        header = "\n".join(doc_title + header_records) + "\n"

        adjacency_matrix = self._structure.adjacency_matrix
        bonds, angles, dihedrals = adjacency_to_dof(adjacency_matrix)

        atom_head = "{:10d} !NATOM".format(len(self.atoms))
        atom_block = "".join([atom.to_psf_line() for atom in self.atoms])

        bond_head = "{:10d} !NBOND: bonds".format(len(bonds))
        bond_block = format_psf_block(bonds, n_per_line=4)

        angle_head = "{:10d} !NTHETA: angles".format(len(angles))
        angle_block = format_psf_block(angles, n_per_line=3)

        dihedral_head = "{:10d} !NPHI: dihedrals".format(len(dihedrals))
        dihedral_block = "\n"
        if dihedrals:
            dihedral_block = format_psf_block(dihedrals, n_per_line=2)

        improper_head = "{:10d} !NIMPHI: impropers".format(len(self.impropers))
        improper_block = "\n"
        if self.impropers:
            improper_block = format_psf_block(self.impropers, n_per_line=2)

        ndon_block = "         0 !NDON: donors\n\n"
        nacc_block = "         0 !NACC: acceptors\n\n"

        nnb_head = "         0 !NNB\n"
        nnb_block = []
        for index in range(len(self.atoms)):
            if index == 0:
                continue
            if not index % 8:
                nnb_block.append("{:10d}".format(0) * 8)
        remaining_atoms = len(self.atoms) % 8
        if remaining_atoms:
            nnb_block.append("{:10d}".format(0) * remaining_atoms)
        nnb_block = "\n".join(nnb_block) + "\n"

        cmap_head = "{:10d} !NCRTERM: cross-terms".format(len(self.cross_maps))
        cmap_block = "\n"
        if self.cross_maps:
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
            ndon_block,
            nacc_block,
            nnb_head,
            nnb_block,
            cmap_head,
            cmap_block,
        ]
        return "\n".join(psf_blocks) + "\n"

    def to_psf_file(self, psf_path) -> None:
        """
        Write the molecule to a PSF file.

        :param psf_path: Path to PSF file.

        """
        with open(psf_path, "w") as psf_file:
            psf_file.write(self.to_psf_block())

    def check_parameters(self, parameter_set: CHARMMParameterFile) -> bool:
        """
        Check if parameters needed by the molecule are present in a
        parameter set. Returns True if all parameters are present,
        otherwise false.

        :param parameter_set: ``CHARMMParameterFile``
        :return: ``True`` if all parameters are matched

        """
        return parameter_set.check_parameters(self)

    @classmethod
    def from_sequence(
        cls,
        sequence: List[str],
        topology: CHARMMResidueTopologyFile,
        patches:  (Dict[str, List[int]], None) = None,
        name: (str, None) = None,
        segid: (str, None) = None
    ):
        """
        Create the Molecule instance from a list of residues, a
        CHARMM topology file and a list containing patches to be
        applied.

        :param sequence: a list of CHARMM residues
        :param topology: a ``CHARMMResidueTopologyFile``
        :param patches: a dict of ``CHARMMPatchResidueDefinition``\
        names from ``topology`` to the indices of the residue the \
        patch is to be applied to.
        :param name: the name of the ``Molecule``
        :param segid: the segment for the PDB/PDB

        """
        if name is None:
            name = []
            for charmm_name in sequence:
                try:
                    name.append(CHARMM_TO_ACID_CODE[charmm_name])
                except KeyError:
                    name.append("-{}-".format(charmm_name))
            name = ''.join(name)
        instance = cls(
            name, patches=patches, topology=topology, segment=segid
        )

        for residue_index, residue_name in enumerate(sequence):
            residue_definition = topology.residues[residue_name]
            residue = residue_definition.to_residue(residue_index)
            instance.residues.append(residue)
        instance.finalize()
        return instance
