"""
Module which stores classes intended to deal with CHARMM forcefield
data.
"""
# from __future__ import annotations
import re
import os
from typing import Tuple, List, Dict, Set, Union

from pygen_structures.mol_containers.atom import Atom
from pygen_structures.mol_containers.structure import Structure
from rdkit import Chem
from pygen_structures._functions_const import iter_nwise, adjacency_to_dof


class CHARMMResidueData:
    """
    Base class for CHARMM residue data. Holds information present in a
    residue topology file (RTF) section.

    """
    def __init__(
        self,
        name: str = "Unset",
        atoms: (list, None) = None,
        bonds: (list, None) = None,
        impropers: (list, None) = None,
        cross_maps: (list, None) = None,
        ics: (list, None) = None,
        rtf_file_name: (str, None) = None,
    ):
        self.name = name
        self.atoms = atoms if atoms else list()
        self.bonds = bonds if bonds else list()
        self.impropers = impropers if impropers else list()
        self.cross_maps = cross_maps if cross_maps else list()
        self.ics = ics if ics else list()
        self.rtf_file_name = rtf_file_name


class CHARMMResidueDefinition(CHARMMResidueData):
    """
    Definition of a CHARMM residue. This is a general representation,
    and doesn't have an associated index.

    :param name: residue name, from RESI record
    :param atoms: list of (atom_name, atom_type, partial_charge),\
    from ATOM records
    :param bonds: lists of tuples of 2 ``atom_name``, from BOND/DOUB\
    records. Bonds in DOUB records appear twice.
    :param impropers:  lists of tuples of 4 ``atom_name`` from IMPR\
    records
    :param cross_maps: lists of tuples of 8 ``atom_name`` from CMAP\
    records
    :param ics: list of information present in IC records. The first\
    4 items are ``atom_name``, the last 5 are floats containing bond\
    length/angle information. Where angle i-j-k-l is an improper,\
    the third ``atom_name`` is prefaced with "*".\
    The floats are:\
    the bond length, i-j\
    the angle, i-j-k\
    the dihedral (or improper) i-j-k-l\
    the angle j-k-l\
    the bond length k-l
    :param first: default patch if residue is first in the chain. If\
    ``None``, this defaults to the default first patch of the\
    ``CHARMMResidueTopologyFile`` this residue definition is in.
    :param last: default patch if residue is last in the chain. Treated\
    much the same way as ``first``.

    """

    def __init__(
        self,
        name: str = "Unset",
        atoms: (list, None) = None,
        bonds: (list, None) = None,
        impropers: (list, None) = None,
        cross_maps: (list, None) = None,
        ics: (list, None) = None,
        rtf_file_name: (str, None) = None,
        first: (str, None) = None,
        last: (str, None) = None,
    ):
        super().__init__(
            name,
            atoms,
            bonds,
            impropers,
            cross_maps,
            ics,
            rtf_file_name
        )
        self.first = first
        self.last = last

    @classmethod
    def from_block(cls, block: str):
        """
        Generate the class from a RESI block in an RTF file. This is
        probably the most useful constructor.

        :param block: a multiline string containing the RESI block.
        :return: an instance of the class

        """
        lines = block.split("\n")
        name = lines.pop(0).split()[1]

        atoms, bonds, impropers, cross_maps = [], [], [], []
        ics = []
        first, last = None, None

        for line in lines:
            if "!" in line:
                line = line[: line.index("!")]
            line = line.rstrip()
            if not line:
                continue

            line = line.split()
            record, data = line[0], line[1:]
            record = record.upper()
            if not data:
                continue

            if record == "ATOM":
                atom_name, atom_type, charge = data[:3]
                atoms.append((atom_name, atom_type, float(charge)))
            elif record == "BOND":
                bonds.extend(list(iter_nwise(data, 2)))
            elif record[:4] == "DOUB":
                bonds.extend(list(iter_nwise(data, 2)))
                bonds.extend(list(iter_nwise(data, 2)))
            elif record == "IMPR":
                impropers.extend(list(iter_nwise(data, 4)))
            elif record == "CMAP":
                cross_maps.extend(list(iter_nwise(data, 8)))
            elif record == "IC":
                data[4:] = [float(x) for x in data[4:]]
                ics.append(data)
            elif record[:4] == "PATC":
                patch_name = data[1]
                if data[0][:4] == "FIRS":
                    first = patch_name
                elif data[0] == "LAST":
                    last = patch_name
                try:
                    patch_name = data[3]
                except IndexError:
                    pass
                else:
                    if data[2][:4] == "FIRS":
                        first = patch_name
                    elif data[2] == "LAST":
                        last = patch_name
        return cls(
            name,
            atoms,
            bonds,
            impropers,
            cross_maps,
            ics,
            None,
            first,
            last
        )

    def to_residue(self, index: int):  # -> CHARMMResidue:
        """
        Generate a ``CHARMMResidue`` from the residue definition.

        :param index: the ``residue_index`` of the created residue.
        :return: the created residue
        :rtype: CHARMMResidue

        """
        return CHARMMResidue.from_residue_definition(self, index)

    def to_structure(self) -> Structure:
        """
        Returns the residue as a ``Structure`` object. These can be used
        for pattern matching of residues.

        :return: the generated ``Structure``

        """
        atoms, bonds, improper_ics = [], [], {}
        atom_name_to_index = {}
        for atom_index, atom_data in enumerate(self.atoms):
            atom_name, atom_type, charge = atom_data
            atom_name_to_index[atom_name] = atom_index
            atom = Atom(
                atom_serial=atom_index + 1,
                atom_name=atom_name,
                residue_name=self.name,
                residue_number=0,
                segment_id=self.name,
                atom_type=atom_type,
                partial_charge=charge,
            )
            atoms.append(atom)
        for (atom_name_x, atom_name_y) in self.bonds:
            try:
                x_idx = atom_name_to_index[atom_name_x]
                y_idx = atom_name_to_index[atom_name_y]
            except KeyError:
                continue
            bonds.append((x_idx, y_idx))
        for ic in self.ics:
            atom_names = ic[:4].copy()
            if "*" not in atom_names[2]:
                continue
            atom_names[2] = atom_names[2].replace("*", "")
            try:
                i, j, k, l = [atom_name_to_index[name] for name in atom_names]
            except KeyError:  # Contains atom in another residue
                continue
            improper = ic[6]
            improper_ics[(i, j, k, l)] = improper
        return Structure(atoms, bonds, improper_ics, set_charges=False)

    def to_fragment_mol(self) -> Chem.RWMol:
        """
        Generate a ``Structure and return the RDKit ``Mol`` generated.
        :return: RDKit ``RWMol``

        """
        structure = self.to_structure()
        return structure.mol

    def to_smarts(self):
        """
        Generate a ``Structure``, grab the RDKit ``Mol``, and call
        ``MolToSMARTS`` after removing hydrogen atoms.

        :return: the SMARTS string.

        """
        return Chem.MolToSmarts(Chem.RemoveHs(self.to_fragment_mol()))


class CHARMMResidue(CHARMMResidueData):
    """
    An actual CHARMM residue, in a molecule. This has an associated
    index. The atoms are unchanged from the residue definition
    (see ``CHARMMResidueDefinition``), but rather than referring to
    ``atom_name``, the ``bonds``, ``impropers``, ``cross_maps`` and ``ics`` now
    refer to ``atom_ids``, which are unique atoms rather than general
    abstract atoms definitions.

    An ``atom_id`` takes the following form:
     - ``(residue_index: int, atom_name: str)``

    :param first: default patch if residue is first in chain. If\
    None, this defaults to the default first patch of the\
    ``CHARMMResidueTopologyFile`` the residue definition is in.
    :param last: default patch if residue is last in the chain.\
    Treated much the same way as ``first``.
    :param index: the ``residue_index``.

    """

    def __init__(
        self,
        name: str = "Unset",
        atoms: (list, None) = None,
        bonds: (list, None) = None,
        impropers: (list, None) = None,
        cross_maps: (list, None) = None,
        ics: (list, None) = None,
        rtf_file_name: (str, None) = None,
        first: (str, None) = None,
        last: (str, None) = None,
        index: int = 0,
    ):
        super().__init__(
            name, atoms, bonds, impropers, cross_maps, ics, rtf_file_name
        )
        self.first = first
        self.last = last
        self.index = index
        self.atom_names = set(atom[0] for atom in atoms)

    @classmethod
    def from_residue_definition(
            cls,
            residue_definition: CHARMMResidueDefinition,
            index: int,
            last_index: (int, None) = None,
            next_index: (int, None) = None,
    ):
        """
        Generate a ``CHARMMResidue`` from the residue definition.\
        This is probably the most useful constructor.

        :param residue_definition: ``CHARMMResidueDefinition`` to clone.
        :param index: the ``residue_index`` of the created residue.
        :param last_index: the ``residue_index`` of the previous residue\
        in the chain. Used for connections where atoms start with '-'.\
        Default: index - 1
        :param last_index: the ``residue_index`` of the next residue in\
        the chain. Used for connections where atoms start with '+'.\
        Default: index - 1
        :return: the created residue.

        """
        if last_index is None:
            last_index = index - 1
        if next_index is None:
            next_index = index + 1

        def name_to_id(raw_atom_name) -> (int, str):
            first_char = raw_atom_name[0]
            re_add_star = False
            if first_char == "*":
                re_add_star = True
                first_char = raw_atom_name[1]
                raw_atom_name = raw_atom_name[1:]

            if first_char not in ('-', '+'):
                residue_index = index
            elif first_char == '-':
                residue_index = last_index
                raw_atom_name = raw_atom_name[1:]
            else:
                residue_index = next_index
                raw_atom_name = raw_atom_name[1:]

            if re_add_star:
                raw_atom_name = "*" + raw_atom_name
            return residue_index, raw_atom_name

        instance = cls(
            name=residue_definition.name,
            atoms=residue_definition.atoms.copy(),
            rtf_file_name=residue_definition.rtf_file_name,
            first=residue_definition.first,
            last=residue_definition.last,
            index=index,
        )
        for bond_definition in residue_definition.bonds:
            bond: List[Tuple[int, str], ...] = list()
            for atom_name in bond_definition:
                atom_id = name_to_id(atom_name)
                bond.append(atom_id)
            instance.bonds.append(tuple(bond))
        for improper_definition in residue_definition.impropers:
            improper: List[Tuple[int, str], ...] = list()
            for atom_name in improper_definition:
                atom_id = name_to_id(atom_name)
                improper.append(atom_id)
            instance.impropers.append(tuple(improper))
        for cross_map_definition in residue_definition.cross_maps:
            cross_map: List[Tuple[int, str], ...] = list()
            for atom_name in cross_map_definition:
                atom_id = name_to_id(atom_name)
                cross_map.append(atom_id)
            instance.cross_maps.append(tuple(cross_map))

        for ic_definition in residue_definition.ics:
            ic = ic_definition.copy()
            for position, atom_name in enumerate(ic[:4]):
                atom_id = name_to_id(atom_name)
                ic[position] = atom_id
            instance.ics.append(ic)
        instance.atom_names = [atom[0] for atom in instance.atoms]
        return instance


class CHARMMPatchResidueDefinition(CHARMMResidueData):
    """
    Definition of a CHARMM patch residue. As CHARMM patches
    can apply to multiple residues, these are more complicated
    to represent than the residues they are applied to.

    ``atoms``, ``bonds``, ``impropers``, ``cross_maps`` and ``ics`` are
    ``n_residues`` long lists of the atoms, bonds, impropers,
    cross_maps and ics for each residue. This enables these to be
    zipped together with the actual residues the patch is to apply them
    to.

    Atoms, due to their simplicity, are kept as the form they are
    stored in the residue:
    - (``atom_name``, ``atom_type``, ``charge``)

    In the other representations, atoms are stored as atom references.
    These are tuples of ``(residue_index: int, atom_name: str)``, where
    ``residue_index`` is a placeholder for the order given in the patch.
    - e.g. 2SG1 -> ``(1, 'SG1')``
    - "BOND 1SG1 2SG1" -> ``((0, 'SG1'), (1, 'SG1'))``

    This is to account for new bonds/impropers/crossmaps/ics which involve
    atoms in different residues.

    :param name: patch name, from PRES record
    :param atoms: an ``n_residues`` long list of lists containing\
    ``(atom_name, atom_type, partial_charge)``
    :param bonds: an ``n_residues`` long list of lists containing\
    tuples of 2 atom references, from BOND/DOUB records. Bonds\
    in DOUB records appear twice.
    :param impropers: an ``n_residues`` long list of lists containing\
    tuples of 4 atom references from IMPR records
    :param cross_maps:  an ``n_residues`` long list of lists containing\
    tuples of 8 atom references from CMAP records
    :param ics:  an ``n_residues`` long list of lists containing\
    lists of information present in IC records.\
    The first 4 items are atom references, the last 5 are floats\
    containing bond and angle information. Where angle i-j-k-l is an\
    improper, the third atom name is prefaced with "*".\
    \
    The floats are:\
    the bond length, i-j\
    the angle, i-j-k\
    the dihedral (or improper) i-j-k-l\
    the angle j-k-l\
    the bond length k-l
    :param deletions: set of atom names which are to be deleted from \
    the residue

    """
    def __init__(
        self,
        name: str = "Unset",
        atoms: (list, None) = None,
        bonds: (list, None) = None,
        impropers: (list, None) = None,
        cross_maps: (list, None) = None,
        ics: (list, None) = None,
        rtf_file_name: (str, None) = None,
        deletions: (set, None) = None,
        n_residues: int = 1,
    ):
        super().__init__(
            name, atoms, bonds, impropers, cross_maps, ics, rtf_file_name
        )
        self.deletions = deletions if deletions else set()
        self.n_residues = n_residues

    @classmethod
    def from_block(cls, block: str):
        """
        Generate the class from a PRES block in an RTF file. This is
        probably the most useful constructor.

        :param block: a multi-line string containing the PRES block.
        :return: An instance of the class

        """
        n_residues = 1

        def name_to_reference(raw_atom_name: str) -> Tuple[int, str]:
            """
            As patches can be applied to multiple residues, it's
            important to reflect this in the way we store data.

            This function converts an atom_name from a PRES to an atom
            reference.

            e.g. "2SG1" -> ``(1, 'SG1')``

            :param raw_atom_name: atom name
            :return: (residue_index, atom_name)
            """
            nonlocal n_residues

            re_add_star = False
            if raw_atom_name[0] == "*":
                re_add_star = True
                raw_atom_name = raw_atom_name[1:]

            try:
                residue_serial = int(raw_atom_name[0])
                raw_atom_name = raw_atom_name[1:]
            except ValueError:
                residue_serial = 1
            patch_residue_index = residue_serial - 1
            n_residues = max(patch_residue_index + 1, n_residues)
            if re_add_star:
                raw_atom_name = "*" + raw_atom_name
            return patch_residue_index, raw_atom_name

        lines = block.split("\n")
        name = lines.pop(0).split()[1]

        atoms, bonds, impropers, cross_maps = [], [], [], []
        ics, deletions = [], set()

        for line in lines:
            if "!" in line:
                line = line[: line.index("!")]
            line = line.rstrip()
            if not line:
                continue

            line = line.split()
            record, data = line[0], line[1:]
            record = record.upper()
            if not data:
                continue

            if record == "ATOM":
                atom_name, atom_type, charge = data[:3]
                atom_reference = name_to_reference(atom_name)
                atom = (atom_reference, atom_type, float(charge))
                atoms.append(atom)
            elif record == "IC":
                data[:4] = [name_to_reference(x) for x in data[:4]]
                data[4:] = [float(x) for x in data[4:]]
                ics.append(data)
            else:
                if record[:4] == "DELE":
                    atom_name = data[1]
                    deletions.add(name_to_reference(atom_name))
                elif record == "BOND":
                    data = [name_to_reference(x) for x in data]
                    bonds.extend(list(iter_nwise(data, 2)))
                elif record[:4] == "DOUB":
                    data = [name_to_reference(x) for x in data]
                    bonds.extend(list(iter_nwise(data, 2)))
                    bonds.extend(list(iter_nwise(data, 2)))
                elif record == "IMPR":
                    data = [name_to_reference(x) for x in data]
                    impropers.extend(list(iter_nwise(data, 4)))
                elif record == "CMAP":
                    data = [name_to_reference(x) for x in data]
                    cross_maps.extend(list(iter_nwise(data, 8)))

        split_atoms = [[] for _ in range(n_residues)]
        split_bonds = [[] for _ in range(n_residues)]
        split_impropers = [[] for _ in range(n_residues)]
        split_cross_maps = [[] for _ in range(n_residues)]
        split_ics = [[] for _ in range(n_residues)]
        for atom in atoms:
            (placeholder_index, atom_name), atom_type, charge = atom
            atom = (atom_name, atom_type, charge)
            split_atoms[placeholder_index].append(atom)
        for bond in bonds:
            (placeholder_index, atom_name) = bond[0]
            split_bonds[placeholder_index].append(bond)
        for improper in impropers:
            (placeholder_index, atom_name) = improper[1]
            split_impropers[placeholder_index].append(improper)
        for cross_map in cross_maps:
            (placeholder_index, atom_name) = cross_map[1]
            split_cross_maps[placeholder_index].append(cross_map)
        for ic in ics:
            (placeholder_index, atom_name) = ic[1]
            split_ics[placeholder_index].append(ic)

        instance = cls(
            name,
            split_atoms,
            split_bonds,
            split_impropers,
            split_cross_maps,
            split_ics,
            None,
            deletions,
            n_residues,
        )
        return instance

    def apply(self, *residues):  # *residues: CHARMMResidue):
        """
        Applies the CHARMM patch to an actual residue (or collection of
        actual residues).

        Atoms to be deleted are removed from the applicable residues,
        with now-dangling bonds/impropers/cross-maps/ICs to these
        residues also removed.

        New atoms, bonds, impropers, cross_maps and ics are then
        appended to the relevant residue.

        Bonds are added to the residue which contains the first atom,
        all other multi-atom data is added to the residue which
        contains the second atom. This is to account for the fact that
        the other connections can feature atoms in the residue before
        and after in the first and last positions - in bonds, these
        tend do be in the first position.
        """

        if len(residues) != self.n_residues:
            n_residues = self.n_residues
            raise ValueError("Patch applies to {} residues".format(n_residues))
        for residue in residues:
            if not isinstance(residue, CHARMMResidue):
                raise ValueError("Patch must be applied to CHARMMResidue")

        # Replacing placeholder tuples with formatted atom names.
        replacements = [residue.index for residue in residues]

        # Handling items to be deleted first. Deletions is a set
        # containing atom_ids, split_deletions is a list
        # containing raw atom_names for each residue.
        deletions = set()
        split_deletions = [set() for _ in range(self.n_residues)]
        for (placeholder_index, atom_name) in self.deletions:
            residue_index = replacements[placeholder_index]
            deletions.add((residue_index, atom_name))
            split_deletions[placeholder_index].add(atom_name)

        for deletion_names, residue in zip(split_deletions, residues):
            residue_atom_names = set(a[0] for a in residue.atoms)
            if deletion_names - residue_atom_names:
                err = "Patch {} not applicable to residue {}"
                raise ValueError(err.format(self.name, residue.name))

        for deletion_names, residue in zip(split_deletions, residues):
            atoms = [a for a in residue.atoms if a[0] not in deletion_names]

            bonds, impropers, cross_maps, ics = [], [], [], []
            for bond in residue.bonds:
                for atom_id in bond:
                    if atom_id in deletions:
                        break
                else:
                    bonds.append(bond)
            for improper in residue.impropers:
                for atom_id in improper:
                    if atom_id in deletions:
                        break
                else:
                    impropers.append(improper)
            for cross_map in residue.cross_maps:
                for atom_id in cross_map:
                    if atom_id in deletions:
                        break
                else:
                    cross_maps.append(cross_map)
            for ic in residue.ics:
                for residue_index, atom_name in ic[:4]:
                    if atom_name.replace("*", "") in deletion_names:
                        break
                else:
                    ics.append(ic)
            residue.atoms = atoms
            residue.bonds = bonds
            residue.impropers = impropers
            residue.cross_maps = cross_maps
            residue.ics = ics

        zipped = zip(
            self.atoms, self.bonds, self.impropers, self.cross_maps, self.ics, residues
        )
        for atoms, bonds, impropers, cross_maps, ics, residue in zipped:
            formatted_bonds = []
            for bond in bonds:
                bond = [(replacements[x], a_name) for (x, a_name) in bond]
                formatted_bonds.append(tuple(bond))
            bonds = formatted_bonds
            formatted_impropers = []
            for improper in impropers:
                improper = [(replacements[x], a_name) for (x, a_name) in improper]
                formatted_impropers.append(tuple(improper))
            impropers = formatted_impropers
            formatted_cross_maps = []
            for cross_map in cross_maps:
                cross_map = [(replacements[x], a_name) for (x, a_name) in cross_map]
                formatted_cross_maps.append(tuple(cross_map))
            cross_maps = formatted_cross_maps

            for atom in atoms:
                atom_name = atom[0]
                for res_atom_index, res_atom in enumerate(residue.atoms):
                    res_atom_name = res_atom[0]
                    if res_atom_name == atom_name:
                        residue.atoms[res_atom_index] = atom
                        break
                else:
                    residue.atoms.append(atom)

            # Because the number of bonds determines the bond order, we
            # need to remove overlapping bonds and re-add them. This is
            # to ensure changes in bond order are captured
            bond_intersection = set(bonds).intersection(residue.bonds)
            for bond in bond_intersection:
                while bond in residue.bonds:
                    residue.bonds.remove(bond)
            for bond in bonds:
                residue.bonds.append(bond)

            # This should not be important for impropers and cross-maps
            # If existing impropers and cross maps have been specified,
            # this is the fault of the RTF file. In any case, sane
            # simulation packages should ignore multiple specifications
            for improper in impropers:
                if improper not in residue.impropers:
                    residue.impropers.append(improper)
            for cross_map in cross_maps:
                if cross_map not in residue.cross_maps:
                    residue.cross_maps.append(cross_map)
            # Similarly to bonds, IC angles may have changed for the
            # same atoms. Working out the intersections is slightly
            # complicated by the extra information present.
            for ic in ics:
                ic = ic.copy()
                for idx, (placeholder_index, atom_name) in enumerate(ic[:4]):
                    residue_index = replacements[placeholder_index]
                    atom_id = (residue_index, atom_name)
                    ic[idx] = atom_id
                while True:
                    for idx, residue_ic in enumerate(residue.ics):
                        if ic[:4] == residue_ic[:4]:
                            deletion_index = idx
                            break
                    else:
                        break
                    del residue.ics[deletion_index]
                residue.ics.append(ic)

    def is_applicable_to(self, residue) -> Union[None, List[int]]:
        """
        Work out the positions where the patch can be applied to the
        residue. If the patch can be applied, return a list of
        positions where the patch can be applied, otherwise return
        ``None``.

        :param residue: The ``CHARMMResidue`` the patch is to be\
        applied to
        :return: A list of indices where the residue could be supplied\
        to the patch, if applicable. Otherwise ``None``.

        """
        deletions = set()
        split_deletions = [set() for _ in range(self.n_residues)]
        for (placeholder_index, atom_name) in self.deletions:
            deletions.add((residue.index, atom_name))
            split_deletions[placeholder_index].add(atom_name)

        applicable_positions = []
        for position, deletion_names in enumerate(split_deletions):
            residue_atom_names = set([a[0] for a in residue.atoms])
            if deletion_names - residue_atom_names:
                continue
            else:
                applicable_positions.append(position)
        if applicable_positions:
            return applicable_positions
        else:
            return


class CHARMMResidueTopologyFile:
    """
    A class which reads and stores patch residues and residues from a
    CHARMM RTF or STR file.

    :param file_name: the name of the RTF file.
    :param residues: a mapping of residue names to the\
    ``CHARMMResidueDefinition`` which represents them.
    :param patches: a mapping of patch residue names to the\
    ``CHARMMPatchResidueDefinition`` which represents them.
    :param masses: a mapping of ``atom_type`` to mass.
    :param first: the name of the patch applied to the first residue\
    in a chain.
    :param last: the name of the patch applied to the last residue\
    in a chain.

    """
    def __init__(
        self,
        file_name: (str, None) = None,
        residues: (dict, None) = None,
        patches: (dict, None) = None,
        masses: (dict, None) = None,
        first: (str, None) = None,
        last: (str, None) = None,
    ):
        self.file_name = file_name
        self.residues = residues if residues else {}
        self.patches = patches if patches else {}
        self.masses = masses if masses else {}
        self.first = first
        self.last = last

    def add_residue_definition(
            self,
            residue: CHARMMResidueDefinition
        ) -> None:
        """
        Add a ``CHARMMResidueDefinition`` to the residues

        :param residue: a ``CHARMMResidueDefinition``

        """
        if not residue.first:
            residue.first = self.first
        if not residue.last:
            residue.last = self.last
        name = residue.name
        self.residues[name] = residue

    def add_patch_definition(
            self,
            patch: CHARMMPatchResidueDefinition
        ) -> None:
        """
        Add a ``CHARMMPatchResidueDefinition`` to the patches

        :param patch: a ``CHARMMPatchResidueDefinition``

        """
        name = patch.name
        self.patches[name] = patch

    def read_file(
            self,
            rtf_path: str,
            update_default_patches=False
        ) -> None:
        """
        Read in a CHARMM RTF (or STR) file and add the patches and
        residues, updating the default patches if requested.

        :param rtf_path: path to the .rtf/.str file
        :param update_default_patches: bool flag, updates ``first`` \
        and ``last`` if True.

        """
        with open(rtf_path) as rtf_file:
            contents = rtf_file.read()
        if os.path.splitext(rtf_path)[-1] == ".str":
            try:
                contents = re.findall("(?smi)(^READ RTF.*?^END)", contents)[0]
            except IndexError:
                return
        rtf_file_name = os.path.split(rtf_path)[1]
        # There is a semantic difference between 'NONE' and None.
        # Whereas None suggests the information is missing and should
        # be filled by the topology file, 'NONE' states explicitly that no
        # patches are to be applied.
        first, last = "NONE", "NONE"
        default_regex = "(?im)^DEFAU?L?T? +(FIRST? +[A-Z0-9]+)? ?(LAST +[A-Z0-9]+)"
        matches = re.findall(default_regex, contents)
        if matches:
            for match in matches[0]:
                position, patch_name = match.split()
                if position[:4].upper() == "FIRS":
                    first = patch_name
                if position.upper() == "LAST":
                    last = patch_name
            if update_default_patches:
                self.first = first
                self.last = last

        mass_matches = re.findall(r"(?m)^MASS\s+-1\s+([A-Z0-9]+)\s+([0-9.]+)", contents)
        for atom_type, mass in mass_matches:
            self.masses[atom_type] = float(mass)

        residue_regex = "(?ms)(^RESI.*?)\n^(?=(?:PRES|END|RESI))"
        residue_blocks = re.findall(residue_regex, contents)
        for residue_block in residue_blocks:
            residue = CHARMMResidueDefinition.from_block(residue_block)
            residue.rtf_file_name = rtf_file_name
            if not residue.first:
                residue.first = first
            if not residue.last:
                residue.last = last
            self.add_residue_definition(residue)

        patch_regex = "(?ms)(^PRES.*?)\n^(?=(?:PRES|END|RESI))"
        patch_blocks = re.findall(patch_regex, contents)
        for patch_block in patch_blocks:
            patch = CHARMMPatchResidueDefinition.from_block(patch_block)
            patch.rtf_file_name = rtf_file_name
            self.add_patch_definition(patch)

    @classmethod
    def from_file(cls, rtf_path: str):
        """
        Instantiate the class from a CHARMM RTF file.

        :param rtf_path: path to .rtf/.str file

        """
        rtf_file_name = os.path.split(rtf_path)[1]
        instance = cls(file_name=rtf_file_name)
        instance.read_file(rtf_path, update_default_patches=True)
        return instance

    def __add__(self, other):
        if not isinstance(other, self.__class__):
            raise NotImplementedError
        residues = self.residues.copy()
        residues.update(other.residues)
        patches = self.patches.copy()
        patches.update(other.patches)
        masses = self.masses.copy()
        masses.update(other.masses)
        first, last = self.first, self.last
        if self.first != other.first:
            first = "NONE"
        if self.last != other.last:
            last = "NONE"
        return self.__class__(residues, patches, masses, first, last)


class CHARMMParameterFile:
    """
    A parser for CHARMM PRM files. Used for checking that parameters
    exist for generated structures.

    This class does not contain the values for the parameters, but
    does ensure that the parameters themselves are present.

    :param bonds: set of tuples containing 2 ``atom_type``
    :param angles: set of tuples containing 3 ``atom_type``
    :param dihedrals: set of tuples containing 4 ``atom_type``
    :param impropers: set of tuples containing 4 ``atom_type``
    :param cross_maps: set of tuples containing 8 ``atom_type``
    """
    def __init__(
        self,
        bonds: (set, None) = None,
        angles: (set, None) = None,
        dihedrals: (set, None) = None,
        impropers: (set, None) = None,
        cross_maps: (set, None) = None,
    ):
        self.bonds = bonds if bonds else set()
        self.angles = angles if angles else set()
        self.dihedrals = dihedrals if dihedrals else set()
        self.impropers = impropers if impropers else set()
        self.cross_maps = cross_maps if cross_maps else set()

    def read_file(self, prm_path: str) -> None:
        """
        Read in data from a CHARMM PRM file.

        :param prm_path: path to .prm/.str file

        """
        with open(prm_path, 'r') as prm_file:
            contents = prm_file.read()
        if os.path.splitext(prm_path)[-1] == ".str":
            try:
                contents = re.findall("(?smi)(^READ PARA.*?^END)", contents)[0]
            except IndexError:
                return
        contents = contents.split('\n')
        headings = ("BOND", "ANGL", "DIHE", "IMPR", "CMAP", "NONB", "CUTN")
        section_heading = None

        for line in contents:
            line = line.split('!')[0].rstrip().upper()
            if not line:
                continue
            if line[:4] in headings:
                section_heading = line[:4]
                continue
            elif not section_heading:
                continue
            line = line.split()
            # Easier to wrap whole thing in try/except than wrap each
            # statement. If a line doesn't match one of the records,
            # we're not interested in it.
            try:
                # All of these degrees can be present in backwards form
                # bar CMAPs, which should be ordered.
                if section_heading == "BOND":
                    x, y = line[:2]
                    self.bonds.add((x, y))
                    self.bonds.add((y, x))
                elif section_heading == "ANGL":
                    x, y, z = line[:3]
                    self.angles.add((x, y, z))
                    self.angles.add((z, y, x))
                elif section_heading == "DIHE":
                    w, x, y, z = line[:4]
                    self.dihedrals.add((w, x, y, z))
                    self.dihedrals.add((z, y, x, w))
                elif section_heading == "IMPR":
                    w, x, y, z = line[:4]
                    self.impropers.add((w, x, y, z))
                    self.impropers.add((z, y, x, w))
                elif section_heading == "CMAP":
                    cmap = line[:8]
                    if not len(cmap) == 8:
                        continue
                    self.cross_maps.add(tuple(cmap))
            except ValueError:
                continue

    @classmethod
    def from_file(cls, prm_path: str):
        """
        Instantiate the class from a CHARMM PRM file.

        :param prm_path: path to the .prm/.str file.

        """
        instance = cls()
        instance.read_file(prm_path)
        return instance

    def get_unmatched(self, molecule) -> Dict[str, Set[Tuple[str]]]:
        """
        Check a finalized ``Molecule`` to verify that all the parameters
        in the ``CHARMMParameterFile`` exist for that molecule, returning
        necessary parameters which are missing.

        :param molecule: the ``Molecule``
        :return: dict of parameter type (i.e. 'bonds'/'angles') to a \
        set of tuples containing unmatched ``atom_type`` tuples.

        """
        atoms = molecule.atoms
        adjacency_matrix = molecule.get_adjacency_matrix()
        bonds, angles, dihedrals = adjacency_to_dof(adjacency_matrix)
        impropers = molecule.impropers
        cross_maps = molecule.cross_maps
        bonds = [tuple(atoms[x].atom_type for x in b) for b in bonds]
        angles = [tuple(atoms[x].atom_type for x in a) for a in angles]
        dihedrals = [tuple(atoms[x].atom_type for x in d) for d in dihedrals]
        impropers = [tuple(atoms[x].atom_type for x in i) for i in impropers]
        cross_maps = [tuple(atoms[x].atom_type for x in c) for c in cross_maps]

        unmatched_bonds = set(bonds) - self.bonds
        unmatched_angles = set(angles) - self.angles
        unmatched_cross_maps = set(cross_maps) - self.cross_maps

        dihedrals = set(dihedrals) - self.dihedrals
        impropers = set(impropers) - self.impropers

        unmatched_dihedrals, unmatched_impropers = set(), set()
        # Dihedrals and impropers can contain wildcards. These aren't
        # permitted for the other parameters.
        for dihedral in dihedrals:
            w, x, y, z = dihedral
            if ("X", x, y, "X") in self.dihedrals:
                break
            else:
                unmatched_dihedrals.add(dihedral)
        for improper in impropers:
            w, x, y, z = improper
            potential_impropers = (
                (w, "X", "X", z),
                ("X", x, y, z),
                (w, x, y, "X"),
                ("X", y, z, "X"),
                (x, y, "X", "X"),
                ("X", "X", z, w),
            )
            for potential_improper in potential_impropers:
                if potential_improper in self.impropers:
                    break
            else:
                unmatched_impropers.add(improper)
        unmatched = {
            'bonds': unmatched_bonds,
            'angles': unmatched_angles,
            'dihedrals': unmatched_dihedrals,
            'impropers': unmatched_impropers,
            'cross_maps': unmatched_cross_maps,
        }
        return unmatched

    def check_parameters(self, molecule) -> bool:
        """
        Check a finalized ``Molecule`` to verify that all the parameters
        in the ``CHARMMParameterFile`` exist for that molecule, returning
        ``True`` if all parameters are present, otherwise ``False``
        """
        for record_type, unmatched in self.get_unmatched(molecule).items():
            if unmatched:
                break
        else:
            return True
        return False

    def __add__(self, other):
        if not isinstance(other, self.__class__):
            raise NotImplementedError("Must be CHARMMParameterFile")
        bonds = self.bonds | other.bonds
        angles = self.angles | other.angles
        dihedrals = self.dihedrals | other.dihedrals
        impropers = self.impropers | other.impropers
        cross_maps = self.cross_maps.copy() | other.cross_maps
        return self.__class__(
            bonds,
            angles,
            dihedrals,
            impropers,
            cross_maps,
        )
