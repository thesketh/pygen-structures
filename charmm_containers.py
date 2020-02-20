"""
Module which stores classes intended to deal with CHARMM forcefield
data.
"""

from __future__ import annotations
import re
import os
from typing import Tuple, List, Dict, Set
from _functions_const import iter_nwise


def load_charmm_directory(
    directory_root: str
) -> (CHARMMResidueTopologyFile, CHARMMParameterFile):
    """
    A convenience function which scans a CHARMM toppar folder and
    reads in the rtf, prm and str files to a CHARMMResidueTopologyFile
    and a CHARMMParameterFile.
    """
    rtf, prm = CHARMMResidueTopologyFile(), CHARMMParameterFile()
    for root, dirs, file_names in os.walk(directory_root, topdown=True):
        dirs[:] = [d for d in dirs if d != "non_charmm"]
        for file_name in file_names:
            extension = os.path.splitext(file_name)[1]
            file_name = os.path.join(root, file_name)
            if extension == ".rtf":
                rtf.read_file(file_name)
            elif extension == ".prm":
                prm.read_file(file_name)
            elif extension == ".str":
                rtf.read_file(file_name)
                prm.read_file(file_name)
            else:
                continue
    return rtf, prm


class CHARMMResidueData:
    """
    Base class for CHARMM residue data. Holds information present in a
    residue topology file (RTF):
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
    and doesn"t have an associated index.

    Contains:
    name:       residue name, from RESI record
    atoms:      list of (atom_name, atom_type, partial_charge),
                from ATOM records
    bonds:      lists of tuples of 2 `atom_name`s, from BOND/DOUB
                records. Bonds in DOUB records appear twice
    impropers:  lists of tuples of 4 `atom_name`s from IMPR records
    cross_maps: lists of tuples of 8 `atom_name`s from CMAP records
    ics:        tuple of information present in IC records. The
                first 4 items are `atom_name`s, the last 5 are
                floats containing bond length/angle information.
                Where angle i-j-k-l is an improper, the third
                `atom_name` is prefaced with '*'.
                The floats are:
                    the bond length, i-j
                    the angle, i-j-k
                    the dihedral (or improper) i-j-k-l
                    the angle j-k-l
                    the bond length k-l
    first:      default patch if residue is first in chain. If
                None, this defaults to the default first patch
                of the CHARMMResidueTopologyFile this residue
                definition is in.
    last:       default patch if residue is last in the chain.
                Treated much the same way as `first`.
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
        :return:      an instance of the class
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

    def to_residue(self, index: int) -> CHARMMResidue:
        """
        Generate a CHARMMResidue from the residue definition.
        :param index: the index of the created residue.
        :return:      the created residue
        """
        return CHARMMResidue.from_residue_definition(self, index)


class CHARMMResidue(CHARMMResidueData):
    """
    An actual CHARMM residue, in a molecule. This has an associated
    index. The atoms are unchanged from the residue definition
    (see CHARMMResidueDefinition), but rather than referring to
    `atom_name`s, the bonds, impropers, cross_maps and ics now refer
    to atom_ids, which are unique atoms rather than general abstract
    atoms definitions.

    An atom_id takes the following form:
        `residue_index:atom_name`

    first:      default patch if residue is first in chain. If
                None, this defaults to the default first patch
                of the CHARMMResidueTopologyFile this residue
                definition is in.
    last:       default patch if residue is last in the chain.
                Treated much the same way as `first`.
    index:      the residue's index.
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

    @classmethod
    def from_residue_definition(
            cls,
            residue_definition: CHARMMResidueDefinition,
            index: int
        ):
        """
        Generate a CHARMMResidue from the residue definition.
        This is probably the most useful constructor.
        :param index:              the index of the created residue.
        :param residue_definition: CHARMMResidueDefinition to clone.
        :return:                   the created residue.
        """

        def format_atom_id(atom_name: str):
            """Converts an atom name to an atom id """
            if atom_name.startswith("+"):
                res_idx = index + 1
                atom_name = atom_name[1:]
            elif atom_name.startswith("-"):
                res_idx = index - 1
                atom_name = atom_name[1:]
            else:
                res_idx = index
            atom_id = "{}:{}".format(res_idx, atom_name)
            return atom_id

        instance = cls(
            name=residue_definition.name,
            atoms=residue_definition.atoms.copy(),
            rtf_file_name=residue_definition.rtf_file_name,
            first=residue_definition.first,
            last=residue_definition.last,
            index=index,
        )
        for bond in residue_definition.bonds:
            instance.bonds.append(tuple([format_atom_id(a) for a in bond]))
        for improper in residue_definition.impropers:
            instance.impropers.append(tuple([format_atom_id(a) for a in improper]))
        for cross_map in residue_definition.cross_maps:
            instance.cross_maps.append(tuple([format_atom_id(a) for a in cross_map]))
        for ic in residue_definition.ics:
            ic = [format_atom_id(a) for a in ic[:4]] + ic[4:]
            instance.ics.append(ic)
        return instance


class CHARMMPatchResidueDefinition(CHARMMResidueData):
    """
    Definition of a CHARMM patch residue. As CHARMM patches
    can apply to multiple residues, these are more complicated
    to represent than the residues they are applied to.

    `atoms`, `bonds`, `impropers`, `cross_maps` and `ics` are
    `n_residues` long lists of the atoms, bonds, impropers, cross_maps
    and ics for each residue. This enables these to be zipped together
    with the actual residues the patch is to apply them to.

    Atoms, due to their simplicity, are kept as the form they are
    stored in the residue:
        (atom_name, atom_type, charge)
    In the other representations, atoms are referenced by tuples of
    (residue_index, atom_name), where residue name is the
    order given in the patch.
        e.g. 2SG1 -> (1, 'SG1')
             BOND 1SG1 2SG1 -> ((0, 'SG1'), (1, 'SG1'))
    This is to account for new bonds/impropers/crossmaps/ics
    which involve atoms in different residues.

    Contains:
    name: patch name, from PRES record
    lists for each residue:
        atoms:      list of (atom_name, atom_type, partial_charge)
        bonds:      list of tuples of 2 atom references, from BOND/DOUB
                    records. Bonds in DOUB records appear twice
        impropers:  list of tuples of 4 atom references from IMPR
                    records
        cross_maps: list of tuples of 8 atom references from CMAP
                    records
        ics:        list of tuples of information present in IC
                    records. The first 4 items are atom references, the
                    last 5 are floats containing bond and angle
                    information. Where angle i-j-k-l is an improper,
                    the third atom name is prefaced with *.
                    The floats are:
                        the bond length, i-j
                        the angle, i-j-k
                        the dihedral (or improper) i-j-k-l
                        the angle j-k-l
                        the bond length k-l
        deletions:  set of atom names which are to be deleted from the
                    residue
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
        :param block: a multiline string containing the PRES block.
        :return: An instance of the class
        """
        n_residues = 1

        def name_to_reference(atom_name: str) -> Tuple[int, str]:
            """
            As patches can be applied to multiple residues, it's
            important to reflect this in the way we store data.

            This function converts an atom_name from a PRES to an atom
            reference.

            e.g. 2SG1 -> (1, 'SG1')

            :param a_name: atom name
            :return: (residue_index, atom_name)
            """
            nonlocal n_residues
            re_add_star = False
            if atom_name[0] == "*":
                re_add_star = True
                atom_name = atom_name[1:]

            try:
                residue_serial = int(atom_name[0])
                atom_name = atom_name[1:]
            except ValueError:
                residue_serial = 1
            residue_index = residue_serial - 1
            n_residues = max(n_residues, residue_serial)
            if re_add_star:
                atom_name = "*" + atom_name
            return residue_index, atom_name

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
                residue_index, atom_name = name_to_reference(atom_name)
                atom = ((residue_index, atom_name), atom_type, float(charge))
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
            (residue_index, atom_name), atom_type, charge = atom
            atom = (atom_name, atom_type, charge)
            split_atoms[residue_index].append(atom)
        for bond in bonds:
            residue_index = bond[0][0]
            split_bonds[residue_index].append(bond)
        for improper in impropers:
            residue_index = improper[1][0]
            split_impropers[residue_index].append(improper)
        for cross_map in cross_maps:
            residue_index = cross_map[1][0]
            split_cross_maps[residue_index].append(cross_map)
        for ic in ics:
            residue_index = ic[1][0]
            split_ics[residue_index].append(ic)

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

    def apply(self, *residues: CHARMMResidue):
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
        the other connections feature can atoms in the residue before
        and after in the first and last positions - in bonds, these
        are only in the first position.
        """

        if len(residues) != self.n_residues:
            n_residues = self.n_residues
            raise ValueError("Patch applies to {} residues".format(n_residues))
        for residue in residues:
            if not isinstance(residue, CHARMMResidue):
                raise ValueError("Patch must be applied to CHARMMResidue")

        # Replacing placeholder tuples with formatted atom names.
        replacements = [residue.index for residue in residues]

        def format_ids(
            iterable: List[Tuple[Tuple[int, str]]]
        ) -> List[Tuple[str]]:
            """
            Formatting lists of tuples of atom references to replace
            them with atom ids of the actual atoms they refer to.
            """
            formatted_iterable = []
            for item in iterable:
                formatted_item = []
                for (placeholder_index, atom_name) in item:
                    residue_index = replacements[placeholder_index]
                    if atom_name.startswith("-"):
                        residue_index -= 1
                        atom_name = atom_name[1:]
                    elif atom_name.startswith("+"):
                        residue_index += 1
                        atom_name = atom_name[1:]
                    atom_id = "{}:{}".format(residue_index, atom_name)
                    formatted_item.append(atom_id)
                formatted_iterable.append(tuple(formatted_item))
            return formatted_iterable

        # Handling items to be deleted first. Deletions is a set
        # containing atom_ids, split_deletions is a list
        # containing raw atom_names for each residue.
        deletions = set()
        split_deletions = [set() for _ in range(self.n_residues)]
        for (placeholder_index, atom_name) in self.deletions:
            residue_index = replacements[placeholder_index]
            deletions.add("{}:{}".format(residue_index, atom_name))
            split_deletions[placeholder_index].add(atom_name)

        for deletion_names, residue in zip(split_deletions, residues):
            residue_atom_names = set([a[0] for a in residue.atoms])
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
                for atom_id in ic[:4]:
                    residue_index, atom_name = atom_id.split(":")
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
            bonds = format_ids(bonds)
            impropers = format_ids(impropers)
            cross_maps = format_ids(cross_maps)
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
            # to ensure that if the bond order has changed, we capture
            # that information.
            bond_intersection = set(bonds).intersection(residue.bonds)
            for bond in bond_intersection:
                while bond in residue.bonds:
                    residue.bonds.remove(bond)
            for bond in bonds:
                residue.bonds.append(bond)

            # This should not be important for impropers and cross-maps
            # If existing impropers and cross maps have been specified,
            # this is the fault of the RTF file. In any case, most sane
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
                    atom_id = "{}:{}".format(residue_index, atom_name)
                    ic[idx] = atom_id
                while True:
                    deletion_index = None
                    for idx, residue_ic in enumerate(residue.ics):
                        if ic[:4] == residue_ic[:4]:
                            deletion_index = idx
                            break
                    if deletion_index:
                        del residue.ics[deletion_index]
                    else:
                        break
                residue.ics.append(ic)


class CHARMMResidueTopologyFile:
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

    def add_residue_definition(self, residue):
        if not isinstance(residue, CHARMMResidueDefinition):
            err = "{} not a CHARMMResidueDefinition"
            raise ValueError(err.format(residue.__class__))
        name = residue.name
        self.residues[name] = residue

    def add_patch_definition(self, patch):
        if not isinstance(patch, CHARMMPatchResidueDefinition):
            err = "{!r} not a CHARMMPatchResidueDefinition"
            raise ValueError(err.format(patch))
        name = patch.name
        self.patches[name] = patch

    def read_file(self, rtf_path, update_default_patches=False):
        with open(rtf_path) as rtf_file:
            contents = rtf_file.read()
        if os.path.split(rtf_path)[-1] == ".str":
            try:
                contents = re.findall("(?si)(.+)READ PARAM.*")[0]
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
    def from_file(cls, rtf_path):
        """
        Instantiate the class from a CHARMM RTF file.
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

    Contains:
    bonds:      set of tuples containing 2 `atom_type`s
    angles:     set of tuples containing 3 `atom_type`s
    dihedrals:  set of tuples containing 4 `atom_type`s
    impropers:  set of tuples containing 4 `atom_type`s
    cross_maps: set of tuples containing 8 `atom_type`s
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
        Read in data from a CHARMM PRM file
        """
        with open(prm_path, 'r') as prm_file:
            contents = prm_file.read()
        if os.path.splitext(prm_path)[-1] == ".str":
            try:
                contents = re.findall("(?si)READ PARAM.*?\n(.+)", contents)[0]
            except IndexError:
                return
        contents = contents.split('\n')
        headings = ("BOND", "ANGL", "DIHE", "IMPR", "CMAP", "NONB", "CUTN")
        section_heading = None
        bonds, angles, dihedrals,  = set(), set(), set()
        impropers, cross_maps = set(), set()

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
    def from_file(cls, prm_path: str) -> CHARMMParameterFile:
        instance = cls()
        instance.read_file(prm_path)
        return instance

    def get_unmatched(self, molecule: Molecule) -> Dict[str, Set[Tuple[str]]]:
        """
        Check a finalized Molecule to verify that all the parameters
        in the CHARMMParameterFile exist for that molecule, returning
        parameters in the Molecule which are unmatched.
        """
        if not isinstance(molecule, Molecule):
            mol_class = molecule.__class__.__name__
            err = "Expected Molecule instance, got {}".format(mol_class)
            raise ValueError(err)
        if not molecule._finalized:
            raise ValueError('Molecule has not been finalized.')
        atoms = molecule.atoms
        bonds, angles, dihedrals = adjacency_to_dof(molecule.adjacency_matrix)
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

    def check_parameters(self, molecule: Molecule) -> bool:
        """
        Check a finalized Molecule to verify that all the parameters
        in the CHARMMParameterFile exist for that molecule, returning
        True if all parameters are present, otherwise False.
        """
        for record_type, unmatched in self.get_unmatched(molecule).items():
            if unmatched:
                break
        else:
            return True
        return False

    def __add__(self, other: CHARMMParameterFile) -> CHARMMParameterFile:
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
