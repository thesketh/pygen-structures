"""
This module contains convenience functions which are likely to be the
main way users interact with the code.
"""
import os
from typing import Dict, Tuple, List
import warnings

from pygen_structures._functions_const import (
    TOPPAR_DIRECTORY, ACID_CODE_TO_CHARMM
)
from pygen_structures.mol_containers import Atom, Molecule
from pygen_structures.charmm_containers import (
    CHARMMParameterFile, CHARMMResidueTopologyFile
)


def load_charmm_dir(
    directory_root: str = TOPPAR_DIRECTORY
) -> (CHARMMResidueTopologyFile, CHARMMParameterFile):
    """
    Scans a CHARMM toppar folder and reads in the rtf, prm and str
    files to a CHARMMResidueTopologyFile and a CHARMMParameterFile.
    Takes the highest version of the same file.

    :param str: the path to the toppar folder. Defaults to builtin.
    :return: ``CHARMMResidueTopologyFile``, ``CHARMMParameterFile``

    """
    rtf, prm = CHARMMResidueTopologyFile(), CHARMMParameterFile()
    rtf_files, prm_files = {}, {}
    for root, dirs, file_names in os.walk(directory_root, topdown=True):
        dirs[:] = [d for d in dirs if d != "non_charmm"]
        for file_name in file_names:
            extension = os.path.splitext(file_name)[1]
            split_name = file_name.split('_')
            try:
                version = split_name[1]
            except IndexError:
                continue

            file_type = "_".join(split_name[2:])
            if not file_type:
                continue
            file_data = (root, file_name, version)

            if extension in (".rtf", ".str"):
                if file_type not in rtf_files:
                    rtf_files[file_type] = []
                rtf_files[file_type].append(file_data)
            if extension in (".prm", ".str"):
                if file_type not in prm_files:
                    prm_files[file_type] = []
                prm_files[file_type].append(file_data)

    for file_type, versions in rtf_files.items():
        newest_version = max(versions, key=lambda x: x[-1])
        root, file_name, version = newest_version
        file_path = os.path.join(root, file_name)
        rtf.read_file(file_path)
    for file_type, versions in prm_files.items():
        newest_version = max(versions, key=lambda x: x[-1])
        root, file_name, version = newest_version
        file_path = os.path.join(root, file_name)
        prm.read_file(file_path)
    return rtf, prm


def pdb_to_mol(
        pdb_path: str,
        rtf: CHARMMResidueTopologyFile,
        patches: (Dict[str, Tuple[(int, str)]], None) = None,
        default_histidine: str = "HSE",
        kept_chain: str = "A",
        segment_name: (str, None) = None,
) -> Molecule:
    """
    Generate a Molecule from a PDB file. PDB files are assumed to
    contain only a single molecule.

    Works by reading the list of residues, assigning the coordinates
    based on the atom names, and then generating the positions for
    the missing atoms. At present, only works for short and very
    simple sequences, and does not try to perceive missing residues.
    Patches must be manually applied manually.

    :param pdb_path: the path to a PDB file
    :param rtf: a ``CHARMMResidueTopologyFile`` to supply patch and\
    residue definitions.
    :param patches: a dict mapping patch names to a list of the\
    indices of the residues they are to be applied to.
    :param default_histidine: default histidine version to use.
    :param kept_chain: PDB chain to keep (atoms with no chain are)\
    always kept.
    :param segment_name: segment ID to be used in the ``Molecule``
    :return: a ``Molecule`` generated from the PDB file.
    """
    warnings.warn(
        "In the current version, this method only works for short and "
        "very simple sequences, and does not perceive residues beyond "
        "simply looking at the residue names and atom names. Further "
        "functionality is in the roadmap."
    )
    atoms, residue_names, name = [], [], "UNTITLED"
    last_residue, residue_serial = None, 0
    old_res_to_new, coordinates = {}, {}
    rtf.residues["HIS"] = rtf.residues[default_histidine]
    with open(pdb_path) as pdb_file:
        for line in pdb_file:
            if line[:6] == "COMPND":
                name = line.split()[1]
            if line[:6] in ("ATOM  ", "HETATM"):
                atom = Atom.from_pdb_line(line)
                if atom.chain not in (kept_chain, " ", ""):
                    continue
                elif atom.residue_name not in rtf.residues:
                    err = "Unrecognised residue {}".format(atom.residue_name)
                    raise ValueError(err)
                if atom.residue_number != last_residue:
                    residue_serial += 1
                    last_residue = atom.residue_number
                    old_res_to_new[atom.residue_number] = residue_serial
                    residue_names.append(atom.residue_name)
                if not segment_name and atom.segment_id:
                    segment_name = atom.segment_id
                atom.residue_number = old_res_to_new[atom.residue_number]
                atom_id = (residue_serial - 1, atom.atom_name)
                coordinates[atom_id] = (atom.x, atom.y, atom.z)
                atoms.append(atom)

    residues = []
    for residue_index, residue_name in enumerate(residue_names):
        if residue_name == "HIS":
            residue_name = default_histidine
        residue_definition = rtf.residues[residue_name]
        residue = residue_definition.to_residue(residue_index)
        residues.append(residue)
    molecule = Molecule(
        name=name,
        residues=residues,
        topology=rtf,
        patches=patches,
        segment=segment_name,
        fixed_atoms=coordinates,
        use_etkdg=False
    )
    molecule.finalize()
    return molecule


def sequence_to_mol(
    sequence: List[str],
    rtf: CHARMMResidueTopologyFile,
    patches: (Dict[str, List[int]], None) = None,
    name: (str, None) = None,
    segid: (str, None) = None,
) -> Molecule:
    """
    A function for the lazy. Calls ``Molecule.from_sequence`` with the
    supplied args.

    :param sequence: a list of CHARMM residue names
    :param rtf: a ``CHARMMResidueTopologyFile`` to supply patch and\
    residue definitions.
    :param patches: a dict mapping patch names to a list of the indices\
    of the residues they are to be applied to.
    :param name: name of the molecule for the PSF/PDB
    :return: a ``Molecule`` generated from the given sequence.

    """
    return Molecule.from_sequence(
        sequence, rtf, patches, name, segid
    )


def code_to_mol(
    sequence: str,
    rtf: CHARMMResidueTopologyFile,
    patches: (Dict[str, List[int]], None) = None,
    name: (str, None) = None,
    segid: (str, None) = None,
    default_histidine="HSE"
) -> Molecule:
    """
    Generate a ``Molecule`` from one letter protein code.

    First and last patches can be supplied in the sequence string by
    providing patch names separated by dashes:

    - e.g. ``'NNEU-AFK-CT2'``

    :param sequence: a sequence of one letter protein codes
    :param rtf: a ``CHARMMResidueTopologyFile`` to supply patch and\
    residue definitions.
    :param patches: a dict mapping patch names to a list of the indices\
    of the residues they are to be applied to.
    :param default_histidine: default histidine version to use.
    :param name: name of the molecule for the PSF/PDB
    :return: A ``Molecule`` with the given sequence

    """
    try:
        first, sequence, last = sequence.split('-')
        if patches is None:
            patches = {}
        patches[first] = ["FIRST"]
        patches[last] = ["LAST"]
    except ValueError:
        pass
    residue_names = []

    next_d = False
    for character in sequence:
        if character == 'd':
            next_d = True
            continue

        if next_d:
            character = 'd' + character
            next_d = False

        charmm_res = ACID_CODE_TO_CHARMM[character]
        if character == 'H':
            charmm_res = default_histidine
        elif character == 'dH':
            charmm_res = 'D' + default_histidine
        residue_names.append(charmm_res)
    return Molecule.from_sequence(
        residue_names, rtf, patches, name, segid
    )
