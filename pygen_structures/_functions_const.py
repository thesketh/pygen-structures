from typing import Iterable, Tuple, List
import pkg_resources

from rdkit import Chem
import numpy as np

IndexCollection = List[Tuple[int, ...]]

TOPPAR_DIRECTORY = pkg_resources.resource_filename(
    'pygen_structures', 'toppar'
)


def adjacency_to_dof(
        adjacency_matrix: np.ndarray,
        sort: bool = True
) -> Tuple[IndexCollection, IndexCollection, IndexCollection]:
    """ Given an adjacency matrix, return the internal degrees of
        freedom of the system. This should work on lists of lists
        or numpy arrays.
    """
    bonds, angles, dihedrals = set(), set(), set()
    for w_index, x_row in enumerate(adjacency_matrix):
        x_indices = np.flatnonzero(x_row)
        for x_index in x_indices:
            bond = (w_index, x_index)
            if bond[::-1] not in bonds:
                bonds.add(bond)
            y_indices = np.flatnonzero(adjacency_matrix[x_index])
            y_indices = [y for y in y_indices if y not in bond]
            for y_index in y_indices:
                angle = (w_index, x_index, y_index)
                if angle[::-1] not in angles:
                    angles.add(angle)
                z_indices = np.flatnonzero(adjacency_matrix[y_index])
                z_indices = [z for z in z_indices if z not in angle]
                for z_index in z_indices:
                    dihedral = (w_index, x_index, y_index, z_index)
                    if dihedral[::-1] not in dihedrals:
                        dihedrals.add(dihedral)
    if sort:
        return sorted(bonds), sorted(angles), sorted(dihedrals)
    else:
        return list(bonds), list(angles), list(dihedrals)


def bonds_to_adjacency(bonds: IndexCollection) -> np.ndarray:
    """
    Build an adjacency matrix from a list of bonds.

    This is an array of integers which indicate bond order.

    E.g. self.adjacency_matrix[x, y] == 2 would indicate a double
    bond between atoms with index x and y.
    """
    n_atoms = 0
    for bond in bonds:
        n_atoms = max(n_atoms, *bond)
    n_atoms += 1

    adjacency_matrix = np.zeros((n_atoms, n_atoms), dtype=int)
    for idx_x, idx_y in bonds:
        adjacency_matrix[idx_x, idx_y] += 1
        adjacency_matrix[idx_y, idx_x] += 1
    return adjacency_matrix


def iter_nwise(iterable: Iterable, n: int = 2) -> Iterable:
    """
    Returns a generator which iterates over an iterable n items at a
    time.
    :param iterable: Sequence which can have `iter` called on it
    :param n:        Number of items the generator should return
                     at one time
    :return:         A generator which iterates over `iterable` `n`
                     items at a time
    """

    iterable = iter(iterable)
    list_of_iterables = [iterable] * n
    return zip(*list_of_iterables)


PERIODIC_TABLE = Chem.GetPeriodicTable()


def get_expected_valence(symbol: str) -> int:
    """
    Return expected valence for a given element symbol.
    """
    return PERIODIC_TABLE.GetDefaultValence(symbol)


BOND_ORDER_TO_TYPE = {
    1: Chem.BondType.SINGLE,
    2: Chem.BondType.DOUBLE,
    3: Chem.BondType.TRIPLE,
    4: Chem.BondType.QUADRUPLE,
}

ACID_CODE_TO_CHARMM = {
    'A': 'ALA',
    'R': 'ARG',
    'N': 'ASN',
    'D': 'ASP',
    'C': 'CYS',
    'E': 'GLU',
    'Q': 'GLN',
    'G': 'GLY',
    'H': 'HSD',
    'I': 'ILE',
    'L': 'LEU',
    'K': 'LYS',
    'M': 'MET',
    'F': 'PHE',
    'P': 'PRO',
    'S': 'SER',
    'T': 'THR',
    'W': 'TRP',
    'Y': 'TYR',
    'V': 'VAL'
}
for code, resi_name in list(ACID_CODE_TO_CHARMM.items()):
    if code != 'G':  # D-glycine isn't real, it can't hurt you.
        ACID_CODE_TO_CHARMM['d' + code] = 'D' + resi_name

CHARMM_TO_ACID_CODE = {
    resi_name: code for code, resi_name in ACID_CODE_TO_CHARMM.items()
}
CHARMM_TO_ACID_CODE["HSE"] = "H"
CHARMM_TO_ACID_CODE["DHSE"] = "dH"
CHARMM_TO_ACID_CODE["HSP"] = "H[+]"
CHARMM_TO_ACID_CODE["DHSP"] = "dH[+]"
