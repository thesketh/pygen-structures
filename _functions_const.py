from rdkit import Chem
import numpy as np
from typing import Iterable, Tuple, List


def adjacency_to_dof(adjacency_matrix, sort=True) -> Tuple[List[int]]:
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
