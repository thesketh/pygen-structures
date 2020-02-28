"""
More specific `typing` types for the data structures used in this
package.

These are not currently in use, but are in the roadmap.
"""
from typing import Tuple, Union

AtomName = str
AtomType = str
AtomPartialCharge = float
AtomIndex = int
ResidueIndex = int
ResiduePlaceholderIndex = int

AtomID = Tuple[ResidueIndex, AtomName]
AtomReference = Tuple[ResiduePlaceholderIndex, AtomName]
AtomData = [AtomName, AtomType, AtomPartialCharge]

Position = Tuple[float, float, float]

# Molecular connections
Bond = Tuple[AtomID, AtomID]
BondPlaceholder = Tuple[AtomReference, AtomReference]
BondDefinition = Tuple[AtomName, AtomName]
IndexBond = Tuple[AtomIndex, AtomIndex]

Angle = Tuple[AtomID, AtomID, AtomID]
AnglePlaceholder = Tuple[AtomReference, AtomReference, AtomReference]
AngleDefinition = Tuple[AtomName, AtomName, AtomName]
IndexAngle = Tuple[AtomIndex, AtomIndex, AtomIndex]

Dihedral = Tuple[AtomID, AtomID, AtomID, AtomID]
DihedralPlaceholder = Tuple[
    AtomReference, AtomReference, AtomReference, AtomReference
]
DihedralDefinition = Tuple[AtomName, AtomName, AtomName, AtomName]
IndexDihedral = Tuple[AtomIndex, AtomIndex, AtomIndex, AtomIndex]

Improper = Tuple[AtomID, AtomID, AtomID, AtomID]
ImproperPlaceholder = Tuple[
    AtomReference, AtomReference, AtomReference, AtomReference
]
ImproperDefinition = Tuple[AtomName, AtomName, AtomName, AtomName]
IndexImproper = Tuple[AtomIndex, AtomIndex, AtomIndex, AtomIndex]

CrossMap = Tuple[
    AtomID, AtomID, AtomID, AtomID, AtomID, AtomID, AtomID, AtomID
]
CrossMapPlaceholder = Tuple[
    AtomReference,
    AtomReference,
    AtomReference,
    AtomReference,
    AtomReference,
    AtomReference,
    AtomReference,
    AtomReference
]
CrossMapDefinition = Tuple[
    AtomName,
    AtomName,
    AtomName,
    AtomName,
    AtomName,
    AtomName,
    AtomName,
    AtomName
]
IndexCrossMap = Tuple[
    AtomIndex,
    AtomIndex,
    AtomIndex,
    AtomIndex,
    AtomIndex,
    AtomIndex,
    AtomIndex,
    AtomIndex,
]
# Generic connections
Connection = Union[
    Bond,
    Angle,
    Dihedral,
    Improper,
    CrossMap
]
ConnectionPlaceholder = Union[
    BondPlaceholder,
    AnglePlaceholder,
    DihedralPlaceholder,
    ImproperPlaceholder,
    CrossMapPlaceholder
]
ConnectionDefinition = Union[
    BondDefinition,
    AngleDefinition,
    DihedralDefinition,
    ImproperDefinition,
    CrossMapDefinition
]
IndexConnection = Union[
    IndexBond,
    IndexAngle,
    IndexDihedral,
    IndexImproper,
    IndexCrossMap
]
