"""
More specific `typing` types for the data structures used in this
package.

These are not currently in use, but are in the roadmap.
"""
AtomName = str
AtomType = str
AtomPartialCharge = float
AtomIndex = int
ResidueIndex = int
ResiduePlaceholderIndex = int

AtomID = tuple[ResidueIndex, AtomName]
AtomReference = tuple[ResiduePlaceholderIndex, AtomName]
AtomData = [AtomName, AtomType, AtomPartialCharge]

Position = tuple[float, float, float]

# Molecular connections
Bond = tuple[AtomID, AtomID]
BondPlaceholder = tuple[AtomReference, AtomReference]
BondDefinition = tuple[AtomName, AtomName]
IndexBond = tuple[AtomIndex, AtomIndex]

Angle = tuple[AtomID, AtomID, AtomID]
AnglePlaceholder = tuple[AtomReference, AtomReference, AtomReference]
AngleDefinition = tuple[AtomName, AtomName, AtomName]
IndexAngle = tuple[AtomIndex, AtomIndex, AtomIndex]

Dihedral = tuple[AtomID, AtomID, AtomID, AtomID]
DihedralPlaceholder = tuple[
    AtomReference, AtomReference, AtomReference, AtomReference
]
DihedralDefinition = tuple[AtomName, AtomName, AtomName, AtomName]
IndexDihedral = tuple[AtomIndex, AtomIndex, AtomIndex, AtomIndex]

Improper = tuple[AtomID, AtomID, AtomID, AtomID]
ImproperPlaceholder = tuple[
    AtomReference, AtomReference, AtomReference, AtomReference
]
ImproperDefinition = tuple[AtomName, AtomName, AtomName, AtomName]
IndexImproper = tuple[AtomIndex, AtomIndex, AtomIndex, AtomIndex]

CrossMap = tuple[
    AtomID, AtomID, AtomID, AtomID, AtomID, AtomID, AtomID, AtomID
]
CrossMapPlaceholder = tuple[
    AtomReference,
    AtomReference,
    AtomReference,
    AtomReference,
    AtomReference,
    AtomReference,
    AtomReference,
    AtomReference
]
CrossMapDefinition = tuple[
    AtomName,
    AtomName,
    AtomName,
    AtomName,
    AtomName,
    AtomName,
    AtomName,
    AtomName
]
IndexCrossMap = tuple[
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
Connection = Bond | Angle | Dihedral | Improper | CrossMap
ConnectionPlaceholder = BondPlaceholder | AnglePlaceholder | DihedralPlaceholder | ImproperPlaceholder | CrossMapPlaceholder
ConnectionDefinition = BondDefinition | AngleDefinition | DihedralDefinition | ImproperDefinition | CrossMapDefinition
IndexConnection = IndexBond | IndexAngle | IndexDihedral | IndexImproper | IndexCrossMap
