"""
pygen-structures (pigeon structures) is a Python utility which allows
for the generation of 3 dimensional molecular structures which can be
used in molecular dynamics or Monte Carlo simulations.

Molecules are generated from a list of residues and patches in the
format of the CHARMM forcefield, and can be written out as valid PSF
and PDB files.

The project contains classes to work with CHARMM topology and parameter
files, classes which represent CHARMM residues and patch residues,
and convenience functions to load CHARMM data and generate molecules
with arbitrary sequences. These molecules can be written out to PSF and
PDB to run simulations.

"""
from pygen_structures.mol_containers import Atom, Structure, Molecule
from pygen_structures.charmm_containers import (
    CHARMMResidueTopologyFile,
    CHARMMParameterFile,
    CHARMMResidueDefinition,
    CHARMMResidue,
    CHARMMPatchResidueDefinition
)
from pygen_structures.convenience_functions import (
    load_charmm_dir,
    sequence_to_mol,
    code_to_mol
)
