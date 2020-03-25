Library Usage
==================

Representations
---------------

pygen-structures contains classes to represent CHARMM data, molecular (and
molecular-like) structures and atoms.

Representing CHARMM File Types
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The CHARMM forcefield views molecules as collections of _residues_, with
_patches_ usually applied to define linkages and modifications (though for some
classes of residue, such as the amino acids, linkages between residues are
implicit).

The residue topology file (``.rtf``) contains the definition of these
patches and residues - the information required to build the structure of the
molecule. This includes the types of the atoms, the  necessary bonds, the
chirality of the tetrahedral centres (through the internal coordinate table),
the charges, and which internal coordinates are necessary for the simulation
(i.e. which improper torsions/cross-map terms are required). The masses
and charges are also present, as these are required to build the ``.psf`` file.

The parameter file (``.prm``) contains the definitions of all the parameterised
interactions (i.e. which bonds, angles, dihedrals, impropers, cross maps are
defined, and what the radii of the atoms are for the non-bonded terms).

For some residue and parameter definitions, particularly extensions to existing
files, the topology and parameter files are combined into a stream file
(``.str``) which contains a topology and residue section concatenated together.

A good overview of the actual layout of the CHARMM forcefield file types (both
`residue topology files`__ and `parameter files`__) is given in the NAMD
tutorial.

__ https://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/node24.html
__ https://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/node25.html

As pygen-structures is concerned with building structures and not with running
simulations, the representation of the ``.rtf`` files is complete,
but the representation of ``.prm`` files is limited to verifying that
the parameters are present rather than probing their values.

**Residue Files**

CHARMM ``.rtf`` files are represented using the
:class:`~pygen_structures.charmm_containers.CHARMMResidueTopologyFile` class,
which can be instantiated from a ``.rtf`` or ``.str`` file using its
:meth:`~pygen_structures.charmm_containers.CHARMMResidueTopologyFile.from_file`
method, and additional files can be added using the
:meth:`~pygen_structures.charmm_containers.CHARMMResidueTopologyFile.read_file`
method. Residue and patch definitions within the file are stored as instances of
:class:`~pygen_structures.charmm_containers.CHARMMResidueDefinition` and
:class:`~pygen_structures.charmm_containers.CHARMMPatchResidueDefinition`
respectively. Additional residue definitions and patch
definitions can be added using
:meth:`~pygen_structures.charmm_containers.CHARMMResidueTopologyFile.add_residue_definition`
and
:meth:`~pygen_structures.charmm_containers.CHARMMResidueTopologyFile.add_patch_definition`.

The residue definitions are stored in the ``dict`` attribute
:attr:`~pygen_structures.charmm_containers.CHARMMResidueTopologyFile.residues`,
and the patch definitions in the ``dict`` attribute
:attr:`~pygen_structures.charmm_containers.CHARMMResidueTopologyFile.patches`,
where the keys are the relevant residue/patch name. The masses of the atoms in
Daltons are stored in ``dict`` attribute
:attr:`~pygen_structures.charmm_containers.CHARMMResidueTopologyFile.masses`,
where keys are the atom types.

**Parameter Files**

CHARMM ``.prm`` files are represented using the
:class:`~pygen_structures.charmm_containers.CHARMMParameterFile` class.
Like
:class:`~pygen_structures.charmm_containers.CHARMMResidueTopologyFile`,
this class can be instantiated using its
:meth:`~pygen_structures.charmm_containers.CHARMMParameterFile.from_file`
method and additional files can be read using the
:meth:`~pygen_structures.charmm_containers.CHARMMParameterFile.read_file`
method. Methods
:meth:`~pygen_structures.charmm_containers.CHARMMParameterFile.get_unmatched`
and
:meth:`~pygen_structures.charmm_containers.CHARMMParameterFile.check_parameters`
can be used to assert that all parameters of a Molecule exist in the parameter
file, with the former returning a ``dict`` of internal coordinate type
to a ``list`` of unmatched tuples and the latter
a ``bool`` (``True`` if all parameters exist, ``False`` otherwise).

This class contains five ``set`` objects:
:attr:`~pygen_structures.charmm_containers.CHARMMResidueTopologyFile.bonds`,
:attr:`~pygen_structures.charmm_containers.CHARMMResidueTopologyFile.angles`,
:attr:`~pygen_structures.charmm_containers.CHARMMResidueTopologyFile.dihedrals`,
:attr:`~pygen_structures.charmm_containers.CHARMMResidueTopologyFile.impropers`,
and
:attr:`~pygen_structures.charmm_containers.CHARMMResidueTopologyFile.cross_maps`,
which contain tuples of atom types from the forcefield. In the first three sets,
the tuples are also stored in reverse, as these degrees of freedom can be
specified in reverse order without altering the relevant parameters.


Representing Residues and Patches
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**The Difference Between Definitions and Residues**

CHARMM residues as presented in ``.rtf`` files are an abstract residue
definition, rather than physical residues. As such, these lack an associated
residue index, and their degrees of freedom (the bonds, impropers and cross
maps) refer to the atom names from the residue definition. These are represented
using the :class:`~pygen_structures.charmm_containers.CHARMMResidueDefinition`
class.

CHARMM residues in an actual structure are not abstract, and refer to real
residues in a physical molecule. As such, these _do_ have an associated index
and their degrees of freedom refer to atom IDs - tuples of
``(residue_index: int, atom_name: str)``. This is to allow for easier removal of
dangling bonds and easier application of patches. These residues are represented
using the :class:`~pygen_structures.charmm_containers.CHARMMResidue` class.
Patch definitions can only be applied to real residues.

**CHARMM Residue Definitions**

**CHARMM Residues**

Representing Molecules
^^^^^^^^^^^^^^^^^^^^^^

pygen-structures uses three main classes for molecular representation.

**Molecules**

**Molecule-like Structures**

**Atoms**

Convenience Functions
---------------------

Convenience functions represent the easiest way to interact with the code as
a library.

**Loading CHARMM Data**

**Creating Molecules**

**Creating Molecules from PDB Files**

Basic Library Usage
-------------------

With the convenience functions, library usage is almost as simple as command
line usage. The default toppar directory can be loaded using
:func:`~pygen_structures.convenience_functions.load_charmm_dir`, and
:class:`~pygen_structures.mol_containers.molecule.Molecule` objects can be
generated using :func:`~pygen_structures.convenience_functions.code_to_mol`
(from a ``str`` of single amino acid codes) or
:func:`~pygen_structures.convenience_functions.sequence_to_mol` (from a ``list``
of generic CHARMM residue names). ``Molecule`` objects can also be loaded from
PDB files using :func:`~pygen_structures.convenience_functions.pdb_to_mol`,
though at present this only works for small molecules and makes no attempt to
pattern match the missing residues.

Examples
--------

All examples should be prefaced with

.. code-block:: python

    >>> import pygen_structures as ps
    >>> # Using the supplied toppar directory
    >>> rtf, prm = ps.load_charmm_dir()

To make ``HIS-GLU-TYR``, as outlined in the :doc:`command_line_usage`:

.. code-block:: python

    >>> mol = ps.code_to_mol('HEY', rtf)
    >>> mol.to_pdb_file('HEY.pdb')
    >>> mol.to_psf_file('HEY.psf')

And should we want to use the protonated form of histidine:

.. code-block:: python

    >>> sequence = ['HSP', 'GLU', 'TYR']
    >>> mol = ps.sequence_to_mol(sequence, rtf)
    >>> # or
    >>> mol = ps.code_to_mol('HEY', rtf, default_histidine='HSP')
    >>> mol.to_pdb_file('HEY.pdb')
    >>> mol.to_psf_file('HEY.psf')

Back to non-protein examples, we could create the trisaccharide raffinose.
This requires the use of a patch. The default segment ID, ``PROT``, is less
applicable here, so we can specify a new segment, ``RAFF``. While we're at it,
we can set the name for the ``COMPND`` header to ``Raffinose``:

.. code-block:: python

    >>> sequence = ['AGLC', 'BFRU', 'AGAL']
    >>> patches = {'RAFF': [0, 1, 2]}
    >>> mol = ps.sequence_to_mol(sequence, rtf, patches, 'Raffinose', 'RAFF')
    >>> mol.to_pdb_file('RAFF.pdb')
    >>> mol.to_psf_file('RAFF.psf')

Verification that parameters exist for created structures is simple:

.. code-block:: python

    >>> mol = ps.code_to_mol('AdP', rtf)
    >>> # This will pass verification in 0.2.3
    >>> mol.check_parameters(prm)
    False
    >>> unmatched = prm.get_unmatched(mol)
    >>> unmatched
    {'bonds': {('CPD1', 'CC')}, 'angles': set(), 'dihedrals': set(), 'impropers': set(), 'cross_maps': set()}
    >>> mol2 = ps.code_to_mol('AP', rtf)
    >>> prm.check_parameters(mol2)
    True
    >>> unmatched = prm.get_unmatched(mol2)
    >>> unmatched
    {'bonds': set(), 'angles': set(), 'dihedrals': set(), 'impropers': set(), 'cross_maps': set()}

Molecules can be created manually, though it is far easier to use the
convenience functions.