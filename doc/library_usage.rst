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


Residue Files
"""""""""""""

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


Parameter Files
"""""""""""""""

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

The Difference Between Definitions and Residues
"""""""""""""""""""""""""""""""""""""""""""""""

CHARMM residues as presented in ``.rtf`` files are an abstract residue
definition, rather than physical residues. As such, these lack an associated
residue index, and their degrees of freedom (the bonds, impropers and cross
maps) refer to the atom names from the residue definition. These are represented
using the :class:`~pygen_structures.charmm_containers.CHARMMResidueDefinition`
class.

CHARMM residues in an actual structure are not abstract, and refer to real
residues in a physical molecule. As such, these _do_ have an associated index
and their degrees of freedom refer to *atom IDs* - tuples of
``(residue_index: int, atom_name: str)``. This is to allow for easier removal of
dangling bonds and easier application of patches. These residues are represented
using the :class:`~pygen_structures.charmm_containers.CHARMMResidue` class.
Patch definitions can only be applied to real residues.

The patch definitions, instead of the atom IDs used in the residues, use
*atom references*. These take the same format as atom IDs, but the residue index
is replaced by an integer referring to the order the residues are provided to
the patch. For the ``RAFF`` patch, which creates the sugar raffinose from
``[AGLC, BFRU, AGAL]``, atoms in ``AGLC`` would take the form
``(0, atom_id: str)`` (``1`` for ``BFRU``, ``2`` for ``AGAL``). Patch
definitions are represented using the
:class:`~pygen_structures.charmm_containers.CHARMMPatchResidueDefinition`
class.

CHARMM Residue Definitions
""""""""""""""""""""""""""

CHARMM residue definitions
(:class:`~pygen_structures.charmm_containers.CHARMMResidueDefinition`)
represent all of the data from a ``RESI`` block in a ``.rtf`` file.
This includes:

-  the ``name``
-  a ``list`` of atoms from the ATOM records (represented as
   ``(atom_name: str, atom_type: str, partial_charge: float)``)
-  a ``list`` of bonds - ``tuple`` of 2 atom_names - from BOND/DOUB records.
   Bonds in DOUB records appear twice.
-  a ``list`` of impropers – ``tuple`` of 4 atom_name from IMPR records.
-  a ``list`` of cross_maps – ``tuple`` of 8 atom_name from CMAP records
-  a list of IC records. The first 4 items are atom names,
   the last 5 are floats containing bond length/angle information.
   Where angle i-j-k-l is an improper, the third atom_name is prefaced with
   ``*``. The floats are:

   -  the bond length, i-j
   -  the angle, i-j-k
   -  the dihedral (or improper) i-j-k-l
   -  the angle j-k-l
   -  the bond length k-l

-  the default patch if residue is first in the chain. If ``None``, this
   defaults to the default first patch of the
   :class:`~pygen_structures.charmm_containers.CHARMMResidueTopologyFile` the
   residue definition is in. There is a semantic difference between ``None``
   and ``'NONE'`` - the former represents missing information, and the latter
   is an explicit "no first patch".
-  the default patch if residue is last in the chain. Treated much the same way
   as the first.

Residue definitions are usually created using the class method
:meth:`~pygen_structures.charmm_containers.CHARMMResidueDefinition.from_block`
when a :class:`~pygen_structures.charmm_containers.CHARMMResidueTopologyFile`
is created. Residue definitions can be converted to residues using the
:meth:`~pygen_structures.charmm_containers.CHARMMResidueDefinition.to_residue`.

Residue definitions can be represented as
:class:`~pygen_structures.mol_containers.structure.Structure` objects,
as RDKit ``Mol`` objects, and as SMARTS patterns using methods
:meth:`~pygen_structures.charmm_containers.CHARMMResidueDefinition.to_structure`,
:meth:`~pygen_structures.charmm_containers.CHARMMResidueDefinition.to_fragment_mol`,
and :meth:`~pygen_structures.charmm_containers.CHARMMResidueDefinition.to_smarts`.
Not much is done with this at present, but in principle this could be used to
match missing residues in a ``.pdb`` file.

CHARMM Residues
"""""""""""""""

CHARMM residues (:class:`~pygen_structures.charmm_containers.CHARMMResidue`)
contain the information from the residue definitions, but in terms
of atom IDs instead of atom names. The representation of the atoms themselves
is unchanged. Residues have an associated index, and may be altered using
patches. As such, they may contain new atoms and bonds which are not present in
the definition alone.

CHARMM Patches
""""""""""""""

CHARMM patch residue definitions
(:class:`~pygen_structures.charmm_containers.CHARMMPatchResidueDefinition`)
contain the information from a ``PRES`` block in a ``.rtf`` file.

In general, the representation is similar to
:class:`~pygen_structures.charmm_containers.CHARMMResidueDefinition`, but
rather than a ```list`` of ``tuple`` for each degree of freedom, these are
stored as ``n_residues`` long lists of these lists. For a patch which is to be
applied to three residues, ``bonds`` might look something like this:

.. code-block:: python

    bonds = [
        [((0, 'C'), (1, 'N')), ...],
        [((1, 'CA'), (1, 'CB')), ...],
        [((2, 'CA'), (2, 'CB')), ...]
    ]

When the patch is applied, ``0``, ``1`` and ``2`` are replaced with the indices
of the residues the patch is applied to, and these new bonds are added to those
residues. Patches also contain ``deletions``, atom names to be deleted in each
residue. When atoms are deleted, all degrees of freedom containing these atoms
are removed.

Patches have two important methods:

-  :meth:`~pygen_structures.charmm_containers.CHARMMPatchResidueDefinition.apply`,
   the method which applies the patch to a list of residues, and
-  :meth:`~pygen_structures.charmm_containers.CHARMMResidueDefinition.is_applicable_to`,
   which can be used to test what positions a residue can be supplied to the patch in.
   This is based on the presence of atoms in the ``deletions`` list for each
   residue position.


What's Special About FIRST and LAST?
""""""""""""""""""""""""""""""""""""

If patches are specified as being applied to ``'FIRST'`` or ``'LAST'``
rather than ``0`` or ``-1``, the default first/last patch is not applied.
This is to distinguish between cases where a different terminal patch is
being applied (such as if protein patch ``CT2`` was being applied instead of
``CTER``) and cases where patches just happen to affect the first or last
residue.


Representing Molecules
^^^^^^^^^^^^^^^^^^^^^^

pygen-structures uses three main classes for molecular representation.


Molecules
"""""""""

Complete molecules from the perspective of the CHARMM forcefield are stored as
:class:`~pygen_structures.mol_containers.molecule.Molecule` classes. These
are directly representable as ``.pdb`` or ``.psf`` files.

There are not many cases where it is advantageous to create a molecule directly.
In most cases, it makes more sense to use the convenience functions to create
molecules.

Important methods for the molecule class include:

-  :meth:`~pygen_structures.mol_containers.molecule.Molecule.finalize`,
   which removes dangling bonds, applies patches, and instantiates the
   :class:`~pygen_structures.mol_containers.structure.Structure` used to
   generate coordinates.
-  :meth:`~pygen_structures.mol_containers.molecule.Molecule.to_pdb_block`/
   :meth:`~pygen_structures.mol_containers.molecule.Molecule.to_pdb_file`, used
   to create a PDB file from the molecule.
-  :meth:`~pygen_structures.mol_containers.molecule.Molecule.to_psf_block`/
   :meth:`~pygen_structures.mol_containers.molecule.Molecule.to_psf_file`, used
   to create a PSF file from the molecule.
-  :meth:`~pygen_structures.mol_containers.molecule.Molecule.check_parameters`,
   used to check parameters against a
   :class:`~pygen_structures.charmm_containers.CHARMMParameterFile`
-  :meth:`~pygen_structures.mol_containers.molecule.Molecule.to_mol`, used to
   obtain the RDKit ``Mol`` generated by the underlying ``Structure``.


Molecule-like Structures
""""""""""""""""""""""""

Molecular-like structures can be whole molecules or simply fragments of
molecules. These are represented using the
:class:`~pygen_structures.mol_containers.structure.Structure` class.

This class generates molecular connection tables and instantiates RDKit ``Mol``
objects based on the provided structures. This enables 3D coordinate generation.
Chirality is set using information from the IC table. Atoms are reordered so
that hydrogen atoms come after their parent heavy atom. This is in line with
the psfgen order.

Structures can be created from a ``list`` of
:class:`~pygen_structures.mol_containers.atom.Atom` and a ``list`` of bond
``tuple`` containing two 0-based atom indices. Generated 3D coordinates are set
in the RDKit ``Mol``, stored in
:attr:`~pygen_structures.mol_containers.structure.Structure.mol`
and in the ``x``, ``y``, and ``z`` attributes of the ``Atom`` objects.

The RDKit ``Mol`` can also be retrieved using the
:meth:`~pygen_structures.mol_containers.structure.Structure.to_mol` method.


Atoms
"""""

Atoms, represented using the :class:`~pygen_structures.mol_containers.atom.Atom`
class, are mostly based on the fields stored in a ``.pdb`` file, with some extra
``.psf`` specific fields. These fields appear in PDB order, with PSF fields
attached on the end.

Atoms can be instantiated from PDB lines using the
:meth:`~pygen_structures.mol_containers.atom.Atom.from_pdb_line` class method,
though this will not contain the PSF data.

Atoms can be represented as PDB atom lines using
:meth:`~pygen_structures.mol_containers.atom.Atom.to_pdb_line` and as PSF lines
using :meth:`~pygen_structures.mol_containers.atom.Atom.to_psf_line`.

Where the atom serial number exceeds 99999, PDB atom serials are switched to
hexadecimal encoding. This is unlikely to happen with the current functionality
(which is limited to around 15 residues). PSF lines are in the extended format.


Convenience Functions
---------------------

Convenience functions represent the easiest way to interact with the code as
a library.


Loading CHARMM Data
^^^^^^^^^^^^^^^^^^^

The CHARMM forcefield topology and parameters can be loaded using
:func:`~pygen_structures.convenience_functions.load_charmm_dir`. This function,
called with no arguments, will load the `most recent (July 2019) CHARMM
distribution (with some modifications to the D-amino acid parameters)`__
and return it as a
:class:`~pygen_structures.charmm_containers.CHARMMResidueTopologyFile` and a
:class:`~pygen_structures.charmm_containers.CHARMMParameterFile`.

By specifying a path as a positional argument, that directory will be loaded
instead. The function will pick the latest versions of the parameter and
topology files (36 over 27, 36m over 36), so if you plan on using an older
version of the forcefield (this is not recommended) you will have to remove the
newer versions and change the file extensions to match the current conventions
(``.rtf`` for topology files and, ``.prm`` for parameter files).

__ https://github.com/thesketh/pygen-structures/blob/master/pygen_structures/toppar/README.md


Creating Molecules
^^^^^^^^^^^^^^^^^^

:class:`~pygen_structures.mol_containers.molecule.Molecule` objects can be
generated using :func:`~pygen_structures.convenience_functions.code_to_mol`
(from a ``str`` of single amino acid codes) or
:func:`~pygen_structures.convenience_functions.sequence_to_mol` (from a ``list``
of generic CHARMM residue names).

These functions accept the same arguments, with the exception that
:func:`~pygen_structures.convenience_functions.code_to_mol` expects a string
protein sequence as the first arguemnt and
:func:`~pygen_structures.convenience_functions.sequence_to_mol` expects a
list of names. :func:`~pygen_structures.convenience_functions.code_to_mol`
accepts an additional keyword argument, ``default_histidine``, which determines
the residue definition which is used for ``H``.

Both functions require as a minimum the sequence/code to create and a
:class:`~pygen_structures.charmm_containers.CHARMMResidueTopologyFile`.


Creating Molecules from PDB Files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:class:`~pygen_structures.mol_containers.molecule.Molecule` objects can also be
loaded from PDB files using
:func:`~pygen_structures.convenience_functions.pdb_to_mol`, though at present
this only works for small molecules and makes no attempt to pattern match the
missing residues.


Examples
--------

All examples should be prefaced with

.. code-block:: python

    >>> import pygen_structures as ps
    >>> # Using the supplied toppar directory
    >>> rtf, prm = ps.load_charmm_dir()


Basic Usage
^^^^^^^^^^^

To make ``HIS-GLU-TYR``, as outlined in the :doc:`command_line_usage`, writing
the structure out to ``HEY.psf`` and ``HEY.pdb``:

.. code-block:: python

    >>> mol = ps.code_to_mol('HEY', rtf)
    >>> mol.to_pdb_file('HEY.pdb')
    >>> mol.to_psf_file('HEY.psf')

And should we want to use the protonated form of histidine:

.. code-block:: python

    >>> mol = ps.code_to_mol('HEY', rtf, default_histidine='HSP')
    >>> # or use the sequence
    >>> sequence = ['HSP', 'GLU', 'TYR']
    >>> mol = ps.sequence_to_mol(sequence, rtf)
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

Loading raffinose back from PDB currently requires us to specify the patches:

.. code-block:: python

    >>> patches = {'RAFF': [0, 1, 2]}
    >>> mol = ps.convenience_functions.pdb_to_mol('RAFF.pdb', rtf, patches)

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

For linkages between protein and sugars, we need to use a patch. We
can make the tripeptide ``ALA-ASN-ALA``, where glucose is linked to asparagine.
This requires use of the ``NGLA`` patch and is a good example of one of the
quirks of the CHARMM forcefield definitions. If we specify ``AGLC`` between
``ASN`` and ``ALA``, the protein chain will be broken - bonds between residues
in the peptide are defined using links to the previous (``-C``) and next
(``+N``) residues and these don't exist where there is a sugar in the way.
Unfortunately, this means there is a "correct" way to specify these residues.

.. code-block:: python

    >>> sequence = ['ALA', 'ASN', 'ALA', 'AGLC']
    >>> patches = {'CTER': [-2], 'NGLA': [1, -1]}
    >>> mol = ps.sequence_to_mol(sequence, rtf, patches)
    >>> # This order wouldn't link ASN and ALA because they aren't adjacent,
    >>> # leading to a free-floating ALA residue.
    >>> sequence = ['ALA', 'ASN', 'AGLC', 'ALA']


Manual Molecule Creation
^^^^^^^^^^^^^^^^^^^^^^^^

If when generating a glycopeptide you **absolutely** require the sugar to be
in the middle of the protein chain, it is possible to obtain this by creating
the residues manually, and specifying the ``last_index`` and ``next_index``.

.. code-block:: python

    >>> def make_residue(residue_name, index, last_index=None, next_index=None):
    ...     return ps.CHARMMResidue.from_residue_definition(
    ...         rtf.residues[residue_name],
    ...         index,
    ...         last_index,
    ...         next_index
    ...     )
    ...
    >>> first = make_residue('ALA', 0)
    >>> second = make_residue('ASN', 1, next_index=3)
    >>> third = make_residue('AGLC', 2)
    >>> last = make_residue('ALA', 3, last_index=1)
    >>> residues = [first, second, third, last]
    >>> patches = {'NGLA': [1, 2]}
    >>> mol = ps.Molecule(residues=residues, topology=rtf, patches=patches)
    >>> mol.finalize()

For the most part, though, creating molecules manually is not much more complex
than using the convenience functions, just slightly more verbose.

.. code-block:: python

    >>> sequence = ['HSE', 'GLU', 'TYR']
    >>> residues = []
    >>> for index, residue_name in enumerate(sequence):
    >>>     residue = rtf.residues[residue_name].to_residue(index)
    >>>     residues.append(residue)
    >>> mol = ps.Molecule('HEY', residues, rtf)
    >>> mol.finalize()
