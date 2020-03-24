Library Usage
==================

Basic Usage
-----------

With the convenience functions, library usage is almost as simple as
command line usage. The default toppar directory can be loaded Using
``load_charmm_dir``, and ``Molecule``s can be generated using ``code_to_mol`` (
using the single amino acid code) or ``sequence_to_mol`` (using generic CHARMM
residues).

To make ``HIS-GLU-TYR``, as outlined in the :doc:`command_line_usage`:

.. code-block:: python

    >>> import pygen_structures as ps
    >>> # Using the supplied toppar directory
    >>> rtf, prm = ps.load_charmm_dir()
    >>> mol = ps.code_to_mol('HEY', rtf)

``Molecule``s can also be loaded from PDB files, though at present this only
works for short residues.
