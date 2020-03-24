Command Line Usage
==================

Basic Usage
-----------

Command line usage for peptides is simple, and takes the following form:

.. code-block:: bash

    pygen-structures SEQUENCE -o OUTPUT_PREFIX

Sequences are specified using the one letter protein code by default, and
terminal patches can be supplied by using hyphens as delimiters (e.g.
``NNEU-AFK-CT2``, note that both termini must be supplied). D-amino acids need
only be preceded by a lowercase 'd' (e.g. ``dA`` for D-alanine).

``OUTPUT_PREFIX.psf`` and ``OUTPUT_PREFIX.pdb`` are created. If ``-o`` is not
specified, the content of the ``.pdb`` file is written to stdout and no ``.psf``
file is generated.

The histidine form used can be set using ``--histidine``, and defaults to
``HSE`` (neutral form with proton on N).

Patches can be supplied using ``--patches``, the name of the patch, and the
0-based indices the patch is to be applied to (or strings 'FIRST'/'LAST'
to refer to the first and last residue).

To generate more complex structures, such as sugars, the residue names should
be supplied (hyphen delimited) and the ``-u``/``--use-charmm-names`` option
selected.

``--name`` and ``--segid`` control the names given in the ``COMPND`` record and
the segment ID (4 character max) respectively.

Examples
--------
To produce a simple peptide sequence, the one letter code can be used. To
produce the peptide ``HIS-GLU-TYR``, creating ``HEY.psf`` and ``HEY.pdb``:

.. code-block:: bash

    pygen-structures HEY -o HEY

Supposing we think that histidine should be protonated, we can change the
protonation state of histidine by specifying a different histidine form:

.. code-block:: bash

    pygen-structures HEY -o HEY --histidine HSP

Or simply use the three letter residue codes by using the ``-u`` flag:

.. code-block:: bash

    pygen-structures -u HSP-GLU-TYR -o HEY

Looking at non-protein examples, we could create the trisaccharide raffinose.
This requires the use of the residue codes and a patch. The default segment ID,
``PROT``, is less applicable here, so we can specify that with the ``--segid``
option, and set the name in the ``COMPND`` header using ``--name``. The
following command produces ``RAFF.psf`` and ``RAFF.pdb``:

.. code-block:: bash

    pygen-structures -u AGLC-BFRU-AGAL --patches RAFF 0 1 2 --segid RAFF --name Raffinose -o RAFF

We can also make glycopeptides. To link alpha-glucose to an arginine residue
(in this case, from an ALA-ASN-ALA peptide), we can use the NGLA patch. Note
that because the protein residue is not the last in the chain, we have to apply
the C-terminus patch manually:

.. code-block:: bash

    pygen-structures -u ALA-ASN-ALA-AGLC --patches CTER -2 NGLA 1 -1 -o ANA-NAGLC

By default, if parameters are missing then the files are not created and the
missing parameters are written to stdout. Using the -v flag will disable
verification:

.. code-block:: bash

   $ # Note that this is fixed in v0.2.3, and will now pass verification
    $ pygen-structures AdP -o AdP
    Missing parameters:
    bonds {('CPD1', 'CC')}
    $ pygen-structures -v AdP -o AdP
    $

A different CHARMM distribution can be loaded using the -t option, with the path
to the folder. pygen-structures ships with the latest CHARMM distribution (July
2019) at the time of writing, with some modifications to correct the
D-amino acid parameters
(`these modifications are highlighted in the toppar README`__). The function
which parses the folder will pick the latest versions of the parameter and
topology files (36 over 27, 36m over 36), so if you plan on using an older
version of the forcefield (this is not recommended) you will have to remove the
newer versions and change the file extensions to match the current conventions
(``.rtf`` for topology files and, ``.prm`` for parameter files).


__ https://github.com/thesketh/pygen-structures/blob/master/pygen_structures/toppar/README.md
