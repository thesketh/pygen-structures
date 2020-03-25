.. pygen-structures documentation master file, created by
   sphinx-quickstart on Thu Feb 27 13:59:41 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pygen-structures Documentation
==============================

|rtd| |travis| |codecov| |joss|

.. |rtd| image:: https://readthedocs.org/projects/pygen-structures/badge/?version=latest
    :target: https://pygen-structures.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. |travis| image:: https://travis-ci.org/thesketh/pygen-structures.svg?branch=master
    :target: https://travis-ci.org/thesketh/pygen-structures
    :alt: Build Status

.. |codecov| image:: https://codecov.io/gh/thesketh/pygen-structures/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/thesketh/pygen-structures
    :alt: Code Coverage

.. |joss| image:: https://joss.theoj.org/papers/57dce0c14dd34c111f89077cd367dbd7/status.svg
    :target: https://joss.theoj.org/papers/57dce0c14dd34c111f89077cd367dbd7
    :alt: JOSS Paper Status

``pygen-structures`` (pigeon structures) is a Python utility which allows for
the generation of 3 dimensional molecular structures which can be used in
molecular dynamics or Monte Carlo simulations. Molecules are generated from a
list of residues and patches in the format of the CHARMM forcefield, and can be
written out as valid PSF and PDB files. The package can be used as a command
line utility, or as a Python library.

pygen-structures can be installed using pip
(``pip install pygen-structures``), but relies upon `RDKit`__, which is
not pip-installable. `RDKit can be installed in many ways`__, but the easiest
way is to use the `conda package manager`__. `numpy`__ is an additional
dependency. For full installation instructions, see :doc:`installation`. Python
3.6 and 3.7 are supported. To run the tests, `pytest`__ and `OpenMM`__ are
required.

In essence, pygen-structures aims (eventually) to be a complete `psfgen`__
replacement with more autonomous functionality. At present, structures for small
molecules can be generated. This should make it significantly easier to perform
combinatorial searches on particular sequence lengths and linkages. This
requires no manual intervention provided the molecules of interest are
reasonably small (small enough that embedding coordinates is possible - fewer
than 15 residues - and that the secondary structure is not vitally important)
and the residue/patch definitions already exist in the forcefield.

**A quick example**

.. code-block:: bash

    $ pygen-structures -u AGLC-BFRU-AGAL -p RAFF 0 1 2 -n Raffinose -s RAFF -o RAFF
    $ head RAFF.pdb
    COMPND    Raffinose
    AUTHOR    pygen-structures v0.2.3
    REMARK  42
    REMARK  42 TOPOLOGY FILES USED
    REMARK  42     top_all36_carb.rtf
    ATOM      1  C1  AGLC    1      -1.674  -0.714  -0.235  1.00  0.00      RAFF C
    ATOM      2  H1  AGLC    1      -1.556   0.136   0.443  1.00  0.00      RAFF H
    ATOM      3  O1  AGLC    1      -2.855  -0.603  -0.910  1.00  0.00      RAFF O
    ATOM      4  C5  AGLC    1       0.623  -0.935  -0.607  1.00  0.00      RAFF C
    ATOM      5  H5  AGLC    1       1.188  -1.649  -1.276  1.00  0.00      RAFF H


__ https://www.rdkit.org/
__ https://www.rdkit.org/docs/Install.html
__ https://docs.conda.io/projects/conda/en/latest/
__ https://numpy.org/
__ https://docs.pytest.org/en/latest/
__ http://openmm.org/
__ https://www.ks.uiuc.edu/Research/vmd/plugins/psfgen/

.. toctree::
    :maxdepth: 3
    :caption: Contents:

    Home <self>
    Installation <installation>
    command_line_usage
    library_usage
    source/pygen_structures
