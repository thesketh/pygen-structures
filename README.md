# pygen-structures

`pygen-structures` (pigeon structures) is a Python utility which allows for the generation of 3 dimensional molecular structures which can be used in molecular dynamics or Monte Carlo simulations. Molecules are generated from a list of residues and patches in the format of the CHARMM forcefield, and can be written out as valid PSF and PDB files.

The package can be installed using pip (`pip install pygen-structures`), but relies upon _RDKit_, which must be installed using the [conda package manager](https://docs.conda.io/projects/conda/en/latest/) (`conda install -c rdkit rdkit`). _numpy_ is an additional dependency, and both _pytest_ and _OpenMM_ are required to run all the tests.

In essence, `pygen-structures` aims (eventually) to be a complete [psfgen](https://www.ks.uiuc.edu/Research/vmd/plugins/psfgen/) replacement with more autonomous functionality and without the necessity for an initial 3D structure. This should make it significantly easier to perform combinatorial searches on particular sequence lengths and linkages, requiring no manual intervention provided the molecules of interest are reasonably small and the residue/patch definitions already exist in the forcefield.

The package can be used as a command line utility, or as a Python library.

## Command line usage

Command line usage for peptides is incredibly simple, and takes the following form:

```
pygen-structures SEQUENCE -o OUTPUT_PREFIX
```

Sequences are specified using the one letter protein code, and terminal patches can be supplied by using hyphens as delimiters (e.g. `NNEU-AFK-CT2`, note that both termini must be supplied). OUTPUT_PREFIX.psf and OUTPUT_PREFIX.pdb are created. If the output prefix is not specified, the PDB file is written to stdout and no PSF file is generated. D-amino acids need only be preceded by a lowercase d.

To generate more complex structures, such as sugars, the residue names should be supplied (hyphen delimited) and the `-u`/`--use-charmm-names` option selected. `--name` and `--segid` control the names given in the COMPND record and the segment id respectively. For raffinose:

```
pygen-structures --use-charmm-names AGLC-BFRU-AGAL --patches RAFF 0 1 2 --segid RAFF --name Raffinose -o RAFF
```

## Library usage

Information about classes and functions for usage as a library can be found on [the project's ReadtheDocs page](https://pygen-structures.readthedocs.io/en/latest/). At present this is limited to an API reference, but there are plans for expansion.

## Running tests

To run the tests using pytest, run the following command in the shell: `pytest --pyargs pygen_structures`. `python -m pytest --pyargs pygen_structures` should also work.
When output is not specified, the PDB file is written directly to stdout.
