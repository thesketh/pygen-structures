# pygen-structures

[![Documentation Status](https://readthedocs.org/projects/pygen-structures/badge/?version=latest)](https://pygen-structures.readthedocs.io/en/latest/?badge=latest)[![Build Status](https://travis-ci.org/thesketh/pygen-structures.svg?branch=master)](https://travis-ci.org/thesketh/pygen-structures)[![codecov](https://codecov.io/gh/thesketh/pygen-structures/branch/master/graph/badge.svg)](https://codecov.io/gh/thesketh/pygen-structures)

`pygen-structures` (pigeon structures) is a Python utility which allows for the generation of 3 dimensional molecular structures which can be used in molecular dynamics or Monte Carlo simulations. Molecules are generated from a list of residues and patches in the format of the CHARMM forcefield, and can be written out as valid PSF and PDB files. The package can be used as a command line utility, or as a Python library.

`pygen-structures` can be installed using pip (`pip install pygen-structures`), but relies upon _RDKit_, which must be installed using the [conda package manager](https://docs.conda.io/projects/conda/en/latest/). For installation instructions, see the 'Installation' section of the readme. Python 3.6 and 3.7 are supported.

In essence, `pygen-structures` aims (eventually) to be a complete [psfgen](https://www.ks.uiuc.edu/Research/vmd/plugins/psfgen/) replacement with more autonomous functionality and without the necessity for an initial 3D structure. This should make it significantly easier to perform combinatorial searches on particular sequence lengths and linkages, requiring no manual intervention provided the molecules of interest are reasonably small and the residue/patch definitions already exist in the forcefield.

## Installation

1. [Install the `conda` package manager](https://docs.anaconda.com/anaconda/install/).
2. Set up a `conda` environment with the relevant dependencies (or install them in your base distribuion). This can be done with the following command: `conda create -n pygen-structures -c rdkit -c omnia 'python>=3.6' 'rdkit>=2018.3' numpy 'openmm>=7.4' pytest`.
3. Activate the `conda` environment: `conda activate pygen-structures`
4. Use `pip` to install `pygen-structures` in this environment: `pip install pygen-structures`.
5. Installation complete! You will have to activate this environment using `conda activate pygen-structures` each time you want to use it.
6. Test the installation using `pytest --pyargs pygen_structures`

## Command line usage

Command line usage for peptides is simple, and takes the following form:

```
pygen-structures SEQUENCE -o OUTPUT_PREFIX
```

Sequences are specified using the one letter protein code, and terminal patches can be supplied by using hyphens as delimiters (e.g. `NNEU-AFK-CT2`, note that both termini must be supplied). D-amino acids need only be preceded by a lowercase 'd'.

OUTPUT_PREFIX.psf and OUTPUT_PREFIX.pdb are created. If `-o` is not specified, the PDB file is written to stdout and no PSF file is generated.

The histidine form used can be set using `--histidine`, and defaults to HSE.

Patches can be supplied using `--patches`, the name of the patch, and the 0-based indices the patch is to be applied to (or 'FIRST'/'LAST').

To generate more complex structures, such as sugars, the residue names should be supplied (hyphen delimited) and the `-u`/`--use-charmm-names` option selected. Some example usage is given below.

`--name` and `--segid` control the names given in the COMPND record and the segment id respectively.

### Examples

To produce a simple peptide sequence, the one letter code can be used. To produce the peptide HIS-GLU-TYR, creating `HEY.psf` and `HEY.pdb`:

```
pygen-structures HEY -o HEY
```

Supposing we think histidine should be protonated, we can change the protonation state of histidine by specifying a different histidine form:

```
pygen-structures HEY -o HEY --histidine HSP
```

Or we could use the three letter codes:

```
pygen-structures -u HSP-GLU-TYR -o HEY
```

For the trisaccharide raffinose, we must use the residue codes. The default segid is PROT, so we can specify a more specific segid using the `--segid` flag, and set the name given in the `COMPND` header using `--name`. The the following command produces `RAFF.psf` and `RAFF.pdb`:

```
pygen-structures -u AGLC-BFRU-AGAL --patches RAFF 0 1 2 --segid RAFF --name Raffinose -o RAFF
```

We can also make glycopeptides. To link alpha-glucose to an arginine residue (in this case, from an ALA-ASN-ALA peptide), we can use the NGLA patch. Note that because the protein residue is not the last in the chain, we have to apply the C-terminus patch manually.

```
pygen-structures -u ALA-ASN-ALA-AGLC --patches CTER -2 NGLA 1 -1 -o ANA-NAGLC
```

By default, if parameters are missing then the files are not created and the missing parameters are written to stdout. Using the `-v` flag will disable verification.

```
$ pygen-structures AdP -o AdP
Missing parameters:
bonds {('CPD1', 'CC')}
$ pygen-structures -v AdP -o AdP
$
```

A different CHARMM distribution can be loaded using the `-t` option, with the path to the folder. `pygen-structures` ships with the latest CHARMM distribution (July 2019) at the time of writing. The function which parses the folder will pick the latest versions of the parameter and topology files (36 over 27, 36m over 36), so if you plan on using an older version of the forcefield (this is not recommended) you will have to remove the newer versions and change the extensions to match the current conventions (.rtf, .prm).

## Library usage

Information about classes and functions for usage as a library can be found on [the project's ReadtheDocs page](https://pygen-structures.readthedocs.io/en/latest/).
