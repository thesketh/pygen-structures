# pygen-structures

`pygen-structures` (pigeon structures) is a Python utility which allows for the generation of 3 dimensional molecular structures which can be used in molecular dynamics or Monte Carlo simulations.

Molecules are generated from a list of residues and patches in the format of the CHARMM forcefield, and can be written out as valid PSF and PDB files.

This should make it significantly easier to perform combinatorial searches on particular sequence lengths and linkages, requiring no manual intervention provided the molecules of interest are reasonably small and the residue/patch definitions already exist in the forcefield.

In essence, `pygen-structures` aims (eventually) to be a complete [psfgen](https://www.ks.uiuc.edu/Research/vmd/plugins/psfgen/) replacement with more autonomous functions and without the necessity for an initial 3D structure.

`pygen-structures` depends upon _rdkit_ and _numpy_ and can be installed via pip:

```
pip install pygen-structures
```

Tests can be run with _pytest_, and the only integration test example at present requires `openmm`:

```
pytest --pyargs pygen_structures
```

The package can be used as a command line utility, or as a Python library. Command line usage is simple:

```
pygen-structures --use-charmm-names AGLC-BFRU-AGAL --patches RAFF 0 1 2 --segid RAFF --name Raffinose -o RAFF
```

will produce the structure for the sugar raffinose, saving the PDB file to RAFF.pdb and the PSF file to RAFF.psf. For pure protein structures, the one letter code can be used (the histidine residue type can be set using `--histidine`, but defaults to HSE):

```
pygen-structures YEET -o YEET
```

D-amino acids need only be preceded with a lowercase d. In one letter code mode, first and last patches can also be specified using hyphens as delimiters, and the patch name to be applied in the relevant location:

```
pygen-structures NNEU-AdFK-CTER -o AdFK
```

When output is not specified, the PDB file is written directly to stdout.
