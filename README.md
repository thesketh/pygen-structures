# pygen-structures

pygen-structures (pigeon structures) is a Python utility which allows for the generation of 3 dimensional molecular structures which can be used in molecular dynamics or Monte Carlo simulations.

Molecules are generated from a list of residues and patches in the format of the CHARMM forcefield, and can be written out as valid PSF and PDB files.

This should make it significantly easier to perform combinatorial searches on particular sequence lengths and linkages, requiring no manual intervention provided the molecules of interest are reasonably small and the residue/patch definitions already exist in the forcefield.

In essence, pygen-structures aims (eventually) to be a complete [psfgen](https://www.ks.uiuc.edu/Research/vmd/plugins/psfgen/) replacement with more autonomous functions and without the necessity for an initial 3D structure.

