# CHARMM Distribution

This is the July 2019 CHARMM distribution, with all outdated (pre CHARMM35/36) and non-CHARMM format files removed.

## Changes from the Standard Distribution

The D-amino acid parameter files have been updated to address the following:

  - Missing parameters, primarily with D-proline. These were taken from the L-amino parameters.
  - First patch for the D-amino acids. The normal protein N-terminal patches were changing the atom type of the first acid back to CT1. This leads to the wrong CMAP terms being used where there are D-amino acids in the first position.
  - Incorrect IC tables for D-threonine (the existing table referred to D-allothreonine) and D-isoleucine (the existing table referred to D-alloisoleucine). The improper IC for CB was inverted.

These changes have not yet been reflected in the CHARMM distribution.

## Updates

Updated files can be downloaded from the [MacKerell lab website](http://mackerell.umaryland.edu/charmm_ff.shtml) when they are made available.

## Citations

These files are [in the public domain](https://www.charmm.org/ubbthreads/ubbthreads.php?ubb=showflat&Number=33804), but the normal academic conventions apply: please make sure to cite the relevant CHARMM papers if you use the forcefields defined in these files in your work. This helps to ensure further funding.

Any molecules produced by `pygen-structures` will note the topology files which were used in the PSF and PDB headers. The reference for the paper corresponding to each set of topology files can be found on the MacKerell website, as noted above.
