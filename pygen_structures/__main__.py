#!/usr/bin/env python3
from argparse import ArgumentParser
from pygen_structures import load_charmm_dir, sequence_to_mol, code_to_mol
from pygen_structures.version import version, year

def main():
    parser = ArgumentParser(
        description="Generate PSF and PDB files for a given sequence.",
        epilog=(
            "v{}, (C) Travis Hesketh {}. Available under the terms of "
            "the 3-clause BSD license. "
            "https://github.com/thesketh/pygen-structures"
        ).format(version, year)
    )
    parser.add_argument(
        "sequence",
        help=("A sequence of amino acid codes (the default) or of CHARMM"
              "residue names (- delimited)"),
        type=str
    )
    parser.add_argument(
        "-n", "--name",
        help="Name of the molecule, for the PDB/PSF file",
        type=str,
        default=None
    )
    parser.add_argument(
        "-s", "--segid",
        help="Segment ID used in the PDB/PSF file",
        type=str,
        default="PROT"
    )
    parser.add_argument(
        "-p", "--patches",
        help=("Patches to be applied to the molecule. The name of the patch,"
              "followed by the indices of the residues the patch is to be "
              "applied to (or FIRST/LAST)"),
        nargs="+",
    )
    parser.add_argument(
        "-t", "--toppar",
        help="Path to a toppar folder.",
        type=str,
        default=None
    )
    parser.add_argument(
        "-v", "--verify",
        help="Verify that the params exist. Default: True.",
        action="store_false"
    )
    parser.add_argument(
        "-o", "--output",
        help="The output prefix for the PSF/PDB files.",
    )
    parser.add_argument(
        "--histidine",
        help="The histidine residue to use for 'H' [HSE]",
        default="HSE"
    )
    parser.add_argument(
        "--use-charmm-names", "-u",
        help=(
            "Use actual residue names instead of the amino acid code. "
            "(- delimited). Default: False."
        ),
        action="store_true"
    )
    args = parser.parse_args()

    if args.toppar:
        rtf, prm = load_charmm_dir(args.toppar)
    else:
        rtf, prm = load_charmm_dir()

    patch_list, patches = args.patches, {}
    if patch_list:
        while patch_list:
            patch_name = patch_list.pop(0)
            patch = rtf.patches[patch_name]
            residues = []
            for _ in range(patch.n_residues):
                residue_index = patch_list.pop(0)
                if residue_index not in ("FIRST", "LAST"):
                    residue_index = int(residue_index)
                residues.append(residue_index)
            patches[patch_name] = residues
    else:
        patches = None

    name, segid = args.name, args.segid
    if args.use_charmm_names:
        sequence = args.sequence.split('-')
        mol = sequence_to_mol(sequence, rtf, patches, name, segid)
    else:
        histidine = args.histidine
        sequence = args.sequence
        mol = code_to_mol(sequence, rtf, patches, name, segid, histidine)

    if args.verify:
        missing_params = prm.get_unmatched(mol)
        missing = {}
        for key, value in missing_params.items():
            if value:
                missing[key] = value
        if missing:
            print("Missing parameters:")
            for key, value in missing.items():
                print(key, value)
            exit()

    output_prefix = args.output
    if output_prefix:
        mol.to_pdb_file('{}.pdb'.format(output_prefix))
        mol.to_psf_file('{}.psf'.format(output_prefix))
    else:
        print(mol.to_pdb_block())

if __name__ == "__main__":
    main()
