from pygen_structures.charmm_containers import *

ALANINE_BLOCK = """RESI ALA          0.00
GROUP
ATOM N    NH1    -0.47  !     |
ATOM HN   H       0.31  !  HN-N
ATOM CA   CT1     0.07  !     |     HB1
ATOM HA   HB1     0.09  !     |    /
GROUP                   !  HA-CA--CB-HB2
ATOM CB   CT3    -0.27  !     |    \\
ATOM HB1  HA3     0.09  !     |     HB3
ATOM HB2  HA3     0.09  !   O=C
ATOM HB3  HA3     0.09  !     |
GROUP                   !
ATOM C    C       0.51
ATOM O    O      -0.51
BOND CB CA  N  HN  N  CA
BOND C  CA  C  +N  CA HA  CB HB1  CB HB2  CB HB3
DOUBLE O  C
IMPR N -C CA HN  C CA +N O
CMAP -C  N  CA  C   N  CA  C  +N
DONOR HN N
ACCEPTOR O C
IC -C   CA   *N   HN    1.3551 126.4900  180.0000 115.4200  0.9996
IC -C   N    CA   C     1.3551 126.4900  180.0000 114.4400  1.5390
IC N    CA   C    +N    1.4592 114.4400  180.0000 116.8400  1.3558
IC +N   CA   *C   O     1.3558 116.8400  180.0000 122.5200  1.2297
IC CA   C    +N   +CA   1.5390 116.8400  180.0000 126.7700  1.4613
IC N    C    *CA  CB    1.4592 114.4400  123.2300 111.0900  1.5461
IC N    C    *CA  HA    1.4592 114.4400 -120.4500 106.3900  1.0840
IC C    CA   CB   HB1   1.5390 111.0900  177.2500 109.6000  1.1109
IC HB1  CA   *CB  HB2   1.1109 109.6000  119.1300 111.0500  1.1119
IC HB1  CA   *CB  HB3   1.1109 109.6000 -119.5800 111.6100  1.1114
PATCHING FIRST NNEU LAST CTER
"""

NNEU_BLOCK = """PRES NNEU         0.00 ! neutral N-terminus; charges from LSN
GROUP                  ! use in generate statement
ATOM N    NH2    -0.96 !
ATOM HT1  H       0.34 !         HT1
ATOM HT2  H       0.34 !        /
                       ! --CA--N--HT2
ATOM CA   CT1     0.19 !   |    ! change to CT2 for neutral N terminal glycine
ATOM HA   HB1     0.09 !  HA    ! change to HA1 and HB2 and add HA2 atom for N terminal glycine
DELETE ATOM HN
BOND HT1 N HT2 N
DONOR HT1 N
DONOR HT2 N
IC HT1  N    CA   C     0.0000  0.0000  180.0000  0.0000  0.0000
IC HT2  CA   *N   HT1   0.0000  0.0000  120.0000  0.0000  0.0000
"""

def test_residue_definition_parser():
    """This tests the RESI block parser."""
    res_def = CHARMMResidueDefinition.from_block(ALANINE_BLOCK)
    assert res_def.name == "ALA"
    assert res_def.first == "NNEU"
    assert res_def.last == "CTER"

    atoms = [
        ("N", "NH1", -0.47),
        ("HN", "H", 0.31),
        ("CA", "CT1", 0.07),
        ("HA", "HB1", 0.09),
        ("CB", "CT3", -0.27),
        ("HB1", "HA3", 0.09),
        ("HB2", "HA3", 0.09),
        ("HB3", "HA3", 0.09),
        ("C", "C", 0.51),
        ("O", "O", -0.51)
    ]
    assert res_def.atoms == atoms

    bonds = [
        ("CB", "CA"), ("N", "HN"), ("N", "CA"), ("C", "CA"), ("C", "+N"),
        ("CA", "HA"), ("CB", "HB1"), ("CB", "HB2"), ("CB", "HB3"),
        ("O", "C"), ("O", "C")
    ]
    assert res_def.bonds == bonds

    impropers = [
        ("N", "-C", "CA", "HN"), ("C", "CA", "+N", "O")
    ]
    assert res_def.impropers == impropers

    cross_maps = [
        ("-C", "N", "CA", "C", "N", "CA", "C", "+N")
    ]
    assert res_def.cross_maps == cross_maps

def test_residue_definition_smarts():
    """
    This tests that the IC table is loaded correctly. Otherwise,
    the chirality may be wrong.
    """
    res_def = CHARMMResidueDefinition.from_block(ALANINE_BLOCK)
    assert res_def.to_smarts() == "[#7]-[#6@@H](-[#6])-[#6]=[#8]"

def test_residue_creation():
    """
    Tests that residue definitions can be converted to residues
    correctly.
    """
    res_def = CHARMMResidueDefinition.from_block(ALANINE_BLOCK)
    residue = res_def.to_residue(index=0)
    assert residue.atoms == res_def.atoms

    bonds = [
        ((0, "CB"), (0, "CA")), ((0, "N"), (0, "HN")),
        ((0, "N"), (0, "CA")), ((0, "C"), (0, "CA")),
        ((0, "C"), (1, "N")), ((0, "CA"), (0, "HA")),
        ((0, "CB"), (0, "HB1")), ((0, "CB"), (0, "HB2")),
        ((0, "CB"), (0, "HB3")),
        ((0, "O"), (0, "C")), ((0, "O"), (0, "C"))
    ]
    assert residue.bonds == bonds

    impropers = [
        ((0, "N"), (-1, "C"), (0, "CA"), (0, "HN")),
        ((0, "C"), (0, "CA"), (1, "N"), (0, "O"))
    ]
    assert residue.impropers == impropers

    cross_maps = [
        (
            (-1, "C"),
            (0, "N"),
            (0, "CA"),
            (0, "C"),
            (0, "N"),
            (0, "CA"),
            (0, "C"),
            (1, "N")
        )
    ]
    assert residue.cross_maps == cross_maps

def test_patch_residue_definition_parser():
    """Tests the PRES block parser."""
    pres_def = CHARMMPatchResidueDefinition.from_block(NNEU_BLOCK)
    assert pres_def.n_residues == 1

    atoms = [
        ("N", "NH2", -0.96),
        ("HT1", "H", 0.34),
        ("HT2", "H", 0.34),
        ("CA", "CT1", 0.19),
        ("HA", "HB1", 0.09)
    ]
    assert pres_def.atoms[0] == atoms

    deletions = {(0, "HN"), }
    assert pres_def.deletions == deletions

    bonds = [
        ((0, "HT1"), (0, "N")), ((0, "HT2"), (0, "N"))
    ]
    assert pres_def.bonds[0] == bonds

    first_ic = [
        (0, "HT1"), (0, "N"), (0, "CA"), (0, "C"), 0.0, 0.0, 180.0, 0.0, 0.0
    ]
    assert pres_def.ics[0][0] == first_ic

def test_patch_apply():
    nneu_patch = CHARMMPatchResidueDefinition.from_block(NNEU_BLOCK)
    ala_res_def = CHARMMResidueDefinition.from_block(ALANINE_BLOCK)
    ala_res = ala_res_def.to_residue(index=0)
    nneu_patch.apply(ala_res)

    HN_exists = False
    for (atom_name, atom_type, partial_charge) in ala_res.atoms:
        if atom_name == "HN":
            HN_exists = True
            break
    assert HN_exists == False

    assert ((0, "HT1"), (0, "N")) in ala_res.bonds
    assert ((0, "HT2"), (0, "N")) in ala_res.bonds
