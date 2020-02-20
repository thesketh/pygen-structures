"""
Module containing the Atom class, used to store PDB/PSF data and
represent it directly as a PSF/PDB atom line.
"""

class Atom:
    """
    A class to hold the information present in a PDB/PSF atom.
    Atom objects are directly representable as an atom line from a PDB
    or PSF file.

    record:         "ATOM" or "HETATM" (PDB field)
    atom_serial:    the atom number (PDB and PSF field)
    atom_name:      the name of the atom in the residue
                    (PDB & PSF field)
    alt_locator:    alternate location indicator (PDB field)
    residue_name:   name of the residue the atom is in
                    (PDB & PSF field)
    chain:          identifier for chain (PDB field)
    residue_number: residue sequence number (PDB & PSF field)
    insertion_code: code for insertion of residues. This allows for
                    branches, to a limited degree. e.g. 52 -> 52A -> 53
                    Not widely used (PDB field)
    x:              x position of atom /Angstroms (PDB field)
    y:              y position of atom /Angstroms (PDB field)
    z:              z position of atom /Angstroms (PDB field)
    occupancy:      occupancy of position from crystal structure.
                    Generally used to indicate presence of multiple
                    conformations (PDB field)
    beta_factor:    temperature factor. Indicates how much the position
                    of the atom fluctuates with temperature (PDB field)
    segment_id:     segment identifier. Used by MD simulation programs
                    to keep track of multiple molecules/groups of
                    molecules. (PDB & PSF field)
    element_symbol: chemical element in capitals (PDB field)
    formal_charge:  the formal charge on the atom (PDB field)
    atom_type:      the type of the atom. To be used by the forcefield
                    when calculating bond/NB terms (PSF field)
    partial_charge: the charge of the atom. To be used by the
                    forcefield when calculating electrostatics
                    (PSF field)
    mass:           mass of the atom. To be used by the forcefield to
                    calculated acceleration from given forces
                    (PSF field)
    """

    def __init__(
        self,
        record: str = "ATOM",
        atom_serial: int = 0,
        atom_name: str = "OH2",
        alt_locator: str = "",
        residue_name: str = "TIP3",
        chain: str = "",
        residue_number: int = 0,
        insertion_code: str = "",
        x: float = 0.0,
        y: float = 0.0,
        z: float = 0.0,
        occupancy: float = 1.0,
        beta_factor: float = 0.0,
        segment_id: str = "0",
        element_symbol: (str, None) = None,
        formal_charge: int = 0,
        atom_type: str = "",
        partial_charge: float = 0.0,
        mass: float = 0.0,
    ):
        self.record = record
        self.atom_serial = atom_serial
        self.atom_name = atom_name
        self.alt_locator = alt_locator
        self.residue_name = residue_name
        self.chain = chain
        self.residue_number = residue_number
        self.insertion_code = insertion_code
        self.x = x
        self.y = y
        self.z = z
        self.occupancy = occupancy
        self.beta_factor = beta_factor
        self.segment_id = segment_id
        if element_symbol:
            self.element_symbol = element_symbol
        else:
            self.element_symbol = self.atom_name[0]
        self.formal_charge = formal_charge
        self.atom_type = atom_type
        self.partial_charge = partial_charge
        self.mass = mass

    def to_pdb_line(self) -> str:
        """
        Create a PDB formatted atom record from the Atom. The fields
        are described in the class docstring.

        Where atom serials exceed the 5 character limit, these are
        converted to hexadecimal.

        Where residue numbers exceed the 4 character limit, these are
        wrapped around to 0.
        """
        format_string = (
            "{:<6.6s}{:5s} {:^5.4s}{:.1s}{:<4.4s} {:.1s}{:>4d} {:.1s}   "
            "{: 8.3f}{: 8.3f}{: 8.3f}{: 6.2f}{: 6.2f}      {:<4.4s}{:>2.2s}"
            "{:>2.2s}\n"
        )

        atom_serial_format = "{: 5d}"
        if self.atom_serial >= 100000:
            atom_serial_format = "{: 5X}"

        formal_charge = ""
        if self.formal_charge:
            formal_charge = "{:+2d}".format(self.formal_charge)[::-1]

        atom_string = format_string.format(
            self.record,
            atom_serial_format.format(self.atom_serial),
            self.atom_name,
            self.alt_locator,
            self.residue_name,
            self.chain,
            self.residue_number % 10000,
            self.insertion_code,
            self.x,
            self.y,
            self.z,
            self.occupancy,
            self.beta_factor,
            self.segment_id,
            self.element_symbol,
            formal_charge,
        )
        return atom_string

    def to_psf_line(self) -> str:
        """
        Create a PDB formatted atom record from the Atom. The fields
        are described in the class docstring.

        Where residue numbers exceed the 4 character limit, these are
        wrapped around to 0.
        """
        format_string = (
            " {: 7d} {:<4.4s} {:<4d} {:<4.4s} {:<4.4s} {:<4.4s}  {: 8.6f}"
            "      {:8.4f}           0\n"
        )

        atom_string = format_string.format(
            self.atom_serial,
            self.segment_id,
            self.residue_number % 10000,
            self.residue_name,
            self.atom_name,
            self.atom_type,
            self.partial_charge,
            self.mass,
        )
        return atom_string

    def __str__(self) -> str: return self.to_pdb_line()
