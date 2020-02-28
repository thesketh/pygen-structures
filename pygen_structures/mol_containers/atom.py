"""
Module containing the Atom class, used to store PDB/PSF data and
represent it directly as a PSF/PDB atom line.

This is in a separate, dependency free class to encourage re-use.
"""


class Atom:
    """
    A class to hold the information present in a PDB/PSF atom.
    Atom objects are directly representable as an atom line from a PDB
    or PSF file.

    :param record: "ATOM" or "HETATM" (PDB field)
    :param atom_serial: the atom number (PDB and PSF field)
    :param atom_name: the name of the atom in the residue (PDB & PSF field)
    :param alt_locator: alternate location indicator (PDB field)
    :param residue_name: name of the residue the atom is in (PDB & PSF field)
    :param chain: identifier for chain (PDB field)
    :param residue_number: residue sequence number (PDB & PSF field)
    :param insertion_code: code for insertion of residues. This allows for\
    branches, to a limited degree. e.g. 52 -> 52A -> 53. Not widely used\
    (PDB field)
    :param x: x position of atom /Angstroms (PDB field)
    :param y: y position of atom /Angstroms (PDB field)
    :param z: z position of atom /Angstroms (PDB field)
    :param occupancy: occupancy of position from crystal structure.\
    Generally used to indicate presence of multiple conformations (PDB field)
    :param beta_factor: temperature factor. Indicates how much the position\
    of the atom fluctuates with temperature (PDB field)
    :param segment_id: segment identifier. Used by MD simulation programs\
    to keep track of multiple molecules/groups of molecules. (PDB & PSF field)
    :param element_symbol: chemical element in capitals (PDB field)
    :param formal_charge: the formal charge on the atom (PDB field)
    :param atom_type: the type of the atom. To be used by the forcefield\
    when calculating bond/NB terms (PSF field)
    :param partial_charge: the charge of the atom. To be used by the\
    forcefield when calculating electrostatics (PSF field)
    :param mass: mass of the atom. To be used by the forcefield to\
    calculated acceleration from given force (PSF field)

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
        Create a PDB formatted atom record from the ``Atom``. The fields
        are described in the class docstring.

        Where atom serials exceed the 5 character limit, these are
        converted to hexadecimal.

        Where residue numbers exceed the 4 character limit, these are
        wrapped around to 0.
        """
        format_string = (
            "{:<6.6s}{:5s} {:^4.4s}{:^1.1s}{:<4.4s}{:^1.1s}{:>4d} {:1.1s}  "
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
        Create a PDB formatted atom record from the ``Atom``. The fields
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

    @classmethod
    def from_pdb_line(cls, line: str):
        """
        Instantiate the class from a line in a PDB file.
        This will not have the atom_type, partial_charge
        or mass.

        :param line: an ATOM/HETATM record from a pdb file.

        """
        record = line[:6].strip()
        try:
            atom_serial = int(line[6:11])
        except ValueError:
            atom_serial = int(line[6:11], 16)
        atom_name = line[12:16].strip()
        alt_locator = line[16]
        residue_name = line[17:21].strip()
        chain = line[21]
        residue_number = int(line[22:26])
        insertion_code = line[27]
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        occupancy = float(line[54:60])
        beta_factor = float(line[60:66])
        segment_id = line[72:76].strip()
        element_symbol = line[76:78].strip()
        charge = 0
        try:
            symbol = line[79]
            if symbol == "-":
                charge = 0 - int(line[78])
            else:
                try:
                    charge = int(line[78])
                except ValueError:
                    pass
        except IndexError:
            pass

        return cls(
            record,
            atom_serial,
            atom_name,
            alt_locator,
            residue_name,
            chain,
            residue_number,
            insertion_code,
            x,
            y,
            z,
            occupancy,
            beta_factor,
            segment_id,
            element_symbol,
            charge
         )
