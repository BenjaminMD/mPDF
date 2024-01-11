"""
Author: Benjamin Fahl

ASE CIF parser wrapper
"""
from ase.io import read  # type: ignore # pylint: disable=import-error


class LatticeParser:
    """
    quick parser wrap for CIF files
    """

    def __init__(self, cif_path):
        self.atoms = read(cif_path)
        self._extract_lattice_info()

    def _extract_lattice_info(self):
        cell = self.atoms.cell
        lattice_parameters = cell.lengths()
        lattice_angles = cell.angles()
        self.a, self.b, self.c = lattice_parameters
        self.alpha, self.beta, self.gamma = lattice_angles

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        params = [
            f"Formula: {self.atoms.get_chemical_formula()}",
            "Lattice parameters:\n",
            f"\ta: {self.a}",
            f"\tb: {self.b}",
            f"\tc: {self.c}",
            "Lattice angles:\n",
            f"\talpha: {self.alpha}",
            f"\tbeta: {self.beta}",
            f"\tgamma: {self.gamma}",
        ]
        return "\n".join(params)
