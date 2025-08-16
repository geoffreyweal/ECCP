"""
centre_molecules.py, Geoffrey Weal, 20/4/22

This script will centre the molecules as close to the origin unit cell.
"""

from SUMELF import centre_molecule_in_cell

def centre_molecules(molecules):
	"""
	This method will centre the molecules as close to the origin unit cell.

	Parameters
	----------
	molecules : list of ase.Atoms
		These are the molecules in the crystal
	"""
	for molecule_name, molecule in molecules.items():
		centre_molecule_in_cell(molecule, molecule.get_cell())