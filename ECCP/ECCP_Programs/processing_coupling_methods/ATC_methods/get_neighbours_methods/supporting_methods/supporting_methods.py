"""
supporting_methods.py, Geoffrey Weal, 5/4/22

This script includes methods for making neighbouring pairs of molecules. These methods are used to compare measure to determine if they are exactly the same neighbouring pair or not spatially.
"""

from SUMELF import get_distance, get_cell_corner_points

def make_neighbouring_pair(molecules, index1, index2, displacement):
	"""
	This method is designed to create the neighbouring pairs from the molecules list and the given displacement.

	Parameters
	----------
	molecules : list of ase.Atoms objects
		These are all the individual molecules identified in the crystal that you want to determine neighbouring pairs for.
	index1 : int
		This is the index of the first molecule in the molecules list to take. 
	index2 : int
		This is the index of the second molecule in the molecules list to take. 
	displacement : numpy.array
		This is the displacement to move the second molecule by.

	Returns
	-------
	neighbouring_pair_ase : ase.Atoms
		The neighbouring pair of molecules
	"""

	# First, make a copy of the molecules to be used in the neighbouring pair. 
	molecule1 = molecules[index1].copy()
	molecule2 = molecules[index2].copy()

	# Second, move the second molecule by the displacement.
	molecule2.set_positions(molecule2.get_positions() + displacement)

	# Third, create the neighbouring pair of molecules.
	neighbouring_pair_ase = molecule1 + molecule2
	
	# Fourth, return the neighbouring pair of molecules.
	return neighbouring_pair_ase

# --------------------------------------------------------------------------------------------------------

def is_molecule_already_recorded(new_molecule, new_molecules, crystal_cell_lattice, super_cell_reach=1):
	"""
	This method will determine if you have already recorded this molecule in the crystal. 

	Parameters
	----------
	new_molecule : ase.Atoms
		This is the new molecule that has just been obtained from a crystal symmetry operation.
	new_molecules : list of ase.Atoms
		This is a list of all the unique ase.Atoms objects in the crystal thathave currently been recorded.
	super_cell_reach : int
		This is the cell point coverage you want to check the new_molecule across.

	Returns
	-------
	True if new_molecule is already in new_molecules, False if not.
	"""

	# First, get all the displacements that surround the cell around the molecule.
	cell_points = get_cell_corner_points(crystal_cell_lattice, super_cell_reach=super_cell_reach)

	# Second, look through all the obtain molecules to see if any are tehr same as new_molecule
	for obtained_molecule in new_molecules:

		# 2.1: If new_molecule does not have the same chemical formula as obtained_molecule, move on
		if not (sorted(new_molecule.get_chemical_symbols()) == sorted(obtained_molecule.get_chemical_symbols())):
			continue

		# 2.2: The obtained_molecule can be in any translated position given by a cell_point in cell_points. Need to try them all out.
		for cell_point in cell_points:

			# 2.2.1: Copy the obtained_molecule
			obtained_molecule_copy = obtained_molecule.copy()

			# 2.2.2: Place obtained_molecule in the translated position
			obtained_molecule_copy.set_positions(obtained_molecule_copy.get_positions() + cell_point)

			# 2.2.3: Are new_molecule and obtained_molecule_copy the same molecule in the same place.
			if same_molecules(new_molecule, obtained_molecule_copy):
				# Found a molecule in new_molecules that is the same as new_molecule, so return True
				return True

	# Third, if new_molecule is not the same as any obtained_molecule in new_molecules, even when translated by any of the cell_point.
	return False

max_diff_distance = 0.0001
def same_molecules(new_molecule, obtained_molecule):
	"""
	This method will determine if these two molecules are the same molecules, particularly spatial so that two molecules are not on top of each other.

	Parameters
	----------
	new_molecule : ase.atoms
		This is the new molecule that has just been obtained from a crystal symmetry operation.
	obtained_molecule : ase.Atoms
		This is one of the already obtained molecules that is currently in the new_molecules list

	Returns
	-------
	True if the molecules are the same, False if not.
	"""

	# First, get the elements of the molecules
	new_molecule_elements      = new_molecule.get_chemical_symbols()
	obtained_molecule_elements = obtained_molecule.get_chemical_symbols()

	# Second, get the positions of the molecules
	new_molecule_positions      = new_molecule.get_positions().tolist()
	obtained_molecule_positions = obtained_molecule.get_positions().tolist()

	# Third, check that the two molecules do not overlap in elements and positions.
	for new_element, new_position in zip(new_molecule_elements, new_molecule_positions):
		for index in range(len(obtained_molecule_elements)):
			obtained_element  = obtained_molecule_elements[index]
			obtained_position = obtained_molecule_positions[index]
			if (obtained_element == new_element) and (get_distance(obtained_position, new_position) < max_diff_distance):
				del obtained_molecule_elements[index]
				del obtained_molecule_positions[index]
				break
		else:
			# If any one atom in obtained_molecule does not fit with the new_molecule, then the two molecules are different
			return False

	# Fourth, if you got to hear, you were able to match up every atom in new_molecule with every atom in obtained_molecule. 
	# Therefore, new_molecule and obtained_molecule are teh same.
	return True

# --------------------------------------------------------------------------------------------------------

