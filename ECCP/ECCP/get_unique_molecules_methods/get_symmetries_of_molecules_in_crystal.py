"""
get_symmetries_of_molecules_in_crystal.py, Geoffrey Weal, 20/3/22

This script will determine which molecules in the crystal are symmetrically identical.
"""
from numpy import eye, matmul, array

from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import PointGroupAnalyzer

from SUMELF import get_cell_corner_points, get_distance

def get_symmetries_of_molecules_in_crystal(molecules, crystal): 
	"""
	This method is designed to determine all the molecules in the crystal that are identical due to the point symmetry of the crystal. 

	Parameters
	----------
	molecules : list of ase.Atoms
		This is a list of all the molecules in your crystal.
	crystal : ase.Atoms
		This is the ase.Atoms object of your crystal.

	Returns
	-------
	symmetric_molecules : list
		A list of the indices of symmetric molecules in the molecules list. 
	"""

	# First, get the names of the molecules 
	mol_names = sorted(molecules.keys())

	# Second, get details of the molecule.
	molecules_symbols   = {mol_name: molecule.get_chemical_symbols() for mol_name, molecule in molecules.items()}
	molecules_positions = {mol_name: molecule.get_positions()        for mol_name, molecule in molecules.items()}
	
	# Third, get details of the crystal. 
	crystal_cell_lattice = crystal.get_cell()
	cell_points = get_cell_corner_points(crystal_cell_lattice,super_cell_reach=1,bottom_of_range=0)
	get_middle_of_cell = (sum(cell_points)/float(len(cell_points)))

	# Fourth, get all the symmetry operators that apply to this crystal.
	#         * The crystal has been wrapped to make sure that PointGroupAnalyzer is 
	#           able to analyse the crystal within the unit cell
	print('Getting point group symmetries within the crystal')
	wrapped_crystal = crystal.copy(); wrapped_crystal.wrap()
	mol = AseAtomsAdaptor.get_molecule(wrapped_crystal)
	pga = PointGroupAnalyzer(mol)
	all_symmetry_operations = pga.get_symmetry_operations()

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Fifth, determine all the symmetrically identical molecules in the crystal.
	print('Checking point group symmetries of molecules within the crystal')
	symmetric_molecule_pairs = {}
	for symmetry_operation in all_symmetry_operations:

		# 5.1: Obtain the rotation matrix for this symmetry operation.
		rotation_matrix = symmetry_operation.rotation_matrix

		# 5.2: Check that the rotation matrix is a 3x3 matrix. 
		if rotation_matrix.shape != (3, 3):
			to_string  = f'Error: Rotation Matrix is not a 3x3 matrix.\n'
			to_string += f'Shape of Rotation Matrix: {rotation_matrix.shape}\n'
			to_string += f'Rotation Matrix (see below):\n'
			to_string += f'{rotation_matrix.shape}\n'
			to_string += 'Check this'
			raise Exception(to_string)

		# 5.3: If the rotation matrix is the identity, continue to the next symmetry operation. 
		if (rotation_matrix == eye(rotation_matrix.shape[0])).all():
			continue

		# 5.4: Obtain the translation matrix for this symmetry operation.
		translation_vector = symmetry_operation.translation_vector

		# CHECK: I want to check if translation_vector should be being used. 
		if not (translation_vector == array([0., 0., 0.])).all():
			import pdb; pdb.set_trace()
			raise Exception('Check symmetry translation. Included or not.')

		# 5.5: Compare all the molecule from the moelcules list together and see if any 
		#      are related to each other by the current symmetry operator (symmetry_operation). 
		for index1 in range(len(mol_names)):

			# 5.5.1: Get the name for the first molecule. 
			mol_name1               = mol_names[index1]

			# 5.5.2: Get mol_name1 molecule elements, and translate it's position.
			#        * Note 1: You want to move the origin to the centre of the cell, as point group 
			#                  symmetry operations are performed about the middle of the cell.
			#        * Note 2: You then add get_middle_of_cell to move the molecule back to where it 
			#                  would have been if rotation_matrix was the identity matrix. 
			molecule1_symbols       = molecules_symbols[mol_name1]
			molecule1_positions_sym = matmul(rotation_matrix,(molecules_positions[mol_name1] - get_middle_of_cell).T).T + get_middle_of_cell

			# 5.5.3: For each other molecule in the mol_names list (based on index1):
			for index2 in range(index1+1,len(mol_names)):

				# 5.5.4: Get the name for the second molecule. 
				mol_name2           = mol_names[index2]

				# 5.5.5: Get mol_name2 molecule elements, and translate it's position
				molecule2_symbols   = molecules_symbols[mol_name2]
				molecule2_positions = molecules_positions[mol_name2]

				# 5.5.6: Determine if the two molecules are symmetruc
				are_molecules_equal, m1_to_m2_indices = are_molecules_structurally_equivalent_due_to_crystal_symmetry(molecule1_symbols, molecule2_symbols, molecule1_positions_sym, molecule2_positions)
				
				# 5.5.7: If the two moelcules are structually equivalent due to the crystal symmetry, add this to symmetric_molecule_pairs
				if are_molecules_equal:
					symmetric_molecule_pairs[(mol_name1,mol_name2)] = m1_to_m2_indices

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Sixth, return information about symmetric molecules based on their positions in the crystal structure.
	return symmetric_molecule_pairs

# ------------------------------------------------------------------------------------------------------------------------

difference_in_distance_tolerance = 0.01
def are_molecules_structurally_equivalent_due_to_crystal_symmetry(molecule1_symbols, molecule2_symbols, molecule1_positions, molecule2_positions):
	"""
	This method is designed to determine if two molecules are pretty much the same.

	This molecule will not consider hydrogens.

	Parameters
	----------
	molecule1_symbols : list of str.
		This is the list of elements in molecule 1.
	molecule2_symbols : list of str.
		This is the list of elements in molecule 2.
	molecule1_positions : numpy.array
		This is the list of positions of atoms in molecule 1.
	molecule2_positions : numpy.array
		This is the list of positions of atoms in molecule 2.

	Attributes
	----------
	difference_in_distance_tolerance : float
		This is the difference in distances between the molecules when placed onto each other based on the symmetry rules of the crystal to be considered equivalent (i.e. same bends and twists and well as structure).

	Returns
	-------
	are_molecules_equal : bool.
		This variable indicates if the two molecule are the same or not. True if they are the same, False if not. 
	m1_to_m2_indices : list of tuple
		This list indices have the indices in molecule 1 relate to indices in molecule 2.
	"""

	# First, if each molecule doesnt have the same types and numbers of atoms, then they are definitely different.
	if not (sorted(molecule1_symbols) == sorted(molecule2_symbols)):
		return False, None

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# Second, obtain the indices for each element in each molecule
	def get_indices_of_symbols(molecule_symbols):
		all_symbols = {}
		for index_symbol in range(len(molecule_symbols)):
			symbol = molecule_symbols[index_symbol]
			all_symbols.setdefault(symbol,[]).append(index_symbol)
		return all_symbols
	indices_of_symbols_of_mol_1 = get_indices_of_symbols(molecule1_symbols)
	indices_of_symbols_of_mol_2 = get_indices_of_symbols(molecule2_symbols)

	# 2.1: Remove the hydrogens from indices_of_symbols_of_mol_1
	if 'H' in indices_of_symbols_of_mol_1:
		indices_of_hydrogens_of_mol_1 = indices_of_symbols_of_mol_1['H']
		del indices_of_symbols_of_mol_1['H']
	else:
		indices_of_hydrogens_of_mol_1 = []

	# 2.2: Remove the hydrogens from indices_of_symbols_of_mol_2
	if 'H' in indices_of_symbols_of_mol_2:
		indices_of_hydrogens_of_mol_2 = indices_of_symbols_of_mol_2['H']
		del indices_of_symbols_of_mol_2['H']
	else:
		indices_of_hydrogens_of_mol_2 = []

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# Third, assign indices in molecule 1 to molecule 2 based on those indices being the same element in the same position.
	mol1_index_to_mol2_index = {}
	for element in indices_of_symbols_of_mol_1.keys():

		# 3.1: Get the indices of elements in molecules 1 and 2.
		indices_of_element_in_mol1 = indices_of_symbols_of_mol_1[element]
		indices_of_element_in_mol2 = indices_of_symbols_of_mol_2[element]
		
		# 3.2: Get the distance between atoms of the same element between molecules (that should be in the same position if they are the same molecule).
		distance_database = {}
		for index1 in indices_of_element_in_mol1:
			molecule1_position = molecule1_positions[index1]
			for index2 in indices_of_element_in_mol2:
				molecule2_position = molecule2_positions[index2]
				distance = get_distance(molecule1_position,molecule2_position)
				distance_database[(index1,index2)] = distance
		
		# 3.3: Assign indices in molecule 1 to indices in molecule 2.
		mol1_index_to_mol2_index_for_this_element = {}
		for (index_mol1, index_mol2), distance in sorted(distance_database.items(),key=lambda x: x[1]):
			if distance > difference_in_distance_tolerance:
				break
			if (index_mol1 in mol1_index_to_mol2_index_for_this_element.keys()) or (index_mol2 in mol1_index_to_mol2_index_for_this_element.values()):
				continue
			mol1_index_to_mol2_index_for_this_element[index_mol1] = index_mol2
		
		# 3.4: Determine if all indices of elements in molecule 1 have been assigned to unique indices in molecule 2.
		#      If not, we know the two molecules are not the same. 
		if not (sorted(indices_of_element_in_mol1) == sorted(mol1_index_to_mol2_index_for_this_element.keys())):
			return False, None
		if not (sorted(indices_of_element_in_mol2) == sorted(mol1_index_to_mol2_index_for_this_element.values())):
			return False, None

		# 3.5: Double checks to make sure that atom indices have not been double countered in indices_of_element_in_mol1 
		#      and indices_of_element_in_mol2.
		#      * This should never happen. if it does, there is a bug.
		if not len(set(indices_of_element_in_mol1)) == len(indices_of_element_in_mol1):
			raise Exception('Error, bug here. check out.')
		if not len(set(indices_of_element_in_mol2)) == len(indices_of_element_in_mol2):
			raise Exception('Error, bug here. check out.')

		# 3.6: Double checks to make sure that atom indices for this element have not already been added to mol1_index_to_mol2_index
		#      * Currently, mol1_index_to_mol2_index should not contain any atom indices for element until we reach step 3.7. 
		#      * This should never happen. if it does, there is a bug.
		if any((index1 in mol1_index_to_mol2_index.keys()) for index1 in mol1_index_to_mol2_index_for_this_element.keys()):
			raise Exception('Error, bug here. check out.')
		if any((index2 in mol1_index_to_mol2_index.values()) for index2 in mol1_index_to_mol2_index_for_this_element.values()):
			raise Exception('Error, bug here. check out.')

		# 3.7: Update mol1_index_to_mol2_index_for_this_element
		mol1_index_to_mol2_index.update(mol1_index_to_mol2_index_for_this_element)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Fourth, at this point all the indices from each molecules should have been assigned. These checks will make sure that
	#         all the atoms in mol1 and mol2 have been assigned/included in mol1_index_to_mol2_index.
	if not sorted(list(mol1_index_to_mol2_index.keys())   + indices_of_hydrogens_of_mol_1) == list(range(len(molecule1_symbols))):
		raise Exception('Error, bug here. check out.')
	if not sorted(list(mol1_index_to_mol2_index.values()) + indices_of_hydrogens_of_mol_2) == list(range(len(molecule2_symbols))):
		raise Exception('Error, bug here. check out.')

	# Fifth, return how indices are related.
	return True, tuple(sorted(mol1_index_to_mol2_index.items()))

# ------------------------------------------------------------------------------------------------------------------------








