"""
invariance_method.py, Geoffrey Weal, 23/2/22

This script is designed to use the procrustes analysis to determine if molecules are rotationally, translationally, and reflectively invariant.
"""

from ECCP.ECCP.get_unique_molecules_methods.invariance_methods.comprehensive_invariance_method               import remove_equivalent_molecules_comprehensive_invariance_method
from ECCP.ECCP.get_unique_molecules_methods.invariance_methods.minimal_elemental_abundance_invariance_method import remove_equivalent_molecules_minimal_elemental_abundance_invariance_method
from ECCP.ECCP.get_unique_molecules_methods.invariance_methods.combination_invariance_method                 import remove_equivalent_molecules_combination_invariance_method

def determine_equivalent_molecules_with_invariance_method(invariance_method_type, unique_molecule_names, molecules, molecule_graphs, neighbouring_molecules_about_molecules={}, max_distance_disparity=0.01, include_hydrogens_in_uniqueness_analysis=False, no_of_cpus=1):
	"""
	This method will determine which invariance method will be used.

	Parameters
	----------
	invariance_method_type : str.
		This is the type of invariance method you would like to use. Options given in the instruction manual on the Github page.
	unique_molecule_names : list of int
		This is a list of names of molecules to compare with each other to determine if they are invariant or not.
	molecules : list of ase.Atoms
		This is the list of molecules in the crystal
	molecule_graphs : list of networkx.Graph
		This is the list of the undirected graph representations of each molecule.
	neighbouring_molecules_about_molecules : dict.
		This is the information about the molecules that surround (in the vicinity of) each moleule in the crystal.
	max_distance_disparity : float
		This is the maximum disparity between two molecules to be considered invariant. Default: 0.01 A
	include_hydrogens_in_uniqueness_analysis: bool. 
		This tag indicates if you want to include hydrogens when accessing uniquess between molecules and dimers. Default: False
	no_of_cpus : int.
		This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.

	Returns
	-------
	symmetric_molecule_pairs: dict.
		These are the indices to convert the atoms in a molecule from mol1_name --> mol2_name for a system that does contain hydrogens.
	"""

	# First, choose the invarience method you would like to use.
	if invariance_method_type == 'comprehensive':
		perform_invariance_method = remove_equivalent_molecules_comprehensive_invariance_method
	elif invariance_method_type == 'minimal_elemental_abundance':
		perform_invariance_method = remove_equivalent_molecules_minimal_elemental_abundance_invariance_method
	elif invariance_method_type == 'combination':
		perform_invariance_method = remove_equivalent_molecules_combination_invariance_method
	else:
		print('Error: The input method for the invariance method can be either (given as invariance_method_type):')
		print('\t* comprehensive: The Comprehensive Invariance Method')
		print('\t* minimal_elemental_abundance: The Minimal Elemental Abundance Invariance Method')
		print('\t* combination: The Combination Invariance Method')
		print('See https://github.com/geoffreyweal/ECCP for more information')
		exit('This program will finish without completing')

	# Second, perform the invarance method.
	symmetric_non_hydrogen_molecule_pairs = perform_invariance_method(unique_molecule_names, molecules, molecule_graphs, neighbouring_molecules_about_molecules=neighbouring_molecules_about_molecules, max_distance_disparity=max_distance_disparity, include_hydrogens_in_uniqueness_analysis=include_hydrogens_in_uniqueness_analysis, no_of_cpus=no_of_cpus)

	# Third, Convert the pairs from symmetric_molecule_pairs without hydrogen to with hydrogens
	symmetric_molecule_pairs = convert_non_hydrogen_molecule_indices_to_hydrogen_molecule_indices(symmetric_non_hydrogen_molecule_pairs, molecules)

	# Fourth, return information about the indices of the symmetric molecules.
	return symmetric_molecule_pairs

def convert_non_hydrogen_molecule_indices_to_hydrogen_molecule_indices(symmetric_non_hydrogen_molecule_pairs, molecules):
	"""
	This method is designed to convert the indices from the non-hydrogen version of the molecules to the hydrogen versions

	Parameters
	----------
	symmetric_non_hydrogen_molecule_pairs: dict.
		These are the indices to convert the atoms in a molecule from mol1_name --> mol2_name for a system that does not contain hydrogens.
	molecules : list of ase.Atoms
		This is the list of molecules in the crystal.

	Returns
	-------
	symmetric_molecule_pairs: dict.
		These are the indices to convert the atoms in a molecule from mol1_name --> mol2_name for a system that does contain hydrogens.
	"""

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# First, initialise a dictionary to help convert non-hydrogen molecule indices to hydrogen-based molecule indices
	non_hydrogen_to_hydrogen_conversion_for_all_molecules = {}

	# Second, for each molecule in the molecules dictionary: 
	for mol_name, molecule in sorted(molecules.items(), key=lambda x: x[0]):

		# 2.1: Set a counter to record the non-hydrogen index.
		non_hydrogen_index = 0

		# 2.2: Initalise a dictionary for converting non-hydrogen molecule indices to hydrogen-based molecule indices for a molecule.
		non_hydrogen_to_hydrogen_conversion = {}

		# 2.3: Obtain the indices for converting from non-hydrogen to hydrogen molecule format. 
		for atom_index, element in enumerate(molecule.get_chemical_symbols()):
			if element not in ['H', 'D', 'T']:
				non_hydrogen_to_hydrogen_conversion[non_hydrogen_index] = atom_index
				non_hydrogen_index += 1

		# 2.4: Add non_hydrogen_to_hydrogen_conversion to non_hydrogen_to_hydrogen_conversion_for_all_molecules
		#      for the current molecule mol_name.
		non_hydrogen_to_hydrogen_conversion_for_all_molecules[mol_name] = non_hydrogen_to_hydrogen_conversion

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Third, initialise a dictionary for holding the atom indices of hydrogen-based molecules
	symmetric_molecule_pairs = {}

	# Fourth, for each equivalent molecule pair in symmetric_non_hydrogen_molecule_pairs
	for (mol1_name, mol2_name), non_hydrogen_atom_index_comparisons in sorted(symmetric_non_hydrogen_molecule_pairs.items(), key=lambda x: x[0]):

		# 4.1: If mol1_name == mol2_name, there is a problem somewhere
		if mol1_name == mol2_name:
			raise Exception('Error: molecule is being compared to itself: mol1_name = '+str(mol1_name))

		# 4.2: If mol1_name > mol2_name, this is a bit odd
		if mol1_name > mol2_name:
			raise Exception('Error: mol1_name > mol2_name. This should be the other way around: mol1_name = '+str(mol1_name)+'; mol2_name = '+str(mol2_name))

		# 4.3: Obtain the non-hydrogen to hydrogen conversion dictionary for mol1_name
		mol1_non_hydrogen_to_hydrogen_conversion = non_hydrogen_to_hydrogen_conversion_for_all_molecules[mol1_name]

		# 4.4: Obtain the non-hydrogen to hydrogen conversion dictionary for mol1_name
		mol2_non_hydrogen_to_hydrogen_conversion = non_hydrogen_to_hydrogen_conversion_for_all_molecules[mol2_name]

		# 4.5: Initialise a new list for placing atom indices for hydrogen-based molecules
		hydrogen_based_atom_index_comparisons = []

		# 4.6: For each index comparison in non_hydrogen_atom_index_comparisons:
		for mol1_atom_index_non_hydrogen, mol2_atom_index_non_hydrogen in non_hydrogen_atom_index_comparisons:

			# 4.6.1: Convert mol1_atom_index from non-hydrogen to hydrogen-based
			mol1_atom_index_hydrogen = mol1_non_hydrogen_to_hydrogen_conversion[mol1_atom_index_non_hydrogen]

			# 4.6.2: Convert mol2_atom_index from non-hydrogen to hydrogen-based
			mol2_atom_index_hydrogen = mol2_non_hydrogen_to_hydrogen_conversion[mol2_atom_index_non_hydrogen]

			# 4.6.3: Create a tuple for converting between mol1 and mol2 in a hydrogen-based molecule system.
			hydrogen_based_atom_index_pair = (mol1_atom_index_hydrogen, mol2_atom_index_hydrogen)

			# 4.6.4: Append hydrogen_based_atom_index_pair to hydrogen_based_atom_index_comparisons
			hydrogen_based_atom_index_comparisons.append(hydrogen_based_atom_index_pair)

		# 4.7: Add hydrogen_based_atom_index_comparisons to symmetric_molecule_pairs for key (mol1_name, mol2_name)
		symmetric_molecule_pairs[(mol1_name, mol2_name)] = hydrogen_based_atom_index_comparisons

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Fifth, return symmetric_molecule_pairs
	return symmetric_molecule_pairs



