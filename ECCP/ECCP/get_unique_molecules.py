"""
get_unique_molecules.py, Geoffrey Weal, 14/4/22

This script will obtain unique molecules from the list of molecules. 

This script is designed to determine which molecules are structurally equivalent, including essentually identical in bond distances and angles.
"""
from copy import deepcopy
from collections import Counter
from itertools import product, combinations

from ECCP.ECCP.get_unique_molecules_methods.get_symmetries_of_molecules_in_crystal    import get_symmetries_of_molecules_in_crystal
from ECCP.ECCP.get_unique_molecules_methods.get_conformationally_equivalent_molecules import get_conformationally_equivalent_molecules
from ECCP.ECCP.get_unique_molecules_methods.invariance_method                         import determine_equivalent_molecules_with_invariance_method
from ECCP.ECCP.get_unique_molecules_methods.get_unique_utility_methods                import get_equivalent_pair_names, convert_equivalent_groups_to_dict

def get_unique_molecules(molecules, molecule_graphs, crystal, molecule_equivalence_method={'method': 'invariance_method', 'type': 'combination'}, neighbouring_molecules_about_molecules={}, include_hydrogens_in_uniqueness_analysis=False, no_of_cpus=1):
	"""
	This method will determine which molecules are unique and which are non-unique, and will report back only those molecules that are unique.

	Parameters
	----------
	molecules : list of ase.Atoms
		These are all the molecules that were obtained from the crystal. 
	molecule_graphs : list of networkx.Graph
		These are the graphs that correspond to the molecules in the molecules list
	crystal : ase.Atoms
		The crystal in Atomic Simulation Environment format
	molecule_equivalence_method : dict.
		This is the information required for determining which molecules are equivalent. 
	neighbouring_molecules_about_molecules : dict.
		This is the information about the molecules that surround (in the vicinity of) each moleule in the crystal.
	include_hydrogens_in_uniqueness_analysis: bool. 
		This tag indicates if you want to include hydrogens when accessing uniquess between molecules and dimers. Default: False
	no_of_cpus : int.
		This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.

	Returns
	-------
	structurally_unique_molecule_names : list of ase.Atoms
		These are all the structurally unique molecules that were identified from the crystal. 
	structurally_equivalent_molecules : dict.
		This is a dictionary of the structurally unique molecule, with a lit of all the other molecules in the crystal that are equivalent to this unique molecule. 
	conformationally_unique_molecule_names : list of ase.Atoms
		These are all the conformationally unique molecules that were identified from the crystal. 
	conformationally_equivalent_molecules : dict.
		This is a dictionary of the conformationally unique molecule, with a lit of all the other molecules in the crystal that are equivalent to this unique molecule. 
	"""

	# First, write a list that will contain the names of the unique molecules.
	unique_molecule_names = sorted(molecules.keys())

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Second, determine from the symmetry of the crystal structure which molecules are obviously symmetric due to crystal symmetry.
	# This process allows us to very very quickly eliminate some molecules are being equivalent before using a more rigous method for uniqueness.
	symmetric_molecule_pairs_due_to_symmetries_in_crystal = get_symmetries_of_molecules_in_crystal(molecules, crystal)

	# Third, obtain all the groups of equivalent molecules. 
	structurally_equivalent_molecule_groups_due_to_symmetries_in_crystal = get_equivalent_pair_names(list(symmetric_molecule_pairs_due_to_symmetries_in_crystal.keys()), all_individuals_names=unique_molecule_names)

	# Fourth, convert structurally_equivalent_molecule_groups_due_to_symmetries_in_crystal from list for to dict form as {unique molecule: list of equivalent molecules}.
	structurally_equivalent_molecule_groups_due_to_symmetries_in_crystal = convert_equivalent_groups_to_dict(structurally_equivalent_molecule_groups_due_to_symmetries_in_crystal, all_individuals_names=unique_molecule_names) # Also referred to as equivalent_molecules or structurally_equivalent_molecules

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# Fifth, remove the symmetric molecules that were determined from the symmetry of the crystal.
	
	# 5.1: Get the list of non-unique molecules
	unique_molecule_names_check = []
	non_unique_molecule_names = []
	for unique_mol_name, equivalent_mol_names in sorted(structurally_equivalent_molecule_groups_due_to_symmetries_in_crystal.items()):
		unique_molecule_names_check.append(unique_mol_name)
		non_unique_molecule_names += list(equivalent_mol_names)
	
	# 5.2: Check that a unique molecule is not also in the non-unique molecules list
	for mol_name in unique_molecule_names_check:
		if mol_name in non_unique_molecule_names:
			to_string  = 'Error: There is an mol_name from the unique_molecule_names_check list that is also in the non_unique_molecule_names list.\n'
			to_string += 'unique_molecule_names_check = '+str(unique_molecule_names_check)+'\n'
			to_string += 'non_unique_molecule_names = '+str(non_unique_molecule_names)+'\n'
			to_string += 'Check this out'
			raise Exception(to_string)
	
	# 5.3: Remove the names from the non_unique_molecule_names from the unique_molecule_names list.
	for mol_name in sorted(non_unique_molecule_names, reverse=True):
		unique_molecule_names.remove(mol_name)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# Sixth, obtain all the structurally equivalent molecule groups for this crystal. 

	# 6.1: Obtaining the get_structurally_unique_molecule_names and structurally_equivalent_molecules variables. 
	symmetric_molecule_pairs = get_structurally_unique_molecule_names(unique_molecule_names, molecules, molecule_graphs, molecule_equivalence_method, neighbouring_molecules_about_molecules, include_hydrogens_in_uniqueness_analysis=include_hydrogens_in_uniqueness_analysis, no_of_cpus=no_of_cpus)

	# 6.2: Obtain all the groups of equivalent molecules. 
	structurally_equivalent_molecule_groups = get_equivalent_pair_names(list(symmetric_molecule_pairs.keys()), all_individuals_names=unique_molecule_names)

	# 6.3: Convert structurally_equivalent_molecule_groups from list for to dict form as {unique molecule: list of equivalent molecules}
	structurally_equivalent_molecule_groups = convert_equivalent_groups_to_dict(structurally_equivalent_molecule_groups, all_individuals_names=unique_molecule_names) # Also referred to as equivalent_molecules or structurally_equivalent_molecules

	# 6.4: Add the symmetric molecules due to crystal symmetries from structurally_equivalent_molecule_groups_due_to_symmetries_in_crystal to structurally_equivalent_molecule_groups.
	structurally_equivalent_molecule_groups = add_crystal_symmetry_molecules_to_structurally_equivalent_molecule_groups(structurally_equivalent_molecule_groups, structurally_equivalent_molecule_groups_due_to_symmetries_in_crystal)

	# 6.5: Obtain the names of the unique molecules in the equivalence groups. 
	structurally_unique_molecule_names = sorted(structurally_equivalent_molecule_groups.keys())

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# Seventh, obtain all the pairs of symmetric molecules in the crystal.

	# 7.1: Initalise the list to hold all the symmetric molecule pairs.
	symmetric_molecule_pairs_including_crystal_symmetry = []

	# 7.2: Add the symmetric molecule pairs due to the symmetry of the crystal.
	for mol1_name, mol2_name in symmetric_molecule_pairs:

		# 7.2.1: Obtain all the equivalent molecules for molecule 1.
		all_equivalent_mol1_molecules = [mol1_name] + list(structurally_equivalent_molecule_groups_due_to_symmetries_in_crystal[mol1_name])

		# 7.2.2: Obtain all the equivalent molecules for molecule 2.
		all_equivalent_mol2_molecules = [mol2_name] + list(structurally_equivalent_molecule_groups_due_to_symmetries_in_crystal[mol2_name])

		# 7.2.3: For each set of molecule pairs, including considerations due to crystal symmetries
		for m1n, m2n in product(all_equivalent_mol1_molecules, all_equivalent_mol2_molecules):

			# 7.2.4: Check that m1n is not the same as m2n.
			if m1n == m2n:
				to_string  = 'Error: the molecules being added as structurally equivalent molecule pairs are the same molecule.'+'\n'
				to_string += 'm1n = '+str(m1i)+'; m2n = '+str(m2n)+'\n'
				to_string += 'This is probably a programming issue. Check this'
				raise Exception(to_string)

			# 7.2.5: Add this molecule pair to symmetric_molecule_pairs_including_crystal_symmetry.
			symmetric_molecule_pairs_including_crystal_symmetry.append(tuple(sorted([m1n, m2n])))

	# 7.3: Add the pairs from structurally_equivalent_molecule_groups_due_to_symmetries_in_crystal to symmetric_molecule_pairs_including_crystal_symmetry
	for unique_molecule, list_of_equivalent_molecules in structurally_equivalent_molecule_groups_due_to_symmetries_in_crystal.items():

		# 7.3.1: Get the structurally equivalent molecule group due to crystal symmetry.
		structurally_equivalence_molecule_group_due_to_crystal_symmetry = sorted([unique_molecule] + list(list_of_equivalent_molecules))

		# 7.3.2: All all the combinations of structurally equivalence molecule group.
		symmetric_molecule_pairs_including_crystal_symmetry += list(combinations(structurally_equivalence_molecule_group_due_to_crystal_symmetry, 2))

	# 7.4: Check that none of the pairs in symmetric_molecule_pairs_including_crystal_symmetry are duplicates
	if not len(symmetric_molecule_pairs_including_crystal_symmetry) == len(set(symmetric_molecule_pairs_including_crystal_symmetry)):
		to_string  = 'Error: Some of the pairs in symmetric_molecule_pairs_including_crystal_symmetry are duplicates.'+'\n'
		duplicate_symmetric_molecule_pairs = Counter(symmetric_molecule_pairs_including_crystal_symmetry) - Counter(set(symmetric_molecule_pairs_including_crystal_symmetry))
		if any([(count < 0) for count in duplicate_symmetric_molecule_pairs.value()]):
			raise Exception('Error: negative counts have been found in duplicate_symmetric_molecule_pairs.\nduplicate_symmetric_molecule_pairs: '+str(duplicate_symmetric_molecule_pairs)+'\nCheck this.')
		duplicate_symmetric_molecule_pairs = sorted([pair for pair, count in duplicate_symmetric_molecule_pairs.items() if (count > 0)])
		to_string += 'Duplicate pairs in symmetric_molecule_pairs_including_crystal_symmetry: '+str(duplicate_symmetric_molecule_pairs)+'\n'
		to_string += 'symmetric_molecule_pairs_including_crystal_symmetry = '+str(symmetric_molecule_pairs_including_crystal_symmetry)+'\n'
		to_string += 'Check this out'
		raise Exception(to_string)

	# 7.5: Sort symmetric_molecule_pairs_including_crystal_symmetry
	symmetric_molecule_pairs_including_crystal_symmetry.sort()

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# Eighth, obtain all the conformationally equivalent molecule groups for this crystal.  

	# 8.1: Obtaining the conformationally_unique_molecule_names and conformationally_equivalent_molecules variables. 
	conformationally_equivalent_molecule_groups = get_conformationally_equivalent_molecules(structurally_equivalent_molecule_groups, molecules, molecule_graphs)

	# 8.2: Obtain the names of the unique molecules in the equivalence groups. 
	conformationally_unique_molecule_names = sorted(conformationally_equivalent_molecule_groups.keys())

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Ninth, return names of unique molecules in the molecules list.
	return structurally_unique_molecule_names, structurally_equivalent_molecule_groups, conformationally_unique_molecule_names, conformationally_equivalent_molecule_groups, symmetric_molecule_pairs_including_crystal_symmetry

# -------------------------------------------------------------------------------------------------------------------------

def get_structurally_unique_molecule_names(unique_molecule_names_original, molecules, molecule_graphs, molecule_equivalence_method, neighbouring_molecules_about_molecules, include_hydrogens_in_uniqueness_analysis=False, no_of_cpus=1):
	"""
	This method is designed to determine which molecules are unique and which molecules are equivalent.

	Parameters
	----------
	unique_molecule_names_original : list of int
		These are the names of the molecules in the crystal unit cell.
	molecules : list of ase.Atoms
		These are all the molecules in the crystal unit cell.
	molecule_graphs : list of networkx.Graph
		These are all the bonding connectivit graphs of molecules in the crystal.
	molecule_equivalence_method : dict.
		This is the information required for determining which molecules are equivalent. 
	neighbouring_molecules_about_molecules : dict.
		This is the information about the molecules that surround (in the vicinity of) each moleule in the crystal.
	include_hydrogens_in_uniqueness_analysis: bool. 
		This tag indicates if you want to include hydrogens when accessing uniquess between molecules and dimers. Default: False
	no_of_cpus : int.
		This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.

	Returns
	-------
	unique_molecules : list of ase.Atoms
		These are all the unique molecules that were identified from the crystal. 
	equivalent_molecules : dict.
		This is a dictionary of the unique molecule, with a lit of all the other molecules in the crystal that are equivalent to this unique molecule. 
	"""

	# First, make a copy of the list of all the molecules in the crystal unit cell.
	unique_molecule_names = deepcopy(unique_molecule_names_original)

	# Second, use a rigous method for determining which from the remaining unique_molecules list are truely unique.
	molecule_equivalence_method_type = molecule_equivalence_method['method']
	if molecule_equivalence_method_type.lower() == 'none' or (molecule_equivalence_method_type.lower() is None):
		symmetric_molecule_pairs = []
	elif molecule_equivalence_method_type == 'invariance_method':
		invariance_method_type = molecule_equivalence_method.get('type','combination')
		max_distance_disparity = molecule_equivalence_method.get('max_distance_disparity',0.01)
		symmetric_molecule_pairs = determine_equivalent_molecules_with_invariance_method(invariance_method_type, unique_molecule_names, molecules, molecule_graphs, neighbouring_molecules_about_molecules=neighbouring_molecules_about_molecules, max_distance_disparity=max_distance_disparity, include_hydrogens_in_uniqueness_analysis=include_hydrogens_in_uniqueness_analysis, no_of_cpus=no_of_cpus)
	else:
		print('Error: The input method for the equivalence method can be either:')
		print('\t* invariance_method: The Invariance Method')
		print('See https://github.com/geoffreyweal/ECCP for more information')
		exit('This program will finish without completing')

	# Third, return unique_molecules and equivalent_molecules variables.
	return symmetric_molecule_pairs

# -------------------------------------------------------------------------------------------------------------------------

def add_crystal_symmetry_molecules_to_structurally_equivalent_molecule_groups(structurally_equivalent_molecule_groups_without_crystal_symmetry, structurally_equivalent_molecule_groups_due_to_symmetries_in_crystal):
	"""
	This method is designed to add the symmetric molecules due to crystal symmetry (from tructurally_equivalent_molecule_groups_due_to_symmetries_in_crystal) to the structurally equivalent molecule groups.

	Parameters
	----------
	structurally_equivalent_molecule_groups_without_crystal_symmetry : dict.
		This dictionary contains all the structurally equivalent groups from the crystal. This does not include molecules that are symmetric due to the symmetry of the crystal. 
	structurally_equivalent_molecule_groups_due_to_symmetries_in_crystal : dict.
		This dictionary contains all the molecules that are structurally equivalent to each other due to the crystal symmetry. 

	Returns
	-------
	structurally_equivalent_molecule_groups : dict.
		This dictionary contains all the structurally equivalent groups from the crystal, having been updated with the structurally equivalent molecules due to the symmetry of the crystal.
	"""

	# First, make a copy of structurally_equivalent_molecule_groups_without_crystal_symmetry
	structurally_equivalent_molecule_groups = deepcopy(structurally_equivalent_molecule_groups_without_crystal_symmetry)

	# Second, for each list of equivalent molecules in the structurally equivalent molecule groups, add molecules that are structurally the same due to the symmetry of the crystal.
	for structurally_unique_molecule_name, structurally_equivalent_molecule_names in structurally_equivalent_molecule_groups.items():
		for structurally_equivalent_molecule_name in deepcopy(structurally_equivalent_molecule_names):
			structurally_equivalent_molecule_groups[structurally_unique_molecule_name] += list(structurally_equivalent_molecule_groups_due_to_symmetries_in_crystal[structurally_equivalent_molecule_name])

	# Third, for each unique molecule in the structurally equivalent molecule groups, add molecules that are structurally the same due to the symmetry of the crystal.
	for structurally_unique_molecule_name in structurally_equivalent_molecule_groups.keys():
		structurally_equivalent_molecule_groups[structurally_unique_molecule_name] += list(structurally_equivalent_molecule_groups_due_to_symmetries_in_crystal[structurally_unique_molecule_name])

	# Fourth, sort all the lists in structurally_equivalent_molecule_groups
	for structurally_unique_molecule_name in structurally_equivalent_molecule_groups.keys():
		structurally_equivalent_molecule_groups[structurally_unique_molecule_name].sort()

	# Fifth, return structurally_equivalent_molecule_groups
	return structurally_equivalent_molecule_groups

# -------------------------------------------------------------------------------------------------------------------------







