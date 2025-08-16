"""
convert_existing_unique_molecule_data.py, Geoffrey Weal, 10/3/2024

This method is designed to convert the unique molecules information from the ECCP_Information folder for the ECCP program.
"""
from collections import Counter
from ECCP.ECCP.convert_ECCP_Information_data.convert_equivalence_groups_from_list_to_dict import convert_equivalence_groups_from_list_to_dict

def convert_existing_unique_molecule_data_from_ECCP_Information(structurally_equivalent_molecule_groups_list, conformationally_equivalent_molecule_groups_list, no_of_molecules):
	"""
	This method is designed to convert the unique molecules information from the ECCP_Information folder for the ECCP program.

	Parameters
	----------
	structurally_equivalent_molecule_groups_list : list.
		This list contains all the structurally equivalent molecule groups in this crystal in list form.
	conformationally_equivalent_molecule_groups_list : list. 
		This list contains all the conformationally equivalent molecule groups in this crystal in list form.
	no_of_molecules : int
		This is the number of molecules in the crystal.

	Returns
	-------
	structurally_unique_molecule_names : list of int.
		These are the names of all the structurally unique molecules. 
	structurally_equivalent_molecule_groups_dict : dict.
		This dictionary contains all the structurally equivalent molecule groups in this crystal in dictionary form.
	conformationally_unique_molecule_names : list of int.
		These are the names of all the conformationally unique molecules. 
	conformationally_equivalent_molecule_groups_dict : dict. 
		This dictionary contains all the conformationally equivalent molecule groups in this crystal in dictionary form.
	"""

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	# First, perform checks on the structurally equivalent molecule groups. 

	# 1.1: Make sure there are no duplicates of molecules within the structurally equivalent molecule groups. 
	for structurally_equivalent_molecule_group in structurally_equivalent_molecule_groups_list:
		if not len(structurally_equivalent_molecule_group) == len(set(structurally_equivalent_molecule_group)):
			to_string  = 'Error: There are duplicates in one or more structurally_equivalent_molecule_group\n'
			duplicates_in_structurally_equivalent_molecule_group = [mol_name for mol_name, count in Counter(structurally_equivalent_molecule_group).items() if (count >= 2)]
			to_string += f'Duplicates in earliest structurally_equivalent_molecule_group: {duplicates_in_structurally_equivalent_molecule_group}\n'
			to_string += f'Current structurally_equivalent_molecule_group: {structurally_equivalent_molecule_group}\n'
			to_string += f'all structurally_equivalent_molecule_groups_list: {structurally_equivalent_molecule_groups_list}\n'
			raise Exception(to_string)

	# 1.2: Obtain all the molecules from structurally_equivalent_molecule_groups_list
	all_molecule_names_in_structurally_equivalent_molecule_groups_list = [j for sub in structurally_equivalent_molecule_groups_list for j in sub]
	
	# 1.3: Make sure that there are no duplicate molecules across all the structurally equivalent molecule groups. 
	if not len(all_molecule_names_in_structurally_equivalent_molecule_groups_list) == len(set(all_molecule_names_in_structurally_equivalent_molecule_groups_list)):
		to_string  = 'Error: There are duplicates in structurally_equivalent_molecule_groups_list\n'
		duplicates_in_structurally_equivalent_molecule_groups = [mol_name for mol_name, count in Counter(all_molecule_names_in_structurally_equivalent_molecule_groups_list).items() if (count >= 2)]
		to_string += f'Duplicates molecules: {sorted(duplicates_in_structurally_equivalent_molecule_groups)}\n'
		to_string += f'all structurally_equivalent_molecule_groups_list: {sorted(structurally_equivalent_molecule_groups_list)}\n'
		raise Exception(to_string)

	# 1.4: Make sure all the expected molecule names are found in structurally_equivalent_molecule_groups_list
	if not sorted(all_molecule_names_in_structurally_equivalent_molecule_groups_list) == list(range(1,no_of_molecules+1)):
		to_string  = 'Error: There are missing molecules, extra odd molecules, or repeated molecules in structurally_equivalent_molecule_groups_list\n'
		to_string += f'All structurally_equivalent_molecule_groups_list: {sorted(structurally_equivalent_molecule_groups_list)}\n'
		to_string += f'All molecules in structurally_equivalent_molecule_groups: {sorted(all_molecule_names_in_structurally_equivalent_molecule_groups_list)}\n'
		to_string += f'Expected molecules: {list(range(1,no_of_molecules+1))}\n'
		raise Exception(to_string)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	# Second, perform checks on the structurally equivalent molecule groups. 

	# 2.1: Make sure there are no duplicates of molecules within the conformationally equivalent molecule groups. 
	for conformationally_equivalent_molecule_group in conformationally_equivalent_molecule_groups_list:
		if not len(conformationally_equivalent_molecule_group) == len(set(conformationally_equivalent_molecule_group)):
			to_string  = 'Error: There are duplicates in one or more conformationally_equivalent_molecule_group\n'
			duplicates_in_conformationally_equivalent_molecule_group = [mol_name for mol_name, count in Counter(conformationally_equivalent_molecule_group).items() if (count >= 2)]
			to_string += f'Duplicates in earliest conformationally_equivalent_molecule_group: {sorted(duplicates_in_conformationally_equivalent_molecule_group)}\n'
			to_string += f'Current conformationally_equivalent_molecule_group: {sorted(conformationally_equivalent_molecule_group)}\n'
			to_string += f'all conformationally_equivalent_molecule_groups_list: {sorted(conformationally_equivalent_molecule_groups_list)}\n'
			raise Exception(to_string)

	# 2.2: Obtain all the molecules from conformationally_equivalent_molecule_groups_list
	all_molecule_names_in_conformationally_equivalent_molecule_groups_list = [j for sub in conformationally_equivalent_molecule_groups_list for j in sub]
	
	# 2.3: Make sure that there are no duplicate molecules across all the conformationally equivalent molecule groups. 
	if not len(all_molecule_names_in_conformationally_equivalent_molecule_groups_list) == len(set(all_molecule_names_in_conformationally_equivalent_molecule_groups_list)):
		to_string  = 'Error: There are duplicates in conformationally_equivalent_molecule_groups_list\n'
		duplicates_in_conformationally_equivalent_molecule_groups = [mol_name for mol_name, count in Counter(all_molecule_names_in_conformationally_equivalent_molecule_groups_list).items() if (count >= 2)]
		to_string += f'Duplicates molecules: {sorted(duplicates_in_conformationally_equivalent_molecule_groups)}\n'
		to_string += f'all conformationally_equivalent_molecule_groups_list: {sorted(conformationally_equivalent_molecule_groups_list)}\n'
		raise Exception(to_string)

	# 2.4: Make sure all the expected molecule names are found in conformationally_equivalent_molecule_groups_list
	if not sorted(all_molecule_names_in_conformationally_equivalent_molecule_groups_list) == list(range(1,no_of_molecules+1)):
		to_string  = 'Error: There are missing molecules, extra odd molecules, or repeated molecules in conformationally_equivalent_molecule_groups_list\n'
		to_string += f'All conformationally_equivalent_molecule_groups_list: {sorted(conformationally_equivalent_molecule_groups_list)}\n'
		to_string += f'All molecules in conformationally_equivalent_molecule_groups: {sorted(all_molecule_names_in_conformationally_equivalent_molecule_groups_list)}\n'
		to_string += f'Expected molecules: {list(range(1,no_of_molecules+1))}\n'
		raise Exception(to_string)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	# Third, obtain the structurally equivalent molecule groups in dictionary format.
	structurally_equivalent_molecule_groups_dict = convert_equivalence_groups_from_list_to_dict(structurally_equivalent_molecule_groups_list)

	# Fourth, obtain the structurally unqiue molecules.
	structurally_unique_molecule_names = list(structurally_equivalent_molecule_groups_dict.keys())

	# Fifth, obtain the conformationally equivalent molecule groups in dictionary format.
	conformationally_equivalent_molecule_groups_dict = convert_equivalence_groups_from_list_to_dict(conformationally_equivalent_molecule_groups_list)

	# Sixth, obtain the conformationally unqiue molecules.
	conformationally_unique_molecule_names = list(conformationally_equivalent_molecule_groups_dict.keys())

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	# Seventh, return structurally_unique_molecule_names, structurally_equivalent_molecule_groups_dict, conformationally_unique_molecule_names, and conformationally_equivalent_molecule_groups_dict
	return structurally_unique_molecule_names, structurally_equivalent_molecule_groups_dict, conformationally_unique_molecule_names, conformationally_equivalent_molecule_groups_dict


