"""
convert_existing_unique_dimer_data.py, Geoffrey Weal, 10/3/2024

This method is designed to convert the unique dimers information from the ECCP_Information folder for the ECCP program.
"""
from collections import Counter
from ECCP.ECCP.convert_ECCP_Information_data.convert_equivalence_groups_from_list_to_dict import convert_equivalence_groups_from_list_to_dict

def convert_existing_unique_dimer_data_from_ECCP_Information(structurally_equivalent_dimer_groups_list, no_of_dimers):
	"""
	This method is designed to convert the unique dimers information from the ECCP_Information folder for the ECCP program.

	Parameters
	----------
	structurally_equivalent_dimer_groups_list : list.
		This list contains all the structurally equivalent dimer groups in this crystal in list form.
	no_of_dimers : int
		This is the number of dimers from the ECCP run.

	Returns
	-------
	structurally_unique_dimer_names : list of int.
		These are the names of all the structurally unique dimers. 
	structurally_equivalent_dimer_groups_dict : dict.
		This dictionary contains all the structurally equivalent dimer groups in this crystal in dictionary form.
	"""

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	# First, perform checks on the structurally equivalent dimer groups. 

	# 1.1: Make sure there are no duplicates of dimers within the structurally equivalent dimer groups. 
	for structurally_equivalent_dimer_group in structurally_equivalent_dimer_groups_list:
		if not len(structurally_equivalent_dimer_group) == len(set(structurally_equivalent_dimer_group)):
			to_string  = 'Error: There are duplicates in one or more structurally_equivalent_dimer_group\n'
			duplicates_in_structurally_equivalent_dimer_group = [mol_name for mol_name, count in Counter(structurally_equivalent_dimer_group).items() if (count >= 2)]
			to_string += f'Duplicates in earliest structurally_equivalent_dimer_group: {sorted(duplicates_in_structurally_equivalent_dimer_group)}\n'
			to_string += f'Current structurally_equivalent_dimer_group: {sorted(structurally_equivalent_dimer_group)}\n'
			to_string += f'all structurally_equivalent_dimer_groups_list: {sorted(structurally_equivalent_dimer_groups_list)}\n'
			raise Exception(to_string)

	# 1.2: Obtain all the dimers from structurally_equivalent_dimer_groups_list
	all_dimer_names_in_structurally_equivalent_dimer_groups_list = [j for sub in structurally_equivalent_dimer_groups_list for j in sub]
	
	# 1.3: Make sure that there are no duplicate dimers across all the structurally equivalent dimer groups. 
	if not len(all_dimer_names_in_structurally_equivalent_dimer_groups_list) == len(set(all_dimer_names_in_structurally_equivalent_dimer_groups_list)):
		to_string  = 'Error: There are duplicates in structurally_equivalent_dimer_groups_list\n'
		duplicates_in_structurally_equivalent_dimer_groups = [mol_name for mol_name, count in Counter(all_dimer_names_in_structurally_equivalent_dimer_groups_list).items() if (count >= 2)]
		to_string += f'Duplicates dimers: {sorted(duplicates_in_structurally_equivalent_dimer_groups)}\n'
		to_string += f'all structurally_equivalent_dimer_groups_list: {sorted(structurally_equivalent_dimer_groups_list)}\n'
		raise Exception(to_string)

	# 1.4: Make sure all the expected dimer names are found in structurally_equivalent_dimer_groups_list
	if not sorted(all_dimer_names_in_structurally_equivalent_dimer_groups_list) == list(range(1,no_of_dimers+1)):
		to_string  = 'Error: There are missing dimers, extra odd dimers, or repeated dimers in structurally_equivalent_dimer_groups_list\n'
		to_string += f'All structurally_equivalent_dimer_groups_list: {sorted(structurally_equivalent_dimer_groups_list)}\n'
		to_string += f'All dimers in structurally_equivalent_dimer_groups: {sorted(all_dimer_names_in_structurally_equivalent_dimer_groups_list)}\n'
		to_string += f'Expected dimers: {list(range(1,no_of_dimers+1))}\n'
		raise Exception(to_string)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	# Second, obtain the structurally equivalent dimer groups in dictionary format.
	structurally_equivalent_dimer_groups_dict = convert_equivalence_groups_from_list_to_dict(structurally_equivalent_dimer_groups_list)

	# Third, obtain the structurally unqiue dimers.
	structurally_unique_dimer_names = list(structurally_equivalent_dimer_groups_dict.keys())

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	# Fourth, return structurally_unique_dimer_names and structurally_equivalent_dimer_groups_dict
	return structurally_unique_dimer_names, structurally_equivalent_dimer_groups_dict
