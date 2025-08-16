"""
get_unique_utility_methods.py, Geoffrey Weal, 19/3/2024

This script is designed to provide methods for get_unique_molecules.py and get_unique_dimers.py
"""
from copy import deepcopy
from collections import Counter

def get_equivalent_pair_names(original_equivalent_pairs, all_individuals_names):
	"""
	This method is designed to determine which pairs can be assigned to a collective set of individuals that are equivalent to each other. 

	This method is design so that each individual in a group is symmetric to each other individual based on original_equivalent_pairs. 

	If an individual is only equivalent to some of the molecules in a group, this individual will be placed in a new group.

	Parameters
	----------
	original_equivalent_pairs : list of (int,int)
		A list of pairs of names that are equivalent to each other. 
	all_individuals_names : list of ints
		This is the names of the individuals in the list being examined. 

	Returns
	-------
	equivalent_individuals_names : list
		These are the names of the equivalent individuals to not include in the unique individuals list.
	unique_individuals_names : list
		A list of unique individuals.
	"""

	# First, make a copy of the original_equivalent_pairs to use for this method. 
	equivalent_pairs = sorted(deepcopy(original_equivalent_pairs))

	# Second, initialise a list that contains all the groups of equivalent individuals
	equivalent_individuals_groups = []

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# Third, go through each pair in equivalent_pairs and put them into groups
	# * Each individual in a group is equivalent to EVERY OTHER individual in the group.
	# * If this can not be done, the individual that can not fit into a group will be put into a new group containing itself.

	# 3.1: This section is designed to create equivalence groups from individuals in equivalent_pairs.
	#      If an individual from (indiv1, indiv2) can not be placed in a equivalence group, then
	#        * It might be placed in another equivalence group during the second for loop, or
	#        * It may be in a group on it's own, which in this case it will be left for step 3.4.
	#
	#      For each individual in equivalent_pairs
	for indiv1, indiv2 in equivalent_pairs:

		# 3.2: Find if (indiv1, indiv2) can be added to a group in equivalent_individuals_groups. 
		for unique_mol_name, equivalent_individuals_group in enumerate(equivalent_individuals_groups):

			# 3.2.1: Determine if the individuals in this individual pair are already both in equivalent_individuals_group
			is_indiv1_in_equivalent_individuals_group = indiv1 in equivalent_individuals_group
			is_indiv2_in_equivalent_individuals_group = indiv2 in equivalent_individuals_group

			# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
			# 3.2.2: If indiv1 and indiv2 are in equivalent_individuals_group, we can move on as we have figured out 
			#        the equivalency of these two individuals.
			if is_indiv1_in_equivalent_individuals_group and is_indiv2_in_equivalent_individuals_group:
				break

			# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
			# 3.2.3: If indiv1 is in equivalent_individuals_group, check if indiv2 could be added to equivalent_individuals_group.
			if is_indiv1_in_equivalent_individuals_group:

				# 3.2.3.1: If indiv2 is equivalent to all individuals in equivalent_individuals_group, then add it to this equivalence group.
				#          * If not, make a new equivalence set for indiv2.
				if check_indiv_is_equivalent_to_all_enteries_in_equivalent_individuals_group(indiv2, equivalent_individuals_group, equivalent_pairs):
					equivalent_individuals_groups[unique_mol_name].append(indiv2)

				# 3.2.3.2: We have found the equivalence group that indiv1 is in, so either indiv2 has been added to 
				#          equivalent_individuals_groups[unique_mol_name] or not.
				#          * If it has been added, good
				#          * If it has not been added, we will leave it till the end of the loop, where we have determined 
				#            if it is already apart of an equivalance group or a new equivalance group needs to be created
				#            with it in it. 
				break

			# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
			# 3.2.4: If indiv2 is in equivalent_individuals_group, check if indiv2 could be added to equivalent_individuals_group
			if is_indiv2_in_equivalent_individuals_group:

				# 3.2.4.1: If indiv1 is equivalent to all individuals in equivalent_individuals_group, then add it to this equivalence group
				#          * If not, make a new equivalence set for indiv1
				if check_indiv_is_equivalent_to_all_enteries_in_equivalent_individuals_group(indiv1, equivalent_individuals_group, equivalent_pairs):
					equivalent_individuals_groups[unique_mol_name].append(indiv1)

				# 3.2.4.2: We have found the equivalent individual group that indiv2 is in, so either indiv1 has been added to 
				#          equivalent_individuals_groups[unique_mol_name] or not
				#          * If it has been added, good
				#          * If it has not been added, we will leave it till the end of the loop, where we have determined 
				#            if it is already apart of an equivalance group or a new equivalance group needs to be created
				#            with it in it. 
				break

		else:

			# 3.3: We have not found neither indiv1 or indiv2 in equivalent_individuals_groups, so make a new equivalence group for these individuals.
			equivalent_individuals_groups.append([indiv1, indiv2])

	# 3.4: Add any individuals that were not included in equivalent_pairs. This means these individuals are the only individuals in their own equivalency group. 
	for individual_no in all_individuals_names:
		for equivalent_individuals_group in equivalent_individuals_groups:
			if individual_no in equivalent_individuals_group:
				break
		else:
			equivalent_individuals_groups.append([individual_no])

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# Fourth, sort all the lists in equivalent_individuals_groups
	for unique_mol_name in range(len(equivalent_individuals_groups)):
		equivalent_individuals_groups[unique_mol_name].sort()

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# Fifth, perform several checks to make sure that equivalent_individuals_groups is consistent

	# 5.1: Create a list that contains all the individuals in equivalent_individuals_groups
	all_dimer_names_in_equivalent_individuals_groups = [j for sub in equivalent_individuals_groups for j in sub]

	# 5.2: Check that none of the groups have the same individuals in them as other groups. 
	if not len(all_dimer_names_in_equivalent_individuals_groups) == len(set(all_dimer_names_in_equivalent_individuals_groups)):
		to_string  = 'Error: There are repeated names in equivalent_individuals_groups'
		double_countered_dimers = Counter(all_dimer_names_in_equivalent_individuals_groups) - Counter(set(all_dimer_names_in_equivalent_individuals_groups))
		if any([individual_dimer_index for individual_dimer_index, count in double_countered_dimers.items() if (count < 0)]):
			raise Exception('Error: negative counts given in double_countered_dimers')
		double_countered_individual_names = sorted([individual_dimer_index for individual_dimer_index, count in double_countered_dimers.items() if (count > 0)])
		to_string += 'Double-counted individuals names in equivalent_individuals_groups = '+str(double_countered_individual_names)+'\n'
		to_string += 'equivalent_individuals_groups                                     = '+str(equivalent_individuals_groups)+'\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# 5.3: Check that all the individual that are expected based on all_individuals_names are in equivalent_individuals_groups (all_dimer_names_in_equivalent_individuals_groups)
	if not sorted(all_dimer_names_in_equivalent_individuals_groups) == sorted(all_individuals_names):
		to_string  = 'Error: Some individuals expected in equivalent_individuals_groups were not found'+'\n'
		to_string += 'equivalent_individuals_groups = '+str(equivalent_individuals_groups)+'\n'
		to_string += 'dimers in equivalent_individuals_groups          = '+str(sorted(all_dimer_names_in_equivalent_individuals_groups))+'\n'
		to_string += 'dimers expected in equivalent_individuals_groups = '+str(sorted(all_individuals_names))+'\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Sixth, return equivalent_individuals_groups
	return equivalent_individuals_groups

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# ---------------------------------------------------------------------------------------------------------------

def check_indiv_is_equivalent_to_all_enteries_in_equivalent_individuals_group(new_indiv, equivalent_individuals_group, equivalent_pairs):
	"""
	This method will check that new_indiv is equivalent to all molecules in equivalent_individuals_group by 
	checking that (new_indiv, molecule in equivalent_individuals_group) is in equivalent_pairs. 

	Parameters
	----------
	new_indiv : int.
		This is the molecule we want to check if it is equivalent to all molecules in equivalent_individuals_group./
	equivalent_individuals_group : list of int.
		This is the equivalent group we want to check if all molecules in it are equivalent to new_indiv
	equivalent_pairs : list of (int, int). 
		These are all the equivalent pairs of molecules in the crystal. 

	Returns
	-------
		True if new_indiv is equivalent to all molecules in equivalent_individuals_group. False otherwise. 
	"""

	# First, make a list that contains all the equivalent individual pairs that must be found in equivalent_pairs in order to allow us to add new_indiv to equivalent_individuals_group.
	individual_pairs_to_check_for_equivalency = [tuple(sorted([new_indiv, other_indiv])) for other_indiv in equivalent_individuals_group if not (new_indiv == other_indiv)]

	# Second, check that all the individual pairs in individual_pairs_to_check_for_equivalency are found in equivalent_pairs.
	#         * This is the same checking if new_indiv is equivalent to all individuals in equivalent_individuals_group.
	are_all_individual_pairs_equivalent = all([(indiv_pair in equivalent_pairs) for indiv_pair in individual_pairs_to_check_for_equivalency])

	# Third, return are_all_individual_pairs_equivalent
	#        * This will be true if all (new_indiv, a molecule in equivalent_individuals_group) are found in equivalent_pairs.
	return are_all_individual_pairs_equivalent

# ---------------------------------------------------------------------------------------------------------------

def convert_equivalent_groups_to_dict(equivalent_individuals_groups, all_individuals_names=None):
	"""
	This method is designed to convert equivalent_individuals_groups into dictionary form

	Parameters
	----------
	equivalent_individuals_groups : list
		These are the names of the equivalent individuals to not include in the unique individuals list.
	all_individuals_names : list of ints
		This is the number of individuals in the list being examined. 

	Returns
	-------
	equivalent_individuals_groups_dict : dict.
		This is the dictionary form of equivalent_individuals_groups_dict, where the key represents the unique molecules in the groups. 
	"""

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# First, perform several checks to make sure that equivalent_individuals_groups is consistent

	# 1.1: Create a list that contains all the individuals in equivalent_individuals_groups
	all_dimer_names_in_equivalent_individuals_groups = [j for sub in equivalent_individuals_groups for j in sub]

	# 1.2: Check that none of the groups have the same individuals in them as other groups. 
	if not len(all_dimer_names_in_equivalent_individuals_groups) == len(set(all_dimer_names_in_equivalent_individuals_groups)):
		to_string  = 'Error: There are repeated names in equivalent_individuals_groups'
		double_countered_dimers = Counter(all_dimer_names_in_equivalent_individuals_groups) - Counter(set(all_dimer_names_in_equivalent_individuals_groups))
		if any([individual_dimer_index for individual_dimer_index, count in double_countered_dimers.items() if (count < 0)]):
			raise Exception('Error: negative counts given in double_countered_dimers')
		double_countered_individual_names = sorted([individual_dimer_index for individual_dimer_index, count in double_countered_dimers.items() if (count > 0)])
		to_string += 'Double-counted individuals names in equivalent_individuals_groups = '+str(double_countered_individual_names)+'\n'
		to_string += 'equivalent_individuals_groups = '+str(equivalent_individuals_groups)+'\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# 1.3: Check that all the individual that are expected based on all_individuals_names are in equivalent_individuals_groups (all_dimer_names_in_equivalent_individuals_groups)
	if all_individuals_names is not None:
		if not sorted(all_dimer_names_in_equivalent_individuals_groups) == sorted(all_individuals_names):
			to_string  = 'Error: Some individuals expected in equivalent_individuals_groups were not found'+'\n'
			to_string += 'equivalent_individuals_groups = '+str(equivalent_individuals_groups)+'\n'
			to_string += 'dimers in equivalent_individuals_groups          = '+str(sorted(all_dimer_names_in_equivalent_individuals_groups))+'\n'
			to_string += 'dimers expected in equivalent_individuals_groups = '+str(sorted(all_individuals_names))+'\n'
			to_string += 'Check this.'
			raise Exception(to_string)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Second, initialise a dictionary to record unique and equivalent individuals for equivalance groups.
	equivalent_individuals_groups_dict = {}

	# Third, for each equivalnce group.
	for equivalent_individuals_group in equivalent_individuals_groups:

		# 3.1: Obtain the group leader that will act as the unique molecule.
		equivalence_group_leader = min(equivalent_individuals_group)

		# 3.2: Obtain the other individuals. These will act as the equivalent individuals. 
		equivalent_individuals = sorted([individual_name for individual_name in equivalent_individuals_group if (not individual_name == equivalence_group_leader)])

		# 3.3: Add the unique and equivalent molecules from equivalent_individuals_groups into equivalent_individuals_groups_dict.
		equivalent_individuals_groups_dict[equivalence_group_leader] = equivalent_individuals

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# Fourth, perform several checks to make sure that equivalent_individuals_groups_dict is consistent.

	# 4.1: Check that there are no double-counted individuals across the individual groups
	for unique_individual_mol_name, list_of_equivalent_individuals in equivalent_individuals_groups_dict.items():
		if not len(list_of_equivalent_individuals) == len(set(list_of_equivalent_individuals)):
			to_string  = 'Error: There are repeated names in one or more of the equivalence groups.'+'\n'
			to_string += 'equivalent_individuals_groups_dict = '+str(equivalent_individuals_groups_dict)+'\n'
			to_string += 'Check this.'
			raise Exception(to_string)

	# 4.2: Create a list that contains all the individuals in equivalent_individuals_groups_dict.
	all_dimer_names_in_equivalent_individuals_groups_dict = []
	for unique_individual_mol_name, list_of_equivalent_individuals in equivalent_individuals_groups_dict.items(): 
		all_dimer_names_in_equivalent_individuals_groups_dict += [unique_individual_mol_name] + list(list_of_equivalent_individuals)

	# 4.3: Check that none of the groups have the same individuals in them as other groups. 
	if not len(all_dimer_names_in_equivalent_individuals_groups_dict) == len(set(all_dimer_names_in_equivalent_individuals_groups_dict)):
		to_string  = 'Error: There are repeated names in equivalent_individuals_groups_dict'+'\n'
		double_countered_dimers = Counter(all_dimer_names_in_equivalent_individuals_groups_dict) - Counter(set(all_dimer_names_in_equivalent_individuals_groups_dict))
		if any([individual_dimer_index for individual_dimer_index, count in double_countered_dimers.items() if (count < 0)]):
			raise Exception('Error: negative counts given in double_countered_dimers')
		double_countered_individual_names = sorted([individual_dimer_index for individual_dimer_index, count in double_countered_dimers.items() if (count > 0)])
		to_string += 'Double-counted individuals names in equivalent_individuals_groups_dict = '+str(double_countered_individual_names)+'\n'
		to_string += 'equivalent_individuals_groups_dict = '+str(equivalent_individuals_groups_dict)+'\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# 4.4: Check that all the individual that are expected based on all_individuals_names are in equivalent_individuals_groups_dict (all_dimer_names_in_equivalent_individuals_groups_dict).
	if all_individuals_names is not None:
		if not sorted(all_dimer_names_in_equivalent_individuals_groups_dict) == sorted(all_individuals_names):
			to_string  = 'Error: Some individuals expected in equivalent_individuals_groups_dict were not found'+'\n'
			to_string += 'equivalent_individuals_groups_dict = '+str(equivalent_individuals_groups_dict)+'\n'
			to_string += 'dimers in equivalent_individuals_groups_dict          = '+str(sorted(all_dimer_names_in_equivalent_individuals_groups_dict))+'\n'
			to_string += 'dimers expected in equivalent_individuals_groups_dict = '+str(sorted(all_individuals_names))+'\n'
			to_string += 'Check this.'
			raise Exception(to_string)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Third, return equivalent_individuals_groups_dict
	return equivalent_individuals_groups_dict

# ---------------------------------------------------------------------------------------------------------------

