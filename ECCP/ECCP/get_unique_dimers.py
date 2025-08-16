"""
get_unique_dimers.py, Geoffrey Weal, 23/2/22

This script will obtain unique dimers from the list of dimers. 
"""
from collections import Counter
from copy        import deepcopy
from itertools   import chain, permutations

from ECCP.ECCP.get_unique_dimers_methods.atomic_distance_method  import remove_equivalent_dimers_atomic_distance_method
from ECCP.ECCP.get_unique_dimers_methods.averaging_method        import remove_equivalent_dimers_averaging_method
from ECCP.ECCP.get_unique_dimers_methods.invariance_method       import remove_equivalent_dimers_invariance_method

def get_unique_dimers(all_dimers_info, molecules, molecule_graphs, dimer_equivalence_method={'method': 'invariance_method'}, neighbouring_molecules_about_dimers={}, include_hydrogens_in_uniqueness_analysis=False, no_of_cpus=1):
	"""
	This method will obtain unique dimers from a list of dimers.

	Parameters
	----------
	all_dimers_info : list
		A list of dimers, consisting of the two molecules involved in the dimer.
	molecules : list of ase.Atoms
		These are the molecules that can be used to make dimers from.
	molecule_graphs : list of network.Graph
		These are the graphs associated with the molecules in the molecules list that can be used to make dimers from.
	dimer_equivalence_method : dict
		This is information required for the equivalence method you would like to use for determining equivalent dimers. See https://github.com/geoffreyweal/ECCP for more information. 
	neighbouring_molecules_about_dimers : dict.
		This is the information about the molecules that surround (in the vicinity of) each dimer in the crystal.
	include_hydrogens_in_uniqueness_analysis: bool. 
		This tag indicates if you want to include hydrogens when accessing uniquess between molecules and dimers. Default: False
	no_of_cpus : int.
		This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.

	Returns
	-------
	unique_dimers_indices : list
		A list of unique dimers.
	structurally_equivalent_dimer_groups : dict.
		THis is a dictionary of all the structurally equivalent dimers groups. Format given as -> representative unique dimer index: list of structurally equivalent dimer indices. 
	structurally_equivalent_dimer_pairs : list
		This is a list of all the dimers that the structurally equivalent to each other. 
	"""

	# First, obtain the list of indices of dimers in the dimer list that are equivalent. 
	dimer_equivalence_method_type = dimer_equivalence_method['method']
	if dimer_equivalence_method_type.lower() == 'none' or (dimer_equivalence_method_type.lower() is None):
		structurally_equivalent_dimers = []
	elif dimer_equivalence_method_type == 'atomic_distance_method':
		structurally_equivalent_dimers = remove_equivalent_dimers_atomic_distance_method(all_dimers_info, molecules, include_hydrogens_in_uniqueness_analysis=include_hydrogens_in_uniqueness_analysis)
	elif dimer_equivalence_method_type == 'averaging_method':
		structurally_equivalent_dimers = remove_equivalent_dimers_averaging_method(all_dimers_info, molecules, include_hydrogens_in_uniqueness_analysis=include_hydrogens_in_uniqueness_analysis)
	elif dimer_equivalence_method_type == 'invariance_method':
		invariance_method_type = dimer_equivalence_method.get('type','combination')
		max_distance_disparity = dimer_equivalence_method.get('max_distance_disparity',None)
		structurally_equivalent_dimers = remove_equivalent_dimers_invariance_method(invariance_method_type, all_dimers_info, molecules, molecule_graphs, neighbouring_molecules_about_dimers=neighbouring_molecules_about_dimers, max_distance_disparity=max_distance_disparity, no_of_cpus=no_of_cpus)
	else:
		print('Error: The input method for the equivalence method can be either:')
		print('\t* atomic_distance_method: The Atomic Distance Method')
		print('\t* invariance_method: The Invariance Method')
		print('\t* averaging_method: The Averaging Method')
		print('See https://github.com/geoffreyweal/ECCP for more information')
		exit('This program will finish without completing')

	# Second, make a copy of the equivalent_dimers variable that does not change and can be returned. 
	structurally_equivalent_dimer_pairs = deepcopy(structurally_equivalent_dimers)

	# Third, cycle through the list of symmetric dimers, and remove those symmetric dimers
	#   * This section works by remove those dimers that are symmetric with most other dimers first.
	#   * This shouldnt likely matter, but just in case this is written in this way.
	structurally_equivalent_dimer_groups = get_structurally_equivalent_dimers_groups(structurally_equivalent_dimers, all_dimer_names=sorted(all_dimers_info.keys()))

	# Fourth, obtain the unique dimers from the stucturally equivalent dimer groups
	unique_dimers_indices = sorted(structurally_equivalent_dimer_groups.keys())

	# Fifth, return unique_dimers_indices, structurally_equivalent_dimer_groups, and structurally_equivalent_dimer_pairs.
	return unique_dimers_indices, structurally_equivalent_dimer_groups, structurally_equivalent_dimer_pairs

# -------------------------------------------------------------------------------------------------------------------------

def get_structurally_equivalent_dimers_groups(original_equivalent_dimers, all_dimer_names):
	"""
	Third, cycle through the list of symmetric dimers, and remove those symmetric dimers
		* This section works by remove those dimers that are symmetric with most other dimers first.
		* This shouldnt likely matter, but just in case this is written in this way.

	Parameters
	----------
	original_equivalent_dimers : list
		A list of the indices of equivalent dimers in the dimers list. 
	all_dimer_names : list of ints
		These are all the names of the dimers sampled. 

	Returns
	-------
	equivalent_dimers_indices : list
		These are the indices of the equivalent dimers to not include in the unique dimers list.
	unique_dimers_indices : list
		A list of unique dimers.
	"""

	# First, make a copy of the original_equivalent_dimers to use for this method. 
	equivalent_dimers = sorted(deepcopy(original_equivalent_dimers))

	# Second, initialise a list that contains all the groups of equivalent dimers
	structurally_equivalent_dimers_groups_as_list = []

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	# Third, go through each pair of dimers in equivalent_dimers and put them into groups
	# * Each dimer in a group is equivalent to every other dimer in the group
	# * If this can not be done, the dimer that can not fit into a group will be put into a new group containing itself.
	# * This will mean that there could be dimers that are equivalent to each other, but in different equivalence groups 
	#   * due to being slightly too different to the other dimers in that equivalence group.
	#   * To prevent any problems arising, we have made the decision to define an equivalence as one where all dimers have been found to be equivalent to each other. 

	# 3.1: This section is designed to create equivalence groups from dimers in equivalent_dimers.
	#      # If a dimer from (dimer1_name, dimer2_name) can not be placed in a equivalence group, then
	#        * It might be placed in another equivalence group during the second for loop, or
	#        * It may be in a group on it's own, which in this case it will be left for step 3.4.
	#
	#      For each dimer in equivalent_dimers
	for dimer1_name, dimer2_name in equivalent_dimers:

		# 3.2: Find if (dimer1_name, dimer2_name) can be added to a group in structurally_equivalent_dimers_groups_as_list. 
		for index, equivalent_dimers_group in enumerate(structurally_equivalent_dimers_groups_as_list):

			# 3.2.1: Determine if the dimers in this dimer pair are already both in equivalent_dimers_group
			is_dimer1_in_equivalent_dimers_group = dimer1_name in equivalent_dimers_group
			is_dimer2_in_equivalent_dimers_group = dimer2_name in equivalent_dimers_group

			# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
			# 3.2.2: If dimer1 and dimer2 are in equivalent_dimers_group, we can move on as we have figured out 
			#        the equivalency of these two dimers.
			if is_dimer1_in_equivalent_dimers_group and is_dimer2_in_equivalent_dimers_group:
				break

			# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
			# 3.2.3: If dimer1 is in equivalent_dimers_group, check if dimer2 could be added to equivalent_dimers_group
			if is_dimer1_in_equivalent_dimers_group:

				# 3.2.3.1: Make a list that contains all the equivalent dimer pairs that must be found in 
				#          equivalent_dimers in order to allow us to add dimer2 to equivalent_dimers_group
				dimer_pairs_to_check_for_equivalency = [tuple(sorted([dimer2_name, other_dimer])) for other_dimer in equivalent_dimers_group if not (dimer2_name == other_dimer)]

				# 3.2.3.2: Check if all the dimer pairs in dimer_pairs_to_check_for_equivalency are found in equivalent_dimers
				#          * This is the same checking if dimer2 is equivalent to all dimers in equivalent_dimers_group
				are_all_dimer_pairs_equivalent = all([(dimer_pair in equivalent_dimers) for dimer_pair in dimer_pairs_to_check_for_equivalency])

				# 3.2.3.3: If dimer2 is equivalent to all dimers in equivalent_dimers_group, then add it to this equivalence group
				#          * If not, make a new equivalence set for dimer2
				if are_all_dimer_pairs_equivalent:
					structurally_equivalent_dimers_groups_as_list[index].append(dimer2_name)

				# 3.2.3.4: We have found the equivalent dimer group that dimer1 is in, so either dimer2 has been added to 
				#          structurally_equivalent_dimers_groups_as_list[index] or not
				#          * If it has been added, good
				#          * If it has not been added, we will leave it till the end of the loop, where we have determined 
				#            if it is already apart of an equivalance group or a new equivalance group needs to be created
				#            with it in it. 
				break

			# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
			# 3.2.4: If dimer2 is in equivalent_dimers_group, check if dimer2 could be added to equivalent_dimers_group
			if is_dimer2_in_equivalent_dimers_group:

				# 3.2.4.1: Make a list that contains all the equivalent dimer pairs that must be found in 
				#          equivalent_dimers in order to allow us to add dimer1 to equivalent_dimers_group
				dimer_pairs_to_check_for_equivalency = [tuple(sorted([dimer1_name, other_dimer])) for other_dimer in equivalent_dimers_group if not (dimer1_name == other_dimer)]

				# 3.2.4.2: Check if all the dimer pairs in dimer_pairs_to_check_for_equivalency are found in equivalent_dimers
				#          * This is the same checking if dimer1 is equivalent to all dimers in equivalent_dimers_group
				are_all_dimer_pairs_equivalent = all([(dimer_pair in equivalent_dimers) for dimer_pair in dimer_pairs_to_check_for_equivalency])

				# 3.2.4.3: If dimer1 is equivalent to all dimers in equivalent_dimers_group, then add it to this equivalence group
				#          * If not, make a new equivalence set for dimer1
				if are_all_dimer_pairs_equivalent:
					structurally_equivalent_dimers_groups_as_list[index].append(dimer1_name)

				# 3.2.4.4: We have found the equivalent dimer group that dimer2 is in, so either dimer1 has been added to 
				#          structurally_equivalent_dimers_groups_as_list[index] or not
				#          * If it has been added, good
				#          * If it has not been added, we will leave it till the end of the loop, where we have determined 
				#            if it is already apart of an equivalance group or a new equivalance group needs to be created
				#            with it in it. 
				break

		else:

			# 3.3: We have not found neither dimer1 or dimer2 in structurally_equivalent_dimers_groups_as_list, so make a new equivalence group for these dimers.
			structurally_equivalent_dimers_groups_as_list.append([dimer1_name, dimer2_name])

		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	# 3.4: Add any dimers that were not included in equivalent_dimers. This means these dimers are the only dimers in their own equivalency group. 
	for dimer_no in all_dimer_names:
		for equivalent_dimers_group in structurally_equivalent_dimers_groups_as_list:
			if dimer_no in equivalent_dimers_group:
				break
		else:
			structurally_equivalent_dimers_groups_as_list.append([dimer_no])

	# 3.5: Sort all the lists in structurally_equivalent_dimers_groups_as_list
	for index in range(len(structurally_equivalent_dimers_groups_as_list)):
		structurally_equivalent_dimers_groups_as_list[index].sort()

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	# Fourth, perform some checks to make sure that structurally_equivalent_dimers_groups_as_list is consistent. 

	# 4.1: Check that there are no dimers have been double-countered within a structurally equivalent dimer group
	for list_of_equivalent_dimer_indices in structurally_equivalent_dimers_groups_as_list:
		if not len(list_of_equivalent_dimer_indices) == len(set(list_of_equivalent_dimer_indices)):
			to_string  = 'Error: At last one of the dimer groups has double countered dimer indices.\n'
			to_string += 'list_of_equivalent_dimer_indices of interest = '+str(list_of_equivalent_dimer_indices)+'\n'
			to_string += 'structurally_equivalent_dimers_groups_as_list = '+str(structurally_equivalent_dimers_groups_as_list)+'\n'
			to_string += 'Check this.'
			raise Exception(toString)

	# 4.2: Create a list for containing all the dimer indices given in structurally_equivalent_dimers_groups_as_list
	all_dimers_collectively_from_structurally_equivalent_dimers_groups_as_list_groups = []
	for list_of_equivalent_dimer_indices in structurally_equivalent_dimers_groups_as_list: 
		all_dimers_collectively_from_structurally_equivalent_dimers_groups_as_list_groups += list(list_of_equivalent_dimer_indices)

	# 4.3: Check for double-counting of dimer indices across dimer groups in structurally_equivalent_dimers_groups_as_list:
	if not len(all_dimers_collectively_from_structurally_equivalent_dimers_groups_as_list_groups) == len(set(all_dimers_collectively_from_structurally_equivalent_dimers_groups_as_list_groups)):
		to_string  = 'Error: At last one of the dimers has been double countered in structurally_equivalent_dimers_groups_as_list.\n'
		double_counted_dimers = Counter(all_dimers_collectively_from_structurally_equivalent_dimers_groups_as_list_groups) - Counter(set(all_dimers_collectively_from_structurally_equivalent_dimers_groups_as_list_groups))
		if not any([(value < 0) for value in double_counted_dimers.values()]):
			raise Exception('Issue with Counter. This is a programming error')
		double_counted_dimers = sorted([dimer_index for dimer_index, count in double_counted_dimers.items() if (count > 0)])
		to_string += 'Double counted dimers: '+str(double_counted_dimers)+'\n'
		to_string += 'Dimers in structurally_equivalent_dimers_groups_as_list = '+str(sorted(all_dimers_collectively_from_structurally_equivalent_dimers_groups_as_list_groups))+'\n'
		to_string += 'Check this.'
		raise Exception(toString)

	# 4.4: Check that all the dimers in original_equivalent_dimers have been accounted for based on all_dimer_names
	if not sorted(all_dimers_collectively_from_structurally_equivalent_dimers_groups_as_list_groups) == list(all_dimer_names):
		to_string  = 'Error: There may be missing dimers indices in structurally_equivalent_dimers_groups_as_list.\n'
		to_string += 'Dimers in structurally_equivalent_dimers_groups_as_list = '+str(sorted(all_dimers_collectively_from_structurally_equivalent_dimers_groups_as_list_groups))+'\n'
		to_string += 'Dimers expected = '+str(list(all_dimer_names))+'\n'
		to_string += 'Check this.'
		raise Exception(toString)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	# Fifth, obtain the unique dimer to act as a representative for each structurally equivalent dimer group. 

	# 5.1: Initialise a dictionary for recording equivalent dimers
	structurally_equivalent_dimers_groups = {}

	# 5.2: Add each structurally equivalent dimer group to structurally_equivalent_dimers_groups
	for equivalent_dimers_group in structurally_equivalent_dimers_groups_as_list:

		# 5.2.1: Obtain the unique dimer to represent this dimer, which is the lowest index dimer in the group
		unique_dimer_index = min(equivalent_dimers_group)

		# 5.2.2: Remove the unique dimer from equivalent_dimers_group
		equivalent_dimers_group = sorted([dimer_index for dimer_index in equivalent_dimers_group if not (dimer_index == unique_dimer_index)])

		# 5.2.3: Add the equivalent dimer group to structurally_equivalent_dimers_groups
		structurally_equivalent_dimers_groups[unique_dimer_index] = tuple(equivalent_dimers_group)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	# Sixth, perform some checks to make sure that structurally_equivalent_dimers_groups is consistent. 

	# 6.1: Check that there are no dimers have been double-countered within a structurally equivalent dimer group
	for unique_dimer_index, list_of_equivalent_dimer_indices in structurally_equivalent_dimers_groups.items():

		# 6.1.1: Make a list that includes all the dimers in the structurally equivalent dimers group.
		equivalent_dimers_group = [unique_dimer_index] + list(list_of_equivalent_dimer_indices)

		# 6.1.2: Check if there are any double-counting of dimer indices in the structurally equivalent dimers group.
		if not len(equivalent_dimers_group) == len(set(equivalent_dimers_group)):
			to_string  = 'Error: At last one of the dimer groups has double countered dimer indices.\n'
			to_string += 'equivalent_dimers_group of interest = '+str(equivalent_dimers_group)+'\n'
			to_string += 'structurally_equivalent_dimers_groups = '+str(structurally_equivalent_dimers_groups)+'\n'
			to_string += 'Check this.'
			raise Exception(toString)

	# 6.2: Create a list for containing all the dimer indices given in structurally_equivalent_dimers_groups
	all_dimers_collectively_from_structurally_equivalent_dimers_groups = []
	for unique_dimer_index, list_of_equivalent_dimer_indices in structurally_equivalent_dimers_groups.items(): 
		all_dimers_collectively_from_structurally_equivalent_dimers_groups += [unique_dimer_index] + list(list_of_equivalent_dimer_indices)

	# 6.3: Check for double-counting of dimer indices across dimer groups in structurally_equivalent_dimers_groups:
	if not len(all_dimers_collectively_from_structurally_equivalent_dimers_groups) == len(set(all_dimers_collectively_from_structurally_equivalent_dimers_groups)):
		to_string  = 'Error: At last one of the dimers has been double countered in structurally_equivalent_dimers_groups.\n'
		double_counted_dimers = Counter(all_dimers_collectively_from_structurally_equivalent_dimers_groups) - Counter(set(all_dimers_collectively_from_structurally_equivalent_dimers_groups))
		if not any([(value < 0) for value in double_counted_dimers.values()]):
			raise Exception('Issue with Counter. This is a programming error')
		double_counted_dimers = sorted([dimer_index for dimer_index, count in double_counted_dimers.items() if (count > 0)])
		to_string += 'Double counted dimers: '+str(double_counted_dimers)+'\n'
		to_string += 'Dimers in structurally_equivalent_dimers_groups = '+str(sorted(all_dimers_collectively_from_structurally_equivalent_dimers_groups))+'\n'
		to_string += 'Check this.'
		raise Exception(toString)

	# 6.4: Check that all the dimers in original_equivalent_dimers have been accounted for based on all_dimer_names
	if not sorted(all_dimers_collectively_from_structurally_equivalent_dimers_groups) == list(all_dimer_names):
		to_string  = 'Error: There may be missing dimers indices in structurally_equivalent_dimers_groups.\n'
		to_string += 'Dimers in structurally_equivalent_dimers_groups = '+str(sorted(all_dimers_collectively_from_structurally_equivalent_dimers_groups))+'\n'
		to_string += 'Dimers expected = '+str(list(all_dimer_names))+'\n'
		to_string += 'Check this.'
		raise Exception(toString)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	# Seventh, return structurally_equivalent_dimers_groups
	return structurally_equivalent_dimers_groups

# -------------------------------------------------------------------------------------------------------------------------

