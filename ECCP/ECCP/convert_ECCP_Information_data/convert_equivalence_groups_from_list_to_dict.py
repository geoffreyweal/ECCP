"""
convert_equivalence_groups_from_list_to_dict.py, Geoffrey Weal, 26/3/2024

This method is designed to convert a set of equivalence groups from list for to dictionary form.
"""

def convert_equivalence_groups_from_list_to_dict(equivalence_groups_list):
	"""
	This method is designed to convert a set of equivalence groups from list for to dictionary form.

	Parameters
	----------
	equivalence_groups_list : list.
		This contains all the equivalence groups in this crystal in list form.

	Returns
	-------
	equivalence_groups_dict : dict.
		This contains all the equivalence groups in this crystal in dictionary form.
	"""

	# First, initialise the equivalence_groups_dict dictionary.
	equivalence_groups_dict = {}

	# Second, convert the equivalence_groups_list from list format into dictionary format.
	for equivalence_group in equivalence_groups_list:

		# 2.1: Obtain the representative unique individual for this equivalent group. 
		unique_individual = min(equivalence_group)

		# 2.2: Obtain a list of the other individuals in the equivalence group.
		#      * This is the list of equivalent individuals in the equivalence group.
		equivalent_individuals = sorted([individual for individual in equivalence_group if not (individual == unique_individual)])

		# 2.3: Add this equivalence group to equivalence_groups_dict
		equivalence_groups_dict[unique_individual] = equivalent_individuals

	# Third, return equivalence_groups_dict
	return equivalence_groups_dict