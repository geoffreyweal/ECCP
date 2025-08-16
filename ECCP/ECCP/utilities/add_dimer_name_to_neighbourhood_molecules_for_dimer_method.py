"""
add_dimer_name_to_neighbourhood_molecules_for_dimer_method.py, 22/5/24

This method is designed to include the name of the dimer to the dimer entries in add_dimer_name_to_neighbourhood_molecules_for_dimer_method
"""
from copy import deepcopy

def add_dimer_name_to_neighbourhood_molecules_for_dimer_method(neighbourhood_molecules_for_dimer_method_original):
	"""
	This method is designed to include the name of the dimer to the dimer entries in add_dimer_name_to_neighbourhood_molecules_for_dimer_method

	Parameters
	----------
	neighbourhood_molecules_for_dimer_method_original : list.
		This is neighbourhood_molecules_for_dimer_method without dimers included.

	Returns
	-------
	neighbourhood_molecules_for_dimer_method : list.
		This is neighbourhood_molecules_for_dimer_method with dimers included.
	"""

	# First, make a copy of the neighbourhood_molecules_for_dimer_method dictionary.
	neighbourhood_molecules_for_dimer_method = deepcopy(neighbourhood_molecules_for_dimer_method_original)

	# Second, for each entry in neighbourhood_molecules_for_dimer_method
	for dimer_index in range(len(neighbourhood_molecules_for_dimer_method)):

		# 2.1: Get the name of the dimer.
		dimer_name = dimer_index + 1

		# 2.2: Update neighbourhood_molecules_for_dimer_method to include the dimers name.
		neighbourhood_molecules_for_dimer_method[dimer_index] = tuple([dimer_name] + list(neighbourhood_molecules_for_dimer_method[dimer_index]))

	# Third, return neighbourhood_molecules_for_dimer_method.
	return neighbourhood_molecules_for_dimer_method