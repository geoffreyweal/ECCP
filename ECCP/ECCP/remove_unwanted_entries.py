"""
remove_unwanted_entries.py, Geoffrey Weal, 13/3/2024

This method is designed to remove unwanted attributes from the ase.Atoms object if they exist.
"""
import numpy as np

unwanted_attributes = ['initial_magmoms', 'masses', 'momenta', 'tags']

def remove_unwanted_entries(original_molecule):
	"""
	This method is designed to remove unwanted attributes from the ase.Atoms object if they exist.
	
	Parameters
	----------
	original_molecule : ase.Atoms
		This is the original molecule that may contain the unwanted attributes.

	Returns
	-------
	molecule : ase.Atoms
		This is the molecule without any unwanted attributes.
	"""

	# First, make a copy of the original molecule.
	molecule = original_molecule.copy()

	# Second, remove any unwanted attributes from the molcule
	for unwanted_attribute in unwanted_attributes:

		# 2.1: If we do not want to remove this attribute, move on.
		if unwanted_attribute not in molecule.arrays.keys():
			continue

		# 2.2: If:
		#    * unwanted_attribute is 'masses', or 
		#    * does not contain a non-zero value,
		# remove the unwanted attribute from molecule.arrays.
		if unwanted_attribute in ['masses']:
			del molecule.arrays[unwanted_attribute]
		elif not contains_non_zero_values(molecule.arrays[unwanted_attribute]):
			del molecule.arrays[unwanted_attribute]

	# Third, return molecule that does not contain unwanted attributes
	return molecule

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

def contains_non_zero_values(list_of_values):
	"""
	This method is designed to check if the list given contains any non-zero values

	Parameters
	----------
	list_of_values : list, touple, numpy.array
		This is a list of value you want to check the values of

	Returns
	-------
	Returns True if there is at least one non-zero value in the list, otherwise return False. 
	"""

	# First, for each value in the list of values given:
	for value in list_of_values:

		# 1.1: If the value is an interger or float:
		if isinstance(value, (int, float)):

			# If value is non-zero, return True
			if value != 0.0:
				return True

		elif isinstance(value, (list, tuple, np.ndarray)):

			# 1.2: If value is a multi-number value, look through each of the values
			for vv in value:

				 # 1.2.1: If vv is non-zero, return True
				if vv != 0.0:
					return True

	# Second, if you have got here, all the values in the list were zero, so return False.
	return False

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

