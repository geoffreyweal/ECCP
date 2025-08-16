"""
utilities.py, Geoffrey Weal, 9/2/24

This script includes auxiliary methods for the ECCP program to use. 
"""
import numpy as np

def get_permutated_indices_list(comparison): 
	"""
	This method will provide the permutation list of how to reorder the atoms indices in a molecule from a dictionary that tell this program how the indices of one molecule translate into another equivalent molecule. 

	Parameters
	----------
	comparison : dict.
		A dictionary that contains how the indices of one molecule translate into another identical molecule. 

	Returns
	-------
	idx : list
		The permutation list that tells the program how to reorder atoms in a list.
	"""
	permutation1 = [value for key, value in sorted(comparison.items())]
	try:
		idx = np.empty_like(permutation1)
		idx[permutation1] = np.arange(len(permutation1))
	except:
		import pdb; pdb.set_trace()
	return idx