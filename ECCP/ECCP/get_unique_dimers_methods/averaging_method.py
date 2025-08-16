"""
averaging_method.py, Geoffrey Weal, Daniel Packwood, 11/3/22

This script is designed to use remove equivalent dimers based on the average distance between atoms in each molecule in the dimer.
"""
import numpy as np

from ase.geometry import get_distances

from SUMELF import remove_hydrogens

def remove_equivalent_dimers_averaging_method(dimers, molecules):
	"""
	This method will remove symmetric dimers from a list of dimers, returning all unique non-symmetric dimers. This analysis will ignore hydrogen atoms. 

	This method will first obtain the average distance between atoms in one molecule to another molecule in a dimer.

	The dimer is then represented by the average distance.

	The average distance between dimers is then used to determine if that dimer is eqivalent or not.

	Parameters
	----------
	dimers : list
		A list of dimers, consisting of the two molecules involved in the dimer.
	molecules : list of ase.Atoms
		This is the list of molecules that can be used to make dimers

	Returns
	-------
	equivalent_dimers : list
		A list of the indices of equivalent dimers in the dimers list. 
	"""

	print('add environment settings here')
	print('add distance comparison to help speed up')
	import pdb; pdb.set_trace()
	
	# First, get all the molecules without hydrogens.
	non_hydrogen_molecules = {}
	for mol_name in molecules.keys():
		non_hydrogen_molecule = remove_hydrogens(molecules[mol_name], graph=None)
		non_hydrogen_molecules[mol_name] = non_hydrogen_molecule

	# Second, obtain a list of all the distances between atoms in one molecule from another molecule
	average_distances_between_dimers = []
	for dimer_index, mol_name1, mol_name2, unit_cell_displacement, displacement, to_move_com_by in dimers:
		molecule1 = non_hydrogen_molecules[mol_name1].copy()
		molecule2 = non_hydrogen_molecules[mol_name2].copy()
		molecule2.set_positions(molecule2.get_positions() + displacement)
		elements1  = molecule1.get_chemical_symbols()
		elements2  = molecule2.get_chemical_symbols()
		positions1 = molecule1.get_positions()
		positions2 = molecule2.get_positions()
		dd_1_to_2  = get_all_distances_from_1_to_2(positions1,positions2,elements1,elements2)
		dd_2_to_1  = get_all_distances_from_1_to_2(positions2,positions1,elements2,elements1)
		all_distances = [j for sub in dd_1_to_2 for j in sub] + [j for sub in dd_2_to_1 for j in sub]
		average_distance = round(sum(all_distances)/float(len(all_distances)),2)
		average_distances_between_dimers.append(average_distance)

	# Third, determine which dimers are equivalent (symmetric)
	equivalent_dimers = []
	for dimer_index1 in range(len(average_distances_between_dimers)):
		dimer1_average_distance = average_distances_between_dimers[dimer_index1]
		for dimer_index2 in range(dimer_index1+1,len(average_distances_between_dimers)):
			dimer2_average_distance = average_distances_between_dimers[dimer_index2]
			if (dimer1_average_distance == dimer2_average_distance):
				equivalent_dimers.append((dimer_index1,dimer_index2))
	
	# Fourth, return a list of equivalent dimers.
	return equivalent_dimers

# ------------------------------------------------------------------------------------------------------------------------------

def get_all_distances_from_1_to_2(positions1: np.array, positions2: np.array, elements1: list, elements2: list) -> list:
	"""
	This method will give a 2D list of distance between atoms on molecule1 and atoms on molecule2.

	Parameters
	----------
	positions1 : np.array
		The position of atoms on molecule1
	positions2 : np.array
		The position of atoms on molecule2
	elements1 : np.array
		The elements of atoms on molecule1
	elements2 : np.array
		The elements of atoms on molecule2

	Returns
	-------
	dd_1_to_2 : list
		A 2D list of distance between atoms on molecule1 and atoms on molecule2. The elements involved in the bond have been included in this list
	"""
	dd_1_to_2 = []
	for position1, element1 in zip(positions1, elements1):
		distances = get_distances(position1,positions2)[1].tolist()[0]
		dd_1_to_2.append(distances)
	dd_1_to_2.sort()
	return dd_1_to_2

# ------------------------------------------------------------------------------------------------------------------------------





