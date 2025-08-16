"""
utilities.py, Geoffrey Weal, 10/3/24

This script contains minor methods for the ECCP program.
"""
import numpy as np
from SUMELF import get_distance, remove_hydrogens

def convert_dimer_details_to_neighbourhood_molecules_for_dimer_method(dimer_details, molecules, unit_cell_lattice_vectors, make_dimer_method):
	"""
	This method is designed to convert the information about the dimers from ``dimer_details`` format to ``dimer_details_to_neighbourhood_molecules`` format.

	Parameters
	----------
	dimer_details : dict
		This dictionary contains the information about all the dimers in the crystal for the given 'max_dimer_distance' setting in the make_dimer_method dictionary from the ``Run_ECCP.py`` file. The format of this file is {dimer_name: (mol1_name, mol2_name, UCV(i), UCV(j), UCV(k), DV(x), DV(y), DV(z))}
	molecules, unit_cell_lattice_vectors, make_dimer_method

	Returns
	-------
	neighbourhood_molecules_for_dimer_method : list
		This is the information from dimer_details in the format of (mol1_name, mol2_name, (UCV(i), UCV(j), UCV(k)), np.array(DV(x), DV(y), DV(z))))
	"""

	# First, make sure that the dimer name in dimer_details are in consecutive order from 1 to len(dimer_details)
	if not sorted(dimer_details.keys()) == sorted(range(1,len(dimer_details)+1)):
		to_string  = 'Error: The dimers names in dimer_details are not in consecutive order from 1 to len(dimer_details)+1\n'
		to_string += 'dimer names in dimer_details: '+str(list(dimer_details.keys()))+'\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# Second, get a copy of the molecules without hydrogens
	molecules_without_hydrogens = {mol_name: remove_hydrogens(molecule.copy()) for mol_name, molecule in molecules.items()}

	# Third, initialise the dimer_details_to_neighbourhood_molecules list
	dimer_details_to_neighbourhood_molecules = []

	# Fourth, go through each dimer in dimer_details and convert it into the format for dimer_details_to_neighbourhood_molecules
	for dimer_name in sorted(dimer_details.keys()):

		# 4.1: Obtain the dimer data from dimer_details[dimer_name]
		mol1_name, mol2_name, UCV_i, UCV_j, UCV_k, DV_x, DV_y, DV_z, move_COM_x, move_COM_y, move_COM_z = dimer_details[dimer_name]

		# 4.2: Get the unit cell ijk values
		unit_cell_ijk_vector = (UCV_i, UCV_j, UCV_k)

		# 4.3: Get the displacement vector from unit_cell_ijk_vector for molecule 2
		displacement_vector = np.matmul(np.array(unit_cell_ijk_vector), unit_cell_lattice_vectors)
		displacement_vector = np.array([round(value,12) for value in displacement_vector])

		# 4.4: Get the displacement vector for molecule 2 from the `All_Dimer_Information.txt` 
		printed_displacement_vector = np.array([round(DV_x,12), round(DV_y,12), round(DV_z,12)])

		# 4.5: Check that the displacement vector calculated from unit_cell_ijk_vector is the same 
		#      as the printed printed_displacement_vector
		if not np.array_equal(displacement_vector, printed_displacement_vector):
			raise Exception('Error: displacement vector are not the same')

		# 4.6: Get the distance between the molecules in the dimer using make_dimer_method
		distance_between_molecule_in_dimer = get_distance_between_molecules_in_dimer(molecules_without_hydrogens[mol1_name], molecules_without_hydrogens[mol2_name], displacement_vector, make_dimer_method)

		# 4.7: Make an entry for dimer_details_to_neighbourhood_molecules 
		dimers_detail_to_neighbourhood_molecules = (mol1_name, mol2_name, unit_cell_ijk_vector, displacement_vector, distance_between_molecule_in_dimer)

		# 4.8: Add dimers_detail_to_neighbourhood_molecules to dimer_details_to_neighbourhood_molecules
		dimer_details_to_neighbourhood_molecules.append(dimers_detail_to_neighbourhood_molecules)

	# Fifth, return dimer_details_to_neighbourhood_molecules
	return dimer_details_to_neighbourhood_molecules

# ------------------------------------------------------------------------------------------------------------------------------------------------

def get_distance_between_molecules_in_dimer(molecule_1, molecule_2, displacement_vector, make_dimer_method):
	"""
	This method is designed to obtain the distance between the two molecules in the dimer based on make_dimer_method.

	Parameters
	----------
	molecule_1 : ase.Atoms
		This is the ase.Atoms object of the first molecule. 
	molecule_2 : ase.Atoms
		This is the ase.Atoms object of the second molecule. 
	displacement_vector : 3D numpy.array
		This is the spatial amount to displace molecule 2 by for this dimer. 
	make_dimer_method : dict
		This contains the information about which method was used determine how to measure the distance between the twio molecules in the dimer. 

	Returns
	-------
	distance : float
		This is the distance between the two molecules in the dimer based on the make_dimer_method chosen to measure the distance.
	"""

	# First, get the name of the make_dimer_method method used
	make_dimer_method_name = make_dimer_method['method']

	# Second, get the distance between the two molecules in the dimer based on the make_dimer_method chosen to measure the distance.
	if   make_dimer_method_name == 'centre_of_mass':
		distance_between_molecule_in_dimer = get_distance_using_centre_of_mass_method(molecule_1, molecule_2, displacement_vector)
	elif make_dimer_method_name == 'centre_of_molecule':
		distance_between_molecule_in_dimer = get_distance_using_centre_of_molecule_method(molecule_1, molecule_2, displacement_vector)
	elif make_dimer_method_name == 'average_distance_method':
		distance_between_molecule_in_dimer = get_distance_using_average_distance_method_method(molecule_1, molecule_2, displacement_vector)
	elif make_dimer_method_name == 'nearest_atoms_method':
		distance_between_molecule_in_dimer = get_distance_using_nearest_atoms_method_method(molecule_1, molecule_2, displacement_vector)

	# Third, return the measured distance between the two molecules in the dimer. 
	return round(distance_between_molecule_in_dimer, 4)

def get_distance_using_centre_of_mass_method(molecule_1, molecule_2, displacement_vector):
	raise Exception('Check this method works')

def get_distance_using_centre_of_molecule_method(molecule_1, molecule_2, displacement_vector):
	raise Exception('Check this method works')

def get_distance_using_average_distance_method_method(molecule_1, molecule_2, displacement_vector):
	raise Exception('Check this method works')

def get_distance_using_nearest_atoms_method_method(molecule_1, molecule_2, displacement_vector):
	"""
	This method is designed to obtain the distance between the two molecules in the dimer based on nearest_atoms_method method.

	Parameters
	----------
	molecule_1 : ase.Atoms
		This is the ase.Atoms object of the first molecule. 
	molecule_2 : ase.Atoms
		This is the ase.Atoms object of the second molecule. 
	displacement_vector : 3D numpy.array
		This is the spatial amount to displace molecule 2 by for this dimer. 

	Returns
	-------
	distance : float
		This is the distance between the two molecules in the dimer based on the make_dimer_method chosen to measure the distance.
	"""

	# First, obtain the positions of all the atoms in the molecules:
	positions1 = molecule_1.get_positions()
	positions2 = molecule_2.get_positions()

	# Second, initialise the variable to record the distance between the two molecules to.
	shortest_distance = float('inf')

	# Third, determine the distance between the nearest atoms between the two molecules in the dimer. 
	for pos_index1 in range(len(positions1)):
		for pos_index2 in range(len(positions2)):
			distance = round(get_distance(positions1[pos_index1], positions2[pos_index2] + displacement_vector), 4)
			# If we have found a distance that is less than max_distance, then these two molcules are a neighbouring pair.
			if (distance <= shortest_distance):
				shortest_distance = distance

	# Fourth, return shortest_distance
	return shortest_distance

# ------------------------------------------------------------------------------------------------------------------------------------------------







