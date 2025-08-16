"""
convert_dimer_details_to_all_dimers_info_method.py, Geoffrey Weal, 22/5/24

This method will obtain dimers and produce all_dimers_info from the dimer_details list.
"""
import numpy as np
from tqdm import tqdm

def convert_dimer_details_to_all_dimers_info_method(dimer_details, neighbourhood_molecules_for_dimer_method):
	"""
	This method will obtain dimers and produce all_dimers_info from the dimer_details list.

	BACKGROUND
	----------

	When performing the original "get_dimers" method on crystals that had already been run by ECCP and had data stored in the "ECCP_Information" folder, 
	Some dimers were shifted by the "centre_molecule_in_cell" method (in SUMELF) by a cell lattice distance out of what was recorded on file (for dimers 
	in xyz files EET files, etc). The reason for this was because when using the distance_methods in SUMELF-->distance_methods.py, there can be floating 
	point errors arising due to computers working in binary rather than base 10. As many of the dimers had centre of masses that were likely in the centre
	of the cells or in positions that were degenerately close to the centre and origin as other positions, many of the floating point errors were enough that
	rounding was not working as expected (at least this is what I think the issue is). This couldn't be solved, as this "bug" is just an inherit problem 
	in computer programming and for computers in general. When attempting to correct for this bug, and same program would just arrise in another part of
	running the "centre_molecule_in_cell" method (in SUMELF). 
	* For this reason, it is easier and more reliable for ECCP to record the move_COM displacement that was applied to centre a dimer as much as possible 
	  to the centre of the origin unit cell. This allows ECCP to be more reliable at rerunning itself of already run Crystals. This is helpful if you 
	  want to easily rerun ECCP on your crystal for different Quanutm Chemistry settings, such as changing from using Gaussian to ORCA (and vise versa), 
	  changing basis sets, changes functionals, etc. 


	Parameters
	----------
	dimer_details : dict.
		This contains information about the recorded dimers. 
	neighbourhood_molecules_for_dimer_method : list
		This list contains all the information about which molecules neighbour each other in the crystal, based on methods given by the dimer_method dictionary. 
	neighbourhood_molecules_for_environment_method : list
		This list contains all the information about which molecules neighbour each other in the crystal, based on methods given by the environment_method dictionary. 

	Returns
	-------
	all_dimers_info : list
		This is a list of all the dimers identified, given as the tuple (index of molecule 1 in molecules list, index of molecule 1 in molecules list, molecule 1, molecule 2). 
	"""

	# First, make sure that the neighbourhood_molecules_for_dimer_method list is sorted by shortest_distance
	neighbourhood_molecules_for_dimer_method.sort(key=lambda x: (x[4], x[0], x[1], x[2][0], x[2][1], x[2][2]))

	# Second, obtain the shortest distances betweeen molecules in dimers from neighbourhood_molecules_for_dimer_method
	#         * This information will be converted to a dictionary where the shortest didstance can be obtained for 
	#           dimer (mol_1, mol_2, UCV_i, UCV_j, UCV_k)
	shortest_distances = get_shortest_distances_dictionary(neighbourhood_molecules_for_dimer_method)

	# Third, initialise a dictionary for holding all the dimers found in this crystal, given the user inputs. 
	all_dimers_info = {}

	# Fourth, create the progress bar for obtaining dimers information from convert_dimer_details_to_all_dimers_info_method.
	pbar = tqdm(sorted(dimer_details.items(), key=lambda x:x[0]), unit='dimer')

	# Fifth, make all the dimers found using the displacements found, and then centre the dimer as close to the centre
	#         of the unit cell as possible while retaining the original periodic positions of the molecules in the dimer. 
	for dimer_name, (mol_name1, mol_name2, unit_cell_displacement_x, unit_cell_displacement_y, unit_cell_displacement_z, displacement_x, displacement_y, displacement_z, move_COM_x, move_COM_y, move_COM_z) in pbar:

		# 5.1: Create the tuple for holding the unit cell displacement vector for moving molecule 2 by.
		unit_cell_displacement = (unit_cell_displacement_x, unit_cell_displacement_y, unit_cell_displacement_z)

		# 5.2: Obtain the associated displacement vector for the unit cell displacement vector for moving molecule 2 by.
		displacement = np.array([displacement_x, displacement_y, displacement_z])

		# 5.3: Get the vector for moving the centre of mass of the dimer by for the saved dimers on file. 
		move_centre_of_mass_by = np.array([move_COM_x, move_COM_y, move_COM_z])

		# 5.4: Get the shortest distance between the atoms of the two molecules in the dimer
		shortest_distance = shortest_distances[(mol_name1, mol_name2, unit_cell_displacement_x, unit_cell_displacement_y, unit_cell_displacement_z)]

		# 5.5: Make a tuple to hold all the information to create the dimer. 
		dimer_pair = (mol_name1, mol_name2, unit_cell_displacement, displacement, move_centre_of_mass_by, shortest_distance)

		# 5.6: Check that dimer_name is not already in the all_dimers_info dictionary
		if dimer_name in all_dimers_info:
			raise Exception(f'Error: dimer_name is already in all_dimers_info: {all_dimers_info}')

		# 5.7: Append the information about how to make the dimer using all_dimers_info
		all_dimers_info[dimer_name] = dimer_pair

	# Sixth, return the list of dimers.
	return all_dimers_info

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

def get_shortest_distances_dictionary(neighbourhood_molecules_for_dimer_method):
	"""
	This method is designed to create a dictionary to allow the user to obtain the shortest distance between atoms in a dimer
	based on the input key (mol_name1, mol_name2, unit_cell_displacement_x, unit_cell_displacement_y, unit_cell_displacement_z) 
	"""

	# First, initialise the dictionary for obtaining the shortest distance between molecules in the selected dimer. 
	shortest_distances = {}

	# Second, for each dimer in neighbourhood_molecules_for_dimer_method.
	for mol_name1, mol_name2, unit_cell_displacement, displacement, shortest_distance in neighbourhood_molecules_for_dimer_method:

		# 2.1: Obtain the ijk values from the unit_cell_displacement tuple.
		unit_cell_displacement_i, unit_cell_displacement_j, unit_cell_displacement_k = unit_cell_displacement

		# 2.2: Obtain the key for this dimer.
		key = (mol_name1, mol_name2, unit_cell_displacement_i, unit_cell_displacement_j, unit_cell_displacement_k)

		# 2.3: Make sure that the key is not already in the shortest_distances dictionary.
		#      * If there is, there is a duplication. 
		if key in shortest_distances.keys():
			to_string  = f'Error: The key {key} is already in shortest_distances\n'
			to_string += f'shortest_distances = {shortest_distances}\n'
			to_string += 'Check this'
			raise Exception(to_string)

		# 2.4: Record shortest_distance in shortest_distances for the key
		shortest_distances[key] = shortest_distance

	# Third, return shortest_distances
	return shortest_distances

# -----------------------------------------------------------------------------------------------------------------------------------------------------
