"""
get_neighbouring_molecules.py, Geoffrey Weal, 2/3/23

This script will evaluate the molecules surrounding each molecule in the unit cell.
"""

from ECCP.ECCP.get_neighbouring_molecules_methods.centre_of_mass_method          import get_neighbours_centre_of_mass_method
from ECCP.ECCP.get_neighbouring_molecules_methods.centre_of_molecule_method      import get_neighbours_centre_of_molecule_method
from ECCP.ECCP.get_neighbouring_molecules_methods.average_distance_method        import get_neighbours_average_distance_method
from ECCP.ECCP.get_neighbouring_molecules_methods.nearest_atoms_method           import get_neighbours_nearest_atoms_method
from SUMELF import centre_molecule_in_cell

def get_neighbouring_molecules(molecules, molecule_graphs, make_dimer_method={'method': 'nearest_atoms_method', 'max_dimer_distance': 8.0}, environment_settings={'environment_radius': 8.0, 'include_environment_in_molecule_calcs': True, 'include_environment_in_dimer_calcs': True}, include_hydrogens_in_neighbour_analysis=False, no_of_cpus=1):
	"""
	This method is designed to obtain all the molecules that neighbour each molecule in the original unit cell within some distance. 

	Parameters
	----------
	molecules : list of ase.Atoms objects.
		These are all the individual molecules identified in the crystal that you want to determine dimers for.
	molecule_graphs : list of networkx.graphs
		These are all the graphs that describe the bonding network in each molecule in the molecules list. 
	make_dimer_method : dict.
		This is the information used to create dimers from the molecules in the crystal.
	environment_settings : dict.
		This is the information used to determine the molecules that surround each molecules in the crystal. 
	include_hydrogens_in_neighbour_analysis: bool. 
		This tag indicates if you want to include hydrogens when accessing neighbours between molecules. Default: False
	no_of_cpus : int.
		This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.

	Returns
	-------
	neighbourhood_molecules_for_dimer_method : list
		This is a list of all the dimers identified, given as the tuple (index of molecule 1 in molecules list, index of molecule 1 in molecules list, molecule 1, molecule 2). 
	neighbourhood_molecules_for_environment_method : list
		This is a list of all the environmental list of neighbouring olecule pairs identified, given as the tuple (index of molecule 1 in molecules list, index of molecule 1 in molecules list, molecule 1, molecule 2). 
	"""

	# First, write some initial strings for any potential error message.
	toString_error_message  = 'Error: The input method for the dimer method can be either:'+'\n'
	toString_error_message += '\t* CMass_method: The Centre of Mass Method'+'\n'
	toString_error_message += '\t* CMol_method: The Centre of Molecule Method'+'\n'
	toString_error_message += '\t* average_distance_method: The Average Distance Method'+'\n'
	toString_error_message += '\t* nearest_atoms_method: The Nearest Atoms Method'+'\n'
	toString_error_message += 'See https://github.com/geoffreyweal/ECCP for more information'+'\n'
	toString_error_message += 'This program will finish without completing'+'\n'

	# Second, check that a method has been given in make_dimer_method
	if not ('method' in make_dimer_method):
		toString_error_message += '\n'
		toString_error_message += 'You have not given an input for the dimer method'+'\n'
		toString_error_message += 'See https://github.com/geoffreyweal/ECCP for more information'+'\n'
		toString_error_message += 'This program will finish without completing'+'\n'
		raise Exception(toString_error_message)

	# Third, obtain information about how to make dimers from molecules in your crystal.
	make_dimer_method_name = make_dimer_method['method']
	if   make_dimer_method_name == 'centre_of_mass':
		max_dimer_distance = make_dimer_method['max_dimer_distance']
	elif make_dimer_method_name == 'centre_of_molecule':
		max_dimer_distance = make_dimer_method['max_dimer_distance']
	elif make_dimer_method_name == 'average_distance_method':
		max_dimer_distance = make_dimer_method['average_distance_method']
	elif make_dimer_method_name == 'nearest_atoms_method':
		max_dimer_distance = make_dimer_method['max_dimer_distance']
	else:
		toString_error_message += '\n'
		toString_error_message += 'The dimer method name you have given is: '+str()+'\n'
		toString_error_message += 'See https://github.com/geoffreyweal/ECCP for more information'+'\n'
		toString_error_message += 'This program will finish without completing'+'\n'
		raise Exception(toString_error_message)

	# Fourth, obtain information about the environmental settings. 
	include_environment_in_molecule_calcs = ('include_environment_in_molecule_calcs' in environment_settings) and environment_settings['include_environment_in_molecule_calcs']
	include_environment_in_dimer_calcs    = ('include_environment_in_dimer_calcs'    in environment_settings) and environment_settings['include_environment_in_dimer_calcs']
	include_environment = include_environment_in_molecule_calcs or include_environment_in_dimer_calcs
	if include_environment:
		max_environment_distance = environment_settings['max_environment_distance']

	# Fifth, obtain the neighbours surrounding the molecules for the dimer method and the environment settings.
	if   make_dimer_method_name == 'centre_of_mass':
		raise Exception('Check if this method includes multiprocessing')
		neighbourhood_molecules_for_dimer_method       = get_neighbours_centre_of_mass_method(molecules, max_dimer_distance, include_hydrogens_in_neighbour_analysis=include_hydrogens_in_neighbour_analysis)
		neighbourhood_molecules_for_environment_method = get_neighbourhood_molecules_for_environment_method(include_environment, molecules, molecule_graphs, max_environment_distance, include_hydrogens_in_neighbour_analysis=include_hydrogens_in_neighbour_analysis, no_of_cpus=no_of_cpus)
	elif make_dimer_method_name == 'centre_of_molecule':
		raise Exception('Check if this method includes multiprocessing')
		neighbourhood_molecules_for_dimer_method       = get_neighbours_centre_of_molecule_method(molecules, max_dimer_distance, include_hydrogens_in_neighbour_analysis=include_hydrogens_in_neighbour_analysis)
		neighbourhood_molecules_for_environment_method = get_neighbourhood_molecules_for_environment_method(include_environment, molecules, molecule_graphs, max_environment_distance, include_hydrogens_in_neighbour_analysis=include_hydrogens_in_neighbour_analysis, no_of_cpus=no_of_cpus)
	elif make_dimer_method_name == 'average_distance_method':
		raise Exception('Check if this method includes multiprocessing')
		neighbourhood_molecules_for_dimer_method       = get_neighbours_average_distance_method(molecules, max_dimer_distance, include_hydrogens_in_neighbour_analysis=include_hydrogens_in_neighbour_analysis)
		neighbourhood_molecules_for_environment_method = get_neighbourhood_molecules_for_environment_method(include_environment, molecules, molecule_graphs, max_environment_distance, include_hydrogens_in_neighbour_analysis=include_hydrogens_in_neighbour_analysis, no_of_cpus=no_of_cpus)
	elif make_dimer_method_name == 'nearest_atoms_method':
		if include_environment:
			max_distance = max([max_dimer_distance, max_environment_distance])
			neighbourhood_molecules_for_dimer_and_environment_method                                 = get_neighbours_nearest_atoms_method(molecules, molecule_graphs, max_distance, include_hydrogens_in_neighbour_analysis=include_hydrogens_in_neighbour_analysis, no_of_cpus=no_of_cpus)
			neighbourhood_molecules_for_dimer_method, neighbourhood_molecules_for_environment_method = create_lists_for_dimer_and_environment_methods(neighbourhood_molecules_for_dimer_and_environment_method, max_dimer_distance, max_environment_distance)
		else:
			neighbourhood_molecules_for_dimer_method       = get_neighbours_nearest_atoms_method(molecules, molecule_graphs, max_dimer_distance, include_hydrogens_in_neighbour_analysis=include_hydrogens_in_neighbour_analysis, no_of_cpus=no_of_cpus)
			neighbourhood_molecules_for_environment_method = []

	# Sixth, sort the dimers in neighbourhood_molecules_for_dimer_method. 
	#        Sorting format: (distance between molecules, mol1_index, mol2_index, unit_cell_i, unit_cell_j, unit_cell_k)
	neighbourhood_molecules_for_dimer_method.sort(key=lambda x: (x[4], x[0], x[1], x[2][0], x[2][1], x[2][2]))

	# Seventh, return neighbourhood_molecules_for_dimer_method and neighbourhood_molecules_for_environment_method
	return neighbourhood_molecules_for_dimer_method, neighbourhood_molecules_for_environment_method

# ====================================================================================================================================================================================

def get_neighbourhood_molecules_for_environment_method(include_environment, molecules, molecule_graphs, max_environment_distance, include_hydrogens_in_neighbour_analysis, no_of_cpus):
	"""
	This method is designed to provide neighbourhood_molecules for the environment method
	"""
	if include_environment:
		neighbourhood_molecules_for_environment_method = get_neighbours_nearest_atoms_method(molecules, molecule_graphs, max_environment_distance, include_hydrogens_in_neighbour_analysis=include_hydrogens_in_neighbour_analysis, no_of_cpus=no_of_cpus)
	else:
		neighbourhood_molecules_for_environment_method = []
	return neighbourhood_molecules_for_environment_method

def create_lists_for_dimer_and_environment_methods(neighbourhood_molecules_for_dimer_and_environment_method, max_dimer_distance, max_environment_distance):
	"""
	This method is designed to create two lists from neighbourhood_molecules_for_dimer_and_environment_method, one for the dimer method and one for environment settings

	Parameters
	----------
	neighbourhood_molecules_for_dimer_and_environment_method : list
		This is the list of neighbourhood_molecules information. 
	max_dimer_distance : float
		This is the maximum distance that two molecules can be from each other to be considered dimers. 
	max_environment_distance : float
		This is the maximum distance that two molecules are within the environment of each other. 

	Returns
	-------
	neighbourhood_molecules_for_dimer_method : list
		This is the list of neighbourhood_molecules information to be used by a dimer method.
	neighbourhood_molecules_for_environment_method : list
		This is the list of neighbourhood_molecules information to be used for environment settings.
	"""

	# First, initiate the lists for storing neighbourhood_molecules info for dimer and environmental methods
	neighbourhood_molecules_for_dimer_method = []
	neighbourhood_molecules_for_environment_method = []

	# Second, copy the information to neighbourhood_molecules_for_dimer_method and neighbourhood_molecules_for_environment_method
	# based on the values of max_dimer_distance and max_environment_distance
	for index1, index2, unit_cell_displacement, displacement, shortest_distance in sorted(neighbourhood_molecules_for_dimer_and_environment_method,key=lambda x:x[4]):
		neighbourhood_molecules = (index1, index2, unit_cell_displacement, displacement, shortest_distance)
		if shortest_distance <= max_dimer_distance:
			neighbourhood_molecules_for_dimer_method.append(neighbourhood_molecules)
		if shortest_distance <= max_environment_distance:
			neighbourhood_molecules_for_environment_method.append(neighbourhood_molecules)

	# Third, return neighbourhood_molecules_for_dimer_method and neighbourhood_molecules_for_environment_method
	return neighbourhood_molecules_for_dimer_method, neighbourhood_molecules_for_environment_method

# ====================================================================================================================================================================================

