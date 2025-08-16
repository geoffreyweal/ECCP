"""
general_methods.py, Geoffrey Weal, 17/2/22

This script includes methods for obtaining neighbours by looking for molecules with non-hydrogen atoms within the vicinity of each other.
"""
import sys
from tqdm import tqdm

import multiprocessing as mp
from tqdm.contrib.concurrent import process_map

from SUMELF import get_distance, make_folder, remove_folder
from ECCP.ECCP.get_neighbouring_molecules_methods.Neighbourhood_Generator           import Neighbourhood_Generator
from ECCP.ECCP.get_neighbouring_molecules_methods.Neighbourhood_Generator_Multi_CPU import Neighbourhood_Generator_Multi_CPU

def get_neighbours_nearest_atoms_method(molecules, molecule_graphs, max_distance, include_hydrogens_in_neighbour_analysis=False, no_of_cpus=1):
	"""
	This method will obtain the molecules in the neighbourhood of molecules in the origin unit cell. 

	Parameters
	----------
	molecules : list of ase.Atoms objects.
		These are all the individual molecules identified in the crystal that you want to determine neighbours for.
	molecule_graphs : list of networkx.graph
		This is a list of all the networkx graphs that describe the bonding system for each molecule. 
	max_distance : float.
		This is the maximum distance that atoms in two molecules can be within each other for the two molecules to be considered neighbouring. Given in Å. 
	include_hydrogens_in_neighbour_analysis: bool. 
		This tag indicates if you want to include hydrogens when accessing neighbours between molecules. Default: False
	no_of_cpus : int.
		This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.

	Returns
	-------
	neighbourhood_molecules_info : list
		This is a list of all the molecules that neighbour each other within the vicinity given by max_distance.
	"""

	# First, use either the single or multi-cpu method for obtaining neighbouring molecules.
	if no_of_cpus == 1:
		neighbourhood_molecules_info = obtain_neighbours_with_single_cpu(molecules, molecule_graphs, max_distance)
	else:
		neighbourhood_molecules_info = obtain_neighbours_with_multi_cpu (molecules, molecule_graphs, max_distance, no_of_cpus)

	# Second, return the list of the information of neighbouring molecules in this crystal
	return neighbourhood_molecules_info

# ===============================================================================================================
# ===============================================================================================================
# ===============================================================================================================

def obtain_neighbours_with_multi_cpu(molecules, molecule_graphs, max_distance, no_of_cpus=1):
	"""
	This method will obtain the molecules in the neighbourhood of molecules in the origin unit cell. 

	Parameters
	----------
	molecules : list of ase.Atoms objects.
		These are all the individual molecules identified in the crystal that you want to determine neighbours for.
	molecule_graphs : list of networkx.graph
		This is a list of all the networkx graphs that describe the bonding system for each molecule. 
	max_distance : float.
		This is the maximum distance that atoms in two molecules can be within each other for the two molecules to be considered neighbouring. Given in Å. 
	no_of_cpus : int.
		This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.

	Returns
	-------
	neighbourhood_molecules_info : list
		This is a list of all the molecules that neighbour each other within the vicinity given by max_distance.
	"""

	# First, create the manager to save lists to
	with mp.Manager() as manager:

		# Second, create the list to collect information on the neighbouring molecules in the crystal. 
		neighbourhood_molecules_info = manager.list()

		# Third, obtain the number of neighbourhood sets.
		no_of_neighbourhood_sets = int((len(molecules) * (len(molecules) + 1)) / 2)

		# Fourth, write warning message to the user.
		print('Obtaining neighbourhoods between molecules (Please wait until after 100%, as the process will still be running.)', file=sys.stderr)

		# Fifth, obtain the input value generator. 
		input_values = tqdm(get_inputs(molecules, molecule_graphs, max_distance, neighbourhood_molecules_info), total=no_of_neighbourhood_sets, desc='Obtaining neighbouring pairs of molecules', unit='calc')

		# Sixth, run the multiprocessing jobs.
		#process_map(obtain_neighbours_method_for_multi_cpu, get_inputs(molecules, molecule_graphs, max_distance, neighbourhood_molecules_info), total=no_of_neighbourhood_sets, desc='Obtaining neighbouring pairs of molecules', unit='calc', max_workers=no_of_cpus)
		pool = mp.Pool(processes=no_of_cpus)
		pool.map_async(obtain_neighbours_method_for_multi_cpu, input_values)
		pool.close()
		pool.join()

		# Seventh, save the list from a multiprocessing Manager list to a regular list
		neighbourhood_molecules_info = list(neighbourhood_molecules_info)
		
	# Eighth, return neighbourhood_molecules_info
	return neighbourhood_molecules_info

def get_inputs(molecules, molecule_graphs, max_distance, neighbourhood_molecules_info):
	"""
	This method a generator designed to obtain the inputs for obtaining neighbours with multiple CPUs.

	Parameters
	----------
	molecules : list of ase.Atoms objects.
		These are all the individual molecules identified in the crystal that you want to determine neighbours for.
	molecule_graphs : list of networkx.graph
		This is a list of all the networkx graphs that describe the bonding system for each molecule. 
	max_distance : float.
		This is the maximum distance that atoms in two molecules can be within each other for the two molecules to be considered neighbouring. Given in Å. 
	neighbourhood_molecules_info : list
		This is a list of all the molecules that neighbour each other within the vicinity given by max_distance.
	"""

	# First, obtain the names of the molecules.
	mol_names = sorted(molecules.keys())

	# Second, get the first molecule/monomer. 
	for index1 in range(len(mol_names)):

		# Third, obtain the molecule name, the molecule, and the associated graph of molecule mol_name1
		mol_name1       = mol_names[index1]
		molecule1       = molecules[mol_name1]
		molecule_graph1 = molecule_graphs[mol_name1]

		# Fourth, get the second molecule/monomer. This could be the same molecule as index1, but will be displaced to a different position.
		for index2 in range(index1,len(molecules)): 

			# Fifth, obtain the molecule name, the molecule, and the associated graph of molecule mol_name2
			mol_name2       = mol_names[index2]
			molecule2       = molecules[mol_name2]
			molecule_graph2 = molecule_graphs[mol_name2]

			# Sixth, yield the inputs.
			yield (mol_name1, mol_name2, molecule1, molecule2, molecule_graph1, molecule_graph2, max_distance, neighbourhood_molecules_info)

def obtain_neighbours_method_for_multi_cpu(input_variables):
	"""
	This method will obtain neighbourhood information between molecules in the crystal based on the distances the closest atoms in each molecule. 

	Parameters
	----------
	molecules : list of ase.Atoms objects.
		These are all the individual molecules identified in the crystal that you want to determine neighbours for.
	molecule_graphs : list of networkx.graph
		This is a list of all the networkx graphs that describe the bonding system for each molecule. 
	max_distance : float.
		This is the maximum distance that atoms in two molecules can be within each other for the two molecules to be considered neighbouring. Given in Å. 
	no_of_cpus : int.
		This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.
	neighbourhood_molecules_info : list
		This is a list of all the molecules that neighbour each other within the vicinity given by max_distance.
	"""

	# First, obtain the variables for processing from the input_variables list.
	mol_name1, mol_name2, molecule1, molecule2, molecule_graph1, molecule_graph2, max_distance, neighbourhood_molecules_info = input_variables

	# Second, create the neighbourhood generator that will create all the neighbouring pairs between molecules that could exist in the crystal.
	neighbourhood_generator_object = Neighbourhood_Generator_Multi_CPU(mol_name1, mol_name2, molecule1, molecule2, molecule_graph1, molecule_graph2, molecule1.get_cell())
	neighbourhood_generator = neighbourhood_generator_object.generator()

	# Third, obtain the list of molecules that are neighbours. These are molecules that are within max_distance distance of each other. 
	obtain_neighbours(neighbourhood_generator, max_distance, 'Neighbourhood_Generator_Multi_CPU', neighbourhood_molecules_info)
			
# ===============================================================================================================
# ===============================================================================================================
# ===============================================================================================================

def obtain_neighbours_with_single_cpu(molecules, molecule_graphs, max_distance):
	"""
	This method will obtain neighbouring pairs of molecules in the crystal based on the distances the closest atoms in each molecule. 

	Parameters
	----------
	molecules : list of ase.Atoms objects.
		These are all the individual molecules identified in the crystal that you want to determine neighbours for.
	molecule_graphs : list of networkx.graph
		This is a list of all the networkx graphs that describe the bonding system for each molecule. 
	max_distance : float.
		This is the maximum distance that atoms in two molecules can be within each other for the two molecules to be considered neighbouring. Given in Å. 

	Returns
	-------
	neighbourhood_molecules_info : list
		This is a list of all the molecules that neighbour each other within the vicinity given by max_distance.
	"""

	# First, get the cell for a molecule in the crystal
	molecule_cell = molecules[tuple(molecules.keys())[0]].get_cell()

	# Second, create the neighbourhood generator that will create all the neighbouring pairs that could exist in the crystal.
	neighbourhood_generator_object = Neighbourhood_Generator(molecules, molecule_graphs, molecule_cell)
	neighbourhood_generator = neighbourhood_generator_object.generator()

	# Third, initialise a list for holding information about neighbouring molecules. 
	neighbourhood_molecules_info = []

	# Fourth, obtain the list of molecules that are neighbours. These are molecules that are within max_distance distance of each other. 
	obtain_neighbours(neighbourhood_generator, max_distance, 'Neighbourhood_Generator', neighbourhood_molecules_info)

	# Fifth, return the list of neighbouring pairs of molecules.
	return neighbourhood_molecules_info

# ===============================================================================================================
# ===============================================================================================================
# ===============================================================================================================

def obtain_neighbours(neighbourhood_generator, max_distance, generator_type, neighbourhood_molecules_info):
	"""
	This is the main method section that is designed to obtain the molecules that neighbour each other. 

	This means they are within max_distance of each other (excluding hydrogen atoms). 

	Parameters
	----------
	neighbourhood_generator : generator
		This is the generator you want to use for obtaining the molecules you want to check if they are neighbours or not. 
	max_distance : float
		This is the maximum distance between molecules to be cosidered neighbours. 
	generator_type : str.
		This is the name of the generator being used. This should be either 'Neighbourhood_Generator' or 'Neighbourhood_Generator_Multi_CPU'. 
	neighbourhood_molecules_info : list
		This list is for recording information about which molecules neighbour each other in the crystal. 
	"""

	# First, identify neighbouring pairs of molecules based on if there are any atoms in each molecule that are within max_distance distance of each other
	for mol_name1, mol_name2, positions1, positions2, displacement, unit_cell_displacement in neighbourhood_generator:
		
		# Second, if any atom between each molecule is within max_distance, you have a neighbouring pair. 
		are_molecules_within_max_neighbour_distance, shortest_distance = determine_if_molecules_within_max_distance(positions1, positions2, displacement, max_distance)
		if are_molecules_within_max_neighbour_distance:
			neighbourhood_molecules_info.append((mol_name1, mol_name2, unit_cell_displacement, displacement, shortest_distance))

		# Third, send the result of if the neighbouring pair of molecules was accepted or not back to the generator.
		end_of_for_loop_check = neighbourhood_generator.send(are_molecules_within_max_neighbour_distance)
		if not (end_of_for_loop_check == 'Go to get_neighbours method for loop'):
			raise Exception(f'Communication error of {generator_type} generator with this for loop.')

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def determine_if_molecules_within_max_distance(positions1, positions2, displacement, max_distance):
	"""
	This method will determine if two molecules are within max_distance of each other. 

	Parameters
	----------
	positions1 : numpy.array
		These are the positions of atoms in molecule 1.
	positions2 : numpy.array
		These are the positions of atoms in molecule 2.
	displacement : numpy.array
		This is the displacement to move molecule 2 by.
	max_distance : float.
		This is the maximum distance that atoms in two molecules can be within each other for the two molecules to be considered neighbouring. Given in Å. 

	Returns
	-------
	True if these two molecules are within max_distance of each other. False if not. 
	"""

	# First, determine the shortest distance between the two molecules.
	shortest_distance = float('inf')
	for pos_index1 in range(len(positions1)):
		for pos_index2 in range(len(positions2)):
			distance = round(get_distance(positions1[pos_index1], positions2[pos_index2] + displacement), 4)
			# If we have found a distance that is less than max_distance, then these two molcules are a neighbouring pair.
			if (distance <= shortest_distance):
				shortest_distance = distance

	# Second, determine if the distance is less than the maximum distance
	return (shortest_distance <= max_distance), shortest_distance			

# ===============================================================================================================
# ===============================================================================================================
# ===============================================================================================================

