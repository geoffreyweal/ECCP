"""
general_methods.py, Geoffrey Weal, 17/2/22

This script includes methods for obtaining dimers.
"""
from tqdm import trange
from SUMELF import centre_molecule_in_cell

def get_dimers(molecules, molecule_graphs, neighbourhood_molecules_for_dimer_method=[], neighbourhood_molecules_for_environment_method=[], no_of_cpus=1):
	"""
	This method will obtain dimers between molecules in the crystal.

	Parameters
	----------
	molecules : list of ase.Atoms objects.
		These are all the individual molecules identified in the crystal that you want to determine dimers for.
	molecule_graphs : list of networkx.graph
		These are the graph associated with each molecule in the molecules list. This graph contains the bonding information about each molecule in the molecules list. 
	neighbourhood_molecules_for_dimer_method : list
		This list contains all the information about which molecules neighbour each other in the crystal, based on methods given by the dimer_method dictionary. 
	neighbourhood_molecules_for_environment_method : list
		This list contains all the information about which molecules neighbour each other in the crystal, based on methods given by the environment_method dictionary. 
	no_of_cpus : int.
		This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.

	Returns
	-------
	dimer_pairs : list
		This is a list of all the dimers identified, given as the tuple (index of molecule 1 in molecules list, index of molecule 1 in molecules list, molecule 1, molecule 2). 
	"""

	# First, make sure that the neighbourhood_molecules_for_dimer_method list is sorted by shortest_distance
	neighbourhood_molecules_for_dimer_method.sort(key=lambda x: (x[4], x[0], x[1], x[2][0], x[2][1], x[2][2]))

	# Second, get the lowest named molecule from the molecules dictionary.
	first_molecule_name = sorted(molecules.keys())[0]

	# Third, get the cell points of the super cell around the origin unit cell with reach 1 
	crystal_cell_lattice = molecules[first_molecule_name].get_cell()

	# Fourth, initialise a dictionary for holding all the dimers found in this crystal, given the user inputs. 
	dimer_pairs = {}

	# Fifth, make all the dimers found using the displacements found, and then centre the dimer as close to the centre
	#        of the unit cell as possible while retaining the original periodic positions of the molecules in the dimer. 
	for dimer_index in trange(len(neighbourhood_molecules_for_dimer_method), unit='dimer'):

		# 5.1: Get information for constructing a dimer from neighbourhood_molecules_for_dimer_method[dimer_index]
		mol_name1, mol_name2, unit_cell_displacement, displacement, shortest_distance = neighbourhood_molecules_for_dimer_method[dimer_index]

		# 5.2: Make a copy of the two molecules to use in the dimer. 
		molecule1 = molecules[mol_name1].copy()
		molecule2 = molecules[mol_name2].copy()

		# 5.3: Translate molecule 2 by the displacement amount
		molecule2.set_positions(molecule2.get_positions() + displacement)

		# 5.4: Make the ase.Atoms object for the dimer based on the two monomers
		dimer_pair_molecules = molecule1 + molecule2

		# 5.5: Determine the displacement needed to move the centre of the dimer to the origin
		move_centre_of_mass_by = centre_molecule_in_cell(dimer_pair_molecules, crystal_cell_lattice, move_molecule=False, dimer_index=dimer_index)

		# 5.6: Make a tuple to hold all the information to create the dimer. 
		dimer_pair = (mol_name1, mol_name2, unit_cell_displacement, displacement, move_centre_of_mass_by, shortest_distance)

		# 5.7: Get the name of the dimer
		dimer_name = dimer_index + 1

		# 5.8: Check that dimer_name is not already in the dimer_pairs dictionary
		if dimer_name in dimer_pairs:
			raise Exception(f'Error: dimer_name is already in dimer_pairs: {dimer_pairs}')

		# 5.9: Append the information about how to make the dimer using dimer_pairs
		dimer_pairs[dimer_name] = dimer_pair

		# 5.10: Update neighbourhood_molecules_for_dimer_method to include the dimers name:
		#neighbourhood_molecules_for_dimer_method[dimer_index] = tuple([dimer_name] + list(neighbourhood_molecules_for_dimer_method[dimer_index]))

	# Sixth, return the list of dimers.
	return dimer_pairs



























	'''
	# Second, order the dimer_pairs list by the shortest distances between molecules in each dimer.
	zipped_lists = zip(shortest_distances, dimer_pairs)
	sorted_zipped_lists = sorted(zipped_lists)
	dimer_pairs = [dimer_pair for _, dimer_pair in sorted_zipped_lists]
	shortest_distances = [shortest_distance for shortest_distance, _ in sorted_zipped_lists]

	# Third, get the cell points of the super cell around the origin unit cell with reach 1 
	crystal_cell_lattice = molecules[0].get_cell()

	# Fourth, indicate if you want to include the the surrounding environment in your dimer calculations where possible
	if 'include_environment_where_possible' in dimer_method:
		include_environment_where_possible = dimer_method['include_environment_where_possible']
	else:
		include_environment_where_possible = False
	if include_environment_where_possible:
		if 'environment_radius' in dimer_method:
			environment_radius = dimer_method['environment_radius']
		else:
			raise Exception('Error, if you have "include_environment_where_possible" = True in your dimer_method, you also need to include a setting for "environment_radius" in the dimer_method dictionary.')
		dimer_pairs_environment = add_environment_to_dimer_pairs(dimer_pairs, molecules, crystal_cell_lattice, environment_radius)

	# Fifth, make all the dimers found using the displacements found, and then centre the dimer as close
	# to the centre of the unit cell as possible while retaining the original periodic positions of the molecules in the dimer. 
	for dimer_index in range(len(dimer_pairs)):
		mol_name1, mol_name2, unit_cell_displacement, displacement = dimer_pairs[dimer_index]
		molecule1 = molecules[mol_name1].copy()
		molecule2 = molecules[mol_name2].copy()
		molecule2.set_positions(molecule2.get_positions() + displacement)
		dimer_pair_molecules = molecule1 + molecule2
		move_centre_of_mass_by = centre_molecule_in_cell(dimer_pair_molecules, crystal_cell_lattice, move_molecule=False)
		dimer_pairs[dimer_index] = (dimer_index, mol_name1, mol_name2, unit_cell_displacement, displacement, move_centre_of_mass_by)

	# Sixth, return the list of dimers.
	return dimer_pairs
	'''







