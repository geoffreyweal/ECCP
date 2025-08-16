"""
get_neighbouring_molecules_about_system.py, Geoffrey Weal, 6/3/23

This script will provide a list of the molecules that neighbour a given system.
"""

from SUMELF import convert_ijk_to_displacement_vector

def get_neighbouring_molecules_about_molecules(environment_settings, neighbourhood_molecules_for_environment_method):
	"""
	This method will provide a list of the molecules that neighbour each molecule in the crystal.

	Parameters
	----------
	environment_settings : dict.
		This dictionary holds information about if you want to obtain the environment about each molecule, for molecule specific calculations. 
	neighbourhood_molecules_for_environment_method : dict.
		This is the list of all the information about the neighbouring molecules in the crystal. 

	Returns
	-------
	neighbouring_molecules_about_molecules : dict:
		This is a list of all the neighbouring molecules surrounding each origin cell molecule in the crystal. 
	"""

	# First, initialise a dictionary that indicates which molecules neighbour the molecules in the origin unit cell.
	neighbouring_molecules_about_molecules = {}

	# Second, determine if you want to include the environment for molecule calculations.
	include_environment_in_molecule_calcs = ('include_environment_in_molecule_calcs' in environment_settings) and environment_settings['include_environment_in_molecule_calcs']
	if not include_environment_in_molecule_calcs:
		return neighbouring_molecules_about_molecules

	raise Exception('Check this method, as have not used it before.')

	# Third, go through the neighbourhood_molecules_for_environment_method and record which molecules neighbour which origin cell molecules. 
	for mol_name1, mol_name2, unit_cell_displacement, displacement, shortest_distance in neighbourhood_molecules_for_environment_method:
		neighbouring_molecules_about_molecules.setdefault(mol_name1,[]).append((mol_name2, unit_cell_displacement, displacement, shortest_distance))

	# Fourth, return neighbouring_molecules_about_molecules
	return neighbouring_molecules_about_molecules

def get_neighbouring_molecules_about_dimers(environment_settings, neighbourhood_molecules_for_dimer_method, neighbourhood_molecules_for_environment_method):
	"""
	This method will provide a list of the molecules that neighbour each dimer in the crystal.

	Parameters
	----------
	environment_settings : dict.
		This dictionary holds information about if you want to obtain the environment about each molecule, for molecule specific calculations. 
	neighbourhood_molecules_for_dimer_method : dict.
		This is the list of all the information required to construct the dimers in the crystal. 
	neighbourhood_molecules_for_environment_method : dict.
		This is the list of all the information about the neighbouring molecules in the crystal. 
	crystal_cell_lattice : np.array
		This is the cell lattice of the crystal.

	Returns
	-------
	neighbouring_molecules_about_dimers : dict:
		This is a list of all the neighbouring molecules surrounding each origin cell molecule in the crystal. 
	"""

	# First, initialise a dictionary that indicates which molecules neighbour the molecules in the origin unit cell.
	neighbouring_molecules_about_dimers = {}

	# Second, determine if you want to include the environment for dimer calculations.
	include_environment_in_dimer_calcs = ('include_environment_in_dimer_calcs' in environment_settings) and environment_settings['include_environment_in_dimer_calcs']
	if not include_environment_in_dimer_calcs:
		return neighbouring_molecules_about_dimers

	raise Exception('Never needed to use this method, Check this before using.')

	# Second, go through the neighbourhood_molecules_for_environment_method and record which molecules neighbour which origin cell molecules. 
	for dimer_name, mol_name1, mol_name2, dimer_unit_cell_displacement, dimer_displacement, dimer_shortest_distance in neighbourhood_molecules_for_dimer_method:

		# 2.1: Obtain the dimer key that describes the details of the dimer. 
		dimer_key = (mol_name1, mol_name2, dimer_unit_cell_displacement)

		# 2.2: Look through each neighbour in the neighbourhood dictionary. 
		for neighbour_mol_name1, neighbour_mol_name2, neighbour_unit_cell_displacement, neighbour_displacement, neighbour_shortest_distance in neighbourhood_molecules_for_environment_method: 

			# 2.2.1: If mol_name1 == neighbour_mol_name1, then neighbour_mol_name2 surrounds mol_name1
			if (mol_name1 == neighbour_mol_name1):

				# 2.2.1.1: Create the tuple which describes the neighbouring molecule we are exploring.
				neighbour_information = (neighbour_mol_name2, neighbour_unit_cell_displacement)

				if (neighbour_mol_name2 == mol_name2) and all([(v1 == v2) for (v1,v2) in zip(neighbour_unit_cell_displacement, dimer_unit_cell_displacement)]):
					# 2.2.1.2: If this statement is true, neighbour mol 2 is dimer mol 2, so we dont want to include this molecule in the dimer neighbours list. 
					pass
				elif (dimer_key in neighbouring_molecules_about_dimers.keys()) and (neighbour_information in neighbouring_molecules_about_dimers[dimer_key].keys()):
					# 2.2.1.3: If this statement is true, we already have included this molecule in the neighbour list, dont't need to double count it.
					pass
				else:
					# 2.2.1.4: Include this molecule with details given in neighbour_information in neighbouring_molecules_about_dimers for this dimer.
					neighbouring_molecules_about_dimers.setdefault(dimer_key,{})[neighbour_information] = neighbour_displacement

			# 2.2.2: If mol_name2 == neighbour_mol_name1, then neighbour_mol_name2 surrounds mol_name2
			if (mol_name2 == neighbour_mol_name1):

				# 2.2.2.1: Create the tuple which describes the neighbouring molecule we are exploring.
				neighbour_of_dm2_unit_cell_displacement = tuple([v1+v2 for (v1,v2) in zip(dimer_unit_cell_displacement, neighbour_unit_cell_displacement)])
				neighbour_information = (neighbour_mol_name2, neighbour_of_dm2_unit_cell_displacement)

				if (neighbour_mol_name2 == mol_name1) and all([(v1 == v2) for (v1,v2) in zip(neighbour_of_dm2_unit_cell_displacement, (0,0,0))]):
					# 2.2.2.2: If the following is true, neighbour mol 2 is dimer mol 1, so we dont want to include this molecule in the dimer neighbours list. 
					# Note: Because dimer1 is in the origin, neighbour 2 is mol_name 1 if dimer_unit_cell_displacement + neighbour_unit_cell_displacement = (0,0,0)
					pass
				elif (dimer_key in neighbouring_molecules_about_dimers.keys()) and (neighbour_information in neighbouring_molecules_about_dimers[dimer_key].keys()):
					# 2.2.2.3: If this statement is true, we already have included this molecule in the neighbour list, dont't need to double count it.
					pass
				else:
					# 2.2.2.4: Include this molecule with details given in neighbour_information in neighbouring_molecules_about_dimers for this dimer.
					#neighbour_of_dm2_displacement_test = convert_ijk_to_displacement_vector(neighbour_of_dm2_unit_cell_displacement, crystal_cell_lattice)
					neighbour_of_dm2_displacement = dimer_displacement + neighbour_displacement
					neighbouring_molecules_about_dimers.setdefault(dimer_key,{})[neighbour_information] = neighbour_of_dm2_displacement

	import pdb; pdb.set_trace()

	# Third, return neighbouring_molecules_about_dimers
	return neighbouring_molecules_about_dimers

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


