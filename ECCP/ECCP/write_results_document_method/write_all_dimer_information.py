"""
write_all_dimer_information.py, Geoffrey Weal, 22/5/24

Obtain the information of all the molecules, where indices have been updated to names
"""
title_components = ['Dimer No', 'First Molecule', 'Second Molecule', 'UCV(i)', 'UCV(j)', 'UCV(k)', 'DV(x)', 'DV(y)', 'DV(z)', 'move_dimer_COM(x)', 'move_dimer_COM(y)', 'move_dimer_COM(z)']

no_of_decimal_places = 12
def write_all_dimer_information(path_to_eccp_folder, all_dimers_info):
	"""
	This method will write all the positional information about the molecules relationships to each other in each dimer obtained about each molecule type in the crystal. 

	Parameters
	----------
	path_to_eccp_folder : str.
		This is the path to this crystal in the ECCP folder
	all_dimers_info : list 
		This is all the information about all the dimers that have been idenified in the crystal using the settings as given by the user.
	"""

	# First, initialise a list for recording the information about the dimers, and include title_components in it
	dimer_information = [title_components]

	# Second, initialise the maximum length of the title components across all dimers recorded.
	max_len_of_title_components = [len(title_component) for title_component in title_components]

	# Third, for each dimer in all_dimers_info:
	for dimer_no, (mol1_name, mol2_name, unit_cell_displacement, displacement, move_centre_of_mass_by, shortest_distance) in sorted(all_dimers_info.items(), key=lambda x: x[0]):

		# 3.1: Initialise a list for recording the information about this dimer. 
		current_dimer_information = []

		# 3.2: Record the name of the dimer.
		current_dimer_information.append(str(dimer_no))

		# 3.3: Record the name of the molecules that make up this dimer.
		current_dimer_information.append(str(mol1_name))
		current_dimer_information.append(str(mol2_name))

		# 3.4: Record the unit cell ijk values for displacing the second molecule to make dimer.
		current_dimer_information.append(str(round(unit_cell_displacement[0], no_of_decimal_places)))
		current_dimer_information.append(str(round(unit_cell_displacement[1], no_of_decimal_places)))
		current_dimer_information.append(str(round(unit_cell_displacement[2], no_of_decimal_places)))

		# 3.5: Record the corresponding unit cell ijk displacement values (in Angstroms) for displacing the second molecule to make dimer.
		current_dimer_information.append(str(round(displacement[0], no_of_decimal_places)))
		current_dimer_information.append(str(round(displacement[1], no_of_decimal_places)))
		current_dimer_information.append(str(round(displacement[2], no_of_decimal_places)))

		# 3.6: Record the displacement for moving the centre of mass by. 
		current_dimer_information.append(str(round(move_centre_of_mass_by[0], no_of_decimal_places)))
		current_dimer_information.append(str(round(move_centre_of_mass_by[1], no_of_decimal_places)))
		current_dimer_information.append(str(round(move_centre_of_mass_by[2], no_of_decimal_places)))

		# 3.7: Record the information about this dimer in dimer_information
		dimer_information.append(current_dimer_information)

		# 3.8: Update the maximum string length for each title component in max_len_of_title_components
		for index in range(len(list(zip(current_dimer_information, max_len_of_title_components)))):
			if len(current_dimer_information[index]) > max_len_of_title_components[index]:
				max_len_of_title_components[index] = len(current_dimer_information[index])

	# Fourth, obtain the largest length amongst the UCV values
	maximum_UCV = max(max_len_of_title_components[3:6])
	max_len_of_title_components[3] = maximum_UCV
	max_len_of_title_components[4] = maximum_UCV
	max_len_of_title_components[5] = maximum_UCV

	# Fifth, obtain the largest length amongst the UCV values
	maximum_DV = max(max_len_of_title_components[6:9])
	max_len_of_title_components[6] = maximum_DV
	max_len_of_title_components[7] = maximum_DV
	max_len_of_title_components[8] = maximum_DV

	# Sixth, obtain the largest length amongst the UCV values
	maximum_move_COM = max(max_len_of_title_components[9:12])
	max_len_of_title_components[9]  = maximum_move_COM
	max_len_of_title_components[10] = maximum_move_COM
	max_len_of_title_components[11] = maximum_move_COM

	# Seventh, initialise a list for recording the formatted dimer information. 
	dimer_information_to_string = []

	# Eighth, for each dimer in dimer_information:
	for current_dimer_components in dimer_information:

		# 8.1: Initalise a string to record information about the dimer to.
		current_dimer_information_to_string = ''

		# 8.2: For each component of the information given for the dimer. 
		for info_component, max_len_of_title_component in zip(current_dimer_components, max_len_of_title_components):

			# 8.2.1: Determine the number of spaces to add before the string for this component. 
			no_of_spaces = max_len_of_title_component - len(info_component)

			# 8.2.2: Create the string for this component, including spaces. 
			dimer_information_component_to_string = ' '*4 + ' '*(no_of_spaces) + str(info_component)

			# 8.2.3: Add the string for this component to current_dimer_information_to_string
			current_dimer_information_to_string += dimer_information_component_to_string

		# 8.3: Add the string for this dimer to dimer_information_to_string
		dimer_information_to_string.append(current_dimer_information_to_string)

	# Ninth, print all the dimer information to text. 
	data_filename = 'All_Dimer_Information.txt'
	with open(path_to_eccp_folder+'/'+data_filename,'w') as Dimer_Information_TXT:
		Dimer_Information_TXT.write('\n'.join(dimer_information_to_string))

# ---------------------------------------------------------------------------------------------------------------

