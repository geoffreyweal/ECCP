"""
check_dimer_duplication.py, Geoffrey Weal, 10/3/2024

This method will check if there are any dimers that have been entered in twice into the neighbourhood_molecules_for_dimer_method list
"""

def check_dimer_duplication(neighbourhood_molecules_for_dimer_method):
	"""
	This method will check if there are any dimers that have been entered in twice into the neighbourhood_molecules_for_dimer_method list

	Parameters
	----------
	neighbourhood_molecules_for_dimer_method : list
		This is the information from dimers_details in the format of (mol1_name, mol2_name, (UCV(i), UCV(j), UCV(k)), np.array(DV(x), DV(y), DV(z))))
	"""

	# First, initialise a list for recording if there are duplicate dimers in neighbourhood_molecules_for_dimer_method
	duplicate_dimers = []

	# Second, for each dimer in neighbourhood_molecules_for_dimer_method
	for index1 in range(len(neighbourhood_molecules_for_dimer_method)):

		# Third, obtain the components from neighbourhood_molecules_for_dimer_method[index1]
		d1_mol1_index, d1_mol2_index, d1_UC_ijk, d1_displacement_of_m2, d1_distance = neighbourhood_molecules_for_dimer_method[index1]

		# Fourth, for each other dimer in neighbourhood_molecules_for_dimer_method
		for index2 in range(index1+1,len(neighbourhood_molecules_for_dimer_method)):

			# Fifth, obtain the components from neighbourhood_molecules_for_dimer_method[index2]
			d2_mol1_index, d2_mol2_index, d2_UC_ijk, d2_displacement_of_m2, d2_distance = neighbourhood_molecules_for_dimer_method[index2]

			# Sixth, make sure that the components of d1_UC_ijk are all floats
			d1_UC_ijk = tuple([int(value) for value in d1_UC_ijk])
			d2_UC_ijk = tuple([int(value) for value in d2_UC_ijk])

			# Seventh, perform comparisons between both dimers
			same_mol1     = d1_mol1_index        == d2_mol1_index
			same_mol2     = d1_mol2_index        == d2_mol2_index
			same_UC_ijk   = d1_UC_ijk            == d2_UC_ijk
			#same_distance = round(d1_distance,4) == round(d2_distance,4)

			# Eighth, if all of these are the same, we have a duplicate dimers
			if same_mol1 and same_mol2 and same_UC_ijk:
				duplicate_dimers.append((d1_mol1_index, d1_mol2_index, d1_UC_ijk))

	# Ninth, report any dimers that are duplicates
	if len(duplicate_dimers) > 0:
		to_string  = 'Error: There are duplicate dimer in the neighbourhood_molecules_for_dimer_method list.\n'
		to_string += 'Duplicate dimers in neighbourhood_molecules_for_dimer_method: '+str(duplicate_dimers)+'\n'
		to_string += 'Check this.'
		raise Exception(to_string)