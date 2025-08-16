"""
get_symmetric_dimer_pairs_MEA.py, Geoffrey Weal, 11/2/24

This method is designed to determine which dimers are spatially symmetric to each other using the minimal elemental abundance invariance method.
"""
import sys, warnings
from tqdm import tqdm

import multiprocessing as mp

from ECCP.ECCP.get_unique_dimers_methods.invariance_methods.minimal_elemental_abundance_invariance_utility_methods.are_dimers_variant_MEA                import are_dimers_variant_MEA
from ECCP.ECCP.get_unique_dimers_methods.invariance_methods.minimal_elemental_abundance_invariance_utility_methods.are_dimers_variant_from_comprehensive import are_dimers_variant_from_comprehensive

def get_symmetric_dimer_pairs_MEA(dimers, non_hydrogen_molecules_elements, non_hydrogen_molecules_positions, all_no_of_H_on_atoms_in_molecule, max_distance_disparity, neighbouring_molecules_about_dimers, non_hydrogen_molecules, non_hydrogen_graphs, no_of_cpus=1):
	"""
	This method is designed to determine which dimers are spatially symmetric to each other using the minimal elemental abundance invariance method.

	Parameters
	----------
	dimers : dict
	    This is a dict of dimer information, given as {name of dimer: (name of molecule 1 in dimer, name of molecule 2 in dimer, unit cell ijk displacement of molecule 2, unit cell displacement of molecule 2, displacement of dimer COM)}
	non_hydrogen_molecules_elements : dict. of lists of str. 
	    These are the element of all the molecules that can make up the dimers.
	non_hydrogen_molecules_positions :  dict. of 2d numpy.arrays
	    These are the positions of all the molecules that can make up the dimers.
	all_no_of_H_on_atoms_in_molecule :  dict. of lists of int. 
	    These are the number of hydrogens bound to each "heavy" atom in each molecule that can make up the dimers. 
	max_distance_disparity : float
	    This is the maximum difference in relative positions of atoms between dimer 1 and dimer 2 to be considered equivalent/identical
	neighbouring_molecules_about_dimers : dict.
	    This is the information about the molecules that surround (in the vicinity of) each dimer in the crystal.
	non_hydrogen_molecules : dict. of ase.Atoms
	    These are Atoms of objects of the molecules in the crystal. These do not contain hydrogen atoms. This is for debugging use only. 
	non_hydrogen_graphs : dict of networkx.Graphs
		These are the graphs of the associated molecule ase.Atoms objects in non_hydrogen_molecules. The hydrogens in this graph has been removed as nodes and instead appended to "heavy" atoms in the molecule as node variables. This is used for the are_dimers_variant_from_comp method. 
	no_of_cpus : int.
	    This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.

	Returns
	-------
	symmetric_dimer_pairs : list object (either a list or a multiprocessing.Manager.list)
	    This list is designed to store all the symmetric dimers identified from this program.
	"""

	# Prestep: Set the number of cores for each comparison of two molecules
	no_of_cpus_for_comparing_two_dimers = no_of_cpus

	# First, write a message to tell the user what is being done. 
	print('Comparing translational, rotational, and reflective invarience between dimers. '+str(len(dimers))+' dimers to be examined. This can take a while with large and complex dimers.')

	# Second, obtain the total number of dimer comparisons that will be performed.
	nn = int((len(dimers)*(len(dimers)-1))/2)

	# Second, determine which pairs of dimers are symmetric.
	if True: #no_of_cpus == 1:

		# 2.1: Initialise the list to store information about symmetric dimers in. 
		symmetric_dimer_pairs = []

		# 2.2: Create a progress bar for running this task
		pbar = tqdm(get_inputs(dimers, non_hydrogen_molecules_elements, non_hydrogen_molecules_positions, all_no_of_H_on_atoms_in_molecule, max_distance_disparity, neighbouring_molecules_about_dimers, non_hydrogen_molecules, non_hydrogen_graphs, symmetric_dimer_pairs, no_of_cpus_for_comparing_two_dimers), total=nn, unit='dimer pair')

		# 2.3: For each comparison of dimers.
		for input_data in pbar:

			# 2.3.1: Obtain the indices of the two dimers being compared. 
			dimer1_name = input_data[0]
			dimer2_name = input_data[1]

			# 2.3.2: Write a message to tell the user what two dimers are being compared. 
			pbar.set_description('Comparing dimers '+str(dimer1_name)+' and '+str(dimer2_name))

			# 2.3.3: Compare the two dimers. 
			compare_if_two_dimers_are_symmetric_single_process(input_data)

		# 2.4: Close the progress bar. 
		pbar.close()

	else:

		# 2.5: Create the manager to save lists to
		with mp.Manager() as manager:

			# 2.6: Initialise the list to store information about symmetric dimers in. 
			symmetric_dimer_pairs = manager.list()

			# 2.7: Write a message for the user
			print('Comparing dimers (Please wait until after 100%, as the process will still be running.)', file=sys.stderr)

			# 2.8: Perform compare_if_two_dimers_are_symmetric_single_process using multiprocessing
			#process_map(  compare_if_two_dimers_are_symmetric_single_process,      get_inputs(dimers, non_hydrogen_molecules_elements, non_hydrogen_molecules_positions, all_no_of_H_on_atoms_in_molecule, max_distance_disparity, neighbouring_molecules_about_dimers, non_hydrogen_molecules, non_hydrogen_graphs, symmetric_dimer_pairs), total=nn, unit='dimer pair', desc='Comparing dimers', max_workers=no_of_cpus)
			pool = mp.Pool(processes=no_of_cpus)
			pool.map_async(compare_if_two_dimers_are_symmetric_single_process, tqdm(get_inputs(dimers, non_hydrogen_molecules_elements, non_hydrogen_molecules_positions, all_no_of_H_on_atoms_in_molecule, max_distance_disparity, neighbouring_molecules_about_dimers, non_hydrogen_molecules, non_hydrogen_graphs, symmetric_dimer_pairs, no_of_cpus_for_comparing_two_dimers), total=nn, unit='dimer pair', desc='Comparing dimers'))
			pool.close()
			pool.join()

			# 2.9: Convert symmetric_dimer_pairs from a manager.list object to a list object.
			symmetric_dimer_pairs = list(symmetric_dimer_pairs)

	# Third, sort the dimers.
	symmetric_dimer_pairs.sort()

	# Fourth, return symmetric_dimer_pairs
	return symmetric_dimer_pairs

# --------------------------------------------------------------------------------------------------------------

shortest_distance_between_dimers_max_threshold_limit = 0.1
def compare_if_two_dimers_are_symmetric_single_process(input_data):
	"""
	This method is designed to determine if two diemrs are symmetric to each other. 

	Parameters
	----------
	dimer1_name : int.
		This is the name of the first dimer being compared.
	dimer2_name : int.
		This is the name of the second dimer being compared.

	d1_m1_original_elements : list of str.
		These are the elements of the atoms in molecule 1 of the first dimer being compared. 
	d1_m1_original_positions : Nx3 numpy.array
		These are the positions of the atoms in molecule 1 of the first dimer being compared. 
	d1_m1_no_H_attached_to_nonH_atoms : list of int.
		These are the number of hydrogens bound to "heavy" atoms in molecule 1 of the first dimer being compared. 

	d1_m2_original_elements : list of str.
		These are the elements of the atoms in molecule 2 of the first dimer being compared. 
	d1_m2_original_positions : Nx3 numpy.array
		These are the positions of the atoms in molecule 2 of the first dimer being compared. 
	d1_m2_no_H_attached_to_nonH_atoms : list of int.
		These are the number of hydrogens bound to "heavy" atoms in molecule 2 of the first dimer being compared. 

	d2_m1_original_elements : list of str.
		These are the elements of the atoms in molecule 1 of the second dimer being compared. 
	d2_m1_original_positions : Nx3 numpy.array
		These are the positions of the atoms in molecule 1 of the second dimer being compared. 
	d2_m1_no_H_attached_to_nonH_atoms : list of int.
		These are the number of hydrogens bound to "heavy" atoms in molecule 1 of the second dimer being compared. 

	d2_m2_original_elements : list of str.
		These are the elements of the atoms in molecule 2 of the second dimer being compared. 
	d2_m2_original_positions : Nx3 numpy.array
		These are the positions of the atoms in molecule 2 of the second dimer being compared. 
	d2_m2_no_H_attached_to_nonH_atoms : list of int.
		These are the number of hydrogens bound to "heavy" atoms in molecule 2 of the second dimer being compared. 

	d1_m1_name : int
		This is the name of molecule 1 that makes up the first dimer. 
	d1_m2_name : int
		This is the name of molecule 2 that makes up the first dimer. 
	unit_cell_disp1 : (int, int, int)
		This is the ijk unit cell displacement to move molecule 2 by with respect to molecule 1 to create the first dimer. 

	d2_m1_name : int
		This is the name of molecule 1 that makes up the second dimer. 
	d2_m2_name : int
		This is the name of molecule 2 that makes up the second dimer. 
	unit_cell_disp2 : (int, int, int)
		This is the ijk unit cell displacement to move molecule 2 by with respect to molecule 1 to create the second dimer. 

	max_distance_disparity : float
	    This is the maximum difference in relative positions of atoms between dimer 1 and dimer 2 to be considered equivalent/identical
	neighbouring_molecules_about_dimers : dict.
	    This is the information about the molecules that surround (in the vicinity of) each dimer in the crystal.
	non_hydrogen_molecules : dict. of ase.Atoms
	    These are Atoms of objects of the molecules in the crystal. These do not contain hydrogen atoms. This is for debugging use only. 
	non_hydrogen_graphs : dict of networkx.Graphs
		These are the graphs of the associated molecule ase.Atoms objects in non_hydrogen_molecules. The hydrogens in this graph has been removed as nodes and instead appended to "heavy" atoms in the molecule as node variables. This is used for the are_dimers_variant_from_comp method. 

	symmetric_dimer_pairs : list
		This list stores which dimers are symmetric to each other. This list will be populated during this method.

	no_of_cpus_for_comparing_two_dimers : int
		This is the number of cpus reserved for comparing two dimers together.
	"""

	# First, separate the input variables from input_data. 
	dimer1_name, dimer2_name, dimer1_details, dimer2_details, info_about_dimers_being_compared, distances_between_molecules_in_dimers, max_distance_disparity, neighbouring_molecules_about_dimers, non_hydrogen_molecules, non_hydrogen_graphs, symmetric_dimer_pairs, no_of_cpus_for_comparing_two_dimers = input_data

	# Second, obtain the shortest distance between the molecules in dimer 1 and dimer 2 based on the dimer invariance method used. 
	d1_shortest_distance, d2_shortest_distance = distances_between_molecules_in_dimers

	# Third, check if the shortest distance between molecules in the two dimers within the given max threshold limit.
	#        * If they are not, then the two dimers are definitely not the same.
	#        * This is not necessary to the functioning of the method, but speeds things up a lot of it is obvious that two dimers are not equivalent.
	if abs(d1_shortest_distance - d2_shortest_distance) > shortest_distance_between_dimers_max_threshold_limit:
		return

	# Fourth, get the individual variables from dimer2_details. 
	((d1_m1_original_elements, d1_m1_original_positions, d1_m1_no_H_attached_to_nonH_atoms), (d1_m2_original_elements, d1_m2_original_positions, d1_m2_no_H_attached_to_nonH_atoms)) = dimer1_details

	# Fifth, get the individual variables from dimer2_details. 
	((d2_m1_original_elements, d2_m1_original_positions, d2_m1_no_H_attached_to_nonH_atoms), (d2_m2_original_elements, d2_m2_original_positions, d2_m2_no_H_attached_to_nonH_atoms)) = dimer2_details

	# Sixth, get the individual variables from info_about_dimers_being_compared. 
	((d1_m1_name, d1_m2_name, unit_cell_disp1), (d2_m1_name, d2_m2_name, unit_cell_disp2)) = info_about_dimers_being_compared

	# ----------------------------------------------------------------------------------------------------------------------------------
	# Seventh, obtain the molecules of the second dimer in these next two loops.
	#          * Arrangement 1
	dimer2_details_arrangement_1 = ((d2_m1_original_elements, d2_m1_original_positions, d2_m1_no_H_attached_to_nonH_atoms), (d2_m2_original_elements, d2_m2_original_positions, d2_m2_no_H_attached_to_nonH_atoms))
	#          * Arrangement 2
	dimer2_details_arrangement_2 = (dimer2_details_arrangement_1[1], dimer2_details_arrangement_1[0])
	# ----------------------------------------------------------------------------------------------------------------------------------

	# Eighth, determine if the two dimers are structurally variant of each other.
	#         * We determine if we can use the Minimal Elemental Abundance Method based on if there are enough atoms in the molecules to use it. 
	#         * If we do not have enough atoms to use this method, default to using the comprehensive method. 
	if (len(d1_m1_original_elements) >= 2) and (len(d1_m2_original_elements) >= 2) and (len(d2_m1_original_elements) >= 2) and (len(d2_m2_original_elements) >= 2):

		# 8.1: Determine if the dimers are variants using the Minimal Elemental Abundance Method with dimer 2 in arrangement 1.
		positions_are_variant_MEA_Arr1, align_vectors_gave_a_warning_Arrangement1 = are_dimers_variant_MEA(dimer1_details, dimer2_details_arrangement_1, max_distance_disparity, info_about_dimers_being_compared, neighbouring_molecules_about_dimers, non_hydrogen_molecules, no_of_cpus=no_of_cpus_for_comparing_two_dimers)

		# 8.2: Found that these dimers are variants using the Minimal Elemental Abundance Method with dimer 2 in arrangement 1. Append this dimer pair to the symmetric_dimer_pairs list.
		if positions_are_variant_MEA_Arr1:
			symmetric_dimer_pairs.append((dimer1_name, dimer2_name))
			return

		# 8.3: Determine if the dimers are variants using the Minimal Elemental Abundance Method with dimer 2 in arrangement 2.
		positions_are_variant_MEA_Arr2, align_vectors_gave_a_warning_Arrangement2 = are_dimers_variant_MEA(dimer1_details, dimer2_details_arrangement_2, max_distance_disparity, info_about_dimers_being_compared, neighbouring_molecules_about_dimers, non_hydrogen_molecules, no_of_cpus=no_of_cpus_for_comparing_two_dimers)

		# 8.4: Found that these dimers are variants using the Minimal Elemental Abundance Method with dimer 2 in arrangement 2. Append this dimer pair to the symmetric_dimer_pairs list.
		if positions_are_variant_MEA_Arr2:
			symmetric_dimer_pairs.append((dimer1_name, dimer2_name))
			return

		# 8.5: If:
		#        1. The Minimal Elemental Abundance (MEA) method could not find a way that these dimers were variant and, 
		#        2. There was a warning issue given when running the Rotation.align_vectors method in scipy, 
		#      We will perform the comprehensive method just to make sure these two dimers are structurally equivalent or not. 
		if (len(align_vectors_gave_a_warning_Arrangement1) > 0) or (len(align_vectors_gave_a_warning_Arrangement2) > 0):

			# 8.6: Print to the user that you obtain a warning from Rotation.align_vectors method in scipy, and did not find a way to show that these two dimers are structurally equivalent,
			#      so we will make sure these two dimers are structurally equivalent or not using the comprehensive method.
			to_string = f'We received a warning from from scipy.spatial.transform.Rotation when using Rotation.align_vectors in the Minimal Elemental Abundance (MEA) method. Because we did not find structural equivalency between dimers {dimer1_name} and {dimer2_name} using the MEA method, we will double check this using the comprehensive method. '
			#warnings.warn(to_string)

			# 8.7: There are not enough different elements to perform the Minimal Elemental Abundance Method, so fall back to the more 
			#      computationally demanding comprehensive method. 
			if are_dimers_variant_from_comprehensive(dimer1_details, dimer2_details_arrangement_1, non_hydrogen_graphs[d1_m1_name], non_hydrogen_graphs[d1_m2_name], non_hydrogen_graphs[d2_m1_name], non_hydrogen_graphs[d2_m2_name], max_distance_disparity, info_about_dimers_being_compared, neighbouring_molecules_about_dimers, non_hydrogen_molecules, no_of_cpus=no_of_cpus_for_comparing_two_dimers) or are_dimers_variant_from_comprehensive(dimer1_details, dimer2_details_arrangement_2, non_hydrogen_graphs[d1_m1_name], non_hydrogen_graphs[d1_m2_name], non_hydrogen_graphs[d2_m2_name], non_hydrogen_graphs[d2_m1_name], max_distance_disparity, info_about_dimers_being_compared, neighbouring_molecules_about_dimers, non_hydrogen_molecules, no_of_cpus=no_of_cpus_for_comparing_two_dimers):
				
				# 8.8: Found that these dimers are variants using the Comprehensive Method. Append this dimer pair to the symmetric_dimer_pairs list.
				symmetric_dimer_pairs.append((dimer1_name, dimer2_name))
				
	else:

		# 8.11: There are not enough different elements to perform the Minimal Elemental Abundance Method, so fall back to the more 
		#      computationally demanding comprehensive method. 
		if are_dimers_variant_from_comprehensive(dimer1_details, dimer2_details_arrangement_1, non_hydrogen_graphs[d1_m1_name], non_hydrogen_graphs[d1_m2_name], non_hydrogen_graphs[d2_m1_name], non_hydrogen_graphs[d2_m2_name], max_distance_disparity, info_about_dimers_being_compared, neighbouring_molecules_about_dimers, non_hydrogen_molecules, no_of_cpus=no_of_cpus_for_comparing_two_dimers) or are_dimers_variant_from_comprehensive(dimer1_details, dimer2_details_arrangement_2, non_hydrogen_graphs[d1_m1_name], non_hydrogen_graphs[d1_m2_name], non_hydrogen_graphs[d2_m2_name], non_hydrogen_graphs[d2_m1_name], max_distance_disparity, info_about_dimers_being_compared, neighbouring_molecules_about_dimers, non_hydrogen_molecules, no_of_cpus=no_of_cpus_for_comparing_two_dimers):
			
			# 8.12: Found that these dimers are variants using the Comprehensive Method. Append this dimer pair to the symmetric_dimer_pairs list.
			symmetric_dimer_pairs.append((dimer1_name, dimer2_name))

def get_inputs(dimers, non_hydrogen_molecules_elements, non_hydrogen_molecules_positions, all_no_of_H_on_atoms_in_molecule, max_distance_disparity, neighbouring_molecules_about_dimers, non_hydrogen_molecules, non_hydrogen_graphs, symmetric_dimer_pairs, no_of_cpus_for_comparing_two_dimers):
	"""
	This generator is designed to return all the input methods required for the compare_if_two_dimers_are_symmetric_single_process method. 

	Parameters
	----------
	dimers : dict of int: (int, int, (int, int, int), numpy.array, numpy.array)
		This list contains all the information about how to construct the dimers from the molecules in the crystal. This information is:
			m1_name : int
				This is the name of molecule 1 in the dimer. 
			m2_name : int
				This is the name of molecule 2 in the dimer. 
			unit_cell_disp : (int, int, int)
				This is the unit cell to move molecule 2 into relative to molecule 1, given as ijk values
			dist : 1D numpy.array
				This is the unit cell to move molecule 2 into relative to molecule 1, given in Angstroms in the x, y, and z directions. 
			move_com_by : 1D numpy.array
				This is the x, y, and z displacements to move the centre of mass of the dimer by to ... (figure this out), given in Angstroms.
	
	non_hydrogen_molecules_elements : dict.
		This is a dictionary of the elements of the molecules that can make up the dimers. This is given as a list of strings to represent the elements in each molecule. 
	non_hydrogen_molecules_positions : dict.
		This is a dictionary of the positions of the molecules that can make up the dimers. This is given as a list of numpy.array to represent the positions in each molecule. 
	all_no_of_H_on_atoms_in_molecule : dict.
		This is a dictionary of the number of hydrogens bound to each "heavy" atom in the molecules that can make up the dimers. This is given as a list of ints to represent the number of bound hydrogens in each molecule. 

	d1_shortest_distance : float
		This is the shortest distance or general distance between two molecules in dimer 1 based on the choice of dimer method that was used.
	d2_shortest_distance : float
		This is the shortest distance or general distance between two molecules in dimer 2 based on the choice of dimer method that was used.

	max_distance_disparity : float
	    This is the maximum difference in relative positions of atoms between dimer 1 and dimer 2 to be considered equivalent/identical
	neighbouring_molecules_about_dimers : dict.
	    This is the information about the molecules that surround (in the vicinity of) each dimer in the crystal.
	non_hydrogen_molecules : dict. of ase.Atoms
	    These are Atoms of objects of the molecules in the crystal. These do not contain hydrogen atoms. This is for debugging use only. 
	non_hydrogen_graphs : dict of networkx.Graphs
		These are the graphs of the associated molecule ase.Atoms objects in non_hydrogen_molecules. The hydrogens in this graph has been removed as nodes and instead appended to "heavy" atoms in the molecule as node variables. This is used for the are_dimers_variant_from_comp method. 

	symmetric_dimer_pairs : list
		This list stores which dimers are symmetric to each other. 

	no_of_cpus_for_comparing_two_dimers : int
		This is the number of cpus reserved for comparing two dimers together.
	"""

	# First, obtain the names of the dimers to look through. 
	dimer_names = sorted(dimers.keys())

	# Second, obtain the molecules of the first dimer in the first two loops.
	for index1 in range(len(dimer_names)):

		# Third, obtain the name of the first dimer to examine. 
		dimer1_name = dimer_names[index1]

		# Fourth, obtain the indices of the molecules in the dimer, as well as the displacement of molecule 2 in dimer 1, and the centre of mass molecule to move dimer 2 by. 
		d1_m1_name, d1_m2_name, unit_cell_disp1, dist1, move_com_by_1, d1_shortest_distance = dimers[dimer1_name]

		# Fifth, obtain the elements, positions, and number of hydrogens bound to "heavy" atoms in molecule 1 of dimer 1
		d1_m1_original_elements           = non_hydrogen_molecules_elements[d1_m1_name]
		d1_m1_original_positions          = non_hydrogen_molecules_positions[d1_m1_name] + move_com_by_1
		d1_m1_no_H_attached_to_nonH_atoms = all_no_of_H_on_atoms_in_molecule[d1_m1_name]

		# Sixth, obtain the elements, positions, and number of hydrogens bound to "heavy" atoms in molecule 2 of dimer 1
		d1_m2_original_elements           = non_hydrogen_molecules_elements[d1_m2_name]
		d1_m2_original_positions          = non_hydrogen_molecules_positions[d1_m2_name] + move_com_by_1 + dist1
		d1_m2_no_H_attached_to_nonH_atoms = all_no_of_H_on_atoms_in_molecule[d1_m2_name]

		# Seventh, collect the above details into a single tuple to keep the data together.
		dimer1_details = ((d1_m1_original_elements, d1_m1_original_positions, d1_m1_no_H_attached_to_nonH_atoms), (d1_m2_original_elements, d1_m2_original_positions, d1_m2_no_H_attached_to_nonH_atoms))

		# Eighth, for every other dimer in dimers
		for index2 in range(index1+1,len(dimer_names)):

			# Ninth, obtain the name of the first dimer to examine. 
			dimer2_name = dimer_names[index2]

			# Tenth, obtain the indices of the molecules in the dimer, as well as the displacement of molecule 2 in dimer 1, and the centre of mass molecule to move dimer 2 by. 
			d2_m1_name, d2_m2_name, unit_cell_disp2, dist2, move_com_by_2, d2_shortest_distance = dimers[dimer2_name]

			# Eleventh, obtain the elements, positions, and number of hydrogens bound to "heavy" atoms in molecule 1 of dimer 2
			d2_m1_original_elements           = non_hydrogen_molecules_elements[d2_m1_name]
			d2_m1_original_positions          = non_hydrogen_molecules_positions[d2_m1_name] + move_com_by_2
			d2_m1_no_H_attached_to_nonH_atoms = all_no_of_H_on_atoms_in_molecule[d2_m1_name]

			# Twelfth, obtain the elements, positions, and number of hydrogens bound to "heavy" atoms in molecule 1 of dimer 2
			d2_m2_original_elements           = non_hydrogen_molecules_elements[d2_m2_name]
			d2_m2_original_positions          = non_hydrogen_molecules_positions[d2_m2_name] + move_com_by_2 + dist2
			d2_m2_no_H_attached_to_nonH_atoms = all_no_of_H_on_atoms_in_molecule[d2_m2_name]

			# Thirteenth, collect the above details into a single tuple to keep the data together.
			dimer2_details = ((d2_m1_original_elements, d2_m1_original_positions, d2_m1_no_H_attached_to_nonH_atoms), (d2_m2_original_elements, d2_m2_original_positions, d2_m2_no_H_attached_to_nonH_atoms))

			# Fourteenth, collect general details about the two dimers being compared.
			info_about_dimers_being_compared = ((d1_m1_name, d1_m2_name, unit_cell_disp1), (d2_m1_name, d2_m2_name, unit_cell_disp2))

			# Fifteenth, collect the distances between molecules in the two dimers 
			distances_between_molecules_in_dimers = (d1_shortest_distance, d2_shortest_distance)

			# Sixteenth, yield dimer details for the two dimers being compared.
			yield (dimer1_name, dimer2_name, dimer1_details, dimer2_details, info_about_dimers_being_compared, distances_between_molecules_in_dimers, max_distance_disparity, neighbouring_molecules_about_dimers, non_hydrogen_molecules, non_hydrogen_graphs, symmetric_dimer_pairs, no_of_cpus_for_comparing_two_dimers)

# --------------------------------------------------------------------------------------------------------------

