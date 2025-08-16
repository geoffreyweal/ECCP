"""
minimal_elemental_abundance_invariance_method.py, Geoffrey Weal, 5/4/22

This method will check that the elements in each dimer are the same and that the two dimers are rotationally variant.
"""
import numpy as np
import multiprocessing as mp
from SUMELF import are_two_lists_within_eachother

from ECCP.ECCP.invariance_methods.common_minimal_elemental_abundance_invariance_utility_methods.methods_for_MEA_invariance_method.get_positions_of_low_abundant_elements_to_scan import get_positions_of_low_abundant_elements_to_scan
from ECCP.ECCP.invariance_methods.common_minimal_elemental_abundance_invariance_utility_methods.methods_for_MEA_invariance_method.get_points_of_first_dimer                      import get_points_of_first_dimer
from ECCP.ECCP.invariance_methods.common_minimal_elemental_abundance_invariance_utility_methods.methods_for_MEA_invariance_method.obtain_possible_3D_points_in_dimer2            import obtain_possible_3D_points_in_dimer2
from ECCP.ECCP.invariance_methods.common_minimal_elemental_abundance_invariance_utility_methods.methods_for_MEA_invariance_method.get_points_of_second_dimer                     import get_points_of_second_dimer
from ECCP.ECCP.invariance_methods.common_minimal_elemental_abundance_invariance_utility_methods.determine_invariance_MEA                                                         import determine_invariance_MEA

# Set parameters for the are_dimers_variant_MEA ,ethod
length_max_distance_disparity = 0.01
dotproduct_max_distance_disparity = 0.01

def are_dimers_variant_MEA(dimer1_details, dimer2_details, max_distance_disparity, dimer_being_compared_info, neighbouring_molecules_about_dimers, non_hydrogen_molecules, no_of_cpus=1): 
	"""
	This method will check that the elements in each dimer are the same and that the two dimers are rotationally variant.
	
	Parameters
	----------
	dimer1_details : tuple 
		This contains the details of the elements and positions of the molecules in the first dimer.
	dimer2_details : tuple 
		This contains the details of the elements and positions of the molecules in the second dimer.
	max_distance_disparity: float
		This is the maximum that any two "could be equivalent" atoms can be between dimer 1 and dimer 2 for dimers 1 and 2 to be considered variant

	dimer_being_compared_info: list
		These are the information about how to make the dimer from the molecules
	neighbouring_molecules_about_dimers : dict.
		This is the information about the molecules that surround (in the vicinity of) each dimer in the crystal.
	non_hydrogen_molecules : list of ase.Atoms
		This is the list of molecules in the crystal, not including hydrogens.

	no_of_cpus : int.
		This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.

	Attributes
	----------
	length_max_distance_disparity : float
		This is the difference in distances in lengths of two dimers to be considered potentially variant. This allows for a quick determination of if the two dimers are definitely not variant. 
	dotproduct_max_distance_disparity : float
		This is the difference in the angles between lengths of two dimers to be considered potentially variant. This allows for a quick determination of if the two dimers are definitely not variant. 

	Returns
	-------
	dimers_are_variant : bool
		This indicates if the two dimers are the same or not as determined using the Minimal Elemental Abundance (MEA) method. 
	align_vectors_gave_a_warning_all : set
		This set contains all the unique warning that were obtained during the Minimal Elemental Abundance (MEA) method (usually from numpy).
	"""

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# First, obtain the elements, positions, and number of hydrogens bound to each "heavy" atom for the first dimer.
	((d1_m1_original_elements, d1_m1_original_positions, d1_m1_no_H_attached_to_nonH_atoms), (d1_m2_original_elements, d1_m2_original_positions, d1_m2_no_H_attached_to_nonH_atoms)) = dimer1_details
	dimer1_elements        = d1_m1_original_elements + d1_m2_original_elements
	dimer1_positions       = np.concatenate([d1_m1_original_positions, d1_m2_original_positions])
	dimer1_no_of_hydrogens = d1_m1_no_H_attached_to_nonH_atoms + d1_m2_no_H_attached_to_nonH_atoms

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Second, obtain the elements, positions, and number of hydrogens bound to each "heavy" atom for the second dimer.
	((d2_m1_original_elements, d2_m1_original_positions, d2_m1_no_H_attached_to_nonH_atoms), (d2_m2_original_elements, d2_m2_original_positions, d2_m2_no_H_attached_to_nonH_atoms)) = dimer2_details
	dimer2_elements        = d2_m1_original_elements + d2_m2_original_elements
	dimer2_positions       = np.concatenate([d2_m1_original_positions, d2_m2_original_positions])
	dimer2_no_of_hydrogens = d2_m1_no_H_attached_to_nonH_atoms + d2_m2_no_H_attached_to_nonH_atoms

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Third, if the order of atoms are not the same, then these two dimers are not made of the same molecules, so are definitely different dimers!
	if   ((sorted(d1_m1_original_elements) == sorted(d2_m1_original_elements)) and (sorted(d1_m2_original_elements) == sorted(d2_m2_original_elements))):
		# 3.1: This means that d1_m1 could go with d2_m1, and d1_m2 could go with d2_m2.
		pass
	elif ((sorted(d1_m1_original_elements) == sorted(d2_m2_original_elements)) and (sorted(d1_m2_original_elements) == sorted(d2_m1_original_elements))):
		# 3.2: This means that d1_m1 could go with d2_m2, and d1_m2 could go with d2_m1.
		pass
	else:
		return False, set()

	# ----------------------------------------------------------------------------------------------------------------------------------------------------------
	# Fourth, get the positions of alignment dimers.

	# 4.1: Get the positions of the two molecules in dimer 1 for alignment.
	d1_m1_positions_of_lowest_elements = get_positions_of_low_abundant_elements_to_scan(d1_m1_original_elements, d1_m1_original_positions, max_number_of_atoms=2)
	d1_m2_positions_of_lowest_elements = get_positions_of_low_abundant_elements_to_scan(d1_m2_original_elements, d1_m2_original_positions, max_number_of_atoms='all')

	# 4.2: Get the positions of the two molecules in dimer 2 for alignment.
	d2_m1_positions_of_lowest_elements = get_positions_of_low_abundant_elements_to_scan(d2_m1_original_elements, d2_m1_original_positions, max_number_of_atoms=2)
	d2_m2_positions_of_lowest_elements = get_positions_of_low_abundant_elements_to_scan(d2_m2_original_elements, d2_m2_original_positions, max_number_of_atoms='all')
	
	# ----------------------------------------------------------------------------------------------------------------------------------------------------------
	# Fifth, obtain the directions, vector lengths, and dotproducts to focus on in dimer 1.

	# 5.1: Obtain the directions, vector lengths, and dotproducts to focus on in dimer 1.
	d1_direction1, d1_direction2, d1_direction3, d1_centre_position, lengths_d1, dotproducts_d1_sorted, indices_of_points_d1, elements_in_points_d1 = get_points_of_first_dimer(d1_m1_positions_of_lowest_elements, d1_m2_positions_of_lowest_elements, len(d1_m1_original_elements), d1_m2_original_elements=d1_m2_original_elements, d1_m2_original_positions=d1_m2_original_positions)
	
	# 5.2: Get the indices of atom point in dimer 1.
	d1_centre_index, d1_point1_index, d1_point2_index, d1_point3_index = indices_of_points_d1
	
	# 5.3: Get the elements of atom point in dimer 1.
	d1_centre_element, d1_point1_element, d1_point2_element, d1_point3_element = elements_in_points_d1
	
	# ----------------------------------------------------------------------------------------------------------------------------------------------------------
	# Sixth, initialise the return values.

	# 6.1: Initialise a boolean to indicate if two molecules have been found to be variants of each other.
	is_variant = False

	# 6.2: Initialise align_vectors_gave_a_warning as an empty set.
	#      * align_vectors_gave_a_warning: Will contain a warning message object if determine_invariance_MEA gives a warning. 
	align_vectors_gave_a_warning_all = set()

	# ----------------------------------------------------------------------------------------------------------------------------------------------------------
	# Seventh, setup the input generator for performing the compare_two_dimers_in_two_index_configurations_single_process method.

	# 7.1: Initialise the direction points of dimer 1.
	d1_directions = (d1_direction1, d1_direction2, d1_direction3)

	# 7.2: Initialise the indices of the atoms for dimer 1.
	d1_indices_in_directions = (d1_point1_index, d1_point2_index, d1_point3_index)

	# 7.3: Create the generator that provides the inputs for the compare_two_molecules_in_two_index_configurations_single_process method.
	input_generator = get_inputs(d2_m1_positions_of_lowest_elements, d2_m2_positions_of_lowest_elements, d1_centre_element, d1_point1_element, d1_point2_element, d1_point3_element, d1_directions, d1_centre_index, d1_indices_in_directions, d2_m1_original_elements, lengths_d1, dotproducts_d1_sorted, length_max_distance_disparity, dotproduct_max_distance_disparity, dimer1_elements, dimer1_positions, dimer1_no_of_hydrogens, dimer2_elements, dimer2_positions, dimer2_no_of_hydrogens, max_distance_disparity, dimer_being_compared_info, neighbouring_molecules_about_dimers, non_hydrogen_molecules)

	# ----------------------------------------------------------------------------------------------------------------------------------------------------------

	# Eighth, perform the compare_two_dimers_in_two_index_configurations_single_process method on all inputs depending on if you are performing the task with one cpu or with multiple cores.
	if no_of_cpus == 1: # If the user only wants to use 1 cpu, perform tasks without using multiprocessing

		# 8.1: Analyse the spatial positions of dimer 2, for all the different possible index arrangements of dimer 2.
		for input_data in input_generator:

			# 8.2: Determine if the two molecules are varient in the given index arrangement of dimer 2.
			is_variant, align_vectors_gave_a_warning_single = compare_two_dimers_in_two_index_configurations_single_process(input_data)

			# 8.3: Update align_vectors_gave_a_warning_all with align_vectors_gave_a_warning_single
			align_vectors_gave_a_warning_all.update(align_vectors_gave_a_warning_single)

			# 8.4: If the dimer is variant, break out of the loop. 
			if is_variant:
				break
			else:
				align_vectors_gave_a_warning_all.update(align_vectors_gave_a_warning_single)

	else:

		# 8.5: Set up the multiprocessing pool.
		with mp.Pool(processes=no_of_cpus) as pool:

			# 8.6: For the output of compare_two_dimers_in_two_index_configurations_single_process for each of the numerous index arrangements of molecule 2.
			for is_variant, align_vectors_gave_a_warning_single in pool.imap(compare_two_dimers_in_two_index_configurations_single_process, input_generator):

				# 8.7: Update align_vectors_gave_a_warning_all with align_vectors_gave_a_warning_single
				align_vectors_gave_a_warning_all.update(align_vectors_gave_a_warning_single)

				# 8.8: If the molecule is variant, break out of the loop. 
				if is_variant:
					break
				else:
					align_vectors_gave_a_warning_all.update(align_vectors_gave_a_warning_single)

	# Ninth, add any warning messages from align_vectors_gave_a_warning_single to align_vectors_gave_a_warning_all
	if is_variant:
		align_vectors_gave_a_warning_all = set()

	# ----------------------------------------------------------------------------------------------------------------------------------------------------------

	# Tenth, return if the two dimer are variants of each other, and any warning messages from numpy. 
	return is_variant, align_vectors_gave_a_warning_all

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

def get_inputs(d2_m1_positions_of_lowest_elements, d2_m2_positions_of_lowest_elements, d1_centre_element, d1_point1_element, d1_point2_element, d1_point3_element, d1_directions, d1_centre_index, d1_indices_in_directions, d2_m1_original_elements, lengths_d1, dotproducts_d1_sorted, length_max_distance_disparity, dotproduct_max_distance_disparity, dimer1_elements, dimer1_positions, dimer1_no_of_hydrogens, dimer2_elements, dimer2_positions, dimer2_no_of_hydrogens, max_distance_disparity, dimer_being_compared_info, neighbouring_molecules_about_dimers, non_hydrogen_molecules):
	"""
	This method is designed to provide the inputs for the ```compare_two_dimers_in_two_index_configurations_single_process``` method

	Parameters
	----------
	d2_m1_positions_of_lowest_elements : list
		This list contains the (element, position, index) of the atoms used to allign the first molecule of dimer 2 with the molecules in dimer 1.
	d2_m2_positions_of_lowest_elements
		This list contains the (element, position, index) of the atoms used to allign the second molecule of dimer 2 with the molecules in dimer 1.
	
	d1_centre_element : str
		This is the element of the centre atom in the MEA analysis of dimer 1. 
	d1_point1_element
		This is the element of the point 1 atom in the MEA analysis of dimer 1. 
	d1_point2_element
		This is the element of the point 2 atom in the MEA analysis of dimer 1. 
	d1_point3_element
		This is the element of the point 3 atom in the MEA analysis of dimer 1. 
	d1_directions : list of numpy.array
		These are the directions from the centre atom to point 1 atom (d1_directions[0]), point 2 atom (d1_directions[1]), and point 3 atom (d1_directions[2]) in dimer 1
	d1_centre_index : int
		This is the index of the centre atom for dimer 1. 
	d1_indices_in_directions : list of ints.
		These are the indices of point 1, 2, and 3 atoms in dimer 1.

	d2_m1_original_elements : list of str.
		These are the elements in the first molecule in dimer 2

	lengths_d1 : list of floats.
		These are the lengths of the direction vectors for dimer 1, given in d1_directions.
	dotproducts_d1_sorted: list of floats.
		These are the dot products (measuring angles between direction vectors) of the direction vectors for dimer 1, given in d1_directions.

	length_max_distance_disparity : float
		This is the maximum difference in the lengths of the direction vectors to be considered the same dimer.
	dotproduct_max_distance_disparity
		This is the maximum difference in the dot products (measuring angles between direction vectors) of the direction vectors to be considered the same dimer.

	dimer1_elements : list of str.
		These are the elements of the atoms in dimer 1.
	dimer1_positions : list of (1x3) numpy array
		These are the positions of the atoms in dimer 1.
	dimer1_no_of_hydrogens : list of ints.
		These are the number of hydrogens attached to each atom in dimer 1.

	dimer2_elements : list of str.
		These are the elements of the atoms in dimer 2.
	dimer2_positions : list of (1x3) numpy array
		These are the positions of the atoms in dimer 2.
	dimer2_no_of_hydrogens : list of ints.
		These are the number of hydrogens attached to each atom in dimer 2.

    max_distance_disparity: float
        This is the maximum that any two "could be equivalent" atoms can be between dimer 1 and dimer 2 for dimer 1 and 2 to be considered variant.

	dimer_being_compared_info : tuple
		This is the information about the two dimers being compared, in the format: dimer info --> (mol1, mol2, unit cell displacement of mol2), 

	neighbouring_molecules_about_dimers : dict.
		This dictionary contains information about the molecules that surround the molecules in the dimer.
	non_hydrogen_molecules : dict of ase.Atoms
		This dictionary contains all the molecules in the crystal without hydrogens.
	
	Returns
	-------
	d2_centre_position : (1x3) numpy.array
		This is the position of the centre atom in dimer 2.
	d2_point1_position : (1x3) numpy.array
		This is the position of the point 1 atom in dimer 2.
	d2_point2_position : (1x3) numpy.array
		This is the position of the point 2 atom in dimer 2.
	d2_point3_position : (1x3) numpy.array
		This is the position of the point 3 atom in dimer 2.

	d1_directions : list of numpy.array
		These are the directions from the centre atom to point 1 atom (d1_directions[0]), point 2 atom (d1_directions[1]), and point 3 atom (d1_directions[2]) in dimer 1
	
	d1_indices_in_directions : list of ints.
		These are the indices of point 1, 2, and 3 atoms in dimer 1.
	d2_indices_in_directions : list of ints.
		These are the indices of point 1, 2, and 3 atoms in dimer 2.

	lengths_d1 : list of floats.
		These are the lengths of the direction vectors for dimer 1, given in d1_directions.
	dotproducts_d1_sorted: list of floats.
		These are the dot products (measuring angles between direction vectors) of the direction vectors for dimer 1, given in d1_directions.

	length_max_distance_disparity : float
		This is the maximum difference in the lengths of the direction vectors to be considered the same dimer.
	dotproduct_max_distance_disparity
		This is the maximum difference in the dot products (measuring angles between direction vectors) of the direction vectors to be considered the same dimer.

	d1_centre_index : int
		This is the index of the centre atom for dimer 1. 
	d2_centre_index : int
		This is the index of the centre atom for dimer 2. 

	dimer1_elements : list of str.
		These are the elements of the atoms in dimer 1.
	dimer1_positions : list of (1x3) numpy array
		These are the positions of the atoms in dimer 1.
	dimer1_no_of_hydrogens : list of ints.
		These are the number of hydrogens attached to each atom in dimer 1.

	dimer2_elements : list of str.
		These are the elements of the atoms in dimer 2.
	dimer2_positions : list of (1x3) numpy array
		These are the positions of the atoms in dimer 2.
	dimer2_no_of_hydrogens : list of ints.
		These are the number of hydrogens attached to each atom in dimer 2.

    max_distance_disparity: float
        This is the maximum that any two "could be equivalent" atoms can be between dimer 1 and dimer 2 for dimer 1 and 2 to be considered variant.

	dimer_being_compared_info : tuple
		This is the information about the two dimers being compared, in the format: dimer info --> (mol1, mol2, unit cell displacement of mol2), 

	neighbouring_molecules_about_dimers : dict.
		This dictionary contains information about the molecules that surround the molecules in the dimer.
	non_hydrogen_molecules : dict of ase.Atoms
		This dictionary contains all the molecules in the crystal without hydrogens.
	"""

	# First, analyse the spatial positions of dimer 2.
	for d2_centre_element, d2_centre_position, d2_centre_index, d2_point1_element, d2_point1_position, d2_point1_index, d2_point2_element, d2_point2_position, d2_point2_index, d2_point3_element, d2_point3_position, d2_point3_index in obtain_possible_3D_points_in_dimer2(d2_m1_positions_of_lowest_elements, d2_m2_positions_of_lowest_elements, d1_centre_element, d1_point1_element, d1_point2_element, d1_point3_element, len(d2_m1_original_elements)):

		# Second, initialise the indices of the atoms for dimer 2.
		d2_indices_in_directions = (d2_point1_index, d2_point2_index, d2_point3_index)

		# Third, return the inputs needed for the ``compare_two_dimers_in_two_index_configurations_single_process`` method
		yield d2_centre_position, d2_point1_position, d2_point2_position, d2_point3_position, d1_directions, d1_indices_in_directions, d2_indices_in_directions, lengths_d1, dotproducts_d1_sorted, length_max_distance_disparity, dotproduct_max_distance_disparity, d1_centre_index, d2_centre_index, dimer1_elements, dimer1_positions, dimer1_no_of_hydrogens, dimer2_elements, dimer2_positions, dimer2_no_of_hydrogens, max_distance_disparity, dimer_being_compared_info, neighbouring_molecules_about_dimers, non_hydrogen_molecules

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

def compare_two_dimers_in_two_index_configurations_single_process(input_data):
	"""
	This method is designed to determine if your two dimers are equivalent to each other via the "Minimal Elemential Abundance" methodology.

	Parameters
	----------
	d2_centre_position : (1x3) numpy.array
		This is the position of the centre atom in dimer 2.
	d2_point1_position : (1x3) numpy.array
		This is the position of the point 1 atom in dimer 2.
	d2_point2_position : (1x3) numpy.array
		This is the position of the point 2 atom in dimer 2.
	d2_point3_position : (1x3) numpy.array
		This is the position of the point 3 atom in dimer 2.

	d1_directions : list of numpy.array
		These are the directions from the centre atom to point 1 atom (d1_directions[0]), point 2 atom (d1_directions[1]), and point 3 atom (d1_directions[2]) in dimer 1
	
	d1_indices_in_directions : list of ints.
		These are the indices of point 1, 2, and 3 atoms in dimer 1.
	d2_indices_in_directions : list of ints.
		These are the indices of point 1, 2, and 3 atoms in dimer 2.

	lengths_d1 : list of floats.
		These are the lengths of the direction vectors for dimer 1, given in d1_directions.
	dotproducts_d1_sorted: list of floats.
		These are the dot products (measuring angles between direction vectors) of the direction vectors for dimer 1, given in d1_directions.

	length_max_distance_disparity : float
		This is the maximum difference in the lengths of the direction vectors to be considered the same dimer.
	dotproduct_max_distance_disparity
		This is the maximum difference in the dot products (measuring angles between direction vectors) of the direction vectors to be considered the same dimer.

	d1_centre_index : int
		This is the index of the centre atom for dimer 1. 
	d2_centre_index : int
		This is the index of the centre atom for dimer 2. 

	dimer1_elements : list of str.
		These are the elements of the atoms in dimer 1.
	dimer1_positions : list of (1x3) numpy array
		These are the positions of the atoms in dimer 1.
	dimer1_no_of_hydrogens : list of ints.
		These are the number of hydrogens attached to each atom in dimer 1.

	dimer2_elements : list of str.
		These are the elements of the atoms in dimer 2.
	dimer2_positions : list of (1x3) numpy array
		These are the positions of the atoms in dimer 2.
	dimer2_no_of_hydrogens : list of ints.
		These are the number of hydrogens attached to each atom in dimer 2.

    max_distance_disparity: float
        This is the maximum that any two "could be equivalent" atoms can be between dimer 1 and dimer 2 for dimer 1 and 2 to be considered variant.

	dimer_being_compared_info : tuple
		This is the information about the two dimers being compared, in the format: dimer info --> (mol1, mol2, unit cell displacement of mol2), 

	neighbouring_molecules_about_dimers : dict.
		This dictionary contains information about the molecules that surround the molecules in the dimer (used for environmental purposes).
	non_hydrogen_molecules : dict of ase.Atoms
		This dictionary contains all the molecules in the crystal without hydrogens.

	Returns
	-------
	dimers_are_variant : bool
		This indicates if the two dimers are the same or not as determined using the Minimal Elemental Abundance (MEA) method. 
	align_vectors_gave_a_warning_single : set
		This set contains all the unique warning that were obtained during the Minimal Elemental Abundance (MEA) method (usually from numpy).
	"""

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	# First: Obtain the inputs from the input data.

	# 1.1: Obtain the inputs from the input data.
	d2_centre_position, d2_point1_position, d2_point2_position, d2_point3_position, d1_directions, d1_indices_in_directions, d2_indices_in_directions, lengths_d1, dotproducts_d1_sorted, length_max_distance_disparity, dotproduct_max_distance_disparity, d1_centre_index, d2_centre_index, dimer1_elements, dimer1_positions, dimer1_no_of_hydrogens, dimer2_elements, dimer2_positions, dimer2_no_of_hydrogens, max_distance_disparity, dimer_being_compared_info, neighbouring_molecules_about_dimers, non_hydrogen_molecules = input_data

	# 1.2: Obtain the direction positions of the points for dimer 1.
	d1_direction1, d1_direction2, d1_direction3 = d1_directions

	# 1.3: Obtain the indices of the points for dimer 1.
	d1_point1_index, d1_point2_index, d1_point3_index = d1_indices_in_directions

	# 1.4: Obtain the indices of the points for dimer 2.
	d2_point1_index, d2_point2_index, d2_point3_index = d2_indices_in_directions

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	# Second, obtain a possible set of directions, vector lengths, and dotproducts to focus on in dimer 2.
	d2_direction1, d2_direction2, d2_direction3, lengths_d2, dotproducts_d2_sorted = get_points_of_second_dimer(d2_centre_position, d2_point1_position, d2_point2_position, d2_point3_position)

	# Third, if the vector lengths of directions and dot-products (angles) between directions are not the same, then move on to another of the possible 3D points in dimer 2.
	if not (are_two_lists_within_eachother(sorted(lengths_d1), sorted(lengths_d2), length_max_distance_disparity) and are_two_lists_within_eachother(dotproducts_d1_sorted, dotproducts_d2_sorted, dotproduct_max_distance_disparity)):
		return False, set()

	# Fourth, get the directions and the indices of directions for dimer 1
	d1_directions = [d1_direction1, d1_direction2, d1_direction3]
	d1_indices_in_directions = [d1_point1_index, d1_point2_index, d1_point3_index]

	# Fifth, get the directions and the indices of directions for dimer 2
	d2_directions = [d2_direction1, d2_direction2, d2_direction3]
	d2_indices_in_directions = [d2_point1_index, d2_point2_index, d2_point3_index]

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	# Sixth, get the matrix to transform dimer 2 upon dimer 1.
	dimers_are_variant, dimer1_m12_to_dimer2_m12_conversion, align_vectors_gave_a_warning_single = determine_invariance_MEA(d1_directions, d1_indices_in_directions, lengths_d1, d2_directions, d2_indices_in_directions, lengths_d2, d1_centre_index, d2_centre_index, dimer1_elements, dimer1_positions, dimer1_no_of_hydrogens, dimer2_elements, dimer2_positions, dimer2_no_of_hydrogens, max_distance_disparity, length_max_distance_disparity, dimer_being_compared_info, neighbouring_molecules_about_dimers, non_hydrogen_molecules)

	# Seventh, return the result from the determine_invariance_MEA method
	return dimers_are_variant, align_vectors_gave_a_warning_single

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

