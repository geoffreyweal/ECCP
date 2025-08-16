"""
are_dimers_variant_from_comprehensive.py, Geoffrey Weal, 9/2/24

This method will check that the elements in each molecule are the same and that the two dimers are rotationally variant. 
"""
import numpy as np
from copy import deepcopy
import multiprocessing as mp

from SUMELF import GraphMatcher

from ECCP.ECCP.invariance_methods.utilities                                                             import get_permutated_indices_list
from ECCP.ECCP.invariance_methods.common_utility_methods_for_all_invariance_methods.are_systems_variant import are_systems_variant

def are_dimers_variant_from_comprehensive(dimer1_details, dimer2_details, D1_M1_non_H_graph, D1_M2_non_H_graph, D2_M1_non_H_graph, D2_M2_non_H_graph, max_distance_disparity, info_about_dimers_being_compared, neighbouring_molecules_about_dimers, non_hydrogen_molecules, no_of_cpus=1): 
	"""
	This method will check that the elements in each molecule are the same and that the two dimers are rotationally variant. 
	
	Parameters
	----------
	dimer1_details : list
		This list contains the elements, positions, and no of hydrogens bound to non-hydrogen atoms in the two molecules that make up dimer 1.
	dimer2_details : list
		This list contains the elements, positions, and no of hydrogens bound to non-hydrogen atoms in the two molecules that make up dimer 2.

	D1_M1_non_H_graph : networkx.graph
		This is the non-hydrogen graph of molecule 1 in dimer 1.
	D1_M2_non_H_graph : networkx.graph
		This is the non-hydrogen graph of molecule 2 in dimer 1.
	D1_M1_non_H_graph : networkx.graph
		This is the non-hydrogen graph of molecule 1 in dimer 2.
	D2_M2_non_H_graph : networkx.graph
		This is the non-hydrogen graph of molecule 2 in dimer 2.

    max_distance_disparity: float
        This is the maximum that any two "could be equivalent" atoms can be between system 1 and system 2 for systems 1 and 2 to be considered variant.
    info_about_dimers_being_compared: list
        These are the names of the molecules and the unit vector displacement that makes up the two dimers being compared. 
    neighbouring_molecules_about_dimers : dict.
        This is the information about the molecules that surround (in the vicinity of) each dimer in the crystal.
    non_hydrogen_molecules : dict. of ase.Atoms
        This is the list of molecules in the crystal, not including hydrogens.

    no_of_cpus : int.
        This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.

	Returns
	-------
	True if the two dimers are variant. False if they are invariant (unique).
	"""

	# First, create the first dimer from d1_m1 and d1_m2 details.
	((d1_m1_original_elements, d1_m1_original_positions, d1_m1_no_H_attached_to_nonH_atoms), (d1_m2_original_elements, d1_m2_original_positions, d1_m2_no_H_attached_to_nonH_atoms)) = dimer1_details

	# Second, obtain the elements and positions of the second dimer
	((d2_m1_original_elements, d2_m1_original_positions, d2_m1_no_H_attached_to_nonH_atoms), (d2_m2_original_elements, d2_m2_original_positions, d2_m2_no_H_attached_to_nonH_atoms)) = dimer2_details

	# Third, if the order of atoms are not the same, then these two dimers are not made of the same molecules, so are definitely different dimers!
	if   ((sorted(d1_m1_original_elements) == sorted(d2_m1_original_elements)) and (sorted(d1_m2_original_elements) == sorted(d2_m2_original_elements))):
		# 3.1: This means that d1_m1 could go with d2_m1, and d1_m2 could go with d2_m2.
		pass
	elif ((sorted(d1_m1_original_elements) == sorted(d2_m2_original_elements)) and (sorted(d1_m2_original_elements) == sorted(d2_m1_original_elements))):
		# 3.2: This means that d1_m1 could go with d2_m2, and d1_m2 could go with d2_m1.
		pass
	else:
		return False

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Fourth, determine how the indices of d2_m1 can be related to d1_m1
	#         * This is set as GraphMatcher(D2_M1_non_H_graph, D1_M1_non_H_graph) rather than GraphMatcher(D1_M1_non_H_graph, D2_M1_non_H_graph)
	#           because we want to convert D2_M1_non_H_graph into D1_M1_non_H_graph
	GM_M1 = GraphMatcher(D2_M1_non_H_graph, D1_M1_non_H_graph)
	em_indices_d2_m1_to_d1_m1 = GM_M1.get_all_unique_isomorphic_graphs()

	# Fifth, determine how the indices of d2_m2 can be related to d2_m1
	#        * This is set as GraphMatcher(D2_M2_non_H_graph, D1_M2_non_H_graph) rather than GraphMatcher(D1_M2_non_H_graph, D2_M2_non_H_graph)
	#          because we want to convert D2_M2_non_H_graph into D1_M2_non_H_graph
	GM_M2 = GraphMatcher(D2_M2_non_H_graph, D1_M2_non_H_graph)
	em_indices_d2_m2_to_d1_m2 = GM_M2.get_all_unique_isomorphic_graphs()

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# Sixth, go through each variation of em_indices_d2_m1_to_d1_m1 and em_indices_d2_m2_to_d1_m2 and check if the two dimers are varients of each other
	#        using a comprehensive scheme.

	# 6.1: Initalise the boolean variable for recording if two dimers are variants of each other or not. 
	dimers_are_variants = False

	# 6.2: Create the input generator for running the ```are_systems_variant``` method.
	input_generator = get_inputs(em_indices_d2_m1_to_d1_m1, em_indices_d2_m2_to_d1_m2,  d1_m1_original_elements, d1_m1_original_positions, d1_m1_no_H_attached_to_nonH_atoms,  d1_m2_original_elements, d1_m2_original_positions, d1_m2_no_H_attached_to_nonH_atoms,  d2_m1_original_elements, d2_m1_original_positions, d2_m1_no_H_attached_to_nonH_atoms,  d2_m2_original_elements, d2_m2_original_positions, d2_m2_no_H_attached_to_nonH_atoms,  max_distance_disparity, info_about_dimers_being_compared, neighbouring_molecules_about_dimers, non_hydrogen_molecules)

	# 6.3: compare the two dimers using the ```are_systems_variant``` method.
	if no_of_cpus == 1: # If the user only wants to use 1 cpu, perform tasks without using multiprocessing

		# 6.3.1: Analyse the inputs variables required to compare dimers 1 and 2 with each other using the comprehensive ```are_systems_variant``` method.
		for input_data in input_generator:

			# 6.3.2: Determine if these two dimers are variants of each other for a given index arrangement for dimer 2. 
			dimers_are_variants = are_systems_variant(*input_data)

			if isinstance(dimers_are_variants, Exception):
				raise dimers_are_variants

			# 6.3.3: If dimers_are_variants is True, break out of this for loop
			if dimers_are_variants:
				break

	else:

		# 6.3.4: Set up the multiprocessing pool.
		with mp.Pool(processes=no_of_cpus) as pool:

			# 6.3.5: For the output of ```are_systems_variant``` for each of the numerous index arrangements of dimer 2.
			for dimers_are_variants in pool.imap(are_systems_variant_parallel, input_generator):

				if isinstance(dimers_are_variants, Exception):
					raise dimers_are_variants

				# 6.3.6: If dimers_are_variants is True, break out of this for loop
				if dimers_are_variants:
					break

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Seventh, if no variance could be found between these two dimers, return False, meaning the dimers are invariant. 
	return dimers_are_variants

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def get_inputs(em_indices_d2_m1_to_d1_m1, em_indices_d2_m2_to_d1_m2,  d1_m1_original_elements, d1_m1_original_positions, d1_m1_no_H_attached_to_nonH_atoms,  d1_m2_original_elements, d1_m2_original_positions, d1_m2_no_H_attached_to_nonH_atoms,  d2_m1_original_elements, d2_m1_original_positions, d2_m1_no_H_attached_to_nonH_atoms,  d2_m2_original_elements, d2_m2_original_positions, d2_m2_no_H_attached_to_nonH_atoms,  max_distance_disparity, info_about_dimers_being_compared, neighbouring_molecules_about_dimers, non_hydrogen_molecules):
	"""
	This method will give all the inputs needed for the ```are_systems_variant``` method.

	Parameters
	----------
	em_indices_d2_m1_to_d1_m1 : list of dict.
		This list contains all the way that the atoms (indices) of molecule 1 of dimer 2 can possibly be transposed onto the atoms (indices) of molecule 1(/2) of dimer 1.
	em_indices_d2_m2_to_d1_m2
		This list contains all the way that the atoms (indices) of molecule 2 of dimer 2 can possibly be transposed onto the atoms (indices) of molecule 2(/1) of dimer 1.

	d1_m1_original_elements : list of str.
		These are the elements of the atoms in molecule 1 of dimer 1 (All of these atoms should be non-hydrogens).
	d1_m1_original_positions : list of (1x3) numpy.array
		These are the positions of the atoms in molecule 1 of dimer 1 (All of these atoms should be non-hydrogens).
	d1_m1_no_H_attached_to_nonH_atoms : list of int
		These list contains the number of hydrogens bound to the (non-hydrogen) atoms in molecule 1 of dimer 1.

	d1_m2_original_elements : list of str.
		These are the elements of the atoms in molecule 2 of dimer 1 (All of these atoms should be non-hydrogens).
	d1_m2_original_positions : list of (1x3) numpy.array
		These are the positions of the atoms in molecule 2 of dimer 1 (All of these atoms should be non-hydrogens).
	d1_m2_no_H_attached_to_nonH_atoms : list of int
		These list contains the number of hydrogens bound to the (non-hydrogen) atoms in molecule 2 of dimer 1.

	d2_m1_original_elements : list of str.
		These are the elements of the atoms in molecule 1 of dimer 2 (All of these atoms should be non-hydrogens).
	d2_m1_original_positions : list of (1x3) numpy.array
		These are the positions of the atoms in molecule 1 of dimer 2 (All of these atoms should be non-hydrogens).
	d2_m1_no_H_attached_to_nonH_atoms : list of int
		These list contains the number of hydrogens bound to the (non-hydrogen) atoms in molecule 1 of dimer 2.

	d2_m2_original_elements : list of str.
		These are the elements of the atoms in molecule 2 of dimer 2 (All of these atoms should be non-hydrogens).
	d2_m2_original_positions : list of (1x3) numpy.array
		These are the positions of the atoms in molecule 2 of dimer 2 (All of these atoms should be non-hydrogens).
	d2_m2_no_H_attached_to_nonH_atoms : list of int
		These list contains the number of hydrogens bound to the (non-hydrogen) atoms in molecule 2 of dimer 2.

    max_distance_disparity: float
        This is the maximum that any two "could be equivalent" atoms can be between system 1 and system 2 for systems 1 and 2 to be considered variant.
    info_about_dimers_being_compared: list
        These are the names of the molecules and the unit vector displacement that makes up the two dimers being compared. 
    neighbouring_molecules_about_dimers : dict.
        This is the information about the molecules that surround (in the vicinity of) each dimer in the crystal.
    non_hydrogen_molecules : dict. of ase.Atoms
        This is the list of molecules in the crystal, not including hydrogens.

	Returns
	-------
	d1_elements : list of str.
		These are the elements of the atoms in dimer 1 (All of these atoms should be non-hydrogens).
	d1_positions : list of (1x3) numpy.array
		These are the positions of the atoms in dimer 1 (All of these atoms should be non-hydrogens).
	dimer1_no_of_hydrogens : list of int
		These list contains the number of hydrogens bound to the (non-hydrogen) atoms dimer 1.

	d2_elements_way_1 : list of str.
		These are the elements of the atoms in dimer 2 (created as M1,M2) (All of these atoms should be non-hydrogens).
	d2_positions_way_1 : list of (1x3) numpy.array
		These are the positions of the atoms in dimer 2 (created as M1,M2) (All of these atoms should be non-hydrogens).
	dimer2_no_of_hydrogens_way_1 : list of int
		These list contains the number of hydrogens bound to the (non-hydrogen) atoms dimer 2 (created as M1,M2).

    max_distance_disparity: float
        This is the maximum that any two "could be equivalent" atoms can be between system 1 and system 2 for systems 1 and 2 to be considered variant.
    info_about_dimers_being_compared: list
        These are the names of the molecules and the unit vector displacement that makes up the two dimers being compared. 
    neighbouring_molecules_about_dimers : dict.
        This is the information about the molecules that surround (in the vicinity of) each dimer in the crystal.
    non_hydrogen_molecules : dict. of ase.Atoms
        This is the list of molecules in the crystal, not including hydrogens.
	"""

	# First, obtain the lists of elements, positions, and no of bonded hydrogen for each non-hydrogen atom in dimer 1
	d1_elements                                         = d1_m1_original_elements + d1_m2_original_elements
	d1_positions                                        = np.concatenate([d1_m1_original_positions, d1_m2_original_positions])
	dimer1_no_of_hydrogens                              = d1_m1_no_H_attached_to_nonH_atoms + d1_m2_no_H_attached_to_nonH_atoms

	# Second, for each permutation of indices between d1 and d2 for m1.
	for comparison1 in em_indices_d2_m1_to_d1_m1:

		# Third, rearrange the element and position indices for d2_m1
		idx_d2_m1                                       = get_permutated_indices_list(comparison1)
		d2_m1_reordered_elements                        = [d2_m1_original_elements[index] for index in idx_d2_m1]
		d2_m1_reordered_positions                       = deepcopy(d2_m1_original_positions)[idx_d2_m1, :]
		d2_m1_no_H_attached_to_nonH_atoms_reordered     = [d2_m1_no_H_attached_to_nonH_atoms[index] for index in idx_d2_m1]

		# Fourth, for each permutation of indices between d1 and d2 for m2.
		for comparison2 in em_indices_d2_m2_to_d1_m2:

			# Fifth, rearrange the element and position indices for d2_m2
			idx_d2_m2                                   = get_permutated_indices_list(comparison2)
			d2_m2_reordered_elements                    = [d2_m2_original_elements[index] for index in idx_d2_m2]
			d2_m2_reordered_positions                   = deepcopy(d2_m2_original_positions)[idx_d2_m2, :]
			d2_m2_no_H_attached_to_nonH_atoms_reordered = [d2_m2_no_H_attached_to_nonH_atoms[index] for index in idx_d2_m2]

			# Sixth, get the elements and positions of atoms in d2 for d2 = d2_m1+d2_m2
			d2_elements_way_1                           = d2_m1_reordered_elements + d2_m2_reordered_elements
			d2_positions_way_1                          = np.concatenate([d2_m1_reordered_positions,d2_m2_reordered_positions])
			dimer2_no_of_hydrogens_way_1                = d2_m1_no_H_attached_to_nonH_atoms_reordered + d2_m2_no_H_attached_to_nonH_atoms_reordered

			# Seventh, bound all the inputs needed for the ```are_systems_variant``` method into a single tuple. 
			input_data = (d1_elements, d1_positions, dimer1_no_of_hydrogens, d2_elements_way_1, d2_positions_way_1, dimer2_no_of_hydrogens_way_1, max_distance_disparity, info_about_dimers_being_compared, neighbouring_molecules_about_dimers, non_hydrogen_molecules)

			# Eighth, yield input_data
			yield input_data

def are_systems_variant_parallel(input_data):
	"""
	This method will extracts the data in input_data for the ```are_systems_variant``` method. 

	Parameters
	----------
	d1_elements : list of str.
		These are the elements of the atoms in dimer 1 (All of these atoms should be non-hydrogens).
	d1_positions : list of (1x3) numpy.array
		These are the positions of the atoms in dimer 1 (All of these atoms should be non-hydrogens).
	dimer1_no_of_hydrogens : list of int
		These list contains the number of hydrogens bound to the (non-hydrogen) atoms dimer 1.

	d2_elements_way_1 : list of str.
		These are the elements of the atoms in dimer 2 (created as M1,M2) (All of these atoms should be non-hydrogens).
	d2_positions_way_1 : list of (1x3) numpy.array
		These are the positions of the atoms in dimer 2 (created as M1,M2) (All of these atoms should be non-hydrogens).
	dimer2_no_of_hydrogens_way_1 : list of int
		These list contains the number of hydrogens bound to the (non-hydrogen) atoms dimer 2 (created as M1,M2).

    max_distance_disparity: float
        This is the maximum that any two "could be equivalent" atoms can be between system 1 and system 2 for systems 1 and 2 to be considered variant.
    info_about_dimers_being_compared: list
        These are the names of the molecules and the unit vector displacement that makes up the two dimers being compared. 
    neighbouring_molecules_about_dimers : dict.
        This is the information about the molecules that surround (in the vicinity of) each dimer in the crystal.
    non_hydrogen_molecules : dict. of ase.Atoms
        This is the list of molecules in the crystal, not including hydrogens.

	Returns
	-------
	dimers_are_variants : bool.
		True if the two dimers are variant. False if they are invariant (unique).
	"""

	# First, extract the information from the ``input_data`` tuple. These variables are needed for the ``are_systems_variant`` method. 
	d1_elements, d1_positions, dimer1_no_of_hydrogens, d2_elements_way_1, d2_positions_way_1, dimer2_no_of_hydrogens_way_1, max_distance_disparity, info_about_dimers_being_compared, neighbouring_molecules_about_dimers, non_hydrogen_molecules = input_data

	# Second, run the ``are_systems_variant`` method to determine if the two dimers of interest (for the given index arrangement for dimer 2) are variants of each other. 
	dimers_are_variants = are_systems_variant(d1_elements, d1_positions, dimer1_no_of_hydrogens, d2_elements_way_1, d2_positions_way_1, dimer2_no_of_hydrogens_way_1, max_distance_disparity=max_distance_disparity, names_being_compared=info_about_dimers_being_compared, neighbouring_molecules_about_systems=neighbouring_molecules_about_dimers, non_hydrogen_systems=non_hydrogen_molecules)

	# Third, return if the two dimers are variants of each other or not based on the result from the ``are_systems_variant`` method.
	return dimers_are_variants

# ---------------------------------------------------------------------------------------------------------------------------------------------------

