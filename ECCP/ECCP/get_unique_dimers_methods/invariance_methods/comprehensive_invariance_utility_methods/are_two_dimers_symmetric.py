"""
are_two_dimers_symmetric.py, Geoffrey Weal, 7/2/24

This script is designed to compare two dimers to each other, and determine if they are symmetric to each other. 
"""
import numpy as np
from ctypes import c_bool
from tqdm import trange, tqdm

import multiprocessing as mp
from tqdm.contrib.concurrent import process_map

from copy import deepcopy
from itertools import product

from SUMELF import GraphMatcher, remove_hydrogens

from ECCP.ECCP.invariance_methods.utilities                                                             import get_permutated_indices_list
from ECCP.ECCP.invariance_methods.common_utility_methods_for_all_invariance_methods.are_systems_variant import are_systems_variant

# ----------------------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------

def are_two_dimers_symmetric_way1(em_indices_d2_m1_to_d1_m1, em_indices_d2_m2_to_d1_m2, d2_m1_original_elements, d2_m1_original_positions, d2_m1_no_H_attached_to_nonH_atoms, d2_m2_original_elements, d2_m2_original_positions, d2_m2_no_H_attached_to_nonH_atoms, d1_elements, d1_positions, d1_no_H_attached_to_nonH_atoms, max_distance_disparity, names_of_dimers_being_compared, neighbouring_molecules_about_dimers, non_hydrogen_molecules, no_of_cpus=1):
    """
    This method is designed to compare 
        * molecule 1 in dimer 2 --> WITH --> molecule 1 in dimer 2, and 
        * molecule 2 in dimer 2 --> WITH --> molecule 2 in dimer 2. 

    Parameters
    ----------
    em_indices_d2_m1_to_d1_m1 : dictionary of (name of d2_m1 in molecules dict, name of d1_m1 in molecules dict): list of ways atoms in d2_m1 can map onto atom in d1_m1
        This is a dictionary of the equivalent molecules (em), along with the indices that allow the atoms in molecule 1 of dimer 2 to map onto molecule 1 of dimer 1. 
    em_indices_d2_m2_to_d1_m2 : dictionary of (name of d2_m2 in molecules dict, name of d1_m2 in molecules dict): list of ways atoms in d2_m2 can map onto atom in d1_m2
        This is a dictionary of the equivalent molecules (em), along with the indices that allow the atoms in molecule 2 of dimer 2 to map onto molecule 2 of dimer 1. 

    d2_m1_original_elements : list of str.
        These are the elements that make up molecule 1 of dimer 2. 
    d2_m1_original_positions : 2D np.array
        These are the positions that make up molecule 1 of dimer 2. 
    d2_m1_no_H_attached_to_nonH_atoms : list of int
        These are the number of hydrogen bound to each "heavy" atom in molecule 1 of dimer 2. 

    d2_m2_original_elements : list of str.
        These are the elements that make up molecule 2 of dimer 2. 
    d2_m2_original_positions : 2D np.array
        These are the positions that make up molecule 2 of dimer 2. 
    d2_m2_no_H_attached_to_nonH_atoms : list of int
        These are the number of hydrogen bound to each "heavy" atoms in molecule 2 of dimer 2. 

    d1_elements : list of str.
        These are the elements of the atoms in dimer 1.
    d1_positions : 2D np.array
        These are the positions of the atoms in dimer 1.
    d1_no_H_attached_to_nonH_atoms : list of int
        These are the number of hydrogen bound to each "heavy" atoms in dimer 1.

    max_distance_disparity : float
        This is the maximum difference in relative positions of atoms between dimer 1 and dimer 2 to be considered equivalent/identical
    names_of_dimers_being_compared : tuple of (int, int)
        These are the names of the two dimers being compared in the list of dimers. 
    neighbouring_molecules_about_dimers : dict.
        This is the information about the molecules that surround (in the vicinity of) each dimer in the crystal.
    non_hydrogen_molecules : list of ase.Atoms
        These are Atoms of objects of the molecules in the crystal. These do not contain hydrogen atoms. This is for debugging use only. 

    no_of_cpus : int.
        This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.

    Returns
    -------
    found_dimer_match_way1 : bool.
        This boolean indicates if the two dimers are symmetric to each other or not. 
    """

    # First, inialise the boolean object. This is an object that can be saved to rather than a bool.
    found_dimer_match_way1 = False

    # Second, obtain the generator for the input data.
    input_generator_way1 = compare_two_dimers_generator_way1(em_indices_d2_m1_to_d1_m1, em_indices_d2_m2_to_d1_m2, d2_m1_original_elements, d2_m1_original_positions, d2_m1_no_H_attached_to_nonH_atoms, d2_m2_original_elements, d2_m2_original_positions, d2_m2_no_H_attached_to_nonH_atoms, d1_elements, d1_positions, d1_no_H_attached_to_nonH_atoms, max_distance_disparity, names_of_dimers_being_compared, neighbouring_molecules_about_dimers, non_hydrogen_molecules)

    # Third, perform the compare_two_dimers method on all inputs depending on if you are performing the task with one cpu or with multiple cores.
    if no_of_cpus == 1: # If the user only want to use 1 cpu, dont use multiprocessing tools. 

        # 3.1: For each way that dimer 2 can be mapped onto dimer 1. 
        for input_data in input_generator_way1:
            
            # 3.2: Determine if the two dimers are variants of each other. 
            found_dimer_match_way1 = compare_two_dimers(input_data)

            # 3.3: If are_molecules_variant_comprehensive_single_process throws an exception, break out of the loop
            if isinstance(found_dimer_match_way1, Exception):
                break

            # 3.4: Make sure that at this point, is_variant is a boolean.
            if not isinstance(found_dimer_match_way1, bool):
                break

            # 3.5: If found_dimer_match_way1 == True, we have found the two dimers are the same, so break out of the loop.
            if found_dimer_match_way1 == True:
                break

    else:

        # 3.6: Set up the multiprocessing pool.
        with mp.Pool(processes=no_of_cpus) as pool:

            # 3.7: Compare the two dimers via all the possible comparisons of the atom indices in dimer 2 to the atom indices in dimer 1. 
            for found_dimer_match_way1 in pool.map_async(compare_two_dimers, input_generator_way1):

                # 3.8: If are_molecules_variant_comprehensive_single_process throws an exception, break out of the loop
                if isinstance(found_dimer_match_way1, Exception):
                    break

                # 3.9: Make sure that at this point, is_variant is a boolean.
                if not isinstance(found_dimer_match_way1, bool):
                    break

                # 3.10: If ``found_dimer_match_way1`` is True, break out of the for loop.
                if found_dimer_match_way1 == True:
                    break

    # Fourth, if there was an exception, indicate what it is.
    if isinstance(found_dimer_match_way1, Exception):
        exception_message  = 'Error: There was a problem in the "are_molecules_variant_comprehensive_single_process" method somewhere\n\n'
        exception_message += str(found_dimer_match_way1)
        raise Exception(exception_message)

    # Fifth, if found_dimer_match_way1 is not a boolean at this point, raise this as a problem.
    if not isinstance(found_dimer_match_way1, bool):
        exception_message  = 'Error: The "found_dimer_match_way1" variable is not a boolean.\n'
        exception_message += '"found_dimer_match_way1" = '+str(found_dimer_match_way1)+'\n'
        exception_message += 'Check this.'
        raise Exception(exception_message)

    # Sixth, return found_dimer_match_way1.
    return found_dimer_match_way1

def compare_two_dimers_generator_way1(em_indices_d2_m1_to_d1_m1, em_indices_d2_m2_to_d1_m2, d2_m1_original_elements, d2_m1_original_positions, d2_m1_no_H_attached_to_nonH_atoms, d2_m2_original_elements, d2_m2_original_positions, d2_m2_no_H_attached_to_nonH_atoms, d1_elements, d1_positions, d1_no_H_attached_to_nonH_atoms, max_distance_disparity, names_of_dimers_being_compared, neighbouring_molecules_about_dimers, non_hydrogen_molecules):
    """
    This method is a generator for the compare_two_dimers method

    This generator is specifically designed to allow the user to compare:
        * molecule 1 in dimer 2 --> WITH --> molecule 1 in dimer 1, and 
        * molecule 2 in dimer 2 --> WITH --> molecule 2 in dimer 1. 

    Parameters
    ----------
    em_indices_d2_m1_to_d1_m1 : dictionary of (name of d2_m1 in molecules dict, name of d1_m1 in molecules dict): list of ways atoms in d2_m1 can map onto atom in d1_m1
        This is a dictionary of the equivalent molecules (em), along with the indices that allow the atoms in molecule 1 of dimer 2 to map onto molecule 1 of dimer 1. 
    em_indices_d2_m2_to_d1_m2 : dictionary of (name of d2_m2 in molecules dict, name of d1_m2 in molecules dict): list of ways atoms in d2_m2 can map onto atom in d1_m2
        This is a dictionary of the equivalent molecules (em), along with the indices that allow the atoms in molecule 2 of dimer 2 to map onto molecule 2 of dimer 1. 

    d2_m1_original_elements : list of str.
        These are the elements that make up molecule 1 of dimer 2. 
    d2_m1_original_positions : 2D np.array
        These are the positions that make up molecule 1 of dimer 2. 
    d2_m1_no_H_attached_to_nonH_atoms : list of int
        These are the number of hydrogen bound to each "heavy" atom in molecule 1 of dimer 2. 

    d2_m2_original_elements : list of str.
        These are the elements that make up molecule 2 of dimer 2. 
    d2_m2_original_positions : 2D np.array
        These are the positions that make up molecule 2 of dimer 2. 
    d2_m2_no_H_attached_to_nonH_atoms : list of int
        These are the number of hydrogen bound to each "heavy" atoms in molecule 2 of dimer 2. 

    d1_elements : list of str.
        These are the elements of the atoms in dimer 1.
    d1_positions : 2D np.array
        These are the positions of the atoms in dimer 1.
    d1_no_H_attached_to_nonH_atoms : list of int
        These are the number of hydrogen bound to each "heavy" atoms in dimer 1.

    max_distance_disparity : float
        This is the maximum difference in relative positions of atoms between dimer 1 and dimer 2 to be considered equivalent/identical
    names_of_dimers_being_compared : tuple of (int, int)
        These are the names of the two dimers being compared in the list of dimers. 
    neighbouring_molecules_about_dimers : dict.
        This is the information about the molecules that surround (in the vicinity of) each dimer in the crystal.
    non_hydrogen_molecules : list of ase.Atoms
        These are Atoms of objects of the molecules in the crystal. These do not contain hydrogen atoms. This is for debugging use only. 

    found_dimer_match : Boolean Object
        This object record if the program found that the two dimers are in fact symmetrical or not. 

    Returns
    -------
    d1_elements : list of str.
        These are the elements of the atoms in dimer 1.
    d1_positions : 2D np.array
        These are the positions of the atoms in dimer 1.
    d1_no_H_attached_to_nonH_atoms : list of int
        These are the number of hydrogen bound to each "heavy" atoms in dimer 1.

    d2_elements_way_1 : list of str.
        These are the elements of the atoms in dimer 2.
    d2_positions_way_1 : 2D np.array
        These are the positions of the atoms in dimer 2.
    d2_no_H_attached_to_nonH_atoms_way_1 : list of int
        These are the number of hydrogen bound to each "heavy" atoms in dimer 2.

    max_distance_disparity : float
        This is the maximum difference in relative positions of atoms between dimer 1 and dimer 2 to be considered equivalent/identical
    names_of_dimers_being_compared : tuple of (int, int)
        These are the names of the two dimers being compared in the list of dimers. 
    neighbouring_molecules_about_dimers : dict.
        This is the information about the molecules that surround (in the vicinity of) each dimer in the crystal.
    non_hydrogen_molecules : list of ase.Atoms
        These are Atoms of objects of the molecules in the crystal. These do not contain hydrogen atoms. This is for debugging use only. 
    """

    # First, for each way that molecule 1 of dimer 2 can be mapped onto molecule 1 or dimer 1.
    for comparison1 in em_indices_d2_m1_to_d1_m1:

        # Second, Rearrange the atoms in molecule 1 of dimer 2. This allows the indices of molecule 1 of dimer 2 to be match to molecule 1 of dimer 1.
        idx_d2_m1                                       = get_permutated_indices_list(comparison1)
        d2_m1_reordered_elements                        = [d2_m1_original_elements[index] for index in idx_d2_m1]
        d2_m1_reordered_positions                       = deepcopy(d2_m1_original_positions)[idx_d2_m1, :]
        d2_m1_reordered_no_H_attached_to_nonH_atoms     = [d2_m1_no_H_attached_to_nonH_atoms[index] for index in idx_d2_m1]

        # Third, for each way that molecule 2 of dimer 2 can be mapped onto molecule 2 or dimer 1.
        for comparison2 in em_indices_d2_m2_to_d1_m2:

            # Fourth, rearrange the atoms in molecule 2 of dimer 2. This allows the indices of molecule 2 of dimer 2 to be match to molecule 2 of dimer 1.
            idx_d2_m2                                   = get_permutated_indices_list(comparison2)
            d2_m2_reordered_elements                    = [d2_m2_original_elements[index] for index in idx_d2_m2]
            d2_m2_reordered_positions                   = deepcopy(d2_m2_original_positions)[idx_d2_m2, :]
            d2_m2_reordered_no_H_attached_to_nonH_atoms = [d2_m2_no_H_attached_to_nonH_atoms[index] for index in idx_d2_m2]

            # Fifth, get the elements and positions of atoms in d2 for d2 = d2_m1 + d2_m2 (this is because for this comparison, d2_m1 goes with d1_m1, and d2_m2 goes with d1_m2)
            d2_elements_way_1                    = d2_m1_reordered_elements + d2_m2_reordered_elements
            d2_positions_way_1                   = np.concatenate([d2_m1_reordered_positions,d2_m2_reordered_positions]) 
            d2_no_H_attached_to_nonH_atoms_way_1 = d2_m1_reordered_no_H_attached_to_nonH_atoms + d2_m2_reordered_no_H_attached_to_nonH_atoms

            # Sixth, yield input variables for the compare_two_dimers method. 
            yield (d1_elements, d1_positions, d1_no_H_attached_to_nonH_atoms), (d2_elements_way_1, d2_positions_way_1, d2_no_H_attached_to_nonH_atoms_way_1), max_distance_disparity, names_of_dimers_being_compared, neighbouring_molecules_about_dimers, non_hydrogen_molecules

# ----------------------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------

def are_two_dimers_symmetric_way2(em_indices_d2_m1_to_d1_m2, em_indices_d2_m2_to_d1_m1, d2_m1_original_elements, d2_m1_original_positions, d2_m1_no_H_attached_to_nonH_atoms, d2_m2_original_elements, d2_m2_original_positions, d2_m2_no_H_attached_to_nonH_atoms, d1_elements, d1_positions, d1_no_H_attached_to_nonH_atoms, max_distance_disparity, names_of_dimers_being_compared, neighbouring_molecules_about_dimers, non_hydrogen_molecules, no_of_cpus=1):
    """
    This method is designed to compare 
        * molecule 1 in dimer 2 --> WITH --> molecule 2 in dimer 1, and 
        * molecule 2 in dimer 2 --> WITH --> molecule 1 in dimer 1. 

    Parameters
    ----------
    em_indices_d2_m1_to_d1_m2 : dictionary of (name of d2_m1 in molecules dict, name of d1_m2 in molecules dict): list of ways atoms in d2_m1 can map onto atom in d1_m2
        This is a dictionary of the equivalent molecules (em), along with the indices that allow the atoms in molecule 1 of dimer 2 to map onto molecule 2 of dimer 1. 
    em_indices_d2_m2_to_d1_m1 : dictionary of (name of d2_m2 in molecules dict, name of d1_m1 in molecules dict): list of ways atoms in d2_m2 can map onto atom in d1_m1
        This is a dictionary of the equivalent molecules (em), along with the indices that allow the atoms in molecule 2 of dimer 2 to map onto molecule 1 of dimer 1. 

    d2_m1_original_elements : list of str.
        These are the elements that make up molecule 1 of dimer 2. 
    d2_m1_original_positions : 2D np.array
        These are the positions that make up molecule 1 of dimer 2. 
    d2_m1_no_H_attached_to_nonH_atoms : list of int
        These are the number of hydrogen bound to each "heavy" atom in molecule 1 of dimer 2. 

    d2_m2_original_elements : list of str.
        These are the elements that make up molecule 2 of dimer 2. 
    d2_m2_original_positions : 2D np.array
        These are the positions that make up molecule 2 of dimer 2. 
    d2_m2_no_H_attached_to_nonH_atoms : list of int
        These are the number of hydrogen bound to each "heavy" atoms in molecule 2 of dimer 2. 

    d1_elements : list of str.
        These are the elements of the atoms in dimer 1.
    d1_positions : 2D np.array
        These are the positions of the atoms in dimer 1.
    d1_no_H_attached_to_nonH_atoms : list of int
        These are the number of hydrogen bound to each "heavy" atoms in dimer 1.

    max_distance_disparity : float
        This is the maximum difference in relative positions of atoms between dimer 1 and dimer 2 to be considered equivalent/identical
    names_of_dimers_being_compared : tuple of (int, int)
        These are the names of the two dimers being compared in the list of dimers. 
    neighbouring_molecules_about_dimers : dict.
        This is the information about the molecules that surround (in the vicinity of) each dimer in the crystal.
    non_hydrogen_molecules : list of ase.Atoms
        These are Atoms of objects of the molecules in the crystal. These do not contain hydrogen atoms. This is for debugging use only. 
    no_of_cpus : int.
        This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.

    Returns
    -------
    found_dimer_match_way2 : bool.
        This boolean indicates if the two dimers are symmetric to each other or not. 
    """

    # First, inialise the boolean object. This is an object that can be saved to rather than a bool.
    found_dimer_match_way2 = False

    # Second, obtain the generator for the input data.
    input_generator_way2 = compare_two_dimers_generator_way2(em_indices_d2_m1_to_d1_m2, em_indices_d2_m2_to_d1_m1, d2_m1_original_elements, d2_m1_original_positions, d2_m1_no_H_attached_to_nonH_atoms, d2_m2_original_elements, d2_m2_original_positions, d2_m2_no_H_attached_to_nonH_atoms, d1_elements, d1_positions, d1_no_H_attached_to_nonH_atoms, max_distance_disparity, names_of_dimers_being_compared, neighbouring_molecules_about_dimers, non_hydrogen_molecules)

    # Third, perform the compare_two_dimers method on all inputs depending on if you are performing the task with one cpu or with multiple cores.
    if no_of_cpus == 1: # If the user only want to use 1 cpu, dont use multiprocessing tools. 

        # 3.1: For each way that dimer 2 can be mapped onto dimer 1. 
        for input_data in input_generator_way2:
            
            # 3.2: Determine if the two dimers are variants of each other. 
            found_dimer_match_way2 = compare_two_dimers(input_data)

            # 3.3: If are_molecules_variant_comprehensive_single_process throws an exception, break out of the loop
            if isinstance(found_dimer_match_way2, Exception):
                break

            # 3.4: Make sure that at this point, is_variant is a boolean.
            if not isinstance(found_dimer_match_way2, bool):
                break

            # 3.5: If found_dimer_match == True, we have found the two dimers are the same, so break out of the loop.
            if found_dimer_match_way2 == True:
                break

    else:

        # 3.6: Set up the multiprocessing pool.
        with mp.Pool(processes=no_of_cpus) as pool:

            # 3.7: Compare the two dimers via all the possible comparisons of the atom indices in dimer 2 to the atom indices in dimer 1. 
            for found_dimer_match_way2 in pool.imap_unordered(compare_two_dimers, input_generator_way2):

                # 3.8: If are_molecules_variant_comprehensive_single_process throws an exception, break out of the loop
                if isinstance(found_dimer_match_way2, Exception):
                    break

                # 3.9: Make sure that at this point, is_variant is a boolean.
                if not isinstance(found_dimer_match_way2, bool):
                    break

                # 3.10: If ``found_dimer_match_way2`` is True, break out of the for loop.
                if found_dimer_match_way2 == True:
                    break

    # Fourth, if there was an exception, indicate what it is.
    if isinstance(found_dimer_match_way2, Exception):
        exception_message  = 'Error: There was a problem in the "are_molecules_variant_comprehensive_single_process" method somewhere\n\n'
        exception_message += str(found_dimer_match_way2)
        raise Exception(exception_message)

    # Fifth, if found_dimer_match_way2 is not a boolean at this point, raise this as a problem.
    if not isinstance(found_dimer_match_way2, bool):
        exception_message  = 'Error: The "found_dimer_match_way2" variable is not a boolean.\n'
        exception_message += '"found_dimer_match_way2" = '+str(found_dimer_match_way2)+'\n'
        exception_message += 'Check this.'
        raise Exception(exception_message)

    # Sixth, return found_dimer_match_way2
    return found_dimer_match_way2

def compare_two_dimers_generator_way2(em_indices_d2_m1_to_d1_m2, em_indices_d2_m2_to_d1_m1, d2_m1_original_elements, d2_m1_original_positions, d2_m1_no_H_attached_to_nonH_atoms, d2_m2_original_elements, d2_m2_original_positions, d2_m2_no_H_attached_to_nonH_atoms, d1_elements, d1_positions, d1_no_H_attached_to_nonH_atoms, max_distance_disparity, names_of_dimers_being_compared, neighbouring_molecules_about_dimers, non_hydrogen_molecules):
    """
    This method is a generator for the compare_two_dimers method

    This generator is specifically designed to allow the user to compare:
        * molecule 1 in dimer 2 --> WITH --> molecule 2 in dimer 1, and 
        * molecule 2 in dimer 2 --> WITH --> molecule 1 in dimer 1. 

    Parameters
    ----------
    em_indices_d2_m1_to_d1_m2 : dictionary of (name of d2_m1 in molecules list, name of d1_m2 in molecules list): list of ways atoms in d2_m1 can map onto atom in d1_m2
        This is a dictionary of the equivalent molecules (em), along with the indices that allow the atoms in molecule 1 of dimer 2 to map onto molecule 2 of dimer 1. 
    em_indices_d2_m2_to_d1_m1 : dictionary of (name of d2_m2 in molecules list, name of d1_m1 in molecules list): list of ways atoms in d2_m2 can map onto atom in d1_m1
        This is a dictionary of the equivalent molecules (em), along with the indices that allow the atoms in molecule 2 of dimer 2 to map onto molecule 1 of dimer 1. 

    d2_m1_original_elements : list of str.
        These are the elements that make up molecule 1 of dimer 2. 
    d2_m1_original_positions : 2D np.array
        These are the positions that make up molecule 1 of dimer 2. 
    d2_m1_no_H_attached_to_nonH_atoms : list of int
        These are the number of hydrogen bound to each "heavy" atom in molecule 1 of dimer 2. 

    d2_m2_original_elements : list of str.
        These are the elements that make up molecule 2 of dimer 2. 
    d2_m2_original_positions : 2D np.array
        These are the positions that make up molecule 2 of dimer 2. 
    d2_m2_no_H_attached_to_nonH_atoms : list of int
        These are the number of hydrogen bound to each "heavy" atoms in molecule 2 of dimer 2. 

    d1_elements : list of str.
        These are the elements of the atoms in dimer 1.
    d1_positions : 2D np.array
        These are the positions of the atoms in dimer 1.
    d1_no_H_attached_to_nonH_atoms : list of int
        These are the number of hydrogen bound to each "heavy" atoms in dimer 1.

    max_distance_disparity : float
        This is the maximum difference in relative positions of atoms between dimer 1 and dimer 2 to be considered equivalent/identical
    names_of_dimers_being_compared : tuple of (int, int)
        These are the names of the two dimers being compared in the list of dimers. 
    neighbouring_molecules_about_dimers : dict.
        This is the information about the molecules that surround (in the vicinity of) each dimer in the crystal.
    non_hydrogen_molecules : list of ase.Atoms
        These are Atoms of objects of the molecules in the crystal. These do not contain hydrogen atoms. This is for debugging use only. 

    Returns
    -------
    d1_elements : list of str.
        These are the elements of the atoms in dimer 1.
    d1_positions : 2D np.array
        These are the positions of the atoms in dimer 1.
    d1_no_H_attached_to_nonH_atoms : list of int
        These are the number of hydrogen bound to each "heavy" atoms in dimer 1.

    d2_elements_way_1 : list of str.
        These are the elements of the atoms in dimer 2.
    d2_positions_way_1 : 2D np.array
        These are the positions of the atoms in dimer 2.
    d2_no_H_attached_to_nonH_atoms_way_1 : list of int
        These are the number of hydrogen bound to each "heavy" atoms in dimer 2.

    max_distance_disparity : float
        This is the maximum difference in relative positions of atoms between dimer 1 and dimer 2 to be considered equivalent/identical
    names_of_dimers_being_compared : tuple of (int, int)
        These are the indices of the two dimers being compared in the list of dimers. 
    neighbouring_molecules_about_dimers : dict.
        This is the information about the molecules that surround (in the vicinity of) each dimer in the crystal.
    non_hydrogen_molecules : list of ase.Atoms
        These are Atoms of objects of the molecules in the crystal. These do not contain hydrogen atoms. This is for debugging use only. 
    """

    # First, for each way that molecule 1 of dimer 2 can be mapped onto molecule 2 or dimer 1.
    for comparison1 in em_indices_d2_m1_to_d1_m2:

        # Second, Rearrange the atoms in molecule 1 of dimer 2. This allows the indices of molecule 1 of dimer 2 to be match to molecule 2 of dimer 1.
        idx_d2_m1                                       = get_permutated_indices_list(comparison1)
        d2_m1_reordered_elements                        = [d2_m1_original_elements[index] for index in idx_d2_m1]
        d2_m1_reordered_positions                       = deepcopy(d2_m1_original_positions)[idx_d2_m1, :]
        d2_m1_reordered_no_H_attached_to_nonH_atoms     = [d2_m1_no_H_attached_to_nonH_atoms[index] for index in idx_d2_m1]

        # Third, for each way that molecule 2 of dimer 2 can be mapped onto molecule 1 or dimer 1.
        for comparison2 in em_indices_d2_m2_to_d1_m1:

            # Fourth, rearrange the atoms in molecule 2 of dimer 2. This allows the indices of molecule 2 of dimer 2 to be match to molecule 1 of dimer 1.
            idx_d2_m2                                   = get_permutated_indices_list(comparison2)
            d2_m2_reordered_elements                    = [d2_m2_original_elements[index] for index in idx_d2_m2]
            d2_m2_reordered_positions                   = deepcopy(d2_m2_original_positions)[idx_d2_m2, :]
            d2_m2_reordered_no_H_attached_to_nonH_atoms = [d2_m2_no_H_attached_to_nonH_atoms[index] for index in idx_d2_m2]

            # Fifth, get the elements and positions of atoms in d2 for d2 = d2_m2 + d2_m1 (this is because d2_m2 goes with d1_m1 here, and d2_m1 goes with d1_m2)
            d2_elements_way_2                    = d2_m2_reordered_elements + d2_m1_reordered_elements
            d2_positions_way_2                   = np.concatenate([d2_m2_reordered_positions,d2_m1_reordered_positions]) 
            d2_no_H_attached_to_nonH_atoms_way_2 = d2_m2_reordered_no_H_attached_to_nonH_atoms + d2_m1_reordered_no_H_attached_to_nonH_atoms

            # Sixth, yield input variables for the compare_two_dimers method. 
            yield (d1_elements, d1_positions, d1_no_H_attached_to_nonH_atoms), (d2_elements_way_2, d2_positions_way_2, d2_no_H_attached_to_nonH_atoms_way_2), max_distance_disparity, names_of_dimers_being_compared, neighbouring_molecules_about_dimers, non_hydrogen_molecules

# ----------------------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------

def compare_two_dimers(input_data):
    """
    This method is designed to compare two dimers together to determine if they are invariant or not. 

    Parameters
    ----------
    d1_elements : list of str.
        This is the list of elements in dimer 1
    d1_positions : 2D np.array
        This is the positions of elements in dimer 1
    d1_no_H_attached_to_nonH_atoms : list of ints
        This is the list of hydrogens bound to each "heavy" atom in dimer 1

    d2_elements : list of str.
        This is the list of elements in dimer 2
    d2_positions : 2D np.array
        This is the positions of elements in dimer 2
    d2_no_H_attached_to_nonH_atoms : list of ints
        This is the list of hydrogens bound to each "heavy" atom in dimer 2

    max_distance_disparity : float
        This is the maximum difference in relative positions of atoms between dimer 1 and dimer 2 to be considered equivalent/identical
    names_of_dimers_being_compared : tuple of (int, int)
        These are the names of the two dimers being compared in the list of dimers. 
    neighbouring_molecules_about_dimers : dict.
        This is the information about the molecules that surround (in the vicinity of) each dimer in the crystal.
    non_hydrogen_molecules : list of ase.Atoms
        These are Atoms of objects of the molecules in the crystal. These do not contain hydrogen atoms. This is for debugging use only. 

    Returns
    -------
    are_dimers_variant : bool
        Returns if the two dimers are variants of each other.
    """

    # First, separate all the variables we need from input_data
    (d1_elements, d1_positions, d1_no_H_attached_to_nonH_atoms), (d2_elements, d2_positions, d2_no_H_attached_to_nonH_atoms), max_distance_disparity, names_of_dimers_being_compared, neighbouring_molecules_about_dimers, non_hydrogen_molecules = input_data

    # Second, determine if this configuation of atom indice positions in d2 allow for invarience with d1 (using the comprehensive method).
    are_dimers_variant = are_systems_variant(d1_elements, d1_positions, d1_no_H_attached_to_nonH_atoms, d2_elements, d2_positions, d2_no_H_attached_to_nonH_atoms, max_distance_disparity=max_distance_disparity, names_being_compared=names_of_dimers_being_compared, neighbouring_molecules_about_systems=neighbouring_molecules_about_dimers, non_hydrogen_systems=non_hydrogen_molecules)

    # Third, return are_dimers_variant
    return are_dimers_variant

# ----------------------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------

