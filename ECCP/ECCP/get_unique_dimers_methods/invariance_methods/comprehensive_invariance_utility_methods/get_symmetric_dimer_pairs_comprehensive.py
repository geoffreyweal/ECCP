"""
get_symmetric_dimer_pairs_comprehensive.py, Geoffrey Weal, 7/2/24

This script is designed to determine which dimers in the dimers list are symmetric to each other.

"""
import sys
import numpy as np
from ctypes import c_bool
from tqdm import trange, tqdm

import multiprocessing as mp
from tqdm.contrib.concurrent import process_map

from copy import deepcopy
from itertools import product

from SUMELF import GraphMatcher, remove_hydrogens

from ECCP.ECCP.get_unique_dimers_methods.invariance_methods.comprehensive_invariance_utility_methods.are_two_dimers_symmetric import are_two_dimers_symmetric_way1, are_two_dimers_symmetric_way2

def get_symmetric_dimer_pairs_comprehensive(dimers, non_hydrogen_molecules_elements, non_hydrogen_molecules_positions, all_no_of_H_on_atoms_in_molecule, equivalent_molecule_atom_indices_comparison, max_distance_disparity, neighbouring_molecules_about_dimers, non_hydrogen_molecules, no_of_cpus=1):
    """
    This method is designed to determine which dimers in the dimers list are equivalent to each other.

    Parameters
    ----------
    dimers : list
        This is a list of dimer information, given as (name of dimer, index of molecule 1 in dimer, index of molecule 2 in dimer, unit cell ijk displacement of molecule 2, unit cell displacement of molecule 2, displacement of dimer COM)
    non_hydrogen_molecules_elements : dict. of lists of str. 
        These are the element of all the molecules that can make up the dimers.
    non_hydrogen_molecules_positions :  dict. of 2d numpy.arrays
        These are the positions of all the molecules that can make up the dimers.
    all_no_of_H_on_atoms_in_molecule :  dict. of lists of int. 
        These are the number of hydrogens bound to each "heavy" atom in each molecule that can make up the dimers. 
    equivalent_molecule_atom_indices_comparison : dict.
        This is all the ways that two molecules in non_hydrogen_graphs can map onto each other. 
    max_distance_disparity : float
        This is the maximum difference in relative positions of atoms between dimer 1 and dimer 2 to be considered equivalent/identical
    neighbouring_molecules_about_dimers : dict.
        This is the information about the molecules that surround (in the vicinity of) each dimer in the crystal.
    non_hydrogen_molecules : dict. of ase.Atoms
        These are Atoms of objects of the molecules in the crystal. These do not contain hydrogen atoms. This is for debugging use only. 
    no_of_cpus : int.
        This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.

    Returns
    -------
    symmetric_dimer_pairs : list object (either a list or a multiprocessing.Manager.list)
        This list is designed to store all the symmetric dimers identified from this program.
    """

    # Prestep: Set the number of cores for each comparison of two molecules
    no_of_cpus_for_comparing_two_molecules = no_of_cpus

    # First, tell the user what you are doing. 
    print('Comparing translational, rotational, and reflective invarience between dimers. '+str(len(dimers))+' dimers to be examined. This can take a while with large and complex dimers.')

    # Second, get the number of comparisons that will be performed. 
    nn = int((len(dimers)*(len(dimers)-1))/2)

    if True: # no_of_cpus == 1: # If you are using a single cpu, do not use multiprocessing methods

        # Third, initialise the list to store all the symmetric dimers in.
        symmetric_dimer_pairs = []
    
        # Forth, initialise the tqdm progress bar to show the user the progress
        pbar = tqdm(get_inputs(dimers, non_hydrogen_molecules_elements, non_hydrogen_molecules_positions, all_no_of_H_on_atoms_in_molecule, equivalent_molecule_atom_indices_comparison, max_distance_disparity, neighbouring_molecules_about_dimers, non_hydrogen_molecules, symmetric_dimer_pairs, no_of_cpus), total=nn, unit='dimer pair')

        # Fifth, for each pair of dimers in the dimers list to compare.
        for input_data in pbar:

            # 5.1: Get information about the dimers being compared
            (dimer1_name, d1_m1_name, d1_m2_name) = input_data[0]
            (dimer2_name, d2_m1_name, d2_m2_name) = input_data[1]

            # 5.2: Send a message about what we are doing to the screen for the user
            pbar.set_description('Comparing dimers '+str(dimer1_name)+' (M'+str(d1_m1_name)+',M'+str(d1_m2_name)+') and '+str(dimer2_name)+' (M'+str(d2_m1_name)+',M'+str(d2_m2_name)+')')

            # 5.3: Determin if the two dimers described in input_data are symmetric or not.
            get_symmetric_dimer_pairs_single_process(input_data)

        # Sixth, close the progress bar. 
        pbar.close()

    else:

        # Seventh, create the manager to save lists to
        with mp.Manager() as manager:

            # Eighth, inialise the boolean object. This is an object that can be saved to rather than a bool.
            symmetric_dimer_pairs = manager.list()

            # Ninth, write a message for the user
            print('Comparing dimers (Please wait until after 100%, as the process will still be running.)', file=sys.stderr)

            # Tenth, perform compare_two_dimers on each way that dimer 2 can be mapped onto dimer 1 using multiprocessing. 
            #process_map(get_symmetric_dimer_pairs_single_process, get_inputs(dimers, non_hydrogen_molecules_elements, non_hydrogen_molecules_positions, all_no_of_H_on_atoms_in_molecule, equivalent_molecule_atom_indices_comparison, max_distance_disparity, neighbouring_molecules_about_dimers, non_hydrogen_molecules, symmetric_dimer_pairs), total=nn, unit='dimer pair', desc="Comparing dimers", max_workers=no_of_cpus)
            pool = mp.Pool(processes=no_of_cpus)
            pool.map_async(get_symmetric_dimer_pairs_single_process, tqdm(get_inputs(dimers, non_hydrogen_molecules_elements, non_hydrogen_molecules_positions, all_no_of_H_on_atoms_in_molecule, equivalent_molecule_atom_indices_comparison, max_distance_disparity, neighbouring_molecules_about_dimers, non_hydrogen_molecules, symmetric_dimer_pairs, no_of_cpus), total=nn, unit='dimer pair', desc="Comparing dimers"))
            pool.close()
            pool.join()

            # Eleventh, get the bool. value from found_dimer_match
            symmetric_dimer_pairs = list(symmetric_dimer_pairs)

    # Twelfth, sort the dimers.
    symmetric_dimer_pairs.sort()

    # Return symmetric_dimer_pairs
    return symmetric_dimer_pairs

# ----------------------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------

shortest_distance_between_dimers_max_threshold_limit = 0.1
def get_symmetric_dimer_pairs_single_process(input_data):
    """
    This method is designed to determine which dimers in the dimers list are equivalent to each other.

    Parameters
    ----------
    dimer1_name : int
        This is the name in the dimers list for the first dimer. 
    d1_m1_name : int
        This is the name of molecule 1 in the dimer. This is the molecules name in the molecules lists, (such as non_hydrogen_molecules_elements, non_hydrogen_molecules_positions, all_no_of_H_on_atoms_in_molecule).
    d1_m2_name : int
        This is the name of molecule 1 in the dimer. This is the molecules name in the molecules lists, (such as non_hydrogen_molecules_elements, non_hydrogen_molecules_positions, all_no_of_H_on_atoms_in_molecule).

    dimer2_name : int
        This is the name in the dimers list for the second dimer. 
    d2_m1_name : int
        This is the name of molecule 1 in the dimer. This is the molecules name in the molecules lists, (such as non_hydrogen_molecules_elements, non_hydrogen_molecules_positions, all_no_of_H_on_atoms_in_molecule).
    d2_m2_name : int
        This is the name of molecule 1 in the dimer. This is the molecules name in the molecules lists, (such as non_hydrogen_molecules_elements, non_hydrogen_molecules_positions, all_no_of_H_on_atoms_in_molecule).

    d1_elements : list of str.
        This is the list of elements of atoms in dimer 1.
    d1_positions : list of np.array
        This is the list of positions of atoms in dimer 1.
    d1_no_H_attached_to_nonH_atoms : list of int
        This is the list of hydrogens on each heavy atom in dimer 1.

    d2_m1_original_elements : list of str.
        This is the list of elements of atoms in the first molecule in dimer 2.
    d2_m1_original_positions : list of np.array
        This is the list of positions of atoms in the first molecule in dimer 2.
    d2_m1_no_H_attached_to_nonH_atoms : list of int
        This is the list of hydrogens attached to "heavy" atoms in the first molecule in dimer 2.

    d2_m2_original_elements : list of str.
        This is the list of elements of atoms in the second molecule in dimer 2.
    d2_m2_original_positions : list of np.array
        This is the list of positions of atoms in the second molecule in dimer 2.
    d2_m2_no_H_attached_to_nonH_atoms : list of int
        This is the list of hydrogens attached to "heavy" atoms in the second molecule in dimer 2.

    equivalent_molecule_atom_indices_comparison : dict.
        This is all the ways that two molecules in non_hydrogen_graphs can map onto each other. 
    max_distance_disparity : float
        This is the maximum difference in relative positions of atoms between dimer 1 and dimer 2 to be considered equivalent/identical
    neighbouring_molecules_about_dimers : dict.
        This is the information about the molecules that surround (in the vicinity of) each dimer in the crystal.
    non_hydrogen_molecules : dict of ase.Atoms
        These are Atoms of objects of the molecules in the crystal. These do not contain hydrogen atoms. This is for debugging use only. 

    symmetric_dimer_pairs : list object (either a list or a multiprocessing.Manager.list)
        This list is designed to store all the symmetric dimers identified from this program.

    no_of_cpus : int.
        This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.
    """

    # First, obtain all the details needed for comparing the dimers together
    dimer_1_molecule_indices, dimer_2_molecule_indices, dimer_1_details, dimer_2_molecule_1_details, dimer_2_molecule_2_details, distances_between_molecules_in_dimers, other_information, symmetric_dimer_pairs, no_of_cpus = input_data
    dimer1_name, d1_m1_name, d1_m2_name                                                                                              = dimer_1_molecule_indices
    dimer2_name, d2_m1_name, d2_m2_name                                                                                              = dimer_2_molecule_indices
    d1_elements, d1_positions, d1_no_H_attached_to_nonH_atoms                                                                        = dimer_1_details
    d2_m1_original_elements, d2_m1_original_positions, d2_m1_no_H_attached_to_nonH_atoms                                             = dimer_2_molecule_1_details
    d2_m2_original_elements, d2_m2_original_positions, d2_m2_no_H_attached_to_nonH_atoms                                             = dimer_2_molecule_2_details
    d1_shortest_distance, d2_shortest_distance                                                                                       = distances_between_molecules_in_dimers
    equivalent_molecule_atom_indices_comparison, max_distance_disparity, neighbouring_molecules_about_dimers, non_hydrogen_molecules = other_information

    # Second, check if the shortest distance between molecules in the two dimers within the given max threshold limit.
    #         * If they are not, then the two dimers are definitely not the same.
    #         * This is not necessary to the functioning of the method, but speeds things up a lot of it is obvious that two dimers are not equivalent.
    if abs(d1_shortest_distance - d2_shortest_distance) > shortest_distance_between_dimers_max_threshold_limit:
        return

    # Third, get all the indices that are equivalent to eachother in each molecule in each dimer. 
    em_indices_d2_m1_to_d1_m1 = equivalent_molecule_atom_indices_comparison[(d2_m1_name,d1_m1_name)]
    em_indices_d2_m2_to_d1_m2 = equivalent_molecule_atom_indices_comparison[(d2_m2_name,d1_m2_name)]

    # Fourth, make a tuple to fold the indices to be compared.  
    names_of_dimers_being_compared = (dimer1_name, dimer2_name)

    # Fifth, in each of the next two for loops
    #        1. Get the permutation list that indices how the indices of d2_m1 go with d1_m1
    #        2. permutate the indices in d2_m1. This will not change the molecule but only the indices of each atom in d2_m1 (and likewise for d2_m2) 
    found_dimer_match = are_two_dimers_symmetric_way1(em_indices_d2_m1_to_d1_m1, em_indices_d2_m2_to_d1_m2, d2_m1_original_elements, d2_m1_original_positions, d2_m1_no_H_attached_to_nonH_atoms, d2_m2_original_elements, d2_m2_original_positions, d2_m2_no_H_attached_to_nonH_atoms, d1_elements, d1_positions, d1_no_H_attached_to_nonH_atoms, max_distance_disparity, names_of_dimers_being_compared, neighbouring_molecules_about_dimers, non_hydrogen_molecules, no_of_cpus=1)

    # Sixth, if d1 is found to be varient with d2 at this point, move on to next set of dimer pairs
    if found_dimer_match:
        symmetric_dimer_pairs.append((dimer1_name, dimer2_name))
        return 

    #####################################################################################################################################################################
    # If you get to this point, we may need to reverse the position ordering of d2_m1 and d2_m2 so that d2 = d2_m2+d2_m1 to find varience, rather than d2 = d2_m1+d2_m2 #
    #####################################################################################################################################################################

    # Seventh, get all the indices that are equivalent to eachother in each molecule in each dimer. 
    em_indices_d2_m1_to_d1_m2 = equivalent_molecule_atom_indices_comparison[(d2_m1_name,d1_m2_name)]
    em_indices_d2_m2_to_d1_m1 = equivalent_molecule_atom_indices_comparison[(d2_m2_name,d1_m1_name)]

    # Eighth, in each of the next two for loops.
    #         1. Get the permutation list that indices how the indices of d2_m1 go with d1_m1
    #         2. permutate the indices in d2_m1. This will not change the molecule but only the indices of each atom in d2_m1 (and likewise for d2_m2) 
    found_dimer_match = are_two_dimers_symmetric_way2(em_indices_d2_m1_to_d1_m2, em_indices_d2_m2_to_d1_m1, d2_m1_original_elements, d2_m1_original_positions, d2_m1_no_H_attached_to_nonH_atoms, d2_m2_original_elements, d2_m2_original_positions, d2_m2_no_H_attached_to_nonH_atoms, d1_elements, d1_positions, d1_no_H_attached_to_nonH_atoms, max_distance_disparity, names_of_dimers_being_compared, neighbouring_molecules_about_dimers, non_hydrogen_molecules, no_of_cpus=1)

    # Ninth, these two dimers are varients of each other, record and break out of loop.
    if found_dimer_match:
        symmetric_dimer_pairs.append((dimer1_name, dimer2_name))

# ----------------------------------------------------------------------------------------------------------------------------------------------------------

def get_inputs(dimers, non_hydrogen_molecules_elements, non_hydrogen_molecules_positions, all_no_of_H_on_atoms_in_molecule, equivalent_molecule_atom_indices_comparison, max_distance_disparity, neighbouring_molecules_about_dimers, non_hydrogen_molecules, symmetric_dimer_pairs, no_of_cpus):
    """
    This generator is designed to return all the input methods required for the get_symmetric_dimer_pairs_single_process method. 

    Parameters
    ----------
    dimers : list
        This is a list of dimer information, given as (index/name of dimer, name of molecule 1 in dimer, name of molecule 2 in dimer, unit cell ijk displacement of molecule 2, unit cell displacement of molecule 2, displacement of dimer COM)

    non_hydrogen_molecules_elements : list
        These are the element of all the molecules that can make up the dimers.
    non_hydrogen_molecules_positions : list
        These are the positions of all the molecules that can make up the dimers.
    all_no_of_H_on_atoms_in_molecule : list
        These are the number of hydrogens bound to each "Heavy" atom in each molecule that can make up the dimers. 

    equivalent_molecule_atom_indices_comparison : dict.
        This is all the ways that two molecules in non_hydrogen_graphs can map onto each other. 
    max_distance_disparity : float
        This is the maximum difference in relative positions of atoms between dimer 1 and dimer 2 to be considered equivalent/identical
    neighbouring_molecules_about_dimers : dict.
        This is the information about the molecules that surround (in the vicinity of) each dimer in the crystal.
    non_hydrogen_molecules : list of ase.Atoms
        These are Atoms of objects of the molecules in the crystal. These do not contain hydrogen atoms. This is for debugging use only. 

    symmetric_dimer_pairs : list object (either a list or a multiprocessing.Manager.list)
        This list is designed to store all the symmetric dimers identified from this program.

    no_of_cpus : int.
        This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.
        
    Returns
    -------
    dimer1_name : int
        This is the name in the dimers list for the first dimer. 
    d1_m1_name : int
        This is the name of molecule 1 in the dimer. This is the molecules name in the molecules lists, (such as non_hydrogen_molecules_elements, non_hydrogen_molecules_positions, all_no_of_H_on_atoms_in_molecule).
    d1_m2_name : int
        This is the name of molecule 1 in the dimer. This is the molecules name in the molecules lists, (such as non_hydrogen_molecules_elements, non_hydrogen_molecules_positions, all_no_of_H_on_atoms_in_molecule).

    dimer2_name : int
        This is the name in the dimers list for the second dimer. 
    d2_m1_name : int
        This is the name of molecule 1 in the dimer. This is the molecules name in the molecules lists, (such as non_hydrogen_molecules_elements, non_hydrogen_molecules_positions, all_no_of_H_on_atoms_in_molecule).
    d2_m2_name : int
        This is the name of molecule 1 in the dimer. This is the molecules name in the molecules lists, (such as non_hydrogen_molecules_elements, non_hydrogen_molecules_positions, all_no_of_H_on_atoms_in_molecule).

    d1_elements : list of str.
        This is the list of elements of atoms in dimer 1.
    d1_positions : list of np.array
        This is the list of positions of atoms in dimer 1.
    d1_no_H_attached_to_nonH_atoms : list of int
        This is the list of hydrogens on each heavy atom in dimer 1.

    d2_m1_original_elements : list of str.
        This is the list of elements of atoms in the first molecule in dimer 2.
    d2_m1_original_positions : list of np.array
        This is the list of positions of atoms in the first molecule in dimer 2.
    d2_m1_no_H_attached_to_nonH_atoms : list of int
        This is the list of hydrogens attached to "heavy" atoms in the first molecule in dimer 2.

    d2_m2_original_elements : list of str.
        This is the list of elements of atoms in the second molecule in dimer 2.
    d2_m2_original_positions : list of np.array
        This is the list of positions of atoms in the second molecule in dimer 2.
    d2_m2_no_H_attached_to_nonH_atoms : list of int
        This is the list of hydrogens attached to "heavy" atoms in the second molecule in dimer 2.

    equivalent_molecule_atom_indices_comparison : dict.
        This is all the ways that two molecules in non_hydrogen_graphs can map onto each other. 
    max_distance_disparity : float
        This is the maximum difference in relative positions of atoms between dimer 1 and dimer 2 to be considered equivalent/identical
    neighbouring_molecules_about_dimers : dict.
        This is the information about the molecules that surround (in the vicinity of) each dimer in the crystal.
    non_hydrogen_molecules : dict. of ase.Atoms
        These are Atoms of objects of the molecules in the crystal. These do not contain hydrogen atoms. This is for debugging use only. 

    symmetric_dimer_pairs : list object (either a list or a multiprocessing.Manager.list)
        This list is designed to store all the symmetric dimers identified from this program.

    no_of_cpus : int.
        This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.
    """

    # First, obtain the names of the dimers you have collected.
    dimer_names = sorted(dimers.keys())

    # Second, for each dimer in the dimers list.
    for index1 in range(len(dimer_names)):

        # Third, obtain the name of the first dimer.
        dimer1_name = dimer_names[index1]

        # Fourth, obtain the names of the molecules in the dimer, as well as the displacement of molecule 2 in dimer 1, and the centre of mass molecule to move dimer 2 by. 
        d1_m1_name, d1_m2_name, _, dist1, move_com_by_1, d1_shortest_distance = dimers[dimer1_name]

        # Fifth, obtain the elements, positions, and number of hydrogens bound to "heavy" atoms in molecule 1 of dimer 1
        d1_m1_original_elements           = non_hydrogen_molecules_elements [d1_m1_name]
        d1_m1_original_positions          = non_hydrogen_molecules_positions[d1_m1_name] + move_com_by_1
        d1_m1_no_H_attached_to_nonH_atoms = all_no_of_H_on_atoms_in_molecule[d1_m1_name]

        # Sixth, obtain the elements, positions, and number of hydrogens bound to "heavy" atoms in molecule 2 of dimer 1
        d1_m2_original_elements           = non_hydrogen_molecules_elements [d1_m2_name]
        d1_m2_original_positions          = non_hydrogen_molecules_positions[d1_m2_name] + move_com_by_1 + dist1
        d1_m2_no_H_attached_to_nonH_atoms = all_no_of_H_on_atoms_in_molecule[d1_m2_name]

        # Seventh, obtain the elements, positions, and number of hydrogens bound to "heavy" atoms in dimer 1 (by concaternating molecule 1 and molecule 2 together in dimer 1)
        d1_elements                       = d1_m1_original_elements + d1_m2_original_elements
        d1_positions                      = np.concatenate([d1_m1_original_positions,d1_m2_original_positions])
        d1_no_H_attached_to_nonH_atoms    = d1_m1_no_H_attached_to_nonH_atoms + d1_m2_no_H_attached_to_nonH_atoms

        # Eighth, for every other dimer in the dimers list (that has a higher indice than dimer1_name)
        for index2 in range(index1+1,len(dimer_names)):

            # Ninth, obtain the name of the second dimer.
            dimer2_name = dimer_names[index2]

            # Tenth, obtain the name of the molecules in the dimer, as well as the displacement of molecule 2 in dimer 1, and the centre of mass molecule to move dimer 2 by. 
            d2_m1_name, d2_m2_name, _, dist2, move_com_by_2, d2_shortest_distance = dimers[dimer2_name]

            # Eleventh, obtain the elements, positions, and number of hydrogens bound to "heavy" atoms in molecule 1 of dimer 2
            d2_m1_original_elements           = non_hydrogen_molecules_elements [d2_m1_name]
            d2_m1_original_positions          = non_hydrogen_molecules_positions[d2_m1_name] + move_com_by_2
            d2_m1_no_H_attached_to_nonH_atoms = all_no_of_H_on_atoms_in_molecule[d2_m1_name]

            # Twelfth, obtain the elements, positions, and number of hydrogens bound to "heavy" atoms in molecule 1 of dimer 2
            d2_m2_original_elements           = non_hydrogen_molecules_elements [d2_m2_name]
            d2_m2_original_positions          = non_hydrogen_molecules_positions[d2_m2_name] + move_com_by_2 + dist2
            d2_m2_no_H_attached_to_nonH_atoms = all_no_of_H_on_atoms_in_molecule[d2_m2_name]

            # Finally, yield input data.
            yield (dimer1_name, d1_m1_name, d1_m2_name), (dimer2_name, d2_m1_name, d2_m2_name), (d1_elements, d1_positions, d1_no_H_attached_to_nonH_atoms), (d2_m1_original_elements, d2_m1_original_positions, d2_m1_no_H_attached_to_nonH_atoms), (d2_m2_original_elements, d2_m2_original_positions, d2_m2_no_H_attached_to_nonH_atoms), (d1_shortest_distance, d2_shortest_distance), (equivalent_molecule_atom_indices_comparison, max_distance_disparity, neighbouring_molecules_about_dimers, non_hydrogen_molecules), symmetric_dimer_pairs, no_of_cpus

# ----------------------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------

