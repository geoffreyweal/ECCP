"""
determine_invariance.py, Geoffrey Weal, 3/3/22

This script is designed to determine if the two systems (being either molecules or dimers) are invariant.
"""
import numpy as np
import sys, warnings

from copy import deepcopy
from itertools import permutations
from scipy.spatial.transform import Rotation

from SUMELF import get_distance, get_unit_vector, are_two_lists_within_eachother, get_centre_of_mass, get_reflection_matrix_from_plane

from ECCP.ECCP.invariance_methods.are_environments_equivalent import are_environments_equivalent

# --------------------------------------------------------------------------------------------------------------

# Below are the reflection matricies that will be used to determine if two systems are invariant or not.
non_reflection_matrix = np.eye(3)
reflection_matrix_x   = get_reflection_matrix_from_plane([1,0,0])
reflection_matrix_y   = get_reflection_matrix_from_plane([0,1,0])
reflection_matrix_z   = get_reflection_matrix_from_plane([0,0,1])
reflection_matrix_xy  = reflection_matrix_y @ reflection_matrix_x
reflection_matrix_xz  = reflection_matrix_z @ reflection_matrix_x
reflection_matrix_yz  = reflection_matrix_z @ reflection_matrix_y
reflection_matrix_xyz = reflection_matrix_z @ reflection_matrix_y @ reflection_matrix_x
reflection_matrices = tuple([non_reflection_matrix, reflection_matrix_x, reflection_matrix_y, reflection_matrix_z, reflection_matrix_xy, reflection_matrix_xz, reflection_matrix_yz, reflection_matrix_xyz])

# --------------------------------------------------------------------------------------------------------------

#xyz_distance_tolerances = np.array([0.00000001]*3)
def determine_invariance_MEA(s1_directions, s1_indices_in_directions, lengths_s1, s2_directions, s2_indices_in_directions, lengths_s2, s1_centre_index, s2_centre_index, s1_original_elements, s1_original_positions, no_of_H_on_atoms_in_system1, s2_original_elements, s2_original_positions, no_of_H_on_atoms_in_system2, max_distance_disparity, length_max_distance_disparity, indices_being_compared, neighbouring_molecules_about_systems, non_hydrogen_systems):
    """
    This method is designed to determine if the two chemical systems are invariant based on obtaining a rotation matrix using four atoms, where
        * Each system has a centre atom
        * Each system has three directions that point from the centre atom to another atoms in the system
        
    If the directions of each system can be mapped onto each other using the procrustes analysis, then the two systems are invariant.

    Each direction should not be parallel to each other (if possible). 

    Note: System (chemical system) here refers to either molecules or dimers. 

    Parameters
    ----------
    s1_directions : list of numpy.array
        These are the directions from the centre atom to three atoms in the system that are not in parallel with each other (if possible) in the first system. 
    s1_indices_in_directions
        These are the indices of the atoms pointed to by the directions in the first system.
    lengths_s1
        These are the vector length of the first system's directions

    s2_directions : list of numpy.array
        These are the directions from the centre atom to three atoms in the system that are not in parallel with each other (if possible) in the second system. 
    s2_indices_in_directions
        These are the indices of the atoms pointed to by the directions in the second system.
    lengths_s2
        These are the vector length of the second system's directions
        
    s1_centre_index
        This is the index on the centre atom to draw directions from in the first system. 
    s2_centre_index
        This is the index on the centre atom to draw directions from in the second system. 

    s1_original_elements : list of str.
        These are the elements of the atoms in system 1.
    s1_original_positions : numpy.array
        These are the positions of the atoms in in system 1.
    no_of_H_on_atoms_in_system1 : list
        A list of the number of hydrogens attached to each atoms in system 1.

    s2_original_elements : list of str.
        These are the elements of the atoms in system 2.
    s2_original_positions : numpy.array
        These are the positions of the atoms in in system 2.
    no_of_H_on_atoms_in_system2 : list
        A list of the number of hydrogens attached to each atoms in system 2.

    max_distance_disparity: float
        This is the maximum that any two "could be equivalent" atoms can be between system 1 and system 2 for systems 1 and 2 to be considered variant.
    length_max_distance_disparity : float
        This is the maximum that any two lengths in system 1 and 2 can be to be considered equivalent.

    names_being_compared: list
        These are the names of the systems being compared.
    neighbouring_molecules_about_systems : dict.
        This is the information about the systems that surround (in the vicinity of) each moleule in the crystal.
    non_hydrogen_systems : list of ase.Atoms
        This is the list of systems in the crystal, not including hydrogens.

    Attributes
    ----------
    reflection_matrices : list of numpy.array objects
        These are the matrices for reflecting the second system in various ways in the x, y and z axes to determine if system 2 is equivalent to system 1.

    Returns
    -------
    True if the two systems are variants of eachother, False if the two systems are invariant.
    """

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # First, centre the systems around their centre of masses

    # 1.1: For system 1.
    system1_elements        = deepcopy(s1_original_elements)
    system1_positions       = deepcopy(s1_original_positions)
    system1_com             = get_centre_of_mass(system1_elements, system1_positions)
    system1_positions      -= system1_com
    system1_no_of_hydrogens = deepcopy(no_of_H_on_atoms_in_system1)

    # 1.2: For system 2.
    system2_elements        = deepcopy(s2_original_elements)
    system2_positions       = deepcopy(s2_original_positions)
    system2_com             = get_centre_of_mass(system2_elements, system2_positions)
    system2_positions      -= system2_com
    system2_no_of_hydrogens = deepcopy(no_of_H_on_atoms_in_system2)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    # Second, make sure that the number of elements in system1 is the same as in system2.
    #         * They should be. If they were entered here and are different. This may indicate 
    #           an issue in an earlier part of the code that uses this determine_invariance_MEA
    if not (len(system1_elements) == len(system2_elements)):
        to_string  = 'Error: len(system1_elements) is not equal to len(system2_elements)\n'
        to_string += f'len(system1_elements) = {len(system1_elements)}\n'
        to_string += f'len(system2_elements) = {len(system2_elements)}\n'
        to_string += 'Check this. There may be a programming error.'
        raise Exception(to_string)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    # Third, set a variable to indicate if Rotation.align_vectors gave a warning
    align_vectors_gave_a_warning = False
    warning_messages = []

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Fourth, determine if the two systems are equivalent using the Minimal Elemental Abundance method. 

    # 4.1: For each reflection matrix
    for reflection_matrix in reflection_matrices:

        # 4.2: Obtain the directions and atom positions of system 2 after being reflected by a reflection matrix.
        reflected_s2_directions     = s2_directions @ reflection_matrix.T
        reflected_system2_positions = system2_positions @ reflection_matrix.T

        # 4.3: For each permutation of system 2 lengths
        for idx_s2 in permutations(list(range(3))):

            # 4.4: Get the lists of system 2 sorted by the vector lengths.
            lengths_s2_permutated = list([lengths_s2[ii] for ii in idx_s2])

            # 4.5: Determine if this arrangement of lengths in each system are the same. 
            #      This means we ensure that we have the right permutation of idx_s2 for the second system.
            if not are_two_lists_within_eachother(lengths_s1, lengths_s2_permutated, length_max_distance_disparity):
                continue

            # 4.6: Get the other lists of system 2 sorted by idx_m2.
            reflected_s2_directions_permutated = deepcopy([reflected_s2_directions[ii] for ii in idx_s2])

            # 4.7: Get the rotation matrix.
            with warnings.catch_warnings(record=True) as warning_watcher:
                estimated_rotation, rssd = Rotation.align_vectors(s1_directions, reflected_s2_directions_permutated)
            rotation_matrix = estimated_rotation.as_matrix()

            # 4.8: Determine if Rotation.align_vectors gave a warning. 
            #      If it does, send a boolean back to try other methods if this method couldn't find a way to 
            if len(warning_watcher) > 0:
                align_vectors_gave_a_warning = True
                for warning_message in warning_watcher:
                    warning_messages.append(str(warning_message.message))

            # 4.9: Obtain the position of system 2 when rotated with the rotation matrix.
            rotated_and_reflected_system2_positions = reflected_system2_positions @ rotation_matrix.T

            # 4.10: Get a list of indices in system 1 and system 2 that overlap each other.
            #       i.e. likely to be the same atom in variant systems.
            comparison_s1_s2_indices = get_comparison_s1_s2_indices(system1_elements, system1_positions, system1_no_of_hydrogens, system2_elements, rotated_and_reflected_system2_positions, system2_no_of_hydrogens, max_distance_disparity)

            # 4.11: Sort the indices that relate system 1 to system 2 from lowest distance to highest comparison_s1_s2_indices
            comparison_s1_s2_indices.sort(key=lambda x: x[2])

            # 4.12: Determine which indices in system 1 related to which indices in system 2.
            #       * Each index in system 1 will only be a key   once in the indices_s1_related_to_moved_s2 dictionary, and
            #       * Each index in system 2 will only be a value once in the indices_s1_related_to_moved_s2 dictionary.
            indices_s1_related_to_moved_s2 = {}
            for index_s1, index_s2, distance in sorted(comparison_s1_s2_indices):
                if (index_s1 in indices_s1_related_to_moved_s2.keys()) or (index_s2 in indices_s1_related_to_moved_s2.values()):
                    break
                indices_s1_related_to_moved_s2[index_s1] = index_s2

            # 4.13: If the length of the indices_s1_related_to_moved_s2 is 
            # the same as the number of atoms in system 1 (and in system 2), 
            # then the two systems are equivalent and variants, so return True
            if len(indices_s1_related_to_moved_s2) == len(system1_elements):
                rotation_reflection_matrix = rotation_matrix @ reflection_matrix
                if are_environments_equivalent(rotation_reflection_matrix, system1_com, system2_com, indices_being_compared, neighbouring_molecules_about_systems, non_hydrogen_systems):
                    warning_messages_from_this_comparison = [str(warning_message.message) for warning_message in warning_watcher]
                    return True, indices_s1_related_to_moved_s2, set() # warning_messages_from_this_comparison

            # 4.14: If you get to here, there is a programming error. 
            if len(indices_s1_related_to_moved_s2) > len(system1_elements):
                to_string  = f'Error: len(indices_s1_related_to_moved_s2) > len(system1_elements)\n'
                to_string += f'This indicates that one or more atoms have been double counted in indices_s1_related_to_moved_s2\n'
                to_string += f'len(indices_s1_related_to_moved_s2) = {len(indices_s1_related_to_moved_s2)}\n'
                to_string += f'len(system1_elements) = {len(system1_elements)}\n'
                to_string += f'indices_s1_related_to_moved_s2 = {indices_s1_related_to_moved_s2}\n'
                to_string += f'system1_elements = {system1_elements}\n'
                to_string += f'Check this.\n'
                raise Exception(to_string)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    # Fifth, if there were any unique warnings given from the Rotation.align_vectors method in scipy, print them
    unique_warning_messages = set(warning_messages)
    '''
    for warning_message in unique_warning_messages:
        #warnings.warn(warning_message.message)
        print(warning_message, file=sys.stderr)
    '''

    # Sixth, could not find a combination where the systems were variants, so here the two are invariant.
    return False, None, unique_warning_messages

# --------------------------------------------------------------------------------------------------------------

def get_comparison_s1_s2_indices(system1_elements, system1_positions, system1_no_of_hydrogens, system2_elements, rotated_and_reflected_system2_positions, system2_no_of_hydrogens, max_distance_disparity):
    """
    This method is designed to determine the indices of system 2 that map onto the atoms of system 1, given that the atoms of system 2 have been rotated and reflected about it's centre of mass. 

    Parameters
    ----------
    system1_elements : list of str.
        These are the elements of atoms in system 1. 
    system1_positions : numpy.array
        These are the positions of atoms in system 1. 
    system1_no_of_hydrogens : list of int.
        These are the number of hydrogens bound to the "heavy" atoms in system 1. 

    system2_elements : list of str.
        These are the elements of atoms in system 2. 
    rotated_and_reflected_system2_positions : numpy.array
        These are the positions of atoms in system 2. These are the positions of system 2 after system 2 has been rotated and reflected about it's centre of mass. 
    system2_no_of_hydrogens : list of int.
        These are the number of hydrogens bound to the "heavy" atoms in system 2. 

    max_distance_disparity: float
        This is the maximum that any two "could be equivalent" atoms can be between system 1 and system 2 for systems 1 and 2 to be considered variant
    """

    # First, initialise the list that will be used to identify atoms in system 1 that overlap with the atoms in system 2 after system 2 has been rotated and reflected about it's centre of mass. 
    comparison_s1_s2_indices = []

    # Second, for each atom in system 1
    for index_s1 in range(len(system1_elements)):

        # Third, obtain the element, position, and number of bound hydrogens of atom at index index_s1 in system 1.
        system1_element  = system1_elements[index_s1]
        system1_position = system1_positions[index_s1]
        system1_nH       = system1_no_of_hydrogens[index_s1]

        # Fourth, for each atom in system 2
        for index_s2 in range(len(system2_elements)):

            # Fifth, obtain the element, position, and number of bound hydrogens of atom at index index_s2 in system 2.
            #        * Here, system 2 has been rotated and reflected about it's centre of mass. 
            system2_element  = system2_elements[index_s2]
            system2_position = rotated_and_reflected_system2_positions[index_s2]
            system2_nH       = system2_no_of_hydrogens[index_s2]

            # Sixth, determine if the two atoms are teh same element. 
            if not (system1_element == system2_element):
                continue

            # Seventh, determine if both atoms have the same number of bound hydrogens or not.
            if not (system1_nH == system2_nH):
                continue

            # Eighth, get the distance between these two atoms.
            distance = get_distance(system1_position, system2_position)

            # Ninth, if the two atoms are within max_distance_disparity, we indicate that these two atoms overlap eachother, so are equivalent atoms. 
            if distance <= max_distance_disparity:
                comparison_s1_s2_indices.append((index_s1, index_s2, distance))

    # Tenth, return comparison_s1_s2_indices
    return comparison_s1_s2_indices

# --------------------------------------------------------------------------------------------------------------





