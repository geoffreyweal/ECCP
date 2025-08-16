"""
minimal_elemental_abundance_invariance_method.py, Geoffrey Weal, 5/4/22

This method will check that the elements in each molecule are the same and that the two molecules are rotationally variant.
"""
import sys
from tqdm import tqdm

from SUMELF import are_two_lists_within_eachother

from ECCP.ECCP.invariance_methods.common_minimal_elemental_abundance_invariance_utility_methods.methods_for_MEA_invariance_method.get_positions_of_low_abundant_elements_to_scan import get_positions_of_low_abundant_elements_to_scan
from ECCP.ECCP.invariance_methods.common_minimal_elemental_abundance_invariance_utility_methods.methods_for_MEA_invariance_method.get_points_of_first_molecule                   import get_points_of_first_molecule
from ECCP.ECCP.invariance_methods.common_minimal_elemental_abundance_invariance_utility_methods.methods_for_MEA_invariance_method.obtain_possible_3D_points_in_molecule2         import obtain_possible_3D_points_in_molecule2
from ECCP.ECCP.invariance_methods.common_minimal_elemental_abundance_invariance_utility_methods.methods_for_MEA_invariance_method.get_points_of_second_molecule                  import get_points_of_second_molecule
from ECCP.ECCP.invariance_methods.common_minimal_elemental_abundance_invariance_utility_methods.determine_invariance_MEA                                                         import determine_invariance_MEA

import multiprocessing as mp

max_number_of_atoms = 'all' # This is arbitary number, not sure what best to set this. This was changed from 10 to 'all' on 26/6/23. GRW
length_max_distance_disparity = ( 3.0 * (0.00000001 ** 2.0) ) ** 0.5
dotproduct_max_distance_disparity = 0.01

def are_molecules_variant_MEA(m1_original_elements, m1_original_positions, no_of_H_on_atoms_in_molecule1, m2_original_elements, m2_original_positions, no_of_H_on_atoms_in_molecule2, max_distance_disparity, molecule_names_being_compared, neighbouring_molecules_about_molecules, non_hydrogen_molecules, no_of_cpus=1): 
    """
    This method will check that the elements in each molecule are the same and that the two molecules are rotationally variant.
    
    Parameters
    ----------
    m1_original_elements : list 
        This contains the elements in the first molecules.
    m1_original_positions : numpy.array 
        This contains the positions in the first molecules.
    no_of_H_on_atoms_in_molecule1 : list
        A list of the number of hydrogens attached to each atoms in molecule 1.

    m2_original_elements : list 
        This contains the elements in the second molecules.
    m2_original_positions : numpy.array 
        This contains the positions in the second molecules.
    no_of_H_on_atoms_in_molecule2 : list
        A list of the number of hydrogens attached to each atoms in molecule 2.

    max_distance_disparity: float
        This is the maximum that any two "could be equivalent" atoms can be between molecule 1 and molecule 2 for molecules 1 and 2 to be considered variant.

    molecule_names_being_compared: list
        These are the indices/names of the molecules being compared.
    neighbouring_molecules_about_molecules : dict.
        This is the information about the molecules that surround (in the vicinity of) each moleule in the crystal.
    non_hydrogen_molecules : list of ase.Atoms
        This is the list of molecules in the crystal, not including hydrogens.

    no_of_cpus : int.
        This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.

    Attributes
    ----------
    length_max_distance_disparity : float
        This is the max difference between two molecules lengths to possibly be the same. if two molecules have different lengths, the are definitely different.
    dotproduct_max_distance_disparity : float
        This is the max difference between two molecules angles to possibly be the same. if two molecules have different angles, the are definitely different.

    Returns
    -------
    positions_are_variant : bool
        This indicates if the two molecules are variants of each other.
    mol1_to_mol2_conversion : dict
        This dictionary indicates which atom indices in molecule 1 go with which indices in molecule 2 to give equivalent molecules. 
    align_vectors_gave_a_warning_all_overall : set
        This set contains all the unique warning that were obtained during the Minimal Elemental Abundance (MEA) method (usually from numpy).
    """

    # First, if the the two molecules do not contain the same types and amounts ofelements, then these are different molecules!
    if not sorted(m1_original_elements) == sorted(m2_original_elements):
        return False, None, set()

    # Second, get the position of the elements in the lowest abundances
    m1_positions_of_lowest_elements = get_positions_of_low_abundant_elements_to_scan(m1_original_elements, m1_original_positions, max_number_of_atoms=max_number_of_atoms)
    m2_positions_of_lowest_elements = get_positions_of_low_abundant_elements_to_scan(m2_original_elements, m2_original_positions, max_number_of_atoms=max_number_of_atoms)

    # -----------------------------------------------------------------------------------------------------------------------------------------------------
    # Third, obtain the directions, vector lengths, and dotproducts to focus on in molecule 1.
    
    # 3.1: Obtain the directions, vector lengths, and dotproducts to focus on in molecule 1.
    m1_direction1, m1_direction2, m1_direction3, m1_centre_position, lengths_m1, dotproducts_m1_sorted, indices_of_points_m1, elements_in_points_m1 = get_points_of_first_molecule(m1_positions_of_lowest_elements)
    
    # 3.2: Get the indices of atom point in molecule 1.
    m1_centre_index, m1_point1_index, m1_point2_index, m1_point3_index = indices_of_points_m1
    
    # 3.3: Get the elements of atom point in molecule 1.
    m1_centre_element, m1_point1_element, m1_point2_element, m1_point3_element = elements_in_points_m1

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    # 3.4: Get the directions and the indices of directions for molecule 1
    m1_directions = [m1_direction1, m1_direction2, m1_direction3]
    m1_indices_in_directions = [m1_point1_index, m1_point2_index, m1_point3_index]

    # -----------------------------------------------------------------------------------------------------------------------------------------------------
    # Fourth, initialise the return values.

    # 4.1: Initialise a boolean to indicate if two molecules have been found to be variants of each other.
    is_variant = False

    # 4.2: This list indicates how to convert the indices of molecule 1 to molecule 2 that shows these two molecules as variants of each other. 
    mol1_to_mol2_conversion = None

    # 4.3: Initialise align_vectors_gave_a_warning as an empty set.
    #      * align_vectors_gave_a_warning: Will contain a warning message object if determine_invariance_MEA gives a warning. 
    align_vectors_gave_a_warning_all_overall = set()

    # -----------------------------------------------------------------------------------------------------------------------------------------------------
    
    # Fourth, create the generator that provides the inputs for the compare_two_molecules_in_two_index_configurations_single_process method.
    input_generator = get_inputs(m2_positions_of_lowest_elements, m1_centre_element, m1_point1_element, m1_point2_element, m1_point3_element, m1_directions, m1_centre_index, m1_indices_in_directions, lengths_m1, dotproducts_m1_sorted, dotproduct_max_distance_disparity, m1_original_elements, m1_original_positions, no_of_H_on_atoms_in_molecule1, m2_original_elements, m2_original_positions, no_of_H_on_atoms_in_molecule2, max_distance_disparity, length_max_distance_disparity, molecule_names_being_compared, neighbouring_molecules_about_molecules, non_hydrogen_molecules)

    # Fifth, perform the compare_two_molecules_in_two_index_configurations_single_process method on all inputs depending on if you are performing the task with one cpu or with multiple cores.
    if no_of_cpus == 1: # If the user only wants to use 1 cpu, perform tasks without using multiprocessing

        # 5.4: Analyse the spatial positions of molecule 2, for all the different index arrangements of molecule 2.
        for input_data in input_generator:

            # 5.5: Determine if the two molecules are varient in the given index arrangement of molecule 2.
            is_variant, mol1_to_mol2_conversion, align_vectors_gave_a_warning_all = compare_two_molecules_in_two_index_configurations_single_process(input_data)

            # 5.6: Update align_vectors_gave_a_warning_all_overall with align_vectors_gave_a_warning_all
            align_vectors_gave_a_warning_all_overall.update(align_vectors_gave_a_warning_all)

            # 5.7: If the molecule is variant, break out of the loop. 
            if is_variant:
                break
            else:
                align_vectors_gave_a_warning_all_overall.update(align_vectors_gave_a_warning_all)

    else:

        # 5.8: Set up the for loop for multiprocessing.
        with mp.Pool(processes=no_of_cpus) as pool:

            # 5.9: For the output of compare_two_molecules_in_two_index_configurations_single_process for each of the numerous index arrangements of molecule 2.
            for is_variant, mol1_to_mol2_conversion, align_vectors_gave_a_warning_all in pool.imap(compare_two_molecules_in_two_index_configurations_single_process, input_generator):

                # 5.10: Update align_vectors_gave_a_warning_all_overall with align_vectors_gave_a_warning_all
                align_vectors_gave_a_warning_all_overall.update(align_vectors_gave_a_warning_all)

                # 5.11: If the molecule is variant, break out of the loop. 
                if is_variant:
                    break
                else:
                    align_vectors_gave_a_warning_all_overall.update(align_vectors_gave_a_warning_all)

    # Sixth, if is_variant is True, we don't need to return anything in align_vectors_gave_a_warning_all_overall
    if is_variant:
        align_vectors_gave_a_warning_all_overall = set()

    # -----------------------------------------------------------------------------------------------------------------------------------------------------

    # Seventh, return if the two molecules are variants of each other, as well as index conversions between the molecules (if variant), and any warning messages from numpy. 
    return is_variant, mol1_to_mol2_conversion, align_vectors_gave_a_warning_all_overall

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

def get_inputs(m2_positions_of_lowest_elements, m1_centre_element, m1_point1_element, m1_point2_element, m1_point3_element, m1_directions, m1_centre_index, m1_indices_in_directions, lengths_m1, dotproducts_m1_sorted, dotproduct_max_distance_disparity, m1_original_elements, m1_original_positions, no_of_H_on_atoms_in_molecule1, m2_original_elements, m2_original_positions, no_of_H_on_atoms_in_molecule2, max_distance_disparity, length_max_distance_disparity, molecule_names_being_compared, neighbouring_molecules_about_molecules, non_hydrogen_molecules):
    """
    This generator will provide all the inputs needed for the compare_two_molecules_in_two_index_configurations_single_process method.

    Parameters
    ----------
    m2_positions_of_lowest_elements : list of (str, (numpy.array, int))
        These are the details of the atoms in the molecule given in order of least abundant elements to most abundant elements. The defails for each atom are given as (element, (position, index of atom in molecule))

    m1_centre_element : str.
        This is the element of the centre atom in molecule 1, used in the Minimal Elemental Abundance (MEA) method.
    m1_point1_element : str.
        This is the element of the point 1 atom in molecule 1, used in the Minimal Elemental Abundance (MEA) method.
    m1_point2_element : str.
        This is the element of the point 2 atom in molecule 1, used in the Minimal Elemental Abundance (MEA) method.
    m1_point3_element : str.
        This is the element of the point 3 atom in molecule 1, used in the Minimal Elemental Abundance (MEA) method.

    m1_directions : list of (1x3) numpy.array
        These are the direction vectors from the centre atom to point 1, 2, and 3 atoms in molecule 1, used in the Minimal Elemental Abundance (MEA) method.

    m1_centre_index : int 
        This is the index of the central atom in molecule 1, used in the Minimal Elemental Abundance (MEA) method.
    m1_indices_in_directions
        These are the indices of the point 1, 2 and 3 atoms in molecule 1, used in the Minimal Elemental Abundance (MEA) method.
    lengths_m1 : list of floats
        These are the lengths of the direction vectors in m1_directions.

    dotproducts_m1_sorted : list of floats
        These are the dot products between the direction vectors (angles between the direction vectors) in molecule 1.
    dotproduct_max_distance_disparity : float
        This is the maximum dot product value between two direction vectors to be considered the same (i.e. the maximum angle between two direction vectors to be considered the same).

    m1_original_elements : list of str.
        This is a list of the elements in the first molecule. 
    m1_original_positions : list of (1x3) numpy.array
        This contains all the position of the atoms in the first molecule. 
    no_of_H_on_atoms_in_molecule1 : list of ints.
        This list contains the number of hydrogens bound to each atom in the first molecule. 

    m2_original_elements : list of str.
        This is a list of the elements in the second molecule. 
    m2_original_positions : list of (1x3) numpy.array
        This contains all the position of the atoms in the second molecule. 
    no_of_H_on_atoms_in_molecule2 : list of ints.
        This list contains the number of hydrogens bound to each atom in the second molecule. 

    max_distance_disparity: float
        This is the maximum that any two "could be equivalent" atoms can be between molecule 1 and molecule 2 for molecule 1 and 2 to be considered variant.
    length_max_distance_disparity : float
        This is the maximum difference in the lengths of the direction vectors to be considered the same dimer.

    molecule_names_being_compared : tuple of ints
        These are the name of the molecule that are being compared. 
    neighbouring_molecules_about_molecules : dict 
        This dictionary contains information about the molecules that surrounding the molecules (used for environmental purposes).
    non_hydrogen_molecules : dict of ase.Atoms
        This dictionary contains all the molecules in the crystal without hydrogens.

    Returns
    -------
    molecule2_position_inputs : list of (1x3) numpy.arrays
        These are the positions of the centre and point 1, 2, and 3 atoms selected for processing with the Minimal Elemental Abundance (MEA) method for the second molecule. 
    dotproducts_m1_sorted : list of float.
        This list contains the dot products (angles) between the direction vectors for molecule 1. This has been sorted from lowest to highest value (angle)
    dotproduct_max_distance_disparity : float
        This is the maximum dot product value between two direction vectors to be considered the same (i.e. the maximum angle between two direction vectors to be considered the same).
    m2_indices_in_directions: 
        These are the indices of the point 1, 2, and 3 atoms in molecule 2
    determine_invariance_MEA_inputs : 
        These are most of the inputs values needed for the ``determine_invariance_MEA`` method that is used in the ``compare_two_molecules_in_two_index_configurations_single_process`` method.
    """

    # First, for each possible 3D points we could have for molecule 2:
    for m2_centre_element, m2_centre_position, m2_centre_index,  m2_point1_element, m2_point1_position, m2_point1_index,  m2_point2_element, m2_point2_position, m2_point2_index,  m2_point3_element, m2_point3_position, m2_point3_index in obtain_possible_3D_points_in_molecule2(m2_positions_of_lowest_elements, m1_centre_element, m1_point1_element, m1_point2_element, m1_point3_element):
        
        # Second, package the position inputs for molecule 2. These are used in inputs by the compare_two_molecules_in_two_index_configurations_single_process method.
        molecule2_MEA_position_inputs = m2_centre_position, m2_point1_position, m2_point2_position, m2_point3_position

        # Third, package the molecule 2 point indices. 
        m2_indices_in_directions = (m2_point1_index, m2_point2_index, m2_point3_index)

        # Fourth, package the position inputs used by the determine_invariance_MEA_inputs method (in the compare_two_molecules_in_two_index_configurations_single_process method).
        determine_invariance_MEA_inputs = m1_directions, m1_indices_in_directions, lengths_m1, m1_centre_index, m2_centre_index, m1_original_elements, m1_original_positions, no_of_H_on_atoms_in_molecule1, m2_original_elements, m2_original_positions, no_of_H_on_atoms_in_molecule2, max_distance_disparity, length_max_distance_disparity, molecule_names_being_compared, neighbouring_molecules_about_molecules, non_hydrogen_molecules

        # Fifth, yield the input data of interest to the compare_two_molecules_in_two_index_configurations_single_process method
        yield molecule2_MEA_position_inputs, dotproducts_m1_sorted, dotproduct_max_distance_disparity, m2_indices_in_directions, determine_invariance_MEA_inputs

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

def compare_two_molecules_in_two_index_configurations_single_process(input_data):
    """
    This method will determine if two molecules are variant based on the current index configuration of atoms in the second molecule.

    Parameters
    ----------
    molecule2_position_inputs : list of (1x3) numpy.arrays
        These are the positions of the centre and point 1, 2, and 3 atoms selected for processing with the Minimal Elemental Abundance (MEA) method for the second molecule. 
    dotproducts_m1_sorted : list of float.
        This list contains the dot products (angles) between the direction vectors for molecule 1. This has been sorted from lowest to highest value (angle)
    dotproduct_max_distance_disparity : float
        This is the maximum dot product value between two direction vectors to be considered the same (i.e. the maximum angle between two direction vectors to be considered the same).
    m2_indices_in_directions: 
        These are the indices of the point 1, 2, and 3 atoms in molecule 2
    determine_invariance_MEA_inputs : 
        These are most of the inputs values needed for the ``determine_invariance_MEA`` method that is used in the ``compare_two_molecules_in_two_index_configurations_single_process`` method.

    Attributes
    ----------
    m2_centre_position : (1x3) numpy.array
        This is the position of the central atom for molecule 2 used in the Minimal Elemental Abundance (MEA) method.
    m2_point1_position : (1x3) numpy.array
        This is the position of the point 1 atom for molecule 2 used in the Minimal Elemental Abundance (MEA) method.
    m2_point2_position : (1x3) numpy.array
        This is the position of the point 2 atom for molecule 2 used in the Minimal Elemental Abundance (MEA) method.
    m2_point3_position : (1x3) numpy.array
        This is the position of the point 3 atom for molecule 2 used in the Minimal Elemental Abundance (MEA) method.

    m2_point1_index : (1x3) numpy.array
        This is the index of the point 1 atom for molecule 2 used in the Minimal Elemental Abundance (MEA) method.
    m2_point2_index : (1x3) numpy.array
        This is the index of the point 2 atom for molecule 2 used in the Minimal Elemental Abundance (MEA) method.
    m2_point3_index : (1x3) numpy.array
        This is the index of the point 3 atom for molecule 2 used in the Minimal Elemental Abundance (MEA) method.

    m1_directions : list of (1x3) numpy.arrays
        These are the directions from the central atom to points 1, 2, and 3 atoms in molecule 1. 
    m1_indices_in_directions : list of (1x3) numpy.arrays
        These are the indices of the points 1, 2, and 3 atoms in molecule 1. 

    lengths_m1 : list of floats
        These are the lengths of the direction vectors in m1_directions.
    m1_centre_index : int 
        This is the index of the central atom in molecule 1, used in the Minimal Elemental Abundance (MEA) method.
    m1_centre_index : int 
        This is the index of the central atom in molecule 2, used in the Minimal Elemental Abundance (MEA) method.

    m1_original_elements : list of str.
        This is a list of the elements in the first molecule. 
    m1_original_positions : list of (1x3) numpy.array
        This contains all the position of the atoms in the first molecule. 
    no_of_H_on_atoms_in_molecule1 : list of ints.
        This list contains the number of hydrogens bound to each atom in the first molecule. 

    m2_original_elements : list of str.
        This is a list of the elements in the second molecule. 
    m2_original_positions : list of (1x3) numpy.array
        This contains all the position of the atoms in the second molecule. 
    no_of_H_on_atoms_in_molecule2 : list of ints.
        This list contains the number of hydrogens bound to each atom in the second molecule. 

    max_distance_disparity: float
        This is the maximum that any two "could be equivalent" atoms can be between molecule 1 and molecule 2 for molecule 1 and 2 to be considered variant.
    length_max_distance_disparity : float
        This is the maximum difference in the lengths of the direction vectors to be considered the same dimer.

    molecule_names_being_compared : tuple of ints
        These are the name of the molecule that are being compared. 
    neighbouring_molecules_about_molecules : dict 
        This dictionary contains information about the molecules that surrounding the molecules (used for environmental purposes).
    non_hydrogen_molecules : dict of ase.Atoms
        This dictionary contains all the molecules in the crystal without hydrogens.

    Returns
    -------
    positions_are_variant : bool
        This indicates if the two molecules are variants of each other.
    mol1_to_mol2_conversion : dict
        This dictionary indicates which atom indices in molecule 1 go with which indices in molecule 2 to give equivalent molecules. 
    align_vectors_gave_a_warning_single : set
        This set contains all the unique warning that were obtained during the Minimal Elemental Abundance (MEA) method (usually from numpy).
    """

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # First: Obtain the inputs from the input data.

    # 1.1: Obtain the inputs from the input data.
    molecule2_position_inputs, dotproducts_m1_sorted, dotproduct_max_distance_disparity, m2_indices_in_directions, determine_invariance_MEA_inputs = input_data

    # 1.2: Extract position data about molecule 2 from molecule2_position_inputs.
    m2_centre_position, m2_point1_position, m2_point2_position, m2_point3_position = molecule2_position_inputs

    # 1.3: Obtain the indices of the points of interest in molecule 2.
    m2_point1_index, m2_point2_index, m2_point3_index = m2_indices_in_directions

    # 1.4: Obtain the other inputs required for the determine_invariance_MEA method.
    m1_directions, m1_indices_in_directions, lengths_m1, m1_centre_index, m2_centre_index, m1_original_elements, m1_original_positions, no_of_H_on_atoms_in_molecule1, m2_original_elements, m2_original_positions, no_of_H_on_atoms_in_molecule2, max_distance_disparity, length_max_distance_disparity, molecule_names_being_compared, neighbouring_molecules_about_molecules, non_hydrogen_molecules = determine_invariance_MEA_inputs

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # Second, obtain a possible set of directions, vector lengths, and dotproducts to focus on in molecule 2.
    m2_direction1, m2_direction2, m2_direction3, lengths_m2, dotproducts_m2_sorted = get_points_of_second_molecule(m2_centre_position, m2_point1_position, m2_point2_position, m2_point3_position)

    # Third, if the vector lengths of directions and dot-products (angles) between directions are not the same, then move on to another of the possible 3D points in molecule 2.
    if not (are_two_lists_within_eachother(sorted(lengths_m1), sorted(lengths_m2), length_max_distance_disparity) and are_two_lists_within_eachother(dotproducts_m1_sorted, dotproducts_m2_sorted, dotproduct_max_distance_disparity)):
        return False, None, set()

    # Fourth, get the directions and the indices of directions for molecule 2
    m2_directions = [m2_direction1, m2_direction2, m2_direction3]
    m2_indices_in_directions = [m2_point1_index, m2_point2_index, m2_point3_index]

    # Fifth, get the matrix to transform molecule 2 upon molecule 1.
    positions_are_variant, mol1_to_mol2_conversion, align_vectors_gave_a_warning_single = determine_invariance_MEA(m1_directions, m1_indices_in_directions, lengths_m1, m2_directions, m2_indices_in_directions, lengths_m2, m1_centre_index, m2_centre_index, m1_original_elements, m1_original_positions, no_of_H_on_atoms_in_molecule1, m2_original_elements, m2_original_positions, no_of_H_on_atoms_in_molecule2, max_distance_disparity, length_max_distance_disparity, molecule_names_being_compared, neighbouring_molecules_about_molecules, non_hydrogen_molecules)

    # Sixth, return the result from the determine_invariance_MEA method
    return positions_are_variant, mol1_to_mol2_conversion, align_vectors_gave_a_warning_single

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

