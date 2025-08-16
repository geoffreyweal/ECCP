"""
are_molecules_variant_from_comprehensive.py, Geoffrey Weal, 9/2/24

If one of the molecules in the comparison has less than 4 atoms, then we can not allign them using the methods in the minimal elemental abundance invariance method.

If this happens, we will default to the main component of the comprehensive invarience method for these two molecules.

This code is based various pieces of code from ECCP/ECCP/get_unique_molecules_methods/set_of_invariance_methods/comprehensive_invariance_method.py
"""
from copy import deepcopy
import multiprocessing as mp

from SUMELF import GraphMatcher

from ECCP.ECCP.invariance_methods.utilities                                                             import get_permutated_indices_list
from ECCP.ECCP.invariance_methods.common_utility_methods_for_all_invariance_methods.are_systems_variant import are_systems_variant

def are_molecules_variant_from_comprehensive(m1_original_elements, m1_original_positions, no_of_H_on_atoms_in_molecule1, m2_original_elements, m2_original_positions, no_of_H_on_atoms_in_molecule2, max_distance_disparity, molecule_names_being_compared, neighbouring_molecules_about_molecules, non_hydrogen_molecules, non_hydrogen_graph_m1, non_hydrogen_graph_m2, no_of_cpus=1):
    """
    If one of the molecules in the comparison has less than 4 atoms, then we can not allign them using the methods in the minimal elemental abundance invariance method.

    If this happens, we will default to the main component of the comprehensive invarience method for these two molecules.

    This code is based various pieces of code from ECCP/ECCP/get_unique_molecules_methods/set_of_invariance_methods/comprehensive_invariance_method.py

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

    max_distance_disparity : float
        This is the maximum disparity between two molecules to be considered invariant. If max_distance_disparity is given as None, the default value will be given. Default: 0.01 Å.

    molecule_names_being_compared : (int, int)
        These are the names of the molecules currently being compared. 

    neighbouring_molecules_about_molecules : dict.
        This is the information about the molecules that surround (in the vicinity of) each moleule in the crystal.
    non_hydrogen_molecules : list of ase.Atoms
        This is the list of molecules in the crystal, not including hydrogens.

    non_hydrogen_graph_m1 : networkx.graph
        This is the graph for molecule 1 that does not contain hydrogens
    non_hydrogen_graph_m2 : networkx.graph
        This is the graph for molecule 2 that does not contain hydrogens

    no_of_cpus : int.
        This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.

    Returns
    -------
    is_variant : bool.
        This indicates if the two molecules being compared are variants of each other or not. 
    mol1_to_mol2_conversion_Comp : None or dict.
        If the two molecules are variants of each other, this dictionary indicates how to the indices in molecule 1 (key) relate to those in molecule 2 (value).
    """

    # First, get all the indices that are equivalent to eachother in each molecule in each dimer. 
    #        * Note, this is GraphMatcher(non_hydrogen_graph_m2, non_hydrogen_graph_m1) rather than 
    #          GraphMatcher(non_hydrogen_graph_m1, non_hydrogen_graph_m2) because we want to convert 
    #          graph 2 into graph 1.
    #          otherwise would put this in the code:
    #
    #          ```
    #          GM = GraphMatcher(non_hydrogen_graph_m1, non_hydrogen_graph_m2)
    #          all_unique_matches = GM.get_all_unique_isomorphic_graphs()
    #          em_indices_m2_to_m1 = [ ({value: key for key, value in all_unique_match.items()}) for all_unique_match in all_unique_matches ]
    #          ```
    GM = GraphMatcher(non_hydrogen_graph_m2, non_hydrogen_graph_m1)
    em_indices_m2_to_m1 = GM.get_all_unique_isomorphic_graphs()

    # Second, initialise the ``is_variant`` boolean, which indicates if the two molecules being compared are variant or not. 
    is_variant = False

    # Third, initialise the ``mol1_to_mol2_conversion_Comp`` boolean for recording how the atom indices of mol1 relate to those of mol2. 
    mol1_to_mol2_conversion_Comp = None

    # Fourth, obtain the generator for the input data.
    input_generator = get_inputs(em_indices_m2_to_m1,  m1_original_elements, m1_original_positions, no_of_H_on_atoms_in_molecule1,  m2_original_elements, m2_original_positions, no_of_H_on_atoms_in_molecule2,  max_distance_disparity, molecule_names_being_compared, neighbouring_molecules_about_molecules, non_hydrogen_molecules)

    # Fifth, perform the are_systems_variant method on all inputs depending on if you are performing the task with one cpu or with multiple cores.
    if no_of_cpus == 1: # If the user only wants to use 1 cpu, perform tasks without using multiprocessing

        # 5.1: For each possibility of permutation in em_indices_m2_to_m1
        for input_data in input_generator:

            # 5.2: Determine if the current mol2->mol1 index comparison allows these two molecules to be seen as variant. 
            is_variant, mol1_to_mol2_conversion_Comp = compare_two_molecules_in_two_index_configurations_single_process(input_data)

            # 5.3: If compare_two_molecules_in_two_index_configurations_single_process throws an exception, break out of the loop
            if isinstance(is_variant, Exception):
                break

            # 5.4: Make sure that at this point, is_variant is a boolean.
            if not isinstance(is_variant, bool):
                break

            # 5.5: If ``is_variant`` is True, break out of the for loop.
            if is_variant:
                break

    else:

        # 5.6: Set up the multiprocessing pool.
        with mp.Pool(processes=no_of_cpus) as pool:

            # 5.7: Compare the two molecules via all the possible comparisons between 
            for is_variant, mol1_to_mol2_conversion_Comp in pool.imap(compare_two_molecules_in_two_index_configurations_single_process, input_generator):

                # 5.8: If compare_two_molecules_in_two_index_configurations_single_process throws an exception, break out of the loop
                if isinstance(is_variant, Exception):
                    break

                # 5.9: Make sure that at this point, is_variant is a boolean.
                if not isinstance(is_variant, bool):
                    break

                # 5.10: If ``is_variant`` is True, break out of the for loop.
                if is_variant:
                    break

    # Sixth, return if the two molecules are variants of each other, 
    #        * If the two molecules are variants of each other, mol1_to_mol2_conversion_Comp will also be given.
    #        * mol1_to_mol2_conversion_Comp is a dictionary that gives information about which atoms (indices) 
    #          in molecule 1 are equivalent to which atoms (indices) in molecule 2. 
    return is_variant, mol1_to_mol2_conversion_Comp

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

def get_inputs(em_indices_m2_to_m1,  m1_original_elements, m1_original_positions, no_of_H_on_atoms_in_molecule1,  m2_original_elements, m2_original_positions, no_of_H_on_atoms_in_molecule2,  max_distance_disparity, molecule_names_being_compared, neighbouring_molecules_about_molecules, non_hydrogen_molecules):
    """
    This generator will provide the input data needed for the ```compare_two_molecules_in_two_index_configurations_single_process``` method. 

    Parameters
    ----------
    em_indices_m2_to_m1 : list of dict.
        These are all the ways that the atoms (indices) in molecule 2 map onto the atoms (indices) of molecule 1.

    m1_original_elements : list of str.
        This is the list of elements in molecule 1. 
    m1_original_positions : list of (1x3) numpy.array
        This list contains all the positions of the atoms in molecule 1.
    no_of_H_on_atoms_in_molecule1 : list of int
        This is the number of hydrogens bound to the atoms in molecule 1. 

    m2_original_elements : list of str.
        This is the list of elements in molecule 2.
    m2_original_positions : list of (1x3) numpy.array
        This list contains all the positions of the atoms in molecule 2.
    no_of_H_on_atoms_in_molecule2 : list of int
        This is the number of hydrogens bound to the atoms in molecule 2. 

    max_distance_disparity : float
        This is the maximum disparity between two molecules to be considered invariant. If max_distance_disparity is given as None, the default value will be given. Default: 0.01 Å.

    molecule_names_being_compared : (int, int)
        These are the names of the molecules currently being compared. 

    neighbouring_molecules_about_molecules : dict.
        This is the information about the molecules that surround (in the vicinity of) each moleule in the crystal.
    non_hydrogen_molecules : list of ase.Atoms
        This is the list of molecules in the crystal, not including hydrogens.

    Returns
    -------
    comparison : dict.
        This is one way that the atoms (indices) in molecule 2 map onto the atoms (indices) of molecule 1.

    m1_original_elements : list of str.
        This is the list of elements in molecule 1. 
    m1_original_positions : list of (1x3) numpy.array
        This list contains all the positions of the atoms in molecule 1.
    no_of_H_on_atoms_in_molecule1 : list of int
        This is the number of hydrogens bound to the atoms in molecule 1. 

    m2_original_elements : list of str.
        This is the list of elements in molecule 2.
    m2_original_positions : list of (1x3) numpy.array
        This list contains all the positions of the atoms in molecule 2.
    no_of_H_on_atoms_in_molecule2 : list of int
        This is the number of hydrogens bound to the atoms in molecule 2. 

    max_distance_disparity : float
        This is the maximum disparity between two molecules to be considered invariant. If max_distance_disparity is given as None, the default value will be given. Default: 0.01 Å.

    molecule_names_being_compared : (int, int)
        These are the names of the molecules currently being compared. 

    neighbouring_molecules_about_molecules : dict.
        This is the information about the molecules that surround (in the vicinity of) each moleule in the crystal.
    non_hydrogen_molecules : list of ase.Atoms
        This is the list of molecules in the crystal, not including hydrogens.
    """

    # First, for each possibility of permutation in em_indices_m2_to_m1
    for comparison in em_indices_m2_to_m1: 

        # Second, yield the input data required for running the ```compare_two_molecules_in_two_index_configurations_single_process``` method. 
        yield (comparison,  m1_original_elements, m1_original_positions, no_of_H_on_atoms_in_molecule1,  m2_original_elements, m2_original_positions, no_of_H_on_atoms_in_molecule2,  max_distance_disparity, molecule_names_being_compared, neighbouring_molecules_about_molecules, non_hydrogen_molecules)

def compare_two_molecules_in_two_index_configurations_single_process(input_data):
    """
    This method will run the ``are_systems_variant`` method for the index order of the atoms for molecule 2. 

    Parameters
    ----------
    comparison : dict.
        This is one way that the atoms (indices) in molecule 2 map onto the atoms (indices) of molecule 1.

    m1_original_elements : list of str.
        This is the list of elements in molecule 1. 
    m1_original_positions : list of (1x3) numpy.array
        This list contains all the positions of the atoms in molecule 1.
    no_of_H_on_atoms_in_molecule1 : list of int
        This is the number of hydrogens bound to the atoms in molecule 1. 

    m2_original_elements : list of str.
        This is the list of elements in molecule 2.
    m2_original_positions : list of (1x3) numpy.array
        This list contains all the positions of the atoms in molecule 2.
    no_of_H_on_atoms_in_molecule2 : list of int
        This is the number of hydrogens bound to the atoms in molecule 2. 

    max_distance_disparity : float
        This is the maximum disparity between two molecules to be considered invariant. If max_distance_disparity is given as None, the default value will be given. Default: 0.01 Å.

    molecule_names_being_compared : (int, int)
        These are the names of the molecules currently being compared. 

    neighbouring_molecules_about_molecules : dict.
        This is the information about the molecules that surround (in the vicinity of) each moleule in the crystal.
    non_hydrogen_molecules : list of ase.Atoms
        This is the list of molecules in the crystal, not including hydrogens.

    Returns
    -------
    is_variant : bool.
        This indicates if the two molecules being compared are variants of each other or not. 
    mol1_to_mol2_conversion_Comp : None or dict.
        If the two molecules are variants of each other, this dictionary indicates how to the indices in molecule 1 (key) relate to those in molecule 2 (value).
    """

    # First, get the input data that is needed to run this method.
    comparison,  m1_original_elements, m1_original_positions, no_of_H_on_atoms_in_molecule1,  m2_original_elements, m2_original_positions, no_of_H_on_atoms_in_molecule2,  max_distance_disparity, molecule_names_being_compared, neighbouring_molecules_about_molecules, non_hydrogen_molecules = input_data

    # Second, reorder the positions and elements of molecule 2
    idx_m2                                  = get_permutated_indices_list(comparison)
    molecule2_elements_reordered            = [m2_original_elements[index] for index in idx_m2]
    molecule2_distances_reordered           = deepcopy(m2_original_positions)[idx_m2, :]
    no_of_H_on_atoms_in_molecule2_reordered = [no_of_H_on_atoms_in_molecule2[index] for index in idx_m2]

    # Third, determine if this configuration of molecule 2 means that molecules 1 and 2 are variant.
    if are_systems_variant(m1_original_elements, m1_original_positions, no_of_H_on_atoms_in_molecule1, molecule2_elements_reordered, molecule2_distances_reordered, no_of_H_on_atoms_in_molecule2_reordered, max_distance_disparity, molecule_names_being_compared, neighbouring_molecules_about_molecules, non_hydrogen_molecules):
        
        # 3.1: As these molecules are variants of each other with the indices reordered for idx_m2,
        #      we can obtain the way to convert molecule 1 into molecule 2. 
        mol1_to_mol2_conversion_Comp = {index_mol1: index_mol2 for index_mol1, index_mol2 in enumerate(idx_m2)}

        # 3.2: Return True (we can confirm these molecules are equivalent to each other), and give the 
        #      dictionary for how to convert molecule 1 to molecule 2. 
        return True, mol1_to_mol2_conversion_Comp

    else:

        # 3.3: The two systems are not variants of each other, so return False
        return False, None

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

