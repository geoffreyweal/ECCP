
"""
get_symmetric_molecule_pairs_comprehensive.py, Geoffrey Weal, 11/2/24

This script is designed to determine all the spatially symmetric molecules in the unique_molecules_names list. 
"""
from copy import deepcopy

import multiprocessing as mp

from ECCP.ECCP.invariance_methods.utilities                                                             import get_permutated_indices_list
from ECCP.ECCP.invariance_methods.common_utility_methods_for_all_invariance_methods.are_systems_variant import are_systems_variant

def are_molecules_variant_comprehensive(molecule_1_information, molecule_2_information, em_indices_m2_to_m1, max_distance_disparity, neighbouring_molecules_about_molecules, non_hydrogen_molecules, symmetric_molecule_pairs, no_of_cpus=1):
    """
    If one of the molecules in the comparison has less than 4 atoms, then we can not allign them using the methods in the minimal elemental abundance invariance method.

    If this happens, we will default to the main component of the comprehensive invarience method for these two molecules.

    This code is based various pieces of code from ECCP/ECCP/get_unique_molecules_methods/set_of_invariance_methods/comprehensive_invariance_method.py

    Parameters
    ----------
    molecule_1_information : list
        This list contains the following information about molecules 1.
            mol_name1 : int
                This is the name of the first molecule we want to compare.
            molecule1_elements : list of str.
                This is the list of elements that make up molecule 1.
            molecule1_positions : 2D np.array
                These are the positions of the atoms in molecule 1
            no_of_H_on_atoms_in_molecule1 : list of int.
                These are the number of hydrogen atoms bound to each "heavy" atom in molecule 1. 

    molecule_2_information : list
        This list contains the following information about molecules 2.
            mol_name2 : int
                This is the name of the second molecule we want to compare.
            molecule2_elements : list of str.
                This is the list of elements that make up molecule 2.
            molecule2_positions : 2D np.array
                These are the positions of the atoms in molecule 2.
            no_of_H_on_atoms_in_molecule2 : list of int.
                These are the number of hydrogen atoms bound to each "heavy" atom in molecule 2.

    em_indices_m2_to_m1 : dict
        This dictionary indicates how atom indices in molecule 2 could be mapped onto molecule 1. 
    max_distance_disparity: float
        This is the maximum that any two "could be equivalent" atoms can be between molecule 1 and molecule 2 for molecules 1 and 2 to be considered variant.
    neighbouring_molecules_about_molecules : dict.
        This is the information about the molecules that surround (in the vicinity of) each moleule in the crystal.
    non_hydrogen_molecules : list of ase.Atoms
        This is the list of molecules in the crystal, not including hydrogens.

    symmetric_molecule_pairs : dict
        This dictionary stores which molecules are symmetric to each other, as well as which atom indices in molecule 1 map onto which atom indices in molecule 2. 

    no_of_cpus : int.
        This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.
    
    Returns
    -------
    is_variant : bool.
        This indicates if the two molecules are variants using the comprehensive method. 
    mol1_to_mol2_conversion_Comp : dict. or None
        This dictionary indicates how the atom indices of molecule 1 relate to the atom indices in molecule 2.
    """

    # First, initialise the ``is_variant`` boolean, which indicates if the two molecules being compared are variant or not. 
    is_variant = False

    # Second, initalise a variable to hold how the indices in molecule 1 relate to the indices in molecule 2.
    mol1_to_mol2_conversion_Comp = None

    # Third, obtain the generator for the input data.
    input_generator = get_inputs(molecule_1_information, molecule_2_information, em_indices_m2_to_m1, max_distance_disparity, neighbouring_molecules_about_molecules, non_hydrogen_molecules, symmetric_molecule_pairs)

    # Fourth, perform the are_molecules_variant_comprehensive_single_process method on all inputs depending on if you are performing the task with one cpu or with multiple cores.
    if no_of_cpus == 1: # If the user only wants to use 1 cpu, perform tasks without using multiprocessing

        # 4.1: For each possibility of permutation in em_indices_m2_to_m1
        for input_data in input_generator:

            # 4.2: Determine if the current mol2->mol1 index comparison allows these two molecules to be seen as variant. 
            is_variant, mol1_to_mol2_conversion_Comp = are_molecules_variant_comprehensive_single_process(input_data)

            # 4.3: If are_molecules_variant_comprehensive_single_process throws an exception, break out of the loop
            if isinstance(is_variant, Exception):
                break

            # 4.4: Make sure that at this point, is_variant is a boolean.
            if not isinstance(is_variant, bool):
                break

            # 4.5: If ``is_variant`` is True, break out of the for loop.
            if is_variant == True:
                break

    else:

        # 4.6: Set up the multiprocessing pool.
        with mp.Pool(processes=no_of_cpus) as pool:

            # 4.7: Compare the two molecules via all the possible comparisons of the atom indices in molecule 2 to the atom indices in molecule 1. 
            for is_variant, mol1_to_mol2_conversion_Comp in pool.imap_unordered(are_molecules_variant_comprehensive_single_process, input_generator):

                # 4.8: If are_molecules_variant_comprehensive_single_process throws an exception, break out of the loop
                if isinstance(is_variant, Exception):
                    break

                # 4.9: Make sure that at this point, is_variant is a boolean.
                if not isinstance(is_variant, bool):
                    break

                # 4.10: If ``is_variant`` is True, break out of the for loop.
                if is_variant == True:
                    break

    # Fifth, indicate what the exception is.
    if isinstance(is_variant, Exception):
        exception_message  = 'Error: There was a problem in the "are_molecules_variant_comprehensive_single_process" method somewhere\n\n'
        exception_message += str(is_variant)
        raise Exception(exception_message)

    # Sixth, if is_variant is not a boolean at this point, raise this as a problem.
    if not isinstance(is_variant, bool):
        exception_message  = 'Error: The "is_variant" variable is not a boolean.\n'
        exception_message += '"is_variant" = '+str(is_variant)+'\n'
        exception_message += 'Check this.'
        raise Exception(exception_message)

    # Seventh, return if the two molecules are variants of each other, 
    #          * If the two molecules are variants of each other, mol1_to_mol2_conversion_Comp will also be given.
    #          * mol1_to_mol2_conversion_Comp is a dictionary that gives information about which atoms (indices) 
    #            in molecule 1 are equivalent to which atoms (indices) in molecule 2. 
    return is_variant, mol1_to_mol2_conversion_Comp

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

def get_inputs(molecule_1_information, molecule_2_information, em_indices_m2_to_m1, max_distance_disparity, neighbouring_molecules_about_molecules, non_hydrogen_molecules, symmetric_molecule_pairs):
    """
    This generator is designed to return all the input methods required for the compare_if_two_molecules_are_symmetric_single_process method. 

    Parameters
    ----------
    molecule_1_information : list
        This list contains the following information about molecules 1.
            mol_name1 : int
                This is the name of the first molecule we want to compare.
            molecule1_elements : list of str.
                This is the list of elements that make up molecule 1.
            molecule1_positions : 2D np.array
                These are the positions of the atoms in molecule 1
            no_of_H_on_atoms_in_molecule1 : list of int.
                These are the number of hydrogen atoms bound to each "heavy" atom in molecule 1. 

    molecule_2_information : list
        This list contains the following information about molecules 2.
            mol_name2 : int
                This is the name of the second molecule we want to compare.
            molecule2_elements : list of str.
                This is the list of elements that make up molecule 2.
            molecule2_positions : 2D np.array
                These are the positions of the atoms in molecule 2.
            no_of_H_on_atoms_in_molecule2 : list of int.
                These are the number of hydrogen atoms bound to each "heavy" atom in molecule 2.

    em_indices_m2_to_m1 : dict
        This dictionary indicates how atom indices in molecule 2 could be mapped onto molecule 1. 
    max_distance_disparity: float
        This is the maximum that any two "could be equivalent" atoms can be between molecule 1 and molecule 2 for molecules 1 and 2 to be considered variant.
    neighbouring_molecules_about_molecules : dict.
        This is the information about the molecules that surround (in the vicinity of) each moleule in the crystal.
    non_hydrogen_molecules : list of ase.Atoms
        This is the list of molecules in the crystal, not including hydrogens.

    symmetric_molecule_pairs : dict
        This dictionary stores which molecules are symmetric to each other, as well as which atom indices in molecule 1 map onto which atom indices in molecule 2. 

    Returns
    -------
    names_being_compared : tuple of (int, int)
        This tuple contains the names of the two molecules being compared to each other, being (mol_name1, mol_name2).

    molecule1_elements : list of str.
        This is the list of elements that make up molecule 1.
    molecule1_positions : 2D np.array
        These are the positions of the atoms in molecule 1
    no_of_H_on_atoms_in_molecule1 : list of int.
        These are the number of hydrogen atoms bound to each "heavy" atom in molecule 1. 

    molecule2_elements : list of str.
        This is the list of elements that make up molecule 2.
    molecule2_positions : 2D np.array
        These are the positions of the atoms in molecule 2.
    no_of_H_on_atoms_in_molecule2 : list of int.
        These are the number of hydrogen atoms bound to each "heavy" atom in molecule 2.

    em_indices_m2_to_m1 : dict
        This dictionary indicates how atoms in molecule 2 could be mapped onto molecule 1. 

    max_distance_disparity: float
        This is the maximum that any two "could be equivalent" atoms can be between molecule 1 and molecule 2 for molecules 1 and 2 to be considered variant.
    neighbouring_molecules_about_molecules : dict.
        This is the information about the molecules that surround (in the vicinity of) each moleule in the crystal.
    non_hydrogen_molecules : list of ase.Atoms
        This is the list of molecules in the crystal, not including hydrogens.
    """

    # First, extract the information about molecule 1 from the molecule_1_information list.
    mol_name1, molecule1_elements, molecule1_positions, no_of_H_on_atoms_in_molecule1 = molecule_1_information

    # Second, extract the information about molecule 2 from the molecule_2_information list.
    mol_name2, molecule2_elements, molecule2_positions, no_of_H_on_atoms_in_molecule2 = molecule_2_information

    # Third, for each way that the atom indices of atoms could be ordered in molecule 2 to map onto molecule 1. 
    for comparison in em_indices_m2_to_m1:

        # 3.1: Reorder the atoms in molecule 2 to allign with atoms in molecule 1. 
        idx_m2                                  = get_permutated_indices_list(comparison)
        molecule2_elements_reordered            = [molecule2_elements[index] for index in idx_m2]
        molecule2_positions_reordered           = deepcopy(molecule2_positions)[idx_m2, :]
        no_of_H_on_atoms_in_molecule2_reordered = [no_of_H_on_atoms_in_molecule2[index] for index in idx_m2]

        # 3.2: Make a tuple of the names of the molecules being compared against each other.
        names_being_compared = (mol_name1, mol_name2)

        # 3.3: Yield the information needed to run the "are_molecules_variant_comprehensive_single_process" method. 
        yield (names_being_compared,  molecule1_elements, molecule1_positions, no_of_H_on_atoms_in_molecule1,  molecule2_elements_reordered, molecule2_positions_reordered, no_of_H_on_atoms_in_molecule2_reordered,  idx_m2,  max_distance_disparity, neighbouring_molecules_about_molecules, non_hydrogen_molecules)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

def are_molecules_variant_comprehensive_single_process(input_values):
    """
    This method will return if two molecules are variant based on a set of rearranged molecule 2 atom indices (compared to molecule 1 atom indices). 

    Parameters
    ----------
    mol_name1 : int
        This is the name of the first molecule we want to compare.
    molecule1_elements : list of str.
        This is the list of elements that make up molecule 1.
    molecule1_positions : 2D np.array
        These are the positions of the atoms in molecule 1
    no_of_H_on_atoms_in_molecule1 : list of int.
        These are the number of hydrogen atoms bound to each "heavy" atom in molecule 1. 

    mol_name2 : int
        This is the name of the second molecule we want to compare.
    molecule2_elements : list of str.
        This is the list of elements that make up molecule 2.
    molecule2_positions : 2D np.array
        These are the positions of the atoms in molecule 2.
    no_of_H_on_atoms_in_molecule2 : list of int.
        These are the number of hydrogen atoms bound to each "heavy" atom in molecule 2.

    em_indices_m2_to_m1 : dict
        This dictionary indicates how atoms in molecule 2 could be mapped onto molecule 1. 

    max_distance_disparity: float
        This is the maximum that any two "could be equivalent" atoms can be between molecule 1 and molecule 2 for molecules 1 and 2 to be considered variant.
    neighbouring_molecules_about_molecules : dict.
        This is the information about the molecules that surround (in the vicinity of) each moleule in the crystal.
    non_hydrogen_molecules : list of ase.Atoms
        This is the list of molecules in the crystal, not including hydrogens.
    
    Returns
    -------
    is_variant_comprehensive : bool.
        This indicates if the two molecules are variants using the comprehensive method. 
    mol1_to_mol2_conversion_Comp : dict. or None
        This dictionary indicates how the atom indices of molecule 1 relate to the atom indices in molecule 2.
    """

    # First, obtain the input variables from the "input_values" tuple.
    names_being_compared,  molecule1_elements, molecule1_positions, no_of_H_on_atoms_in_molecule1,  molecule2_elements_reordered, molecule2_positions_reordered, no_of_H_on_atoms_in_molecule2_reordered,  idx_m2,  max_distance_disparity, neighbouring_molecules_about_molecules, non_hydrogen_molecules = input_values

    # Second, determine if the two molecules are variant.
    is_variant_comprehensive = are_systems_variant(molecule1_elements, molecule1_positions, no_of_H_on_atoms_in_molecule1, molecule2_elements_reordered, molecule2_positions_reordered, no_of_H_on_atoms_in_molecule2_reordered, max_distance_disparity=max_distance_disparity, names_being_compared=names_being_compared, neighbouring_molecules_about_systems=neighbouring_molecules_about_molecules, non_hydrogen_systems=non_hydrogen_molecules)

    # Third, if the molecule is variant, give a dictionary that indicates how the atom indices of molecule 1 relate to the atom indices in molecule 2.
    if is_variant_comprehensive:
        #mol1_to_mol2_conversion_Comp = {atomic_index_of_mol_1: atomic_index_of_mol_2 for atomic_index_of_mol_1, atomic_index_of_mol_2 in zip(range(len(molecule1_elements)), idx_m2)}
        mol1_to_mol2_conversion_Comp = list([(atomic_index_of_mol_1, atomic_index_of_mol_2) for atomic_index_of_mol_1, atomic_index_of_mol_2 in zip(range(len(molecule1_elements)), idx_m2)])
    else:
        mol1_to_mol2_conversion_Comp = None

    # Fourth, return if the two molecules are variants of each other or not.
    return is_variant_comprehensive, mol1_to_mol2_conversion_Comp

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

