"""
are_systems_variant.py, Geoffrey Weal, 19/3/22

This script is designed to determine if two chemical systems (such as two molecules or dimers) are rotationally variant.
"""
import numpy as np
from scipy.linalg import orthogonal_procrustes

from SUMELF import get_centre_of_mass
from ECCP.ECCP.invariance_methods.common_utility_methods_for_all_invariance_methods.are_systems_variant_utility_methods.are_environments_equivalent import are_environments_equivalent

def are_systems_variant(system1_elements, system1_positions, no_of_H_on_atoms_in_system1, system2_elements, system2_positions, no_of_H_on_atoms_in_system2, max_distance_disparity, names_being_compared, neighbouring_molecules_about_systems, non_hydrogen_systems): 
    """
    This method is designed to determine if two chemical systems (such as two molecules or dimers) are rotationally variant.

    Note: System (chemical system) is either a molecule or dimer.
    
    Parameters
    ----------
    system1_elements : list
        A list of the elements of atoms in system 1.
    system1_positions : np.array
        A numpy array of the positions of atoms in system 1.
    no_of_H_on_atoms_in_system1 : list
        A list of the number of hydrogens attached to each atoms in system 1.

    system2_elements : list
        A list of the elements of atoms in system 2.
    system2_positions : np.array
        A numpy array of the positions of atoms in system 2.
    no_of_H_on_atoms_in_system2 : list
        A list of the number of hydrogens attached to each atoms in system 2.

    max_distance_disparity: float
        This is the maximum that any two "could be equivalent" atoms can be between system 1 and system 2 for systems 1 and 2 to be considered variant.
    names_being_compared: list
        These are the names of the systems being compared.
    neighbouring_molecules_about_systems : dict.
        This is the information about the molecules that surround (in the vicinity of) each system (molecule or dimer) in the crystal.

    non_hydrogen_systems : list of ase.Atoms
        This is the list of systems in the crystal, not including hydrogens.

    Returns
    -------
    True if the two systems are variant, False if they are invariant (unique)
    """

    # First, if the order of atoms are not the same, something actually may have gone wrong, so check this out
    if not (sorted(system1_elements) == sorted(system2_elements)):
        print('Error in comprehensive_invariance_method.py')
        print('The list of elements for system1 and system2 are not the same. ')
        print('However, they should be the same as the networkx graph nodes have been given information about the element for each atom.')
        print('Therefore, this should have been picked up by GraphMatcher object.')
        print('Check this out')
        import pdb; pdb.set_trace()
        exit('This program will finish without completing.')

    # Second, check that all the lists are the same length, as the system to this point have the same number of atoms. 
    if not (len(system1_elements) == len(system2_elements) == len(system1_positions) == len(system2_positions) == len(no_of_H_on_atoms_in_system1) == len(no_of_H_on_atoms_in_system2)):
        print('Error in comprehensive_invariance_method.py')
        print('Not all the lists between system 1 and system 2 are the same.')
        print('These lists should make be all the same ')
        print('system1_elements = '+str(system1_elements))
        print('system2_elements = '+str(system2_elements))
        print('system1_positions = '+str(system1_positions))
        print('system2_positions = '+str(system2_positions))
        print('no_of_H_on_atoms_in_system1 = '+str(no_of_H_on_atoms_in_system1))
        print('no_of_H_on_atoms_in_system2 = '+str(no_of_H_on_atoms_in_system2))
        print('Check this out')
        import pdb; pdb.set_trace()
        exit('This program will finish without completing.')

    # Third, check that all the hydrogens bound to each atoms of a system are the same between system 1 and system 2. 
    if not all([(nH1 == nH2) for nH1, nH2 in zip(no_of_H_on_atoms_in_system1, no_of_H_on_atoms_in_system2)]):
        return False

    # Fourth, determine if the systems are variant given the particular ordering of atoms in systems 1 and 2.
    positions_are_variant, rotation_reflection_matrix, system1_com, system2_com = determine_if_positions_are_variant(system1_elements, system1_positions, system2_elements, system2_positions, max_distance_disparity=max_distance_disparity)
    
    # Fifth, if the two systems are not varient, return False
    if not positions_are_variant:
        return False 

    # Sixth, determine if the environment around these two systems are the same. If not, return False
    if not are_environments_equivalent(rotation_reflection_matrix, system1_com, system2_com, names_being_compared, neighbouring_molecules_about_systems, non_hydrogen_systems):
        return False

    # Seventh, if got to here, systems and their enivornments are the same, so return True
    return True

# ==================================================================================================================================

def determine_if_positions_are_variant(elements1, positions1, elements2, positions2, max_distance_disparity):
    """
    This method will determine if two systems are translationally, rotationally, and reflectively variant.

    This method have been modified from the _procrustes.py from scipy, see below: https://github.com/scipy/scipy/blob/v1.8.0/scipy/spatial/_procrustes.py#L15-L130

    Note: System (chemical system) is either a molecule or dimer.

    Parameters
    ----------
    elements1 : list
        A list of the elements of atoms in system1
    positions1 : numpy.array
        A array of positions from system1

    elements2 : list
        A list of the elements of atoms in system2
    positions2 : numpy.array
        A array of positions from system2

    max_distance_disparity: float
        This is the maximum that any two "could be equivalent" atoms can be between dimer 1 and dimer 2 for dimers 1 and 2 to be considered variant.

    Returns
    -------
    True if the two dimers are translationally, rotationally, and reflectively variant. False if not. 

    """

    # First, copy the positions of the two systems. 
    mtx1 = np.array(positions1, dtype=np.double, copy=True)
    mtx2 = np.array(positions2, dtype=np.double, copy=True)

    # Second, make sure that the dimensions of the two positions are the same.
    if mtx1.ndim != 2 or mtx2.ndim != 2:
        raise ValueError("Input matrices must be two-dimensional")
    if mtx1.shape != mtx2.shape:
        raise ValueError("Input matrices must be of same shape")
    if mtx1.size == 0:
        raise ValueError("Input matrices must be >0 rows and >0 cols")

    # Third, translate all the data to the origin.
    #mtx1 -= np.mean(mtx1, 0)
    #mtx2 -= np.mean(mtx2, 0)
    # Instead, translate to the centre of mass. 
    mtx1 -= get_centre_of_mass(elements1, mtx1)
    mtx2 -= get_centre_of_mass(elements2, mtx2)

    '''
    # We dont require scaling variance in our analysis, so ignored this part of the scipy code. 
    norm1 = np.linalg.norm(mtx1)
    norm2 = np.linalg.norm(mtx2)

    if norm1 == 0 or norm2 == 0:
        raise ValueError("Input matrices must contain >1 unique points")

    # change scaling of data (in rows) such that trace(mtx*mtx') = 1
    mtx1 /= norm1
    mtx2 /= norm2
    '''

    # Fourth, transform mtx2 to minimize disparity. This method also takes into account reflective variance.
    R, ss = orthogonal_procrustes(mtx1, mtx2) 
    #estimated_rotation, rmsd = Rotation.align_vectors(mtx1, mtx2)
    #R = estimated_rotation.as_matrix()

    # Fifth, obtain the rotated version of mtx2.
    mtx2_rotated = np.dot(mtx2, R.T) ##* ss
    #mtx2_rotated = np.matmul(R,mtx2.T).T 

    # Sixth, get the difference between atom positions between each dimer
    difference_in_atom_positions = np.linalg.norm(mtx1 - mtx2_rotated, axis=1)

    # Seventh, determine if the atoms of each system overlap spatially after rotating mtx2. 
    do_systems_overlap = all([(distance <= max_distance_disparity) for distance in difference_in_atom_positions])

    # Eighth, return the following:
    #         * The result of do_systems_overlap, 
    #         * the rotation matrix used to rotate mtx2 onto mtx1, and
    #         * the numpy array version of positions1 and positions2 which have been translated to move their centre of masses to the origin.
    return do_systems_overlap, R.T, mtx1, mtx2

# ==================================================================================================================================


