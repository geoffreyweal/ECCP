"""
get_points_of_first_molecule.py, Geoffrey Weal, 12/2/24

These methods are designed to give the centre, 1st, 2nd, and 3rd points for the first molecule. 
"""
import numpy as np
from SUMELF import get_unit_vector, are_two_values_within_eachother, get_distance

def get_points_of_first_molecule(m1_positions_of_lowest_elements):
    """
    This method will give the centre, 1st, 2nd, and 3rd points for the first molecule. 

    Parameters
    ----------
    m1_positions_of_lowest_elements : list of (element, position) 
        These are the elements and positions of atoms in molecule 1, from lowest abundance to highest abundance.

    Returns
    -------
    m1_direction1 : numpy.array
        This is the first direction in dimer 1.
    m1_direction2 : numpy.array
        This is the second direction in dimer 1.
    m1_direction3 : numpy.array
        This is the third direction in dimer 1.

    lengths_m1 : list of floats
        These are the lengths of the direction vectors in dimer 1.
    dotproducts_m1 : list of floats
        These are the dot products (angles) between direction vectors in dimer 1.

    indices_of_points_m1 : list of ints
        These are the indices of the points given for dimer 1.
    elements_in_points : list of str.
        These are the element of the points given for dimer 1.
    """

    # First, obtain the centre position.
    m1_centre_element, (m1_centre_position, m1_centre_index) = m1_positions_of_lowest_elements[0]

    # Second, obtain the first direction.
    m1_point1_element, (m1_point1_position, m1_point1_index) = m1_positions_of_lowest_elements[1]
    m1_direction1 = m1_point1_position - m1_centre_position

    # ------------------------------------------------------------------------------------------------------------
    # Third, obtain the second direction.
    for index2 in range(2,len(m1_positions_of_lowest_elements)):
        m1_point2_element, m1_point2_position, m1_point2_index, m1_direction2, dotproducts_m1_set1 = get_vector_details(m1_positions_of_lowest_elements, index2, m1_centre_position, [m1_direction1])
        if all([(not are_two_values_within_eachother(dotproduct_m1_set1,  1.0, 0.1)) for dotproduct_m1_set1 in dotproducts_m1_set1]):
            break
        if all([(not are_two_values_within_eachother(dotproduct_m1_set1, -1.0, 0.1)) for dotproduct_m1_set1 in dotproducts_m1_set1]):
            break
    else:
        raise Exception('MESSAGE from def get_points_of_first_molecule, methods_for_MEA_invariance_method.py in ECCP. All atoms have been found roughly along a line (a linear molecule).'+'\n'+'Check if this is a problem? Not sure of this yet and have not been able to find a good example to test this.')
        index2 = 2
        m1_point2_element, m1_point2_position, m1_point2_index, m1_direction2, dotproducts_m1_set1 = get_vector_details(m1_positions_of_lowest_elements, index2, m1_centre_position, [m1_direction1])
    # ------------------------------------------------------------------------------------------------------------
    # Fourth, obtain the third direction.
    for index3 in range(index2+1,len(m1_positions_of_lowest_elements)):
        m1_point3_element, m1_point3_position, m1_point3_index, m1_direction3, dotproducts_m1_set2 = get_vector_details(m1_positions_of_lowest_elements, index3, m1_centre_position, [m1_direction1, m1_direction2])
        #import pdb; pdb.set_trace()
        if all([(not are_two_values_within_eachother(dotproduct_m1_set2,  1.0, 0.1)) for dotproduct_m1_set2 in dotproducts_m1_set2]):
            break
        if all([(not are_two_values_within_eachother(dotproduct_m1_set2, -1.0, 0.1)) for dotproduct_m1_set2 in dotproducts_m1_set2]):
            break
    else:
        raise Exception('MESSAGE from def get_points_of_first_molecule, methods_for_MEA_invariance_method.py in ECCP. All atoms have been found roughly along a line (a linear molecule).'+'\n'+'Check this at some point. For now if this message comes up, dont use molecule?')
        index3 = index2+1
        m1_point3_element, m1_point3_position, m1_point3_index, m1_direction3, dotproducts_m1_set2 = get_vector_details(m1_positions_of_lowest_elements, index3, m1_centre_position, [m1_direction1, m1_direction2])
    # ------------------------------------------------------------------------------------------------------------

    # Fifth, obtain the lengths and dot products of directions.
    lengths_m1 = [np.linalg.norm(m1_direction) for m1_direction in [m1_direction1, m1_direction2, m1_direction3]]
    dotproducts_m1_sorted = sorted(dotproducts_m1_set1 + dotproducts_m1_set2)

    # Sixth, give a tuple of the indices and elements of the points in the first dimer.
    indices_of_points_m1 = (m1_centre_index, m1_point1_index, m1_point2_index, m1_point3_index)
    elements_in_points_m1 = (m1_centre_element, m1_point1_element, m1_point2_element, m1_point3_element)

    # Seventh, return direction details for dimer 1.
    return m1_direction1, m1_direction2, m1_direction3, m1_centre_position, lengths_m1, dotproducts_m1_sorted, indices_of_points_m1, elements_in_points_m1

# --------------------------------------------------------------------------------------------------------------

def get_vector_details(mol_positions_of_lowest_elements, index, mol_centre_position, mol_other_directions):
    """
    This method will return details about the point given by index in mol_positions_of_lowest_elements

    Parameters
    ----------
    mol_positions_of_lowest_elements : list of (element, position) 
        These are the elements and positions of atoms in the molecule from lowest abundance to highest abundance.
    index : int
        This is the index to get information from in mol_positions_of_lowest_elements.
    mol_centre_position : numpy.array
        This is the centre point in dimer 1.
    mol_other_directions : list of numpy.array
        These are all the other directions that have been obtained to get dot-products of.

    Returns
    -------
    mol_point_element : str.
        This is the element of the point given by mol_positions_of_lowest_elements[index]
    mol_point_position : numpy.array
        This is the position of the point given by mol_positions_of_lowest_elements[index]
    mol_point_index : int
        This is the index of this atom point in the dimer
    mol_direction : numpy.array
        This is the direction from the centre point.
    dotproducts_mol : list of float
        This is a list of the dot products from this direction and all other directions obtained in mol_other_directions.
        
    """

    # First, get the element and position of the point given by the index from mol_positions_of_lowest_elements. 
    mol_point_element, (mol_point_position, mol_point_index) = mol_positions_of_lowest_elements[index]

    # Second, get the direction from this point to the centre point.
    mol_direction = mol_point_position - mol_centre_position
    if get_distance(mol_direction,np.array([0.,0.,0.])) < 0.001:
        raise Exception('Error, mol_direction is close to 0, which it should not be. mol_direction = '+str(mol_direction))

    # Third, get the dot product between this direction and all other directions in mol_other_directions.
    dotproducts_mol = sorted([np.dot(get_unit_vector(mol_direction), get_unit_vector(mol_other_direction)) for mol_other_direction in mol_other_directions])
    
    # Fourth, return quantities.
    return mol_point_element, mol_point_position, mol_point_index, mol_direction, dotproducts_mol

# --------------------------------------------------------------------------------------------------------------