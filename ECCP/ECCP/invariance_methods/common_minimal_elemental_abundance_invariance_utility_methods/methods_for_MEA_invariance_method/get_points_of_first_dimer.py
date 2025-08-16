"""
get_points_of_first_dimer.py, Geoffrey Weal, 12/2/24

These methods are designed to give the centre, 1st, 2nd, and 3rd points for the first dimer. 
"""
import numpy as np
from ase           import Atoms
from ase.visualize import view
from SUMELF        import get_unit_vector, are_two_values_within_eachother

def get_points_of_first_dimer(d1_m1_positions_of_lowest_elements, d1_m2_positions_of_lowest_elements, no_of_atoms_in_d1_m1, d1_m2_original_elements=None, d1_m2_original_positions=None):
    """
    This method will give the centre, 1st, 2nd, and 3rd points for the first dimer. 

    Parameters
    ----------
    d1_m1_positions_of_lowest_elements : list of (element, position) 
        These are the elements and positions of atoms in molecule 1, dimer 1, from lowest abundance to highest abundance.
    d1_m2_positions_of_lowest_elements : list of (element, position) 
        These are the elements and positions of atoms in molecule 2, dimer 1, from lowest abundance to highest abundance.

    no_of_atoms_in_d1_m1 : int
        This is the number of atoms in molecule 1 or dimer 1.

    d1_m2_original_elements : list
        These are the elements of atoms in the second molecule of dimer 1. Used for debugging issues. Default: None. 
    d1_m2_original_positions : numpy.array
        These are the positions of atoms in the second molecule of dimer 1. Used for debugging issues. Default: None. 

    Returns
    -------
    d1_direction1 : numpy.array
        This is the first direction in dimer 1.
    d1_direction2 : numpy.array
        This is the second direction in dimer 1.
    d1_direction3 : numpy.array
        This is the third direction in dimer 1.

    lengths_d1 : list of floats
        These are the lengths of the direction vectors in dimer 1.
    dotproducts_d1 : list of floats
        These are the dot products (angles) between direction vectors in dimer 1.

    indices_of_points_d1 : list of ints
        These are the indices of the points given for dimer 1.
    elements_in_points : list of str.
        These are the element of the points given for dimer 1.
    """

    # First, obtain the centre position.
    d1_centre_element, (d1_centre_position, d1_centre_index) = d1_m1_positions_of_lowest_elements[0]

    # Second, obtain the first direction.
    d1_point1_element, (d1_point1_position, d1_point1_index) = d1_m1_positions_of_lowest_elements[1]
    d1_direction1 = d1_point1_position - d1_centre_position

    # Third, obtain the second direction.
    for index2 in range(len(d1_m2_positions_of_lowest_elements)):
        d1_point2_element, d1_point2_position, d1_point2_index, d1_direction2, dotproducts_d1_set1 = get_vector_details(d1_m2_positions_of_lowest_elements, index2, d1_centre_position, [d1_direction1])
        if all([(not are_two_values_within_eachother(dotproduct_d1_set1,  1.0, 0.1)) for dotproduct_d1_set1 in dotproducts_d1_set1]):
            break
        if all([(not are_two_values_within_eachother(dotproduct_d1_set1, -1.0, 0.1)) for dotproduct_d1_set1 in dotproducts_d1_set1]):
            break
    else:
        # If either these if statements are true, then it possible that every point in molecule 2 can be traced along a line given by:
        # y = (d1_direction1-d1_centre_element)x + d1_centre_element
        # In this case, it doesn't matter which atom to take from d1_m2_positions_of_lowest_elements, so take the first one.
        index2 = 0
        d1_point2_element, d1_point2_position, d1_point2_index, d1_direction2, dotproducts_d1_set1 = get_vector_details(d1_m2_positions_of_lowest_elements, index2, d1_centre_position, [d1_direction1])

    # Fourth, obtain the third direction.
    for index3 in range(index2+1,len(d1_m2_positions_of_lowest_elements)):
        d1_point3_element, d1_point3_position, d1_point3_index, d1_direction3, dotproducts_d1_set2 = get_vector_details(d1_m2_positions_of_lowest_elements, index3, d1_centre_position, [d1_direction1, d1_direction2])
        if all([(not are_two_values_within_eachother(dotproduct_d1_set2,  1.0, 0.1)) for dotproduct_d1_set2 in dotproducts_d1_set2]):
            break
        if all([(not are_two_values_within_eachother(dotproduct_d1_set2, -1.0, 0.1)) for dotproduct_d1_set2 in dotproducts_d1_set2]):
            break
    else:
        # If either these if statements are true, then it possible that every point in molecule 2 can be traced along a plane given by:
        # ax + by + cz = d, where (a,b,c) = n = (d1_direction1-d1_centre_element) cross (d1_direction2-d1_centre_element)
        # In this case, take an index from d1_m2_positions_of_lowest_elements that is not the same as index2
        if not index2+1 >= len(d1_m2_positions_of_lowest_elements):
            index3 = index2+1
        else:
            for index3 in range(len(d1_m2_positions_of_lowest_elements)):
                if not index3 == index2:
                    break
            else:
                if (d1_m2_original_elements is not None) and (d1_m2_original_positions is not None):
                    from ase import Atoms
                    from ase.visualize import view
                    molecule2 = Atoms(symbols=d1_m2_original_elements, positions=d1_m2_original_positions)
                    view(molecule2)
                    raise Exception('Error in def get_points_of_first_dimer, in methods_for_quick_invariance_method.py: Your molecule2 may only has one atom in it?')
        d1_point3_element, d1_point3_position, d1_point3_index, d1_direction3, dotproducts_d1_set2 = get_vector_details(d1_m2_positions_of_lowest_elements, index3, d1_centre_position, [d1_direction1, d1_direction2])

    # Fifth, btain the lengths and dot products of directions.
    lengths_d1 = [np.linalg.norm(d1_direction) for d1_direction in [d1_direction1, d1_direction2, d1_direction3]]
    dotproducts_d1_sorted = sorted(dotproducts_d1_set1 + dotproducts_d1_set2)

    # Sixth, give a tuple of the indices and elements of the points in the first dimer.
    indices_of_points_d1 = (d1_centre_index, d1_point1_index, no_of_atoms_in_d1_m1+d1_point2_index, no_of_atoms_in_d1_m1+d1_point3_index)
    elements_in_points_d1 = (d1_centre_element, d1_point1_element, d1_point2_element, d1_point3_element)

    # Seventh, return direction details for dimer 1.
    return d1_direction1, d1_direction2, d1_direction3, d1_centre_position, lengths_d1, dotproducts_d1_sorted, indices_of_points_d1, elements_in_points_d1

# --------------------------------------------------------------------------------------------------------------

def get_vector_details(d1_m2_positions_of_lowest_elements, index, d1_centre_position, d1_other_directions):
    """
    This method will return details about the point given by index in d1_m2_positions_of_lowest_elements

    Parameters
    ----------
    d1_m2_positions_of_lowest_elements : list of (element, position) 
        These are the elements and positions of atoms in molecule 2, dimer 1, from lowest abundance to highest abundance.
    index : int
        This is the index to get information from in d1_m2_positions_of_lowest_elements.
    d1_centre_position : numpy.array
        This is the centre point in dimer 1.
    d1_other_directions : list of numpy.array
        These are all the other directions that have been obtained to get dot-products of.

    Returns
    -------
    d1_point_element : str.
        This is the element of the point given by d1_m2_positions_of_lowest_elements[index]
    d1_point_position : numpy.array
        This is the position of the point given by d1_m2_positions_of_lowest_elements[index]
    d1_point_index : int
        This is the index of this atom point in the dimer

    d1_direction : numpy.array
        This is the direction from the centre point.
    dotproducts_d1 : list of float
        This is a list of the dot products from this direction and all other directions obtained in d1_other_directions.    
    """

    # First, get the element and position of the point given by the index from d1_m2_positions_of_lowest_elements. 
    d1_point_element, (d1_point_position, d1_point_index) = d1_m2_positions_of_lowest_elements[index]

    # Second, get the direction from this point to the centre point.
    d1_direction = d1_point_position - d1_centre_position

    # Third, get the dot product between this direction and all other directions in d1_other_directions.
    dotproducts_d1 = sorted([np.dot(get_unit_vector(d1_direction), get_unit_vector(d1_other_direction)) for d1_other_direction in d1_other_directions])
    
    # Fourth, return quantities.
    return d1_point_element, d1_point_position, d1_point_index, d1_direction, dotproducts_d1

# --------------------------------------------------------------------------------------------------------------




