"""
get_points_of_second_molecule.py, Geoffrey Weal, 12/2/24

This method will give the centre, 1st, 2nd, and 3rd points for the second molecule. 
"""
import numpy as np
from SUMELF import get_unit_vector

def get_points_of_second_molecule(m2_centre_position, m2_point1_position, m2_point2_position, m2_point3_position):
    """
    This method will give the centre, 1st, 2nd, and 3rd points for the second molecule. 

    Parameters
    ----------
    m2_centre_position : numpy.array
        This is the centre point given for the second molecule. 
    m2_point1_position : numpy.array
        This is the first point given for the second molecule. 
    m2_point2_position : numpy.array
        This is the second point given for the second molecule. 
    m2_point3_position : numpy.array
        This is the third point given for the second molecule. 

    Returns
    -------
    m2_direction1 : numpy.array
        This is the first direction in molecule 2.
    m2_direction2 : numpy.array
        This is the second direction in molecule 2.
    m2_direction3 : numpy.array
        This is the third direction in molecule 2.

    lengths_m2 : list of floats
        These are the lengths of the direction vectors in molecule 2.
    dotproducts_m2 : list of floats
        These are the dot products (angles) between direction vectors in molecule 2.
    """
    
    # First, obtain the directions for the second molecule using the four input points.
    m2_direction1 = m2_point1_position - m2_centre_position
    m2_direction2 = m2_point2_position - m2_centre_position
    m2_direction3 = m2_point3_position - m2_centre_position

    # Second, give the directions as a list.
    m2_directions = (m2_direction1, m2_direction2, m2_direction3)

    # Third, get the vector length of directions in molecule 2.
    lengths_m2 = [np.linalg.norm(m2_direction) for m2_direction in m2_directions]

    # Fourth, get the dot products (angles) between directions for molecule 2.
    dotproducts_m2 = []
    for index1 in range(len(m2_directions)):
        m2_direction_a = m2_directions[index1]
        for index2 in range(index1+1,len(m2_directions)):
            m2_direction_b = m2_directions[index2]
            dotproduct_m2 = np.dot(get_unit_vector(m2_direction_a), get_unit_vector(m2_direction_b))
            dotproducts_m2.append(dotproduct_m2)
    dotproducts_m2_sorted = sorted(dotproducts_m2)

    # Fifth, return direction details for molecule 2.
    return m2_direction1, m2_direction2, m2_direction3, lengths_m2, dotproducts_m2_sorted