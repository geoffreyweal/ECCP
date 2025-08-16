"""
get_points_of_second_dimer.py, Geoffrey Weal, 12/2/24

This method will give the centre, 1st, 2nd, and 3rd points for the second dimer. 
"""
import numpy as np
from SUMELF import get_unit_vector

def get_points_of_second_dimer(d2_centre_position, d2_point1_position, d2_point2_position, d2_point3_position):
    """
    This method will give the centre, 1st, 2nd, and 3rd points for the second dimer. 

    Parameters
    ----------
    d2_centre_position : numpy.array
        This is the centre point given for dimer 1.
    d2_point1_position : numpy.array
        This is the first point given for dimer 1.
    d2_point2_position : numpy.array
        This is the second point given for dimer 1.
    d2_point3_position : numpy.array
        This is the third point given for dimer 1.

    Returns
    -------
    d2_direction1 : numpy.array
        This is the first direction in dimer 2.
    d2_direction2 : numpy.array
        This is the second direction in dimer 2.
    d2_direction3 : numpy.array
        This is the third direction in dimer 2.

    lengths_d2 : list of floats
        These are the lengths of the direction vectors in dimer 2.
    dotproducts_d2 : list of floats
        These are the dot products (angles) between direction vectors in dimer 2.
    """
    
    # First, obtain the directions for the second dimer using the four input points.
    d2_direction1 = d2_point1_position - d2_centre_position
    d2_direction2 = d2_point2_position - d2_centre_position
    d2_direction3 = d2_point3_position - d2_centre_position

    # Second, give the directions as a list.
    d2_directions = (d2_direction1, d2_direction2, d2_direction3)

    # Third, get the vector length of directions in dimer 2.
    lengths_d2 = [np.linalg.norm(d2_direction) for d2_direction in d2_directions]

    # Fourth, get the dot products (angles) between directions for dimer 2.
    dotproducts_d2 = []
    for index1 in range(len(d2_directions)):
        d2_direction_a = d2_directions[index1]
        for index2 in range(index1+1,len(d2_directions)):
            d2_direction_b = d2_directions[index2]
            dotproduct_d2 = np.dot(get_unit_vector(d2_direction_a), get_unit_vector(d2_direction_b))
            dotproducts_d2.append(dotproduct_d2)
    dotproducts_d2_sorted = sorted(dotproducts_d2)

    # Fifth, return direction details for dimer 2.
    return d2_direction1, d2_direction2, d2_direction3, lengths_d2, dotproducts_d2_sorted