"""
obtain_possible_3D_points_in_dimer2.py, Geoffrey Weal, 12/2/24

This method provides the dimer indices of the centre, 1st, 2nd, and 3rd point to obtain spatial information about the dimer
"""

def obtain_possible_3D_points_in_dimer2(d2_m1_positions_of_lowest_elements, d2_m2_positions_of_lowest_elements, d1_centre_element, d1_point1_element, d1_point2_element, d1_point3_element, no_of_atoms_in_d1_m1):
    """
    This method provides the dimer indices of the centre, 1st, 2nd, and 3rd point to obtain spatial information about the dimer
    
    Parameters
    ----------
    d2_m1_positions_of_lowest_elements : list
        The elements and positions of the elements in lowest abunance in the first mlecules of the second dimer, from lowest to highest abundance.
    d2_m2_positions_of_lowest_elements : list
        The elements and positions of the elements in lowest abunance in the second mlecules of the second dimer, from lowest to highest abundance.

    d1_centre_element : str. 
        This is the element of the centre point atom in the first dimer.
    d1_point1_element : str. 
        This is the element of the 1st point atom in the first dimer.
    d1_point2_element : str. 
        This is the element of the 2nd point atom in the first dimer.
    d1_point3_element : str. 
        This is the element of the 3rd point atom in the first dimer.

    no_of_atoms_in_d1_m1 : int
        This is the number of atoms in molecule 1 of dimer 2

    Returns
    -------
    d2_centre_element : str.
        This is the element of the centre atom in dimer 2.
    d2_centre_position : numpy.array
        This is the position of the centre atom in dimer 2.
    d2_centre_index : int
        This is the index of the centre atom in dimer 2.

    d2_point1_element : str.
        This is the element of point 1 in dimer 2.
    d2_point1_position : numpy.array
        This is the position of point 1 in dimer 2.
    d2_point1_index : int
        This is the index of point 1 in dimer 2.

    d2_point2_element : str.
        This is the element of point 2 in dimer 2.
    d2_point2_position : numpy.array
        This is the position of point 2 in dimer 2.
    d2_point2_index : int
        This is the index of point 2 in dimer 2.

    d2_point3_element : str.
        This is the element of point 3 in dimer 2.
    d2_point3_position : numpy.array
        This is the position of point 3 in dimer 2.
    d2_point3_index : int
        This is the index of point 3 in dimer 2.
    """

    # 1. Look through the first molecule
    for d2_m1_index1 in range(len(d2_m1_positions_of_lowest_elements)):
        # obtain the position of the centre point of dimer 2.
        d2_centre_element, (d2_centre_position, d2_centre_index) = d2_m1_positions_of_lowest_elements[d2_m1_index1]
        if not (d2_centre_element == d1_centre_element):
            continue
        for d2_m1_index2 in range(len(d2_m1_positions_of_lowest_elements)):
            # obtain the position of the 1st point of dimer 2.
            if d2_m1_index1 == d2_m1_index2:
                continue
            d2_point1_element, (d2_point1_position, d2_point1_index) = d2_m1_positions_of_lowest_elements[d2_m1_index2]
            if not (d2_point1_element == d1_point1_element):
                continue

            # 2. Look through the first molecule
            for d2_m2_index1 in range(len(d2_m2_positions_of_lowest_elements)):
                # obtain the position of the 2nd point of dimer 2.
                d2_point2_element, (d2_point2_position, d2_point2_index) = d2_m2_positions_of_lowest_elements[d2_m2_index1]
                d2_point2_index += no_of_atoms_in_d1_m1
                if not (d2_point2_element == d1_point2_element):
                    continue
                for d2_m2_index2 in range(len(d2_m2_positions_of_lowest_elements)):
                    # obtain the position of the 3rd point of dimer 2.
                    if d2_m2_index1 == d2_m2_index2:
                        continue
                    d2_point3_element, (d2_point3_position, d2_point3_index) = d2_m2_positions_of_lowest_elements[d2_m2_index2]
                    d2_point3_index += no_of_atoms_in_d1_m1
                    if not (d2_point3_element == d1_point3_element):
                        continue

                    # Have centre, 1st, 2nd, and 3rd points to look at in the second dimer.
                    yield (d2_centre_element, d2_centre_position, d2_centre_index, d2_point1_element, d2_point1_position, d2_point1_index, d2_point2_element, d2_point2_position, d2_point2_index, d2_point3_element, d2_point3_position, d2_point3_index)
