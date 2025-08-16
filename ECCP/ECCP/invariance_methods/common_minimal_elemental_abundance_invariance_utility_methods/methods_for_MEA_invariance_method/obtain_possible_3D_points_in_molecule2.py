"""
obtain_possible_3D_points_in_molecule2.py, Geoffrey Weal, 12/2/24

This method provides the indices of the centre, 1st, 2nd, and 3rd point to obtain spatial information about the second molecule.
"""

def obtain_possible_3D_points_in_molecule2(m2_positions_of_lowest_elements, m1_centre_element, m1_point1_element, m1_point2_element, m1_point3_element):
    """
    This method provides the indices of the centre, 1st, 2nd, and 3rd point to obtain spatial information about the second molecule.
    
    Parameters
    ----------
    m2_positions_of_lowest_elements : list
        The elements and positions of the elements in lowest abunance in the second molecule from lowest to highest abundance.
    m1_centre_element : str. 
        This is the element of the centre point atom in the first molecule.
    m1_point1_element : str. 
        This is the element of the 1st point atom in the first molecule.
    m1_point2_element : str. 
        This is the element of the 2nd point atom in the first molecule.
    m1_point3_element : str. 
        This is the element of the 3rd point atom in the first molecule.

    Returns
    -------
    m2_centre_element : str.
        This is the element of the centre atom in molecule 2.
    m2_centre_position : numpy.array
        This is the position of the centre atom in molecule 2.
    m2_centre_index : int
        This is the index of the centre atom in molecule 2.

    m2_point1_element : str.
        This is the element of point 1 in molecule 2.
    m2_point1_position : numpy.array
        This is the position of point 1 in molecule 2.
    m2_point1_index : int
        This is the index of point 1 in molecule 2.

    m2_point2_element : str.
        This is the element of point 2 in molecule 2.
    m2_point2_position : numpy.array
        This is the position of point 2 in molecule 2.
    m2_point2_index : int
        This is the index of point 2 in molecule 2.

    m2_point3_element : str.
        This is the element of point 3 in molecule 2.
    m2_point3_position : numpy.array
        This is the position of point 3 in molecule 2.
    m2_point3_index : int
        This is the index of point 3 in molecule 2.
    """

    for m2_index1 in range(len(m2_positions_of_lowest_elements)):

        # 1. Obtain the position of the centre point of molecule 2.
        m2_centre_element, (m2_centre_position, m2_centre_index) = m2_positions_of_lowest_elements[m2_index1]
        if not (m2_centre_element == m1_centre_element):
            continue

        for m2_index2 in range(len(m2_positions_of_lowest_elements)):

            # 2. Obtain the position of the 1st point of molecule 2.
            if m2_index1 == m2_index2:
                continue
            m2_point1_element, (m2_point1_position, m2_point1_index) = m2_positions_of_lowest_elements[m2_index2]
            if not (m2_point1_element == m1_point1_element):
                continue

            for m2_index3 in range(len(m2_positions_of_lowest_elements)):

                # 3. Obtain the position of the 2nd point of molecule 2.
                if m2_index3 in [m2_index1, m2_index2]:
                    continue
                m2_point2_element, (m2_point2_position, m2_point2_index) = m2_positions_of_lowest_elements[m2_index3]
                if not (m2_point2_element == m1_point2_element):
                    continue

                for m2_index4 in range(len(m2_positions_of_lowest_elements)):

                    # 4. Obtain the position of the 3rd point of molecule 2.
                    if m2_index4 in [m2_index1, m2_index2, m2_index3]:
                        continue
                    m2_point3_element, (m2_point3_position, m2_point3_index) = m2_positions_of_lowest_elements[m2_index4]
                    if not (m2_point3_element == m1_point3_element):
                        continue

                    # Have centre, 1st, 2nd, and 3rd points to look at in the second dimer.
                    yield (m2_centre_element, m2_centre_position, m2_centre_index, m2_point1_element, m2_point1_position, m2_point1_index, m2_point2_element, m2_point2_position, m2_point2_index, m2_point3_element, m2_point3_position, m2_point3_index)
