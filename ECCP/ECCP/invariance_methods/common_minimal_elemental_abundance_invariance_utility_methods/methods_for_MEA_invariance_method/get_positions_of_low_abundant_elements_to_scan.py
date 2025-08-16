"""
get_positions_of_low_abundant_elements_to_scan.py, Geoffrey Weal, 12/2/24

This method is designed to determine the positions of the elements from least to most abundant, up to max_number_of_atoms.
"""

def get_positions_of_low_abundant_elements_to_scan(elements, positions, max_number_of_atoms):
    """
    This method is designed to determine the positions of the elements from least to most abundant, up to max_number_of_atoms.

    Parameters
    ----------
    elements : list
        The elements in the dimer.
    positions : numpy.array
        The positions in the dimer.
    max_number_of_atoms : int or str. 
        This the maximum number of atoms that you want to include in this list. All positions of elements will be give in the output list. Write "all" if you want to include all the atoms of the dimer in this list.
        
    Returns
    -------
    positions_ordered_by_element_abundance : list of (element, position)
        This is a list of the elements and position of atoms in the list, from lowest element abundance to highest element abundance. 
    """

    # First, just make sure that the dimer has enough atoms for adding to the list given the value of max_number_of_atoms
    if (not max_number_of_atoms == 'all') and (len(elements) < max_number_of_atoms):
        print('Error with def get_positions_of_low_abundant_elements_to_scan in methods_for_MEA_invariance_method.py')
        print('the number of elements given is less than the max_number_of_atoms to examine')
        print('len(elements) = '+str(len(elements)))
        print('max_number_of_atoms = '+str(max_number_of_atoms))
        print('Check this')
        print('(Warning: This error message may not be needed, for further investigation. GRW 26/36/23)')
        import pdb; pdb.set_trace()
        raise Exception('check this out')

    # Second, obtain the positions for elements in the dimer.
    positions_of_elements = {}
    for index in range(len(elements)):
        element  = elements[index]
        position = positions[index]
        positions_of_elements.setdefault(element,[]).append((position,index))

    # Third, give a list of elements and positions of atoms in the dimer from lowest abundance element to highest abundance element.
    positions_ordered_by_element_abundance = []
    for element, positions in sorted(positions_of_elements.items(), key=lambda x:len(x[1])):
        for position, index in positions:
            positions_ordered_by_element_abundance.append((element, (position, index)))
        # if you have obtain enough atoms up to max_number_of_atoms, terminate.
        if (not max_number_of_atoms == 'all') and (len(positions_ordered_by_element_abundance) >= max_number_of_atoms):
            break

    # Fourth, return the positions_ordered_by_element_abundance list
    return positions_ordered_by_element_abundance
