'''
Geoffrey Weal, get_electronic_coupling_data_methods.py, 19/4/22

This script contains methods for get_electronic_coupling_data.py
'''
import os
from numpy import array

def read_dimer_data_from_All_Dimer_Information(path_to_All_Dimer_Information):
    """
    This method is designed to obtain all the unit cell spatial information about the molecules in the crystals that make up all the possible dimers in the crystal.

    Parameters
    ----------
    path_to_All_Dimer_Information : str.
        This is the path to the All_Dimer_Information text file.

    Returns
    -------
    all_dimer_information : dict.
        This dictionary includes all the spatial information about the dimers in the crystal, with respect to the unit cell vectors.
    """

    # First, create the all_dimer_information dictionary.
    all_dimer_information = {}

    # Second, open the All_Dimer_Information.txt file.
    with open(path_to_All_Dimer_Information+'/All_Dimer_Information.txt','r') as All_Dimer_InformationTXT:

        # Third, pass over the first line in All_Dimer_Information.txt
        All_Dimer_InformationTXT.readline()

        # Fourth, for each line in All_Dimer_Information.txt, starting from the second line.
        for line in All_Dimer_InformationTXT:

            # Fifth, strip the line.
            line = line.rstrip().split()

            # Sixth, obtain the dimer number.
            dimer_no = int(line[0])

            # Seventh, obtain the numbers for each molecule used from the crystal.
            molecule1_no = int(line[1])
            molecule2_no = int(line[2])

            # Eighth, obtain the unit_cell_displacement_vector
            unit_cell_displacement_vector = (int(line[3]), int(line[4]), int(line[5]))

            # Ninth, obtain the displacement_vector
            displacement_vector = array([float(line[6]), float(line[7]), float(line[8])])

            # Tenth, make sure that dimer_no is not already in all_dimer_information
            if dimer_no in all_dimer_information:
                raise Exception('huh?')

            # Eleventh, write the dimer data to all_dimer_information
            all_dimer_information[dimer_no] = (molecule1_no, molecule2_no, unit_cell_displacement_vector)

    # Twelfth, return all_dimer_information
    return all_dimer_information

# ========================================================================================================

def read_dimer_data_from_Unique_Dimer_Information(path_to_Unique_Dimer_Information):
    """
    This method is designed to obtain information about which non-unique dimers are the same as the unique dimers in the crystal.

    Parameters
    ----------
    path_to_Unique_Dimer_Information : str.
        This is the path to the Unique_Dimer_Information text file.

    Returns
    -------
    symmetric_to_unique_dimer : dict.
        This dictionary includes the dimers that are the same as a unique dimer, which has had it's EET recorded in the dimers_EET_information dictionary.
    """

    # First, create the unique_dimer_information dictionary.
    symmetric_to_unique_dimer = {}

    # Second, open the Unique_Dimer_Information.txt file.
    path_to_unique_dimers = path_to_Unique_Dimer_Information+'/Unique_Dimer_Information.txt'
    if os.path.exists(path_to_unique_dimers):
        with open(path_to_unique_dimers,'r') as Unique_Dimer_InformationTXT:
            
            # Third, pass over the first line in Unique_Dimer_Information.txt
            Unique_Dimer_InformationTXT.readline()
            
            # Fourth, for each line in Unique_Dimer_Information.txt, starting from the second line.
            for line in Unique_Dimer_InformationTXT:

                # Fifth, strip the line.
                line = line.rstrip().split()

                # Sixth, obtain the dimer number for the symmetric (non-unique) dimer.
                symmetric_dimer_number = int(line[0])

                # Seventh, obtain the dimer number of the unique dimer that is the same as this symmetric (non-unique) dimer.
                unique_dimer_number = int(line[2])

                # Eighth, make sure that dimer_no is not already in symmetric_to_unique_dimer
                if symmetric_dimer_number in symmetric_to_unique_dimer:
                    raise Exception('huh?')

                # Ninth, write the dimer data to symmetric_to_unique_dimer
                symmetric_to_unique_dimer[symmetric_dimer_number] = unique_dimer_number
    
    # Tenth, return symmetric_to_unique_dimer
    return symmetric_to_unique_dimer

# ========================================================================================================

def get_EET_calculation_information(path_to_dimer_information_file, crystal_name, functional_and_basis_set):
    """
    This method is designed to obtain the EET data from the individual EET text files. 

    Parameters
    ----------
    path_to_dimer_information_file : str.
        This is the path to the text file that contains the EET energies for this crystal and for the paricular functional and basis set.
    crystal_name : str.
        This is the name of the crystal.
    functional_and_basis_set : str.
        This is the name of the functional and basis set used in calculations to obtain EET and ATC information. 

    Returns
    -------
    dimers_EET_information : dict.
        This dictionary includes all the energetic EET information about the unique dimers in the crystal.
    """

    # First, create the unique_dimer_information dictionary.
    dimers_EET_information = {}

    # Second, obtain the file name for the file that contains the EET data for this crystal, and for this functional and basis set.
    filename = str(crystal_name)+'_'+str(functional_and_basis_set)+'.txt'

    # Third, open the EET containing text file.
    with open(path_to_dimer_information_file+'/'+filename, 'r') as Dimer_CouplingTXT:

        # Fourth, obtain the energy units for the EET values being read
        energy_units = Dimer_CouplingTXT.readline().rstrip().split()[-1]

        # Fifth, obtain the energy conversion value based on the energy_units obtain previously.
        if   'meV' in energy_units:
            eV_conversion_value = 1000.0
        elif 'eV'  in energy_units:
            eV_conversion_value = 1.0
        else:
            raise Exception('Error: Energy units for EET values need to be given in: '+str(path_to_dimer_information_file+'/'+filename))

        # Sixth, for each line in Dimer_CouplingTXT after the second line
        for line in Dimer_CouplingTXT:

            # Seventh, obtain the dimer name and the EET coupling energy for the dimer
            dimer_name, coupling_energy = line.rstrip().split()

            # Eighth, get info about the molcules and dimer names
            dimer_no, mol1_name, mol2_name = dimer_name.split('_')
            dimer_no = int(dimer_no.replace('Dimer',''))
            mol1_no = int(mol1_name.replace('M',''))
            mol2_no = int(mol2_name.replace('M',''))

            # Ninth, get the EET coupling energy for the dimer.
            coupling_energy = float(coupling_energy) / eV_conversion_value

            # Tenth, set information
            dimers_EET_information[dimer_no] = (mol1_no, mol2_no, coupling_energy)

    # Eleventh, return dimers_EET_information
    return dimers_EET_information

# ========================================================================================================



