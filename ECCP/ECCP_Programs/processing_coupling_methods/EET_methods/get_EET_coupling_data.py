'''
Geoffrey Weal, get_EET_coupling_data.py, 16/3/23

This script is designed to retrieve the EET data for each crystal. 
'''
import os
from ECCP.ECCP_Programs.processing_EET_methods.get_EET_data import get_EET_data
from ECCP.ECCP_Programs.processing_coupling_methods.EET_methods.get_coupling_data_methods import read_dimer_data_from_All_Dimer_Information, read_dimer_data_from_Unique_Dimer_Information

def get_EET_coupling_data(overall_path, log_filename, start_time):
    """
    This method is designed to retrieve the EET data for each crystal. 

    Parameters
    ----------
    overall_path : str.
        This is the path to the current place the program was executed from. 
    log_filename : str.
        This is the name of the output.log file for EET calculations. 
    start_time : float
        This is the start time for the process.

    Returns
    -------
    new_EET_coupling_data : dict.
        This is the dictionary that EET data will be gathered into. 
    """

    # First, obtain the EET coupling data from the Gaussian output.log files. 
    EET_coupling_data, issues = get_EET_data(overall_path+'/'+'Unique_EET_Gaussian_Jobs', log_filename, start_time)

    # Second, report back if there were any issues with the EET calculations. 
    if len(issues) > 0:
        raise Exception('add message here')

    # Third, convert the EET data given above into a format that will be used later on in this method. 
    EET_coupling_data = convert_EET_data(EET_coupling_data)

    # Fourth, obtain all the names of all the crystals being analysed. 
    crystal_names = list(set([crystal_name for crystal_name in EET_coupling_data.keys()]))

    # Fifth, initialise the new EET coupling data dictionary.
    new_EET_coupling_data = {}

    # Sixth, obtain the data for the EET coupling data dictionary.
    for crystal_name in crystal_names:

        # 6.1: Get the spatial information about the molecules that make up the dimer.
        all_dimer_information     = read_dimer_data_from_All_Dimer_Information(overall_path+'/All_Dimer_Information/'+crystal_name)

        # 6.2: Determine which symmetric dimers are associated with which unique dimer that EET calculations were performed upon.
        symmetric_to_unique_dimer = read_dimer_data_from_Unique_Dimer_Information(overall_path+'/All_Dimer_Information/'+crystal_name)

        # 6.3: Obtain the EET coupling data for the crystal of interest. 
        EET_coupling_data_for_crystal = EET_coupling_data[crystal_name]

        # 6.4: Initialise the new EET coupling data dictionary for each functional and basis set. 
        functional_and_basis_set_data = {}

        # 6.5: For each functional and basis set.
        for functional_and_basis_set in EET_coupling_data_for_crystal.keys():

            # 6.5.1: Gather all spatial and EET calculation information together
            EET_coupling_data_for_crystal_and_fb  = combine_dict_information(all_dimer_information, EET_coupling_data_for_crystal[functional_and_basis_set], symmetric_to_unique_dimer)

            # 6.5.2: Record the EET data for this crystal and functional+basis set to functional_and_basis_set_data
            functional_and_basis_set_data[functional_and_basis_set] = EET_coupling_data_for_crystal_and_fb

        # 6.6: Record the EET data for this crystal.
        new_EET_coupling_data[crystal_name] = functional_and_basis_set_data

    # Seventh, return the electronic_coupling_data dictionary.
    return new_EET_coupling_data

# ====================================================================================================================================

def convert_EET_data(EET_coupling_data):
    """
    This method will convert the EET_coupling_data dictionary into a form to be used in this ECCP program.

    Parameters
    ----------
    EET_coupling_data : dict.
        This is the dictionary that EET data is gathered into. 

    Returns
    -------
    new_EET_coupling_data : dict.
        This is the converted dictionary that EET data will be gathered into. 
    """

    # First, initialise the new_EET_coupling_data dictionary.
    new_EET_coupling_data = {}

    # Second, convert the EET_coupling_data dictionary into a new format. 
    for (crystal_name, dimer_name, functional_and_basis_set), (root, EET_data) in EET_coupling_data.items():

        # 2.1: First, obtain information about the dimer from it's name.
        dimer_no, molecule_1, molecule_2 = dimer_name.split('_')
        dimer_no     = int(dimer_no.replace('Dimer',''))
        molecule1_no = int(molecule_1.replace('M',''))
        molecule2_no = int(molecule_2.replace('M','')) 
        
        # 2.2: Obtain the EET coupling value
        EET_coupling_value_eV = float(EET_data[-1])

        # 2.3: Initialise sub-dictionaries in new_EET_coupling_data if they do not exist yet.
        if not crystal_name in new_EET_coupling_data:
            new_EET_coupling_data[crystal_name] = {}
        if not functional_and_basis_set in new_EET_coupling_data[crystal_name]:
            new_EET_coupling_data[crystal_name][functional_and_basis_set] = {}
        if not dimer_no in new_EET_coupling_data[crystal_name][functional_and_basis_set]:
            new_EET_coupling_data[crystal_name][functional_and_basis_set][dimer_no] = {}

        # 2.4: Add EET information to new_EET_coupling_data
        new_EET_coupling_data[crystal_name][functional_and_basis_set][dimer_no] = (molecule1_no, molecule2_no, EET_coupling_value_eV)
    
    # Third, return new_EET_coupling_data
    return new_EET_coupling_data

def combine_dict_information(all_dimer_information, dimers_EET_information, symmetric_to_unique_dimer):
    """
    This method will take the all_dimer_information, dimers_EET_information, and symmetric_to_unique_dimer dictionaries, and combine them together to give all the electronic information about the dimers in the crystal.
    
    Parameters
    ----------
    all_dimer_information : dict.
        This dictionary includes all the spatial information about the dimers in the crystal, with respect to the unit cell vectors.
    dimers_EET_information : dict.    
        This dictionary includes all the energetic EET information about the unique dimers in the crystal.
    symmetric_to_unique_dimer : dict.
        This dictionary includes the dimers that are the same as a unique dimer, which has had it's EET recorded in the dimers_EET_information dictionary.

    Returns
    ------- 
    electronic_coupling_data : dict.
        This is all the energetic energy for each dimer involving each molecule from the unit cell.
    """

    # First, setup the electronic_coupling_data dictionary.
    electronic_coupling_data = {}

    # Second, for each dimer in the all_dimer_information dictionary.
    for dimer_no in all_dimer_information:

        # Third, obtain the molecules and the displacement vector that make up each dimer in all_dimer_information
        molecule1_no, molecule2_no, displacement_vector = all_dimer_information[dimer_no]

        # Fourth, obtain the dimer that contains the EET energy, either
        #   * For the unique dimer itself.
        #   * For non-unique dimer, where it's EET value is the same as it unique dimer counterpart.
        if dimer_no in symmetric_to_unique_dimer.keys():
            dimer_no_for_EET = symmetric_to_unique_dimer[dimer_no]
        else:
            dimer_no_for_EET = dimer_no

        # Fifth, obtain the EET coupling energy for the dimer of interest
        mol1_no, mol2_no, coupling_energy = dimers_EET_information[dimer_no_for_EET]

        # Sixth, make sure that the molecules involved in the non-unique dimer was also the same molecules in the unique dimer.
        if (dimer_no_for_EET == dimer_no) and (not ((molecule1_no == mol1_no) and (molecule2_no == mol2_no))):
            raise Exception('Huh?')

        # Seventh, record the electronic data into the electronic_coupling_data, where the key describes the relative displacement between each molecule in the system.
        electronic_coupling_data[(molecule1_no,molecule2_no,displacement_vector[0],displacement_vector[1],displacement_vector[2])] = coupling_energy #(molecule1_no, molecule2_no, displacement_vector, total_electronic_coupling_energy)

    # Eighth, return the electronic_coupling_data dictionary.
    return electronic_coupling_data

# ====================================================================================================================================







