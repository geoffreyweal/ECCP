'''
Geoffrey Weal, processing_matrix_data.py, 10/6/22

This script is designed to process the matrix eigendata and obtain data from it, just as hole and electron transfer energies. 

'''
import numpy as np

from ECCP.ECCP_Programs.processing_ICT_methods.processing_ICT_data_methods import create_MO_coefficients_matrix

def processing_matrix_data(monomer_data, dimer_orbital_overlap_matrix, dimer_MO_energies, dimer_MO_coefficients_matrix, dimer_MO_orbital_names):
    """
    This method is designed to rocess the matrix eigendata and obtain data from it, just as hole and electron transfer energies. 

    Parameters
    ----------
    monomer_data : list
        This is a list that contains the molecular orbitals (MO) coefficients, HOMO index, and LUMO index for a monomer.
    dimer_orbital_overlap_matrix : numpy.array
        This is the overlap matrix of the dimer.
    dimer_MO_energies : numpy.array
        This is a list of the energies of each MO in the dimer
    dimer_MO_coefficients_matrix
        These are coefficients of the MOs in the dimer
    dimer_MO_orbital_names
        These are the names of the MO orbtials in the dimer.
    
    Returns
    -------
    hole_transfer : float
        This is the energy require for a hole to move from one monomer HOMO to another monomer HOMO.
    electron_charge_transfer : float
        This is the energy require for an electron to move from one monomer LUMO to another monomer LUMO.
    """

    # ------------------------------------------------------------------------------------------

    # First, get the molecular coefficients (MO) of the monomers

    # 1.1: Get the MO_coefficients_data and the MO indices of the HOMO and LUMO for monomer 1.
    mol1_MO_coefficients_data, mol1_HOMO_index, mol1_LUMO_index = monomer_data[0]
    mol1_MO_coefficients_matrix = create_MO_coefficients_matrix(mol1_MO_coefficients_data)

    # 1.2: Get the MO_coefficients_data and the MO indices of the HOMO and LUMO for monomer 2.
    mol2_MO_coefficients_data, mol2_HOMO_index, mol2_LUMO_index = monomer_data[1]
    mol2_MO_coefficients_matrix = create_MO_coefficients_matrix(mol2_MO_coefficients_data)

    # ------------------------------------------------------------------------------------------

    # Second, get the coefficients of the HOMO and LUMO in each monomer in th edimer.

    # 2.1.1: Obtain the MO coefficients for the HOMO of monomer 1. 
    mol1_HOMO_coefficients = mol1_MO_coefficients_matrix[:,mol1_HOMO_index] # Get all the row for column mol1_HOMO_index. This makes a column matrix

    # 2.1.2: Obtain the MO coefficients for the LUMO of monomer 1. 
    mol1_LUMO_coefficients = mol1_MO_coefficients_matrix[:,mol1_LUMO_index] # Get all the row for column mol1_LUMO_index. This makes a column matrix

    # 2.2.1: Obtain the MO coefficients for the HOMO of monomer 2. 
    mol2_HOMO_coefficients = mol2_MO_coefficients_matrix[:,mol2_HOMO_index] # Get all the row for column mol2_HOMO_index. This makes a column matrix

    # 2.2.2: Obtain the MO coefficients for the LUMO of monomer 2. 
    mol2_LUMO_coefficients = mol2_MO_coefficients_matrix[:,mol2_LUMO_index] # Get all the row for column mol2_LUMO_index. This makes a column matrix

    # ------------------------------------------------------------------------------------------

    # Third, extend the lengths of the monomer MO column matrices so they are the same length as the dimer matrices.

    # 3.1.1: Obtain the lengths of the MO coefficients of the HOMO and LUMO of monomer 1 and 2
    original_length_of_mol1_HOMO_coefficients = len(mol1_HOMO_coefficients)
    original_length_of_mol1_LUMO_coefficients = len(mol1_LUMO_coefficients)
    original_length_of_mol2_HOMO_coefficients = len(mol2_HOMO_coefficients)
    original_length_of_mol2_LUMO_coefficients = len(mol2_LUMO_coefficients)

    # 3.2.1 : Add extra zeros to the HOMO of monomer 1 so it has the same length as the dimer (which contains both monomer 1 and 2).
    mol1_HOMO_coefficients = np.matrix(np.hstack([mol1_HOMO_coefficients, np.zeros(original_length_of_mol2_HOMO_coefficients)])).T

    # 3.2.2 : Add extra zeros to the LUMO of monomer 1 so it has the same length as the dimer (which contains both monomer 1 and 2).
    mol1_LUMO_coefficients = np.matrix(np.hstack([mol1_LUMO_coefficients, np.zeros(original_length_of_mol2_LUMO_coefficients)])).T

    # 3.3.1 : Add extra zeros to the HOMO of monomer 2 so it has the same length as the dimer (which contains both monomer 1 and 2).
    mol2_HOMO_coefficients = np.matrix(np.hstack([np.zeros(original_length_of_mol1_HOMO_coefficients), mol2_HOMO_coefficients])).T

    # 3.3.2 : Add extra zeros to the LUMO of monomer 2 so it has the same length as the dimer (which contains both monomer 1 and 2).
    mol2_LUMO_coefficients = np.matrix(np.hstack([np.zeros(original_length_of_mol1_LUMO_coefficients), mol2_LUMO_coefficients])).T

    # ------------------------------------------------------------------------------

    # Fourth, Perform the matrix calculations to obtain hole and electron transportation coupling energy.
    # See https://pubs.rsc.org/en/content/articlepdf/2010/cp/c002337j for more information, Eq 16, 17, and 18.

    # 4.1, create a diagonal matrix from the dimer MO energies.
    dimer_diagonal_energy_matrix = np.diag(dimer_MO_energies.T.tolist()[0])

    # 4.2: Perform the matrix calculations to obtain hole transportation coupling energy.
    mol1_HOMO_proj_on_dimer = mol1_HOMO_coefficients.H @ dimer_orbital_overlap_matrix @ dimer_MO_coefficients_matrix
    mol2_HOMO_proj_on_dimer = mol2_HOMO_coefficients.H @ dimer_orbital_overlap_matrix @ dimer_MO_coefficients_matrix
    hole_transfer = float(mol1_HOMO_proj_on_dimer @ dimer_diagonal_energy_matrix @ mol2_HOMO_proj_on_dimer.T) # energy in eV
    hole_transfer *= 1000.0 # convert this energy value to meV

    # 4.3: Perform the matrix calculations to obtain electron transportation coupling energy.
    mol1_LUMO_proj_on_dimer = mol1_LUMO_coefficients.H @ dimer_orbital_overlap_matrix @ dimer_MO_coefficients_matrix
    mol2_LUMO_proj_on_dimer = mol2_LUMO_coefficients.H @ dimer_orbital_overlap_matrix @ dimer_MO_coefficients_matrix
    electron_charge_transfer = float(mol1_LUMO_proj_on_dimer @ dimer_diagonal_energy_matrix @ mol2_LUMO_proj_on_dimer.T) # energy in eV
    electron_charge_transfer *= 1000.0 # convert this energy value to meV

    # ------------------------------------------------------------------------------

    # Fifth, return the hole and electron tranfer energies in meV

    return hole_transfer, electron_charge_transfer







