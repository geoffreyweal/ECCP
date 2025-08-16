'''
Geoffrey Weal, get_EET_coupling_data.py, 16/3/23

This script is designed to retrieve the ATC data for each crystal. 
'''
import os
from itertools import count
from tqdm import tqdm
from SUMELF import get_cell_corner_points, obtain_graph
from ECCP.ECCP_Programs.processing_coupling_methods.ATC_methods.get_ATC_data                                import get_ATC_data
from ECCP.ECCP_Programs.processing_coupling_methods.ATC_methods.processing_ATC_data_methods                 import get_structurally_equivalent_molecules, obtain_all_molecules_with_ATC_charges_in_crystal
from ECCP.ECCP_Programs.processing_coupling_methods.ATC_methods.get_neighbours_methods.nearest_atoms_method import get_neighbours_nearest_atoms_method

def get_ATC_coupling_data(overall_path, log_filename, start_time):
    """
    This method is designed to retrieve the ATC data for each crystal. 

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
    new_ATC_coupling_data : dict.
        This is the dictionary that EET data will be gathered into. 
    """

    # First, obtain the ATC coupling data from the Gaussian output.log files. 
    ATC_coupling_data, issues = get_ATC_data(overall_path+'/'+'Unique_ATC_Gaussian_Jobs', log_filename, start_time)

    # Second, report back if there were any issues with the ATC calculations. 
    if len(issues) > 0:
        raise Exception('add message here')

    # Third, convert the ATC data given above into a format that will be used later on in this method. 
    ATC_coupling_data = convert_ATC_data(ATC_coupling_data)

    # Fourth, obtain all the names of all the crystals being analysed. 
    crystal_names = list(set([crystal_name for crystal_name in ATC_coupling_data.keys()]))

    # Fifth, initialise the new EET coupling data dictionary.
    new_ATC_coupling_data  = {}
    molecules_in_crystals  = {}
    unit_cells_of_crystals = {}

    # Sixth, obtain the data for the EET coupling data dictionary.
    for crystal_name in crystal_names:

        # 6.1: Obtain the molecules for the current crystal.
        ATC_coupling_data_for_crystal = ATC_coupling_data[crystal_name]

        # 6.2: Obtain the structurally unique molecules from ATC_coupling_data_for_crystal
        molecule_nos = list(ATC_coupling_data_for_crystal[tuple(ATC_coupling_data_for_crystal.keys())[0]].keys())

        # 6.3: Obtain which structurally equivalent molecule are the same as which unique molecule. 
        path_to_txtfile = overall_path+'/ECCP_Information/'+crystal_name
        structurally_equivalent_molecules = get_structurally_equivalent_molecules(path_to_txtfile, molecule_nos)

        # 6.4: Obtain the molecules in the crystal, and have ATC charges assigned to each atom in each molecule. Bonding graphs for these molecules will also be obtained. 
        molecules_in_crystal, molecules_in_crystal_graphs = obtain_all_molecules_with_ATC_charges_in_crystal(path_to_txtfile, ATC_coupling_data_for_crystal, structurally_equivalent_molecules)
        if not (list(sorted(molecules_in_crystal_graphs.keys())) == list(range(1,len(molecules_in_crystal_graphs)+1))):
            raise Exception('huh?')
        molecules_in_crystal_graphs = [molecule for mol_no, molecule in sorted(molecules_in_crystal_graphs.items())]

        # 6.5: Initialise the new_ATC_coupling_data_for_crystal dictionary
        new_ATC_coupling_data_for_crystal = {}

        # 6.6: Obtain the coupling between molecules in the crystal
        for functional_basis_set, molecules_in_crystal_for_fbs in molecules_in_crystal.items():

            # 6.6.1: Make a copy of the molecules in this crystal.
            if not crystal_name in molecules_in_crystals:
                molecules_in_crystal_for_fbs_copy = {}
                for no, mol in sorted(molecules_in_crystal_for_fbs.items()):
                    mol_copy = mol.copy()
                    for atom in mol_copy:
                        atom.charge = 0.0
                    molecules_in_crystal_for_fbs_copy[no] = mol_copy
                molecules_in_crystals[crystal_name] = molecules_in_crystal_for_fbs_copy

            # 6.6.1: Change the format of the molecules_in_crystal_for_fbs into list mode
            molecules_in_crystal_for_fbs_list  = [mol for no, mol in sorted(molecules_in_crystal_for_fbs.items())]

            # 6.6.2: Add the unit cell information for this crystal to unit_cells_of_crystals
            if not (crystal_name in unit_cells_of_crystals):
                unit_cells_of_crystals[crystal_name] = molecules_in_crystal_for_fbs_list[0].get_cell()

            # 6.6.3: Get the ATC coupling values between molecules in the crystal for a certain functional and basis set.
            coupling_values_of_crystal_for_fbs = get_neighbours_nearest_atoms_method(molecules_in_crystal_for_fbs_list, molecules_in_crystal_graphs, no_of_cpus=1)

            # 6.6.4: Record the ATC coupling values between molecules in the crystal for a certain functional and basis set.
            if (functional_basis_set in new_ATC_coupling_data_for_crystal):
                raise Exception('Huh?')
            new_ATC_coupling_data_for_crystal[functional_basis_set] = coupling_values_of_crystal_for_fbs

        # 6.7: Record the ATC coupling values between molecules in the crystal for each functional and basis set.
        if (crystal_name in new_ATC_coupling_data):
            raise Exception('Huh?') 
        new_ATC_coupling_data[crystal_name] = new_ATC_coupling_data_for_crystal

    # Seventh, return the electronic_coupling_data dictionary.
    return new_ATC_coupling_data, molecules_in_crystals, unit_cells_of_crystals

# ====================================================================================================================================

def convert_ATC_data(ATC_coupling_data):
    """
    This method will convert the ATC_coupling_data dictionary into a form to be used in this ECCP program.

    Parameters
    ----------
    ATC_coupling_data : dict.
        This is the dictionary that EET data is gathered into. 

    Returns
    -------
    new_ATC_coupling_data : dict.
        This is the converted dictionary that EET data will be gathered into. 
    """

    # First, initialise the new_ATC_coupling_data dictionary.
    new_ATC_coupling_data = {}

    # Second, convert the EET_coupling_data dictionary into a new format. 
    for (crystal_name, molecule_name, functional_and_basis_set), (root, molecule) in ATC_coupling_data.items():

        # 2.1: First, obtain information about the dimer from it's name.
        molecule_no = int(molecule_name.replace('molecule_',''))

        # 2.3: Initialise sub-dictionaries in new_ATC_coupling_data if they do not exist yet.
        if not crystal_name in new_ATC_coupling_data:
            new_ATC_coupling_data[crystal_name] = {}
        if not functional_and_basis_set in new_ATC_coupling_data[crystal_name]:
            new_ATC_coupling_data[crystal_name][functional_and_basis_set] = {}
        if not molecule_no in new_ATC_coupling_data[crystal_name][functional_and_basis_set]:
            new_ATC_coupling_data[crystal_name][functional_and_basis_set][molecule_no] = {}

        # 2.4: Add EET information to new_ATC_coupling_data
        new_ATC_coupling_data[crystal_name][functional_and_basis_set][molecule_no] = molecule
    
    # Third, return new_EET_coupling_data
    return new_ATC_coupling_data

# ====================================================================================================================================


"""


            # 6.6.1: Convert the molecules_in_crystal_for_fbs dictionary into a sorted list
            molecules_in_crystal_for_fbs = sorted(molecules_in_crystal_for_fbs.items())

            # 6.6.2: Initialise the ATC_coupling_in_crystal_for_fbs dictionary
            ATC_coupling_in_crystal_for_fbs = {}

            # 6.6.3: prepare progress bar for ATC calculations. 
            print('Analysing ATC for '+str(crystal_name)+': '+str(functional_basis_set))
            pbar = tqdm(count(start=1))

            # 6.6.4: For each maximum ijk size
            for super_cell_reach in pbar:

                # 6.6.5: Get the unit cell positions to examine dimers for
                cell_corner_displacements, cell_corner_ijk_units = get_cell_corner_points(unit_cell, super_cell_reach=super_cell_reach, get_corrspeonding_ijk_values=True)

                # 6.6.6: Initialise the highest_coupling_value float.
                highest_coupling_value = 0.0 # eV

                for cell_corner_displacement, cell_corner_ijk_unit in zip(cell_corner_displacements, cell_corner_ijk_units):

                    if not (super_cell_reach in cell_corner_ijk_unit):
                        continue

                    for molecule1_no, molecule1 in molecules_in_crystal_for_fbs:

                        for molecule2_no, molecule2 in molecules_in_crystal_for_fbs:

                            if (molecule1_no == molecule2_no) and (cell_corner_ijk_unit == (0,0,0)):
                                continue

                            dimer_information = (molecule1_no-1, molecule2_no-1, cell_corner_ijk_unit)
                            pbar.set_description(str(dimer_information))

                            coupling_value = get_coulomb_energy(molecule1, molecule2, cell_corner_displacement, relative_permittivity=1.0)

                            ATC_coupling_in_crystal_for_fbs[dimer_information] = coupling_value

                            if abs(coupling_value) > highest_coupling_value:
                                highest_coupling_value = abs(coupling_value)


                if (super_cell_reach > 3) and (highest_coupling_value < 0.1):
                    break

            import pdb; pdb.set_trace()




"""





