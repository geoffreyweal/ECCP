'''
Geoffrey Weal, ECCP_processing_coupling_data.py, 16/3/23

This program is designed to allow the user to identify the coupling between molecules in the crystal with the various techniques.

'''
import os, time
import numpy as np
from datetime import datetime, timedelta
from ase import Atoms

from SUMELF import remove_folder, make_folder, get_centre_of_mass
from ECCP.ECCP_Programs.processing_coupling_methods.EET_methods.get_EET_coupling_data import get_EET_coupling_data
from ECCP.ECCP_Programs.processing_coupling_methods.ATC_methods.get_ATC_coupling_data import get_ATC_coupling_data

# ---------------------------------------------------------------------

class CLICommand:
    """Will process RE Data into text files and an excel file.
    """

    @staticmethod
    def add_arguments(parser):
        pass

    @staticmethod
    def run(args):
        Run_method()

# ---------------------------------------------------------------------

def Run_method():
    """
    This method is the main method for running this program
    """

    # General variables for processing data.
    overall_path = os.getcwd()
    log_filename = 'output.log'
    coupling_data_foldername = 'Coupling_Data/Coupling_Data'
    individual_coupling_data_foldername = 'Coupling_Data/Individual_Coupling_Data'
    for foldername in [coupling_data_foldername, individual_coupling_data_foldername]:
        remove_folder(foldername)
        make_folder(foldername)

    # Second, print starting lines
    print('------------------------------------------------')
    dt_string = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    print('Starting the process_coupling program at ('+str(dt_string)+')')
    start_time = time.time()
    print('------------------------------------------------')

    # Third, obtain EET coupling data for dimers in the crystals of interest. 
    EET_coupling_data                                                = get_EET_coupling_data(overall_path, log_filename, start_time)

    ATC_coupling_data, molecules_in_crystals, unit_cells_of_crystals = get_ATC_coupling_data(overall_path, log_filename, start_time)



    for crystal_name, unit_cell_info in unit_cells_of_crystals.items():

        EET_coupling_data_for_crystal = EET_coupling_data[crystal_name]

        ATC_coupling_data_for_crystal = ATC_coupling_data[crystal_name]

        for functional_and_basis_set in EET_coupling_data_for_crystal.keys():

            EET_coupling_data_for_crystal_fbs = EET_coupling_data_for_crystal[functional_and_basis_set]

            ATC_coupling_data_for_crystal_fbs = ATC_coupling_data_for_crystal[functional_and_basis_set]

            ATC_coupling_data_for_crystal_fbs_dict = {}
            for index1, index2, unit_cell_displacement, displacement, coupling_value in ATC_coupling_data_for_crystal_fbs:
                mol1_no = index1 + 1
                mol2_no = index2 + 1
                ii, jj, kk = unit_cell_displacement
                entry = (mol1_no, mol2_no, ii, jj, kk)
                if not entry in ATC_coupling_data_for_crystal_fbs_dict.keys():
                    ATC_coupling_data_for_crystal_fbs_dict[entry] = coupling_value
            ATC_coupling_data_for_crystal_fbs = ATC_coupling_data_for_crystal_fbs_dict



            all_configs = list(EET_coupling_data_for_crystal_fbs.keys())
            for entry in ATC_coupling_data_for_crystal_fbs.keys():
                if not entry in all_configs:
                    all_configs.append(entry)
            all_configs.sort()

            molecules = molecules_in_crystals[crystal_name]

            for mol_no, molecule in molecules.items():

                coupling_around_molecule = molecule.copy()

                #import pdb; pdb.set_trace()

                coupling_around_molecule_EET_coupling        = [0.0]*len(coupling_around_molecule)
                coupling_around_molecule_ATC_coupling        = [0.0]*len(coupling_around_molecule)
                coupling_around_molecule_EET_coupling_colour = [0.0]*len(coupling_around_molecule)
                coupling_around_molecule_ATC_coupling_colour = [0.0]*len(coupling_around_molecule)

                configs_including_mol_no_as_first_index = [config for config in all_configs if (config[0] == mol_no)]

                for mol1_no, mol2_no, ii, jj, kk in configs_including_mol_no_as_first_index:

                    unit_cell_ijk = (ii, jj, kk)

                    coupled_molecule = molecules[mol2_no]

                    displacement = np.array(unit_cell_ijk) @ unit_cell_info

                    centre_of_mass = get_centre_of_mass(coupled_molecule.get_chemical_symbols(), coupled_molecule.get_positions()) + displacement

                    coupled_molecule = Atoms(symbols='O', positions=[centre_of_mass])

                    EET_coupling_value = EET_coupling_data_for_crystal_fbs.get((mol1_no, mol2_no, ii, jj, kk), 0.0) * 1000.0
                    ATC_coupling_value = ATC_coupling_data_for_crystal_fbs.get((mol1_no, mol2_no, ii, jj, kk), 0.0) * 1000.0 # weird?

                    #import pdb; pdb.set_trace()

                    coupling_around_molecule_EET_coupling        += [EET_coupling_value]     *len(coupled_molecule)
                    coupling_around_molecule_ATC_coupling        += [ATC_coupling_value]     *len(coupled_molecule)
                    coupling_around_molecule_EET_coupling_colour += [abs(EET_coupling_value)]*len(coupled_molecule)
                    coupling_around_molecule_ATC_coupling_colour += [abs(ATC_coupling_value)]*len(coupled_molecule)

                    coupling_around_molecule += coupled_molecule


                coupling_around_molecule.set_array('EET_coupling', np.array(coupling_around_molecule_EET_coupling))
                coupling_around_molecule.set_array('ATC_coupling', np.array(coupling_around_molecule_ATC_coupling))
                coupling_around_molecule.set_array('EET_coupling_colour', np.array(coupling_around_molecule_EET_coupling_colour))
                coupling_around_molecule.set_array('ATC_coupling_colour', np.array(coupling_around_molecule_ATC_coupling_colour))


                from ase.visualize import view
                from ase.io import write
                view(coupling_around_molecule)
                write('test.xyz', coupling_around_molecule)
                import pdb; pdb.set_trace()
                exit()
                    
















