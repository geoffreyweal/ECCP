'''
Geoffrey Weal, ECCP_processing_ICT_data.py, 7/6/22

This program is designed to process the ICT data from Gaussian output files

'''
import os, time
from datetime import datetime, timedelta

from SUMELF import remove_folder, make_folder

from ECCP.ECCP_Programs.processing_Eigendata_methods.process_Eigendata_to_disk import process_Eigendata_to_disk
from ECCP.ECCP_Programs.processing_ICT_methods.processing_ICT_data_methods import get_matrix_from_file, get_MO_orbital_names, get_MO_occupancies, assign_MO_coefficients_with_atoms
from ECCP.ECCP_Programs.processing_ICT_methods.processing_matrix_data import processing_matrix_data
from ECCP.ECCP_Programs.processing_ICT_methods.write_data_to_excel import write_data_to_excel

# ---------------------------------------------------------------------

class CLICommand:
    """Will process ICT Data from your Eigendata into text files and an excel file. 

    Will also process and extract the eigendata from your output.log files and place them in individual txt files.
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
    This method is the main method for running this program. 
    """

    # First, get the general variables for processing data.
    overall_path = os.getcwd()
    log_filename = 'output.log'
    ict_data_foldername = 'ICT_Data'
    individual_ict_data_foldername = 'Individual_ICT_Data'
    for foldername in [ict_data_foldername, individual_ict_data_foldername]:
        remove_folder(foldername)
        make_folder(foldername)

    # Second, start the process and start timing
    print('------------------------------------------------')
    dt_string = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    print('Starting the process_ICT program at ('+str(dt_string)+')')
    start_time = time.time()
    print('------------------------------------------------')

    # Third, obtain the eigendata from the Gaussian output.log files. 
    print('Processing and Extracting Eigen-data from output.log files')
    print('Note: This program may take some time if you have recorded eigendata, such as orbital overlap matrices, as these matrices can be very large depending on the number of atoms in your dimer.')
    eigendata = {}
    issues = []
    for root, dirs, files in os.walk(overall_path):
        dirs.sort()

        # Fourth, check to see if you have these files in your folder. 
        # If you do, then you are in the right place for obtaining eigendata on your dimer and molecules.
        if not (('Dimer' in dirs) and ('Monomer_1' in dirs) and ('Monomer_2' in dirs)):
            continue

        # Fifth, indicate that a Gaussian job has been found
        print('------------------------------------------------')
        print(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+' - Found a Gaussian job: '+str(root))
        dirs[:]  = []
        files[:] = []

        # Sixth, process the eigendata.
        issue = process_Eigendata_to_disk(root, log_filename, start_time)
        if issue is not None:
            issues.append(issue)
            print(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+' - Their was an issue with this job. Will skip processing it further.')
            print('------------------------------------------------')
            continue

        print(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+' - Performing Calculations using Eigendata.')
        dirs[:]  = []
        files[:] = []

        # Seventh, get the names of crystal, dimers, monomers, and functional_and_basis_set
        crystal_name = root.split('/')[-3]
        dimer_name, mon1_name, mon2_name = root.split('/')[-2].split('_')
        dimer_name = '_'.join([dimer_name, mon1_name, mon2_name])
        functional_and_basis_set_name = root.split('/')[-1]

        # Eighth, get the matrix data for the monomers.
        monomer_paths = [('Monomer_1', mon1_name), ('Monomer_2', mon2_name)]
        monomer_data = []
        for monomer_index in range(len(monomer_paths)):
            monomer_foldername, mon_name = monomer_paths[monomer_index]

            # 8.1: Get the path to where matrices are stored
            path_to_matrices = root+'/'+monomer_foldername

            # 8.2: Get the MO_coefficients and break it apart and assign each row to it's corresponding atom
            MO_coefficients        = get_matrix_from_file(path_to_matrices+'/MO_coefficients.txt')
            MO_orbital_names       = get_MO_orbital_names(path_to_matrices+'/MO_orbital_names.txt')
            MO_coefficients_data   = assign_MO_coefficients_with_atoms(MO_coefficients, MO_orbital_names)

            # 8.3: Obtain the occupancies of the mononer orbtials. This will help to determine the HOMO and LUMO for the monomers.
            MO_occupancies         = get_MO_occupancies(path_to_matrices+'/MO_occupancies.txt')

            # 8.4: Determine the index of the HOMO and LUMO
            HOMO_index = MO_occupancies.index('V')-1
            LUMO_index = MO_occupancies.index('V')

            # 8.5: Write the MO data to memory
            monomer_data.append((MO_coefficients_data, HOMO_index, LUMO_index)) #(orbital_overlap_matrix, MO_energies, MO_occupancies, all_MO_coefficients_data)
            
        # ------------------------------------------------------------
        # Ninth, get the matrix data for the dimer.

        # 9.1: Get the path to where matrices are stored
        path_to_matrices = root+'/'+'Dimer'

        # 9.2: Get the orbtial overlap matrix for the dimer
        dimer_orbital_overlap_matrix   = get_matrix_from_file(path_to_matrices+'/orbital_overlap_matrix.txt')

        # 9.3: Get the MO energies for the dimer, and convert it into a diagonal matrix
        dimer_MO_energies              = get_matrix_from_file(path_to_matrices+'/MO_energies.txt')

        # 9.4: Get the MO_coefficients and break it apart and assign each row to it's corresponding atom
        dimer_MO_coefficients_matrix   = get_matrix_from_file(path_to_matrices+'/MO_coefficients.txt')
        dimer_MO_orbital_names         = get_MO_orbital_names(path_to_matrices+'/MO_orbital_names.txt')

        # ------------------------------------------------------------
        # Tenth, construct the HOMO and LUMO for monomer 1 and monomer 2 in the correct atom order as in the dimer
        hole_transfer, electron_charge_transfer = processing_matrix_data(monomer_data, dimer_orbital_overlap_matrix, dimer_MO_energies, dimer_MO_coefficients_matrix, dimer_MO_orbital_names)

        # Eleventh, store the data. 
        eigendata[(crystal_name, dimer_name, functional_and_basis_set_name)] = (root, [hole_transfer, electron_charge_transfer])

        # Twelfth, indicate processing on this dimer has finished.
        print(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+' - ICT Calculations were successfully performed upon '+str(root))
        print('Current program running time (HH:MM:SS): '+str(timedelta(seconds=time.time() - start_time)))
        print('------------------------------------------------')

    # Thirteenth, write the eigendata to an excel file.
    write_data_to_excel(eigendata, ict_data_foldername, individual_ict_data_foldername, start_time)

    # Fourteenth, write any issues to the terminal.
    print('------------------------------------------------')
    if len(issues) > 0:
        print('The following Gaussian jobs could not be processed because they have not finished running or did not complete successfully.')
        for issue in issues:
            print(issue)
        print('------------------------------------------------')
    print(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+' - This ICT calculations program has finished successfully!')
    print('Total running time (HH:MM:SS): '+str(timedelta(seconds=time.time() - start_time)))
    print('------------------------------------------------')

# ------------------------------------------------------------------------------------------------------------------------------------------------------------







