'''
Geoffrey Weal, processing_OPV_Dimer_data.py, 7/3/22

This program is designed to process the data from the Gaussian TD DFT calculations and put them into formats that the user can use for further work. 

'''
import os, time
import numpy as np
from datetime import datetime, timedelta

from SUMELF import remove_folder, make_folder
from ECCP.ECCP_Programs.shared_general_methods.shared_gaussian_methods   import found_a_gaussian_job_that_has_run
from ECCP.ECCP_Programs.shared_general_methods.shared_orca_methods       import found_an_orca_job_that_has_run
from ECCP.ECCP_Programs.processing_RE_methods.processing_RE_data_methods import found_a_re_jobset
from ECCP.ECCP_Programs.processing_RE_methods.obtain_gaussian_RE_data    import obtain_gaussian_RE_data
from ECCP.ECCP_Programs.processing_RE_methods.obtain_orca_RE_data        import obtain_orca_RE_data
from ECCP.ECCP_Programs.processing_RE_methods.write_data_to_excel        import write_data_to_excel

# ---------------------------------------------------------------------

class CLICommand:
    """Will process RE Data into text files and an excel file.
    """

    @staticmethod
    def add_arguments(parser):
        parser.add_argument('--accept_slightly_negative_frequency', nargs='?', help='Will print values with slightly negative frequencies (-100 units). If input value is a number, this is the lowest negative frequency that will be acceptioned. True: -100, (Default: False).')
        parser.add_argument('--analyse_frequencies', nargs='?', help='This indicates if you want to check the frequency calculations to make sure geometric optimisation calculations found a local minimum. (Default: True)')

    @staticmethod
    def run(args):

        # First, determine if the user want to accept slightly negative frequencies
        accept_slightly_negative_frequency = args.accept_slightly_negative_frequency
        if   accept_slightly_negative_frequency is None:
            lower_limit_negative_frequency = 0.0
        elif accept_slightly_negative_frequency.lower() in ['true', 't']:
            lower_limit_negative_frequency = -100.0
        elif accept_slightly_negative_frequency.lower() in ['false', 'f']:
            lower_limit_negative_frequency = 0.0
        elif is_a_number(accept_slightly_negative_frequency):
            lower_limit_negative_frequency = -abs(float(accept_slightly_negative_frequency))
        else:
            raise Exception('Error: accept_slightly_negative_frequency must be either True, False, or a number')

        # Second, determine if the user want to analyse frequencies to determine if the geoemtric optimisations found a local minimum or not
        analyse_frequencies = args.analyse_frequencies
        if analyse_frequencies is None:
            analyse_frequencies = True
        elif analyse_frequencies.lower() in ['true', 't']:
            analyse_frequencies = True
        elif analyse_frequencies.lower() in ['false', 'f']:
            analyse_frequencies = False
        else:
            raise Exception('Error: analyse_frequencies must be either True or False')

        # Third, run method
        Run_method(lower_limit_negative_frequency, analyse_frequencies)

def is_a_number(value):
    try:
        float(value)
    except Exception as exception:
        return False
    return True
# ---------------------------------------------------------------------

ground_structure_foldername  = 'ground_structure'
excited_structure_foldername = 'excited_structure'

def Run_method(lower_limit_negative_frequency, analyse_frequencies):
    """
    This method is the main method for running this program
    """
    # First, set general variables for processing data.
    overall_path = os.getcwd()
    re_data_foldername = 'RE_Data'
    individual_re_data_foldername = 'Individual_RE_Data'
    for foldername in [re_data_foldername, individual_re_data_foldername]:
        remove_folder(foldername)
        make_folder(foldername)

    # Second, print starting message
    print('------------------------------------------------')
    dt_string = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    print('Starting the process_RE program at ('+str(dt_string)+')')
    start_time = time.time()
    print('------------------------------------------------')

    # Third, obtain the reorganisation energy data from the Gaussian output files. 
    print('Gathering Reorganisation Energy Data')
    reorganisation_energy_data = {}
    issues = []
    for root, dirs, files in os.walk(overall_path):
        dirs.sort()

        # 3.1: Determine if their is a Gaussian/ORCA job that has run.
        if found_a_gaussian_job_that_has_run(root, files) or found_an_orca_job_that_has_run(root, files):
            # Found a log file and gjf files before finding the ground_structure and excited_structure folders. 
            # This means we are in a non-reorganisation energy folder, so don't look further into here. 
            dirs[:] = []
            files[:] = []
            continue

        # 3.2: Check to see if we found reorganisation energy data by looking for a folder called ground_structure or excited_structure folders in root
        found_ground_structure_folder, found_excited_structure_folder = found_a_re_jobset(root)

        # 3.3: If "ground_structure" or "excited_structure" folders exist, process reorganisation energy data.
        if found_ground_structure_folder or found_excited_structure_folder:

            # 3.3.1: Check to see if the files in this reorganisation folder are Gaussian files.
            are_gaussian_files_in_ground_state_folder  = found_a_gaussian_job_that_has_run(root+'/'+ground_structure_foldername,  [file for file in os.listdir(root+'/'+ground_structure_foldername)  if os.path.isfile(root+'/'+ground_structure_foldername +'/'+file)])
            are_gaussian_files_in_excited_state_folder = found_a_gaussian_job_that_has_run(root+'/'+excited_structure_foldername, [file for file in os.listdir(root+'/'+excited_structure_foldername) if os.path.isfile(root+'/'+excited_structure_foldername+'/'+file)])

            # 3.3.2: Check to see if the files in this reorganisation folder are ORCA files.
            are_orca_files_in_ground_state_folder      = found_an_orca_job_that_has_run(root+'/'+ground_structure_foldername,  [file for file in os.listdir(root+'/'+ground_structure_foldername)  if os.path.isfile(root+'/'+ground_structure_foldername +'/'+file)])
            are_orca_files_in_excited_state_folder     = found_an_orca_job_that_has_run(root+'/'+excited_structure_foldername, [file for file in os.listdir(root+'/'+excited_structure_foldername) if os.path.isfile(root+'/'+excited_structure_foldername+'/'+file)])

            # 3.3.3: If we are dealing with a Gaussian/ORCA job, obtain information on reorganisation energy. 
            if       (are_gaussian_files_in_ground_state_folder and are_gaussian_files_in_excited_state_folder) and not (are_orca_files_in_ground_state_folder and are_orca_files_in_excited_state_folder):
                obtain_gaussian_RE_data(root, reorganisation_energy_data, start_time, ground_structure_foldername, excited_structure_foldername, lower_limit_negative_frequency, analyse_frequencies, issues)
                dirs[:] = []
                files[:] = []
            elif not (are_gaussian_files_in_ground_state_folder and are_gaussian_files_in_excited_state_folder) and     (are_orca_files_in_ground_state_folder and are_orca_files_in_excited_state_folder):
                obtain_orca_RE_data    (root, reorganisation_energy_data, start_time, ground_structure_foldername, excited_structure_foldername, lower_limit_negative_frequency, analyse_frequencies, issues)
                dirs[:] = []
                files[:] = []
            else:
                print('Note: Some reorganisation energy files in '+str(root)+'have run or running, and some not run yet. Will pass looking at this reorganisation energy dataset for now..')
                continue

    # Fourth, write the EET to an excel file.
    write_data_to_excel(reorganisation_energy_data, re_data_foldername, individual_re_data_foldername, start_time)

    # Fifth, write any issues about reorganisation energy calculations to the terminal.
    print('#'*20)
    if len(issues) > 0:
        print('The following Gaussian jobs could not be processed because they have not finished running or did not complete successfully.')
        print()
        for root, (got_eGS_gGS_energy, got_eGS_gGS_freq, got_eES_gGS_energy, got_eES_gES_energy, got_eES_gES_freq, got_eGS_gES_energy) in issues:
            print(root)
            # ==========================================================================================
            # List any ground state issues with this reorganisation energy calculation.
            print('Ground State Optimisation Calculations')
            if isinstance(got_eGS_gGS_energy,(list,tuple)):
                print('eGS_gGS_energy: Energy - '+str(got_eGS_gGS_energy[0])+' Ha, maximum_force_converged - '+str(got_eGS_gGS_energy[1])+', rms_force_converged - '+str(got_eGS_gGS_energy[2]))
            else:
                got_eGS_gGS_energy = 'Completed' if got_eGS_gGS_energy else 'Not Yet Completed'
                print('eGS_gGS_energy: '+str(got_eGS_gGS_energy))
            if isinstance(got_eGS_gGS_freq,bool):
                got_eGS_gGS_freq = 'Completed' if got_eGS_gGS_freq else 'Not Yet Completed'
            print('eGS_gGS negative frequencies: '+str(got_eGS_gGS_freq))
            if not isinstance(got_eES_gGS_energy,bool):
                print('eES_gGS_energy: Energy - '+str(got_eES_gGS_energy)+' Ha')
            else:
                got_eES_gGS_energy = 'Completed' if got_eES_gGS_energy else 'Not Yet Completed'
                print('eES_gGS_energy: '+str(got_eES_gGS_energy))
            # ==========================================================================================
            # ==========================================================================================
            # List any excited state issues with this reorganisation energy calculation.
            print('Excited State Optimisation Calculations')
            if isinstance(got_eES_gES_energy,(list,tuple)):
                print('eES_gES_energy: Energy - '+str(got_eES_gES_energy[0])+' Ha, maximum_force_converged - '+str(got_eES_gES_energy[1])+', rms_force_converged - '+str(got_eES_gES_energy[2]))
            else:
                got_eES_gES_energy = 'Completed' if got_eES_gES_energy else 'Not Yet Completed'
                print('eGS_gGS_energy: '+str(got_eES_gES_energy))
            if isinstance(got_eES_gES_freq,bool):
                got_eES_gES_freq = 'Completed' if got_eES_gES_freq else 'Not Yet Completed'
            print('eES_gES negative frequencies: '+str(got_eES_gES_freq))
            if not isinstance(got_eGS_gES_energy,bool):
                print('eGS_gES_energy: Energy - '+str(got_eGS_gES_energy)+' Ha')
            else:
                got_eGS_gES_energy = 'Completed' if got_eGS_gES_energy else 'Not Yet Completed'
                print('eGS_gES_energy: '+str(got_eGS_gES_energy))
            # ==========================================================================================
            print('-'*20)
        print('#'*20)
    print(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+' - This RE calculations program has finished successfully!')
    print('Total running time (HH:MM:SS): '+str(timedelta(seconds=time.time() - start_time)))
    print('#'*20)

# ------------------------------------------------------------------------------------------------------------------------------------------------------------






