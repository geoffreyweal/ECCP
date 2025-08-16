'''
Geoffrey Weal, get_EET_data.py, 16/3/23

This program is designed to gather EET data.
'''
import os, time
from datetime import datetime, timedelta
from ECCP.ECCP_Programs.shared_general_methods.shared_gaussian_methods     import found_a_gaussian_job_that_has_run, did_gaussian_job_complete, gaussian_temp_files_to_remove
from ECCP.ECCP_Programs.processing_EET_methods.processing_EET_data_methods import is_this_calc_an_eet_calc, get_electronic_coupling_of_lowest_TD_state

def get_EET_data(overall_path, log_filename, start_time):
    """
    This method will write txt files that contains the coupling energies for the dimers for a crystal.

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
    electronic_coupling_data : dict.
        This is the dictionary that EET data will be gathered into. 
    issues: list
        These are the list of paths of EET data that had an issue for one reason or another. 
    """
    print('Gathering Gaussian EET data')
    electronic_coupling_data = {}
    issues = []
    for root, dirs, files in os.walk(overall_path):
        dirs.sort()
        # Determine if their is a Gaussian job that has run.
        if found_a_gaussian_job_that_has_run(root, files):
            dirs[:] = []
            path_to_log_file = root+'/'+log_filename
            if not is_this_calc_an_eet_calc(path_to_log_file):
                continue
            print('------------------------------------------------')
            print(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+': Found a Gaussian job.')
            # Determine if the job completed or not. 
            if did_gaussian_job_complete(path_to_log_file):
                print(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+': Processing: '+str(root))
                electronic_coupling_datum = get_electronic_coupling_of_lowest_TD_state(path_to_log_file)
                # If error is returned, report this Gaussian job as an issue. 
                if electronic_coupling_datum == 'error':
                    print(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+': Their was an issue with this job. Path: '+str(root))
                    issues.append(root)
                    continue
                # Record data, and remove temp Gaussian files.
                calculation_details = tuple(root.split('/')[-3:]) # This tuple contains (crystal_name, Dimer_name, Functional_and_Basis_Set_name)
                electronic_coupling_data[calculation_details] = (root, electronic_coupling_datum)
                print('Processed EET data in (HH:MM:SS): '+str(timedelta(seconds=time.time() - start_time)))
                # Remove unnecessary files.
                gaussian_temp_files_to_remove(root, files, remove_fort7_file=True)
            else:
                # Report this Gaussian job.
                print(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+': Their was an issue with this job. Path: '+str(root))
                issues.append(root)
    return electronic_coupling_data, issues

