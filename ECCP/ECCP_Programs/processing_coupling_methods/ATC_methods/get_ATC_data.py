'''
Geoffrey Weal, get_EET_data.py, 16/3/23

This program is designed to gather EET data.
'''
import os, time
from datetime import datetime, timedelta
from SUMELF import import_CHG_file
from ECCP.ECCP_Programs.shared_general_methods.shared_gaussian_methods                      import found_a_gaussian_job_that_has_run, did_gaussian_job_complete, gaussian_temp_files_to_remove
from ECCP.ECCP_Programs.processing_coupling_methods.ATC_methods.processing_ATC_data_methods import is_this_calc_an_atc_calc

def get_ATC_data(overall_path, log_filename, start_time):
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
    print('Gathering Gaussian ATC data')
    ATC_coupling_data = {}
    issues = []
    for root, dirs, files in os.walk(overall_path):
        dirs.sort()
        # Determine if their is a Gaussian job that has run.
        if found_a_gaussian_job_that_has_run(root, files):
            dirs[:] = []
            if not is_this_calc_an_atc_calc(root):
                continue
            print('------------------------------------------------')
            print(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+': Found a Gaussian job.')
            # Determine if the job completed or not. 
            path_to_log_file = root+'/'+'output.log'
            if did_gaussian_job_complete(path_to_log_file):
                print(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+': Processing: '+str(root))
                # Import the output.chg file
                CHG_filepath = root + '/'+'output.chg'
                atoms = import_CHG_file(CHG_filepath)
                # Record data, and remove temp Gaussian files.
                calculation_details = tuple(root.split('/')[-3:])
                ATC_coupling_data[calculation_details] = (root, atoms)
                print('Processed EET data in (HH:MM:SS): '+str(timedelta(seconds=time.time() - start_time)))
                # Remove unnecessary files.
                gaussian_temp_files_to_remove(root, files, remove_fort7_file=True)
            else:
                # Report this Gaussian job.
                print(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+': Their was an issue with this job. Path: '+str(root))
                issues.append(root)
    return ATC_coupling_data, issues



    