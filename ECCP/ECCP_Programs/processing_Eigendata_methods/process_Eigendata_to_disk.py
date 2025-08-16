'''
Geoffrey Weal, ECCP_processing_Eigendata_data.py, 10/6/22

This program is designed to process the Eigen data from Gaussian output files

'''
import os
from datetime import datetime

from ECCP.ECCP_Programs.shared_general_methods.shared_gaussian_methods import found_a_gaussian_job_that_has_run, did_gaussian_job_complete, gaussian_temp_files_to_remove
from ECCP.ECCP_Programs.processing_Eigendata_methods.found_data        import found_eigendata_files, should_this_calc_contain_eigendata
from ECCP.ECCP_Programs.processing_Eigendata_methods.get_eigenfiles    import get_eigenfiles

def process_Eigendata_to_disk(root, log_filename, start_time):
    """
    This method will perform the main method for writing eigendata to disk. Include doing any checks.

    Parameters
    ----------
    root : str.
        This is the path to the output.log file to process.
    log_filename : str.
        This is the filename of the ouput.log file.
    start_time
        This is the start time of this program.

    Returns
    -------
    issue : str.
        Returns the root path if their was an issue, return None if their was no problem.

    """

    # First, check to make sure all the files needed to obtain eigendata are located in the folders below
    check_folders = []
    for folder_name in ['Dimer', 'Monomer_1', 'Monomer_2']:

        # 1.1: Make sure that the Gaussian job has finished
        path_to_eigendata = root+'/'+folder_name
        files = [file_name for file_name in os.listdir(path_to_eigendata) if os.path.isfile(path_to_eigendata+'/'+file_name)]

        if not found_a_gaussian_job_that_has_run(path_to_eigendata, files):
            check_folders.append(False)
            break
        if not did_gaussian_job_complete(path_to_eigendata+'/'+log_filename):
            check_folders.append(False)
            break

        # 1.4: Make sure that either the eigendata files are in these folders, or if not the calculations
        if all(found_eigendata_files(path_to_eigendata)) or should_this_calc_contain_eigendata(path_to_eigendata, log_filename):
            check_folders.append(True)
        else:
            check_folders.append(False)
            break

    # Second, if their are any folders that are not ready to be analysed for ICT, then leave this for now and tellthe user. Otherwise, continue on to obtain and use the eigendata.
    if not all(check_folders):
        # report this Gaussian job.
        return root

    # Third, process output files and fort.7 files for eigendata if this has not already been done in these three folders
    for folder_name in ['Dimer', 'Monomer_1', 'Monomer_2']:
        path_to_eigendata = root+'/'+folder_name
        if not all(found_eigendata_files(path_to_eigendata)):
            dt_string = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
            print(str(dt_string)+' - Gathering and writing Eigendata from: '+str(path_to_eigendata+'/'+log_filename))
            get_eigenfiles(path_to_eigendata, log_filename, remove_eigendata_from_outputLOG_file=True, get_MO_data_from_fort7_file=True)
            files = [file_name for file_name in os.listdir(path_to_eigendata) if os.path.isfile(path_to_eigendata+'/'+file_name)]
            gaussian_temp_files_to_remove(path_to_eigendata, files)

    return None

# ----------------------------------------------------------------------------------------------------------------------------------
