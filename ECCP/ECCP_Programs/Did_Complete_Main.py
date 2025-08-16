'''
Did_Complete_Main.py, Geoffrey Weal, 08/03/2019

This program will determine which of your dimers have been successfully calculated in Gaussian.
'''
import os, sys
from tqdm import tqdm

from ECCP.ECCP_Programs.shared_general_methods.shared_general_methods       import reverse_readline
from ECCP.ECCP_Programs.Did_Complete_Main_methods.did_finish_calc_on_system import did_finish_calc_on_system

# -------------------------------------------------------------------------------

def Did_Complete_Main(general_path):
    """
    This method will go through folders in search of output.log files, and will determine from those output.log files if they had finished successfully or not.

    Parameters
    ----------
    general_path : str.
        This is the overall directory to search through for output.log files.

    Returns
    -------
    atc_jobs_results : str.
        These are the results of ATC jobs.
    re_jobs_results : str.
        These are the results of RE jobs.
    fc_jobs_results : str.
        These are the results of FC jobs.
    eet_jobs_results : str.
        These are the results of EET jobs.
    ict_jobs_results : str.
        These are the results of ICT jobs.
    unalligned_jobs : str.
        These are all the unalligned jobs 
    """

    # First, initialise lists to record results to.
    atc_jobs_finished_successfully = []
    atc_jobs_finished_unsuccessfully = []
    atc_jobs_not_begun = []

    re_jobs_finished_successfully = []
    re_jobs_finished_unsuccessfully = []
    re_jobs_not_begun = []

    fc_jobs_finished_successfully = []
    fc_jobs_finished_unsuccessfully = []
    fc_jobs_not_begun = []

    eet_jobs_finished_successfully = []
    eet_jobs_finished_unsuccessfully = []
    eet_jobs_not_begun = []

    eigendata_jobs_finished_successfully = []
    eigendata_jobs_finished_unsuccessfully = []
    eigendata_jobs_not_begun = []

    unalligned_jobs = []
    errored_jobs = []

    # Second, determine all the Gaussian jobs to check. 
    original_path = os.getcwd()
    pbar = tqdm(os.walk(general_path), unit='Jobs')
    for root, dirs, files in pbar:

        # 2.1: Sort dirs and files so that things come out in alphabetical order.
        dirs.sort()
        files.sort()

        # 2.2: Determine if a .gjf or inp. file is found. If so, we have found a Gaussian/ORCA job.
        if   ('eGS_gGS_main_preopt.gjf' in files) or ('eGS_gGS_main_opt.gjf' in files):
            software_type = 'Gaussian'
            input_file_name      = 'eGS_gGS_main_opt.gjf'
            output_file_name     = 'eGS_gGS_main_opt.log'
        elif ('eES_gES_main_preopt.gjf' in files) or ('eES_gES_main_opt.gjf' in files) or ('gaussian_parameters_ES.txt' in files):
            software_type = 'Gaussian'
            input_file_name      = 'eES_gES_main_opt.gjf'
            output_file_name     = 'eES_gES_main_opt.log'
        elif ('eGS_gGS_main_preopt.inp' in files) or ('eGS_gGS_main_opt.inp' in files):
            software_type = 'ORCA'
            input_file_name      = 'eGS_gGS_main_opt.inp'
            output_file_name     = 'eGS_gGS_main_opt.out'
        elif ('eES_gES_main_preopt.inp' in files) or ('eES_gES_main_opt.inp' in files) or ('orca_parameters_ES.txt' in files):
            software_type = 'ORCA'
            input_file_name      = 'eES_gES_main_opt.inp'
            output_file_name     = 'eES_gES_main_opt.out'
        else:
            for file in files:
                # Is this a Gaussian Calculation
                if '.gjf' in file:
                    software_type = 'Gaussian'
                    input_file_name      = file
                    output_file_name     = 'output.log'
                    break
                # Is this an ORCA Calculation
                if '.inp' in file:
                    software_type = 'ORCA'
                    input_file_name      = file
                    output_file_name     = 'output.out'
                    break
            else:
                continue

        # 2.3: Print details of where the program is up to:
        pbar.set_description(root.replace(original_path+'/',''))

        # 2.4: Go through the output.log file to see if the job finished successfully or not.
        try:
            job_type, completion_stage, re_details = did_finish_calc_on_system(root, software_type, input_file_name, output_file_name)
            successfully_analysed_calculations = True
        except Exception as exception_message:
            successfully_analysed_calculations = False
            error_message = exception_message

        # 2.5: Add job path to appropriate list
        if successfully_analysed_calculations:
            if job_type == 'ATC':
                add_to_list(root, completion_stage, atc_jobs_finished_successfully, atc_jobs_finished_unsuccessfully, atc_jobs_not_begun)
            elif job_type == 'RE':
                add_to_list((root, re_details), completion_stage, re_jobs_finished_successfully,  re_jobs_finished_unsuccessfully,  re_jobs_not_begun)
            elif job_type == 'FC':
                add_to_list(root, completion_stage, fc_jobs_finished_successfully,  fc_jobs_finished_unsuccessfully,  fc_jobs_not_begun)
            elif job_type == 'EET':
                add_to_list(root, completion_stage, eet_jobs_finished_successfully, eet_jobs_finished_unsuccessfully, eet_jobs_not_begun)
            else:
                unalligned_jobs.append(root)
        else:
            errored_jobs.append((root,error_message))

        # 2.6: This will prevent the program looking further down the directories.
        dirs[:] = []
        files[:] = []

    # Third, sort each list alphabetically and combine together for easy management of data.
    atc_jobs_results       = (sorted(atc_jobs_finished_successfully), sorted(atc_jobs_finished_unsuccessfully), sorted(atc_jobs_not_begun))
    re_jobs_results        = (sorted(re_jobs_finished_successfully),  sorted(re_jobs_finished_unsuccessfully),  sorted(re_jobs_not_begun))
    fc_jobs_results        = (sorted(fc_jobs_finished_successfully),  sorted(fc_jobs_finished_unsuccessfully),  sorted(fc_jobs_not_begun))
    eet_jobs_results       = (sorted(eet_jobs_finished_successfully), sorted(eet_jobs_finished_unsuccessfully), sorted(eet_jobs_not_begun))
    eigendata_jobs_results = (sorted(eigendata_jobs_finished_successfully), sorted(eigendata_jobs_finished_unsuccessfully), sorted(eigendata_jobs_not_begun))

    # Fourth, return results of ECCP Gaussian jobs. 
    return atc_jobs_results, re_jobs_results, fc_jobs_results, eet_jobs_results, eigendata_jobs_results, unalligned_jobs, errored_jobs

def add_to_list(root_and_stuff, completion_stage, jobs_finished_successfully, jobs_finished_unsuccessfully, jobs_not_begun):
    if completion_stage == 'NBY':
        jobs_not_begun.append(root_and_stuff)
    elif completion_stage == 'NC':
        jobs_finished_unsuccessfully.append(root_and_stuff)
    elif completion_stage == 'NCh':
        jobs_finished_unsuccessfully.append(root_and_stuff)
    elif completion_stage == 'C' or (completion_stage == 'NF'):
        jobs_finished_successfully.append(root_and_stuff)
    else:
        raise Exception('Error: completion_stage must be either "NBY", "NC", "NCh" or "C". completion_stage = '+str(completion_stage))



