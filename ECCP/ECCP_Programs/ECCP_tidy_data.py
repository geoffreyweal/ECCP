'''
Geoffrey Weal, ECCP_tidy_data.py, 18/4/22

This program is designed to tidy up your data folder and get rid of unnecessary files, particularly those very large files. 
'''
import os, shutil

from ECCP.ECCP_Programs.shared_general_methods.shared_gaussian_methods import folder_contains_RE_files
from ECCP.ECCP_Programs.shared_general_methods.shared_gaussian_methods import found_a_gaussian_job_that_has_run
from ECCP.ECCP_Programs.shared_general_methods.shared_gaussian_methods import did_gaussian_job_complete
from ECCP.ECCP_Programs.shared_general_methods.shared_gaussian_methods import gaussian_temp_files_to_remove
from ECCP.ECCP_Programs.shared_general_methods.shared_gaussian_methods import remove_slurm_output_files

class CLICommand:
    """Will tidy up your data folder and get rid of unnecessary files, particularly those very large files
    """

    @staticmethod
    def add_arguments(parser):
        pass

    @staticmethod
    def run(args):
        Run_method()

def Run_method():
    """
    This method will remove all unnecessary files for jobs that have completed.
    """
    
    # First, setup all the initial variables.
    current_path = os.getcwd()
    print('----------------------------------------------')
    print('Tidying Folders from the root path: '+str(current_path))
    print('----------------------------------------------')
    did_find_job = False
    did_not_tidy_jobs = []

    # Second, go through each subdirectory in the parent directory. 
    for root, dirs, files in os.walk(current_path):

        # 2.1: sort the directory just to make tidying happen in alphabetical order.
        dirs.sort()

        # 2.2: What type of calculation type are we dealing with
        contains_RE_files, RE_type = folder_contains_RE_files(root)
        if contains_RE_files:

            # 2.2.1: We are looking at a reorganisation energy calculation
            did_find_job = True
            if RE_type == 'GS':

                # 2.2.1.1: Determine which of these ground state reorganisation energy jobs has finished.
                #has_GS_GS_main_preopt_completed = did_gaussian_job_complete(root+'/GS_GS_main_opt_preopt.log')
                has_GS_GS_main_opt_completed    = did_gaussian_job_complete(root+'/GS_GS_main_opt.log')
                has_GS_GS_freq_completed        = did_gaussian_job_complete(root+'/GS_GS_freq.log')
                has_GS_ES_completed             = did_gaussian_job_complete(root+'/GS_ES.log')

                # 2.2.1.2: Did this ground state reorganisation energy finish?
                if has_GS_GS_main_opt_completed and has_GS_GS_freq_completed and has_GS_ES_completed:
                    gaussian_temp_files_to_remove(root, files, remove_chk_file=True, remove_fort7_file=True) 
                    remove_slurm_output_files(root)
                else:
                    did_not_tidy_jobs.append(root)

            elif RE_type == 'ES':

                # 2.2.2.1: Determine which of these excited state reorganisation energy jobs has finished.
                #has_ES_ES_main_preopt_completed = did_gaussian_job_complete(root+'/ES_ES_main_opt_preopt.log')
                has_ES_ES_main_opt_completed    = did_gaussian_job_complete(root+'/ES_ES_main_opt.log')
                has_ES_ES_freq_completed        = did_gaussian_job_complete(root+'/ES_ES_freq.log')
                has_ES_GS_completed             = did_gaussian_job_complete(root+'/ES_GS.log')

                # 2.2.2.2: Did this excited state reorganisation energy finish?
                if has_ES_ES_main_opt_completed and has_ES_ES_freq_completed and has_ES_GS_completed:
                    gaussian_temp_files_to_remove(root, files, remove_chk_file=True, remove_fort7_file=True) 
                    remove_slurm_output_files(root)
                else:
                    did_not_tidy_jobs.append(root)

            else:
                raise Exception('huh?')

            # 2.2.2: Do not need to move further down the subdirectories anymore, remove all dirs and files lists.
            dirs[:] = []
            files[:] = []

        elif found_a_gaussian_job_that_has_run(root, files): 

            # 2.3.1: We are looking at a non-reorganisation energy ECCP Gaussian calculation.
            #        If the output.log file shows that the program finished successfully, remove all temp files. 
            did_find_job = True
            if did_gaussian_job_complete(root+'/output.log'):
                gaussian_temp_files_to_remove(root, files, remove_chk_file=True, remove_fort7_file=True) # Check this for ICT calcs.
                remove_slurm_output_files(root)
            else:
                did_not_tidy_jobs.append(root)

            # 2.3.2: Do not need to move further down the subdirectories anymore, remove all dirs and files lists.
            dirs[:] = []
            files[:] = []
            
    # Third, print information about this tidying run.
    if not did_find_job:
        print('Finshed, but no temp files were found for tidying.')
    elif len(did_not_tidy_jobs) > 0:
        print('----------------------------------------------')
        print('The following jobs were not tidied. These are either still running or failed.')
        print()
        for root in did_not_tidy_jobs:
            print(root)
    else:
        print('Finished tidying ECCP jobs.')
    print('----------------------------------------------')

