'''
Geoffrey Weal, Run_Adsorber_submitSL_slurm.py, 16/06/2021

This program is designed to submit all sl files called submit.sl to slurm.
'''
import os
from subprocess import Popen, PIPE, TimeoutExpired

from ECCP.ECCP_Programs.ECCP_submit_jobs_to_slurm_settings_methods.settings_methods                     import check_submit_settingsTXT, change_settings, read_submit_settingsTXT_file
from ECCP.ECCP_Programs.ECCP_submit_jobs_to_slurm_methods.determine_quantum_computing_software_type     import determine_quantum_computing_software_type
from ECCP.ECCP_Programs.ECCP_submit_jobs_to_slurm_methods.ECCP_submit_gaussian_jobs_to_slurm            import general_gaussian_submission, RE_GStructure_gaussian_submission, RE_EStructure_gaussian_submission
from ECCP.ECCP_Programs.ECCP_submit_jobs_to_slurm_methods.ECCP_submit_orca_jobs_to_slurm                import general_orca_submission, RE_GStructure_orca_submission, RE_EStructure_orca_submission
from ECCP.ECCP_Programs.ECCP_submit_jobs_to_slurm_methods.check_max_jobs_in_queue_after_next_submission import check_max_jobs_in_queue_after_next_submission
from ECCP.ECCP_Programs.ECCP_submit_jobs_to_slurm_methods.countdown                                     import countdown
from ECCP.ECCP_Programs.ECCP_submit_jobs_to_slurm_methods.wait_for_pending_slurm_job_queue_decrease     import wait_for_pending_slurm_job_queue_decrease

# Get the path to the settings script.
this_scripts_path = os.path.dirname(os.path.abspath(__file__))
submit_settings_name = 'ECCP_submit_jobs_to_slurm_settings_methods/submit_settings.txt'
path_to_settings_txt_file = this_scripts_path+'/'+submit_settings_name

class CLICommand:
    """Submit ECCP jobs to slurm.
    """

    @staticmethod
    def add_arguments(parser):
        parser.add_argument('are_RE_jobs_running_currently', nargs='*', help='Enter in if there are any ECCP jobs currently running on slurm (Default: True).')
        parser.add_argument('--run_solvents', nargs=1, help='If True, dont run solvents.',default='True')

    @staticmethod
    def run(args_submit):

        # First, determine if there are any ECCP jobs are currently running in slurm, as specified by the user.
        are_RE_jobs_running_currently = args_submit.are_RE_jobs_running_currently

        if   len(are_RE_jobs_running_currently) == 0:
            are_RE_jobs_running_currently = True
        elif len(are_RE_jobs_running_currently) >= 2:
            exit('Error: You can only input one variable for are_RE_jobs_running_currently. Given are_RE_jobs_running_currently: '+str(are_RE_jobs_running_currently))
        else:
            are_RE_jobs_running_currently = are_RE_jobs_running_currently[0]
            if   are_RE_jobs_running_currently.lower() in ['t','true']:
                are_RE_jobs_running_currently = True
            elif are_RE_jobs_running_currently.lower() in ['f','false']:
                are_RE_jobs_running_currently = False
            else:
                exit('Error: input for are_RE_jobs_running_currently must be either "True" or "False". Given are_RE_jobs_running_currently: '+str(are_RE_jobs_running_currently))

        # Second, determine if you want to run solvents
        run_solvents = args_submit.run_solvents[0].lower()
        if run_solvents in ['t', 'true']:
            run_solvents = True
        elif run_solvents in ['f', 'false']:
            run_solvents = False
        else:
            raise Exception('Error: --run_solvents must be either True or False. Given input: '+str(run_solvents))

        # Second, use this method to create a settings.txt file if it doesn't already exist, and check that the current settings.txt file can be read without problems.
        check_submit_settingsTXT(path_to_settings_txt_file)

        # Third, run this program.
        Run_method(are_RE_jobs_running_currently, run_solvents)

# =========================================================================================================================================

def Run_method(are_RE_jobs_running_currently, run_solvents):
    '''
    This program is designed to submit all sl files called submit.sl to slurm.

    Parameters
    ----------
    are_RE_jobs_running_currently : bool.
        This tag indicates if any ECCP jobs are currently running on slurm.
    '''

    print('###########################################################################')
    print('###########################################################################')
    print('Run_submitSl_slurm.py')
    print('###########################################################################')
    print('This program is designed to submit all your submit.sl scripts appropriately to slurm.')
    print('###########################################################################')
    print('###########################################################################')

    # Second, read the settings from the settings file. 
    Max_jobs_in_queue_at_any_one_time, Max_jobs_pending_in_queue_from_ECCP_mass_submit, Max_jobs_running_in_queue_from_ECCP_mass_submit, time_to_wait_before_next_submission, time_to_wait_max_queue, time_to_wait_before_next_submission_due_to_temp_submission_issue, number_of_consecutive_error_before_exitting, time_to_wait_before_next_submission_due_to_not_waiting_between_submissions = read_submit_settingsTXT_file(path_to_settings_txt_file)
    wait_between_submissions = False

    if are_RE_jobs_running_currently:
        print('Will assume that RE jobs are currently running')
    else:
        print('Will assume that RE jobs are not currently running, and will also submit frequency calculations and single point calculations if the main optimisation has completed satisfactory.')

    # Third, indicate if there will be a small wait between jobs.
    if wait_between_submissions == True:
        print('This program will wait one minute between submitting jobs.')
    else:
        print('This program will not wait between submitting jobs.')
    print('Will begin to search for submit.sl and other .sl files.')
    print('***************************************************************************')

    # Fourth, time to submit all the GA scripts! Lets get this stuff going!
    if not wait_between_submissions:
        max_consec_counter = 250
        consec_counter = 0
    errors_list = []
    error_counter = 0
    path = os.getcwd()
    for (dirpath, dirnames, filenames) in os.walk(path):
        dirnames.sort()
        filenames.sort()

        # 4.1: Determine if the following submit scripts are in this folder
        is_submitSL_in_filenames                    =  'submit.sl'                  in filenames

        is_eGS_gGS_main_opt_submitSL_in_filenames   =  'eGS_gGS_main_opt_submit.sl' in filenames
        is_eGS_gGS_freq_submitSL_in_filenames       =  'eGS_gGS_freq_submit.sl'     in filenames
        is_eES_gGS_submitSL_in_filenames            =  'eES_gGS_submit.sl'          in filenames

        is_eES_gES_main_opt_submitSL_in_filenames   = ('eES_gES_main_opt_submit.sl' in filenames)
        is_eES_gES_freq_submitSL_in_filenames       =  'eES_gES_freq_submit.sl'     in filenames
        is_eGS_gES_submitSL_in_filenames            =  'eGS_gES_submit.sl'          in filenames

        # 4.2: If there is a submit file to submit, submit it.        
        if any([is_submitSL_in_filenames, is_eGS_gGS_main_opt_submitSL_in_filenames, is_eGS_gGS_freq_submitSL_in_filenames, is_eES_gGS_submitSL_in_filenames, is_eES_gES_main_opt_submitSL_in_filenames, is_eES_gES_freq_submitSL_in_filenames, is_eGS_gES_submitSL_in_filenames]):
            
            # 4.2.1, determine what calculations we are looking at.
            software_type = determine_quantum_computing_software_type(dirpath, filenames)

            # 4.2.2: Figure out which submit files in the folder should be submitted to slurm.
            submission_filenames = []
            if is_submitSL_in_filenames:
                # Submitting either ATC or EET calculation.
                if software_type == 'Gaussian':
                    submission_filenames += general_gaussian_submission(filenames, are_RE_jobs_running_currently)
                elif software_type == 'ORCA':
                    submission_filenames += general_orca_submission(filenames, are_RE_jobs_running_currently)
                else:
                    raise Exception('ERROR: Could not determine what software will be used in this submission file.')

            elif (is_eGS_gGS_main_opt_submitSL_in_filenames or is_eGS_gGS_freq_submitSL_in_filenames or is_eES_gGS_submitSL_in_filenames):
                # Submitting ground structure reorganisation calculations.
                if software_type == 'Gaussian':
                    submission_filenames += RE_GStructure_gaussian_submission(filenames, dirpath, are_RE_jobs_running_currently)
                elif software_type == 'ORCA':
                    submission_filenames += RE_GStructure_orca_submission(filenames, dirpath, are_RE_jobs_running_currently)
                else:
                    raise Exception('ERROR: Could not determine what software will be used in this submission file.')

            elif (is_eES_gES_main_opt_submitSL_in_filenames or is_eES_gES_freq_submitSL_in_filenames or is_eGS_gES_submitSL_in_filenames):
                # Submitting excited structure reorganisation calculations.
                if software_type == 'Gaussian':
                    submission_filenames += RE_EStructure_gaussian_submission(filenames, dirpath, are_RE_jobs_running_currently)
                elif software_type == 'ORCA':
                    submission_filenames += RE_EStructure_orca_submission(filenames, dirpath, are_RE_jobs_running_currently)
                else:
                    raise Exception('ERROR: Could not determine what software will be used in this submission file.')

            # 4.2.3: If we do not have jobs to submit to slurm, move on. 
            if len(submission_filenames) == 0:
                dirnames[:] = []
                filenames[:] = []
                continue

            # ====================================================================================================
            # 4.4: Submit the jobs to slurm
            for submission_filename in submission_filenames:

                # ----------------------------------------------------------------
                # 4.4.1: Determine if it is the right time to submit jobs
                print('*****************************************************************************')
                while True:
                    reached_max_jobs, number_in_queue = check_max_jobs_in_queue_after_next_submission(dirpath, Max_jobs_in_queue_at_any_one_time)
                    if reached_max_jobs:
                        print('-----------------------------------------------------------------------------')
                        print('You can not have any more jobs in the queue before submitting the mass_sub. Will wait a bit of time for some of them to complete')
                        print('Number of Jobs in the queue = '+str(number_in_queue))
                        countdown(time_to_wait_before_next_submission)
                        print('-----------------------------------------------------------------------------')
                    else:
                        print('The number of jobs in the queue currently is: '+str(number_in_queue))
                        break
                
                # 4.4.2: Submit the jobs
                os.chdir(dirpath)
                name = dirpath.replace(path, '').split('/', -1)[1:]
                
                # If you dont want to run solvents, check if this is a solvent molecule. 
                if not run_solvents:
                    if name[-3].endswith('S'): # Prevent solvents from running
                        continue
                
                name = "_".join(str(x) for x in name)
                print("Submitting " + str(name) + " to slurm.")
                print('Submission .sl file found in: '+str(os.getcwd()))
                print('Submission filename: '+str(submission_filename))
                error_counter = 0
                while True:
                    if error_counter == number_of_consecutive_error_before_exitting:
                        break
                    else:
                        submitting_command = ['sbatch', str(submission_filename)]
                        proc = Popen(submitting_command, stdout=PIPE, stderr=PIPE) # shell=True, 
                        try:
                            if not (proc.wait(timeout=(2*60)) == 0): # 120 seconds
                                # A problem occurred during the submission. Report this and wait a bit before trying again.
                                error_counter += 1
                                if error_counter == number_of_consecutive_error_before_exitting:
                                    print('----------------------------------------------')
                                    print('Error in submitting submit script to slurm.')
                                    print('I got '+str(number_of_consecutive_error_before_exitting)+" consecutive errors. Something must not be working right somewhere. I'm going to stop here just in case something is not working.")
                                    print('')
                                    print('The following submit.sl scripts WERE NOT SUBMITTED TO SLURM')
                                    print('')
                                else:
                                    stdout, stderr = proc.communicate()
                                    print('----------------------------------------------')
                                    print('Error in submitting submit script to slurm. This error was:')
                                    print(stderr)
                                    print('Number of consecutive errors: '+str(error_counter))
                                    print('Run_submitSL_slurm.py will retry submitting this job to slurm after '+str(time_to_wait_before_next_submission_due_to_temp_submission_issue)+' seconds of wait time')
                                    print('----------------------------------------------')
                                    countdown(time_to_wait_before_next_submission_due_to_temp_submission_issue)
                            else:
                                # Submission successful, report this and move on.
                                stdout, stderr = proc.communicate()
                                job_number = int(stdout.decode("utf-8").replace('Submitted batch job',''))
                                print("Submitted " + str(name) + " to slurm: "+str(job_number))
                                # Wait until the running and pending queue for this program is available to move on.
                                wait_for_pending_slurm_job_queue_decrease(job_number, Max_jobs_pending_in_queue_from_ECCP_mass_submit, Max_jobs_running_in_queue_from_ECCP_mass_submit)
                                break
                        except TimeoutExpired:
                            # A problem occurred during the submission, sbatch timedout. Report this and wait a bit before trying again.
                            proc.kill()
                            error_counter += 1
                            if error_counter == number_of_consecutive_error_before_exitting:
                                print('----------------------------------------------')
                                print('Error in submitting submit script to slurm. Job timed-out after 2 minutes.')
                                print('I got '+str(number_of_consecutive_error_before_exitting)+" consecutive errors. Something must not be working right somewhere. I'm going to stop here just in case something is not working.")
                                print('')
                                print('The following submit.sl scripts WERE NOT SUBMITTED TO SLURM')
                                print('')
                            else:
                                print('----------------------------------------------')
                                print('Error in submitting submit script to slurm. Job timed-out after 2 minutes.')
                                print('Number of consecutive errors: '+str(error_counter))
                                print('Run_submitSL_slurm.py will retry submitting this job to slurm after '+str(time_to_wait_before_next_submission_due_to_temp_submission_issue)+' seconds of wait time')
                                print('----------------------------------------------')
                                countdown(time_to_wait_before_next_submission_due_to_temp_submission_issue)

                # 4.4.3: Check that there were any errors, and wait until it is possible to submit another job without going over the maximum limit for this user.
                if error_counter == number_of_consecutive_error_before_exitting:
                    print(dirpath)
                    errors_list.append(dirpath)
                else:
                    if wait_between_submissions:
                        reached_max_jobs, number_in_queue = check_max_jobs_in_queue_after_next_submission(dirpath)
                        print('The number of jobs in the queue after submitting job is currently is: '+str(number_in_queue))
                        #print('Will wait for '+str(time_to_wait_max_queue)+' to give time between consecutive submissions')
                        countdown(time_to_wait_max_queue)
                        print('*****************************************************************************')

                # 4.4.4: If you are waiting between 
                dirnames[:] = []
                filenames[:] = []
                if not wait_between_submissions:
                    if consec_counter >= max_consec_counter:
                        print('----------------------------------------------')
                        print('As you are not waiting between consecutive submissions, it is good practise to wait for a minute at some stage')
                        print(str(max_consec_counter) +' have been submitted consecutively. Will not wait for '+str(time_to_wait_before_next_submission_due_to_not_waiting_between_submissions)+' s before continuing')
                        print('----------------------------------------------')
                        countdown(time_to_wait_before_next_submission_due_to_not_waiting_between_submissions)
                        consec_counter = 0
                    else:
                        consec_counter += 1

            # ====================================================================================================

    # Fifth, check out if there were any issues that meant that this program has to finish prematurally. 
    if len(errors_list) > 0:
        print('----------------------------------------------')
        print()
        print('"Run_submitSL_slurm.py" will finish WITHOUT HAVING SUBMITTED ALL JOBS.')
        print()
        print('*****************************************************************************')
        print('The following submit.sl SCRIPTS WERE NOT SUBMITTED SUCCESSFULLY.')
        print()
        for error_dirpath in errors_list:
            print(error_dirpath)
        print('*****************************************************************************') 
    else:
        print('*****************************************************************************')
        print('*****************************************************************************')
        print('*****************************************************************************')
        print('All submit.sl scripts have been submitted to slurm successfully.')
        print('*****************************************************************************')
        print('*****************************************************************************')
        print('*****************************************************************************')

# ------------------------------------------------------------------------------------------------






