"""
update_settings_methods.py, Geoffrey Weal, 3/1/2023

These methods will allow the user to change settings for this program
"""
import os

# =========================================================================================================================================
# These are the default settings for this program.

Max_jobs_in_queue_at_any_one_time_DEFAULT = 5000
Max_jobs_pending_in_queue_from_ECCP_mass_submit_DEFAULT = 10
Max_jobs_running_in_queue_from_ECCP_mass_submit_DEFAULT = Max_jobs_in_queue_at_any_one_time_DEFAULT - Max_jobs_pending_in_queue_from_ECCP_mass_submit_DEFAULT - 1

time_to_wait_before_next_submission_DEFAULT = 20.0
time_to_wait_max_queue_DEFAULT = 60.0

time_to_wait_before_next_submission_due_to_temp_submission_issue_DEFAULT = 10.0
number_of_consecutive_error_before_exitting_DEFAULT = 20

time_to_wait_before_next_submission_due_to_not_waiting_between_submissions_DEFAULT = 60.0
# =========================================================================================================================================

def check_submit_settingsTXT(path_to_settings_txt_file):
    """
    This method will read in the settings for the program. 
    """

    # First, read in the settings file. If this doesnt exist yet.
    # This method will create a new settings file if it doesn't already exist. 
    if not os.path.exists(path_to_settings_txt_file):
        write_submit_settingsTXT_file(path_to_settings_txt_file)
    try:
        settings = read_submit_settingsTXT_file(path_to_settings_txt_file)
    except Exception as ee:
        write_submit_settingsTXT_file(path_to_settings_txt_file)
        settings = read_submit_settingsTXT_file(path_to_settings_txt_file)

    # Second, return the settings. 
    return settings

def read_submit_settingsTXT_file(path_to_settings_txt_file):
    """
    This method will read the settings file.
    """

    # First, read in the settings file.
    variables_found = []
    with open(path_to_settings_txt_file,'r') as submit_settingsTXT:
        for line in submit_settingsTXT:
            if   'Max_jobs_in_queue_at_any_one_time = ' in line:
                line = line.rstrip().replace('Max_jobs_in_queue_at_any_one_time = ','')
                Max_jobs_in_queue_at_any_one_time = int(line)
                variables_found.append('Max_jobs_in_queue_at_any_one_time')
            elif 'Max_jobs_pending_in_queue_from_ECCP_mass_submit = ' in line:
                line = line.rstrip().replace('Max_jobs_pending_in_queue_from_ECCP_mass_submit = ','')
                Max_jobs_pending_in_queue_from_ECCP_mass_submit = float(line)
                variables_found.append('Max_jobs_pending_in_queue_from_ECCP_mass_submit')
            elif 'Max_jobs_running_in_queue_from_ECCP_mass_submit = ' in line:
                line = line.rstrip().replace('Max_jobs_running_in_queue_from_ECCP_mass_submit = ','')
                Max_jobs_running_in_queue_from_ECCP_mass_submit = float(line)
                variables_found.append('Max_jobs_running_in_queue_from_ECCP_mass_submit')
            elif 'time_to_wait_before_next_submission = ' in line:
                line = line.rstrip().replace('time_to_wait_before_next_submission = ','')
                time_to_wait_before_next_submission = float(line)
                variables_found.append('time_to_wait_before_next_submission')
            elif 'time_to_wait_max_queue = ' in line:
                line = line.rstrip().replace('time_to_wait_max_queue = ','')
                time_to_wait_max_queue = float(line)
                variables_found.append('time_to_wait_max_queue')
            elif 'time_to_wait_before_next_submission_due_to_temp_submission_issue = ' in line:
                line = line.rstrip().replace('time_to_wait_before_next_submission_due_to_temp_submission_issue = ','')
                time_to_wait_before_next_submission_due_to_temp_submission_issue = float(line)
                variables_found.append('time_to_wait_before_next_submission_due_to_temp_submission_issue')
            elif 'number_of_consecutive_error_before_exitting = ' in line:
                line = line.rstrip().replace('number_of_consecutive_error_before_exitting = ','')
                number_of_consecutive_error_before_exitting = int(line)
                variables_found.append('number_of_consecutive_error_before_exitting')
            elif 'time_to_wait_before_next_submission_due_to_not_waiting_between_submissions = ' in line:
                line = line.rstrip().replace('time_to_wait_before_next_submission_due_to_not_waiting_between_submissions = ','')
                time_to_wait_before_next_submission_due_to_not_waiting_between_submissions = float(line)
                variables_found.append('time_to_wait_before_next_submission_due_to_not_waiting_between_submissions')

    # Second, check that all the variables have been obtained from the settings file. 
    # 2.1: Determine which variables are contained in the settings file. 
    variables_needed = ['Max_jobs_in_queue_at_any_one_time', 'Max_jobs_pending_in_queue_from_ECCP_mass_submit', 'Max_jobs_running_in_queue_from_ECCP_mass_submit', 'time_to_wait_before_next_submission', 'time_to_wait_max_queue', 'time_to_wait_before_next_submission_due_to_temp_submission_issue', 'number_of_consecutive_error_before_exitting', 'time_to_wait_before_next_submission_due_to_not_waiting_between_submissions']
    variables_you_do_not_have_in_settingsTXT = []
    for variable in variables_needed:
        if not variable in locals():
            variables_you_do_not_have_in_settingsTXT.append(variable)
    if not len(variables_you_do_not_have_in_settingsTXT) == 0:
        print(variables_you_do_not_have_in_settingsTXT)
        import pdb; pdb.set_trace()
        exit('Error')

    # Third, check that no variables have been entered twice:
    if not (len(variables_found) == len(set(variables_found))):
        print(variables_you_do_not_have_in_settingsTXT)
        import pdb; pdb.set_trace()
        exit('Error')

    # Fourth, return all the settings from the settings file.
    return Max_jobs_in_queue_at_any_one_time, Max_jobs_pending_in_queue_from_ECCP_mass_submit, Max_jobs_running_in_queue_from_ECCP_mass_submit, time_to_wait_before_next_submission, time_to_wait_max_queue, time_to_wait_before_next_submission_due_to_temp_submission_issue, number_of_consecutive_error_before_exitting, time_to_wait_before_next_submission_due_to_not_waiting_between_submissions

def write_submit_settingsTXT_file(path_to_settings_txt_file, Max_jobs_in_queue_at_any_one_time=Max_jobs_in_queue_at_any_one_time_DEFAULT, Max_jobs_pending_in_queue_from_ECCP_mass_submit=Max_jobs_pending_in_queue_from_ECCP_mass_submit_DEFAULT, Max_jobs_running_in_queue_from_ECCP_mass_submit=Max_jobs_running_in_queue_from_ECCP_mass_submit_DEFAULT, time_to_wait_before_next_submission=time_to_wait_before_next_submission_DEFAULT, time_to_wait_max_queue=time_to_wait_max_queue_DEFAULT, time_to_wait_before_next_submission_due_to_temp_submission_issue=time_to_wait_before_next_submission_due_to_temp_submission_issue_DEFAULT, number_of_consecutive_error_before_exitting=number_of_consecutive_error_before_exitting_DEFAULT, time_to_wait_before_next_submission_due_to_not_waiting_between_submissions=time_to_wait_before_next_submission_due_to_not_waiting_between_submissions_DEFAULT):
    """
    This method will write the new settings to the settings file.
    """
    with open(path_to_settings_txt_file,'w') as submit_settingsTXT:
        submit_settingsTXT.write('Max_jobs_in_queue_at_any_one_time = '+str(Max_jobs_in_queue_at_any_one_time)+'\n')
        submit_settingsTXT.write('Max_jobs_pending_in_queue_from_ECCP_mass_submit = '+str(Max_jobs_pending_in_queue_from_ECCP_mass_submit)+'\n')
        submit_settingsTXT.write('Max_jobs_running_in_queue_from_ECCP_mass_submit = '+str(Max_jobs_running_in_queue_from_ECCP_mass_submit)+'\n')
        submit_settingsTXT.write('time_to_wait_before_next_submission = '+str(time_to_wait_before_next_submission)+'\n')
        submit_settingsTXT.write('time_to_wait_max_queue = '+str(time_to_wait_max_queue)+'\n')
        submit_settingsTXT.write('time_to_wait_before_next_submission_due_to_temp_submission_issue = '+str(time_to_wait_before_next_submission_due_to_temp_submission_issue)+'\n')
        submit_settingsTXT.write('number_of_consecutive_error_before_exitting = '+str(number_of_consecutive_error_before_exitting)+'\n')
        submit_settingsTXT.write('time_to_wait_before_next_submission_due_to_not_waiting_between_submissions = '+str(time_to_wait_before_next_submission_due_to_not_waiting_between_submissions)+'\n')

# =========================================================================================================================================

def change_settings(args):
    """
    This method will allow the user to change the settings in the settings file. 
    """
    raise Exception('Fix this if you want to run it.')
    continue_running = False
    wait_between_submissions = False
    if len(args) > 2:
        print('Error in changing submit settings:')
        print('You can only enter at least two arguments after "Adsorber slurm"')
        print('Your input arguments are: '+str(args))
        exit('This program is closing without changing any settings')
    if   args[0] == 'max':
        if len(args) == 1:
            print('Setting Max_jobs_in_queue_at_any_one_time to default ('+str(Max_jobs_in_queue_at_any_one_time_DEFAULT)+')')
            Max_jobs_in_queue_at_any_one_time = Max_jobs_in_queue_at_any_one_time_DEFAULT
        elif not instance(args[1],int):
            print('Error in changing submit sertings:')
            print('You need to enter an int if you want to change the settings for Max_jobs_in_queue_at_any_one_time')
            print('Your input value is: '+str(args[1]))
            exit('This program is closing without changing any settings')
        else:
            Max_jobs_in_queue_at_any_one_time = args[1]
        write_submit_settingsTXT_file(path_to_settings_txt_file,Max_jobs_in_queue_at_any_one_time)
    elif args[0] == 'wait':
        if len(args) == 1:
            print('Setting time_to_wait_before_next_submission to default ('+str(time_to_wait_before_next_submission_DEFAULT)+')')
            time_to_wait_before_next_submission = time_to_wait_before_next_submission_DEFAULT
        elif not instance(args[1],float):
            print('Error in changing submit sertings:')
            print('You need to enter a float if you want to change the settings for time_to_wait_before_next_submission')
            print('Your input value is: '+str(args[1]))
            exit('This program is closing without changing any settings')
        else:
            time_to_wait_before_next_submission = args[1]
        write_submit_settingsTXT_file(path_to_settings_txt_file,time_to_wait_before_next_submission)
    elif args[0] == 'wait_max_queue':
        if len(args) == 2:
            print('Setting time_to_wait_max_queue to default ('+str(time_to_wait_max_queue_DEFAULT)+')')
            time_to_wait_max_queue = time_to_wait_max_queue_DEFAULT
        elif not instance(args[1],float):
            print('Error in changing submit sertings:')
            print('You need to enter a float if you want to change the settings for time_to_wait_max_queue')
            print('Your input value is: '+str(args[1]))
            exit('This program is closing without changing any settings')
        else:
            time_to_wait_max_queue = args[1]
        write_submit_settingsTXT_file(path_to_settings_txt_file,time_to_wait_max_queue)
    elif args[0] == 'wait_error':
        if len(args) == 2:
            print('Setting time_to_wait_before_next_submission_due_to_temp_submission_issue to default ('+str(time_to_wait_before_next_submission_due_to_temp_submission_issue_DEFAULT)+')')
            time_to_wait_before_next_submission_due_to_temp_submission_issue = time_to_wait_before_next_submission_due_to_temp_submission_issue_DEFAULT
        elif not instance(args[1],float):
            print('Error in changing submit sertings:')
            print('You need to enter a float if you want to change the settings for time_to_wait_before_next_submission_due_to_temp_submission_issue')
            print('Your input value is: '+str(args[1]))
            exit('This program is closing without changing any settings')
        else:
            time_to_wait_before_next_submission_due_to_temp_submission_issue = args[1]
        write_submit_settingsTXT_file(path_to_settings_txt_file,time_to_wait_before_next_submission_due_to_temp_submission_issue)
    elif args[0] == 'no_conse_errors':
        if len(args) == 2:
            print('Setting number_of_consecutive_error_before_exitting to default ('+str(number_of_consecutive_error_before_exitting_DEFAULT)+')')
            number_of_consecutive_error_before_exitting = number_of_consecutive_error_before_exitting_DEFAULT
        elif not instance(args[1],int):
            print('Error in changing submit sertings:')
            print('You need to enter a int if you want to change the settings for number_of_consecutive_error_before_exitting')
            print('Your input value is: '+str(args[1]))
            exit('This program is closing without changing any settings')
        else:
            number_of_consecutive_error_before_exitting = args[1]
        write_submit_settingsTXT_file(path_to_settings_txt_file,number_of_consecutive_error_before_exitting)
    elif args[0] == 'wait_mass':
        if len(args) == 2:
            print('Setting time_to_wait_before_next_submission_due_to_not_waiting_between_submissions to default ('+str(time_to_wait_before_next_submission_due_to_not_waiting_between_submissions_DEFAULT)+')')
            time_to_wait_before_next_submission_due_to_not_waiting_between_submissions = time_to_wait_before_next_submission_due_to_not_waiting_between_submissions_DEFAULT
        elif not instance(args[1],float):
            print('Error in changing submit sertings:')
            print('You need to enter a float if you want to change the settings for time_to_wait_before_next_submission_due_to_not_waiting_between_submissions')
            print('Your input value is: '+str(args[1]))
            exit('This program is closing without changing any settings')
        else:
            time_to_wait_before_next_submission_due_to_not_waiting_between_submissions = args[1]
        write_submit_settingsTXT_file(path_to_settings_txt_file,time_to_wait_before_next_submission_due_to_not_waiting_between_submissions)
    elif args[0] == 'reset':
        write_submit_settingsTXT_file(path_to_settings_txt_file)
    else:
        print('No changes to seegins has been passed through to the ECCP program.')

# =========================================================================================================================================



'''
elif args[0] in ['t','true']:
    wait_between_submissions = True
elif args[0] in ['f','false']:
    wait_between_submissions = False
'''

