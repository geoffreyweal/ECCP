"""
check_max_jobs_in_queue_after_next_submission.py, Geoffrey Weal, 4/1/2023

This method will check to make sure that the next job can be submitted and not go over the allocated number of maximum jobs.
"""
import getpass
import subprocess

username = getpass.getuser()
command = "squeue -r -u "+str(username)
def check_max_jobs_in_queue_after_next_submission(dirpath, Max_jobs_in_queue_at_any_one_time):
    """
    This method will check to make sure that the next job can be submitted and not go over the allocated number of maximum jobs.

    Parameters
    ----------
    dirpath : str.
        This is the settings? not needed to be programmed yet.
    Max_jobs_in_queue_at_any_one_time : int
        This is the maximum limit of jobs that can be in the user queue

    Returns
    -------

    """
    while True:
        text = myrun(command)
        nlines = len(text.splitlines())-1
        if not (nlines == -1):
            break
        else:
            print('Could not get the number of jobs in the slurm queue. Retrying to get this value.')
    number_of_trials_to_be_submitted = get_number_to_trials_that_will_be_submitted_by_submitSL(dirpath)
    if nlines > Max_jobs_in_queue_at_any_one_time - number_of_trials_to_be_submitted:
        return True, nlines
    else:
        return False, nlines

def get_number_to_trials_that_will_be_submitted_by_submitSL(dirpath):
    return 1

def myrun(cmd):
    """
    This method will run a task in the terminal/bash.

    Parameters
    ----------
    cmd : str. 
        This is the command you would like to run in the terminal/bash.
    """
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout_lines = list(iter(proc.stdout.readline,b''))
    stdout_lines = [line.decode() for line in stdout_lines]
    return ''.join(stdout_lines)

