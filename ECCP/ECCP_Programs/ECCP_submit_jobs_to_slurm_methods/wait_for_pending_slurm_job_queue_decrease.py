"""
wait_for_pending_slurm_job_queue_decrease.py, Geoffrey Weal, 3/1/2023

This method are designed to allow the program to only have a number of jobs in the queue at any one time.

This prevents the queue from being swamped with jobs submitted by this program.
"""
import getpass
from time import sleep
from subprocess import Popen, PIPE

# This list contains the slurm jobs that are in the current queue that were submitted by this program.
running_slurm_jobs_queue = []
pending_slurm_jobs_queue = []

def wait_for_pending_slurm_job_queue_decrease(job_number, Max_jobs_pending_in_queue_from_ECCP_mass_submit, Max_jobs_running_in_queue_from_ECCP_mass_submit):
    """
    This method will use the "check_if_can_submit_next_job" method to wait until the slurm queue is not full.
    
    Parameters
    ----------
    job_number : int
        This is the job number of the just submitted job by this program to slurm.
    Max_jobs_pending_in_queue_from_ECCP_mass_submit : int
        This is the maximum number of jobs that we want in the pending queue that were submitted by this program.
    Max_jobs_running_in_queue_from_ECCP_mass_submit : int
        This is the maximum number of jobs that we want in the running queue that were submitted by this program.
    """
    # First, append the just submitted job to the "pending_slurm_jobs_queue" list.
    pending_slurm_jobs_queue.append(job_number)
    waited = False

    # Second, wait until the number of jobs submitted by this program to slurm that are still pending is less than a given value.
    if not check_if_can_submit_next_job(pending_slurm_jobs_queue, Max_jobs_pending_in_queue_from_ECCP_mass_submit, Max_jobs_running_in_queue_from_ECCP_mass_submit, 'pending'):
        waited = True
        print('-----------------------------------------------------------------------------')
        print('The Pending Slurm Queue is full. Will wait to submit more job until other jobs have turned from pending to running.')
        while not check_if_can_submit_next_job(pending_slurm_jobs_queue, Max_jobs_pending_in_queue_from_ECCP_mass_submit, Max_jobs_running_in_queue_from_ECCP_mass_submit, 'pending'):
            sleep(15)
        print('The Pending Slurm Queue is now NOT full. Will submit this job now and continue submitting other jobs.')

    # Third, wait until the number of jobs submitted by this program to slurm that are still running is less than a given value.
    if not check_if_can_submit_next_job(pending_slurm_jobs_queue, Max_jobs_pending_in_queue_from_ECCP_mass_submit, Max_jobs_running_in_queue_from_ECCP_mass_submit, 'running'):
        if not waited:
            print('-----------------------------------------------------------------------------')
        waited = True
        print('The Running Slurm Queue is full. Will wait to submit this job until jobs have finished running.')
        while not check_if_can_submit_next_job(pending_slurm_jobs_queue, Max_jobs_pending_in_queue_from_ECCP_mass_submit, Max_jobs_running_in_queue_from_ECCP_mass_submit, 'running'):
            sleep(15)
        print('The Running Slurm Queue is now NOT full. Will submit this job now and continue submitting other jobs.')

    if waited:
        print('-----------------------------------------------------------------------------')

def check_if_can_submit_next_job(pending_slurm_jobs_queue, Max_jobs_pending_in_queue_from_ECCP_mass_submit, Max_jobs_running_in_queue_from_ECCP_mass_submit, list_type_to_check):
    """
    This method is designed to check if there are jobs run by this program currently in the slurm queue.

    If there are "Max_jobs_pending_in_queue_from_ECCP_mass_submit" number of these jobs in the slurm queue, this method will wait until some of the jobs have completed.
    
    Parameters
    ----------
    job_number : int
        This is the job number of the just submitted job by this program to slurm.
    Max_jobs_pending_in_queue_from_ECCP_mass_submit : int
        This is the maximum number of jobs that we want in the pending queue that were submitted by this program.
    Max_jobs_running_in_queue_from_ECCP_mass_submit : int
        This is the maximum number of jobs that we want in the running queue that were submitted by this program.

    Returns
    -------
    Are the number of jobs pending (submitted by this program) less than Max_jobs_pending_in_queue_from_ECCP_mass_submit
    """

    # First, check to see what jobs of yours are currently in the slurm queue.
    username = getpass.getuser()
    check_queue_command = ['squeue','-r','-u',username]
    process = Popen(check_queue_command, stdout=PIPE)
    process.wait()
    out, err = process.communicate()
    output = out.decode()

    # Second, determine which of the jobs submitted by this program are still pending or running in the slurm queue.
    live_running_queue = []
    live_pending_queue = []
    for line in output.split('\n')[1:-1:1]:
        #print(line.rstrip(),file=sys.stderr)
        line = line.split()
        jn = line[0]
        if '_' in jn:
            jn = jn.split('_')[0]
        running = line[4]
        #jn, _, _, _, running, _, _, _ = line.split()
        if running == 'R':
            live_running_queue.append(int(jn))
        if running == 'PD':
            live_pending_queue.append(int(jn))

    # Fourth, return the result depending on if you are wanting to analyse your pending or runnning queue. 
    if list_type_to_check == 'pending':
        return check_pending_queue(live_pending_queue, pending_slurm_jobs_queue, running_slurm_jobs_queue, Max_jobs_pending_in_queue_from_ECCP_mass_submit)
    elif list_type_to_check == 'running':
        return check_running_queue(live_running_queue, running_slurm_jobs_queue, Max_jobs_running_in_queue_from_ECCP_mass_submit)
    else:
        raise Exception('Error: list_type_to_check should be either "pending" or "running". list_type_to_check = '+str(list_type_to_check)+'. Check this out.')

def check_pending_queue(live_pending_queue, pending_slurm_jobs_queue, running_slurm_jobs_queue, Max_jobs_pending_in_queue_from_ECCP_mass_submit):
    """
    This method will compare which of the pending jobs in your slurm queue were submitted by this program.

    Parameters
    ----------
    live_pending_queue : list
        This is the pending queue that is currently in your slurm queue.
    pending_slurm_jobs_queue : list
        These are the pending jobs that were submitted by this program.
    running_slurm_jobs_queue : list
        These are the running jobs that were submitted by this program.
    Max_jobs_pending_in_queue_from_ECCP_mass_submit : int
        This is the maximum number of jobs that you want in your live slurm pending queue.

    Returns
    -------
    If the number of pending jobs is below the limit.
    """

    # First, check the pending_slurm_jobs_queue with the live_pending_queue
    for index in range(len(pending_slurm_jobs_queue)-1,-1,-1):
        job_no = pending_slurm_jobs_queue[index]
        if job_no not in live_pending_queue:
            running_slurm_jobs_queue.append(pending_slurm_jobs_queue.pop(index))

    # Second, return if the number of pending jobs is below the limit.
    return len(pending_slurm_jobs_queue) < Max_jobs_pending_in_queue_from_ECCP_mass_submit

def check_running_queue(live_running_queue, running_slurm_jobs_queue, Max_jobs_running_in_queue_from_ECCP_mass_submit):
    """
    This method will compare which of the running jobs in your slurm queue were submitted by this program.

    Parameters
    ----------
    live_running_queue : list
        This is the running queue that is currently in your slurm queue.
    running_slurm_jobs_queue : list
        These are the running jobs that were submitted by this program.
    Max_jobs_running_in_queue_from_ECCP_mass_submit : int
        This is the maximum number of jobs that you want in your live slurm running queue.

    Returns
    -------
    If the number of running jobs is below the limit.
    """

    # First, check the running_slurm_jobs_queue with the live_running_queue
    for index in range(len(running_slurm_jobs_queue)-1,-1,-1):
        job_no = running_slurm_jobs_queue[index]
        if job_no not in live_running_queue:
            del running_slurm_jobs_queue[index]

    # Second, return if the number of pending jobs is below the limit.
    return len(running_slurm_jobs_queue) < Max_jobs_running_in_queue_from_ECCP_mass_submit

# =========================================================================================================================================




