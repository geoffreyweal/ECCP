'''
analyse_EET_output.py, Geoffrey Weal, 29/12/22

This method is designed to check if an EET Gaussian job has completed or not.
'''
import os
from ECCP.ECCP_Programs.shared_general_methods.shared_general_methods import reverse_readline

def analyse_EET_output(software_type, path_to_output):
    """
    This method is designed to check if an EET Gaussian job has completed or not.

    Parameters
    ----------
    path_to_output : str.
        This is the path to the EET log/out file. 

    Returns
    -------
    The results of the job: str.
        * 'NBY': Not begun yet.
        * 'NC' : Not complete.
        * 'C'  : Complete.
    """

    if software_type == "Gaussian":

        # First, if the output.log file does not exist, return this.
        if not os.path.exists(path_to_output):
            return 'NBY' # Not begun yet.

        # Second, check to see if this EET calculation has finished successfully.
        did_terminate_normally = False
        have_EET_details = False
        counter = 0
        for line in reverse_readline(path_to_output):
            # Check if the gaussian file has terminated normally
            if 'Normal termination of Gaussian' in line:
                did_terminate_normally = True
            # If not found the termination signal after 20 lines from the end of file, 
            # The job probably did not terminate properly
            if counter >= 20 and not did_terminate_normally:
                break
            # Indicate if EET header is found, indicate true
            if 'Electronic Coupling for Excitation Energy Tranfer' in line:
                have_EET_details = True
            # If this is found, prbably should have found the EET by then, so dont bother anymore
            if 'initial guesses have been made' in line:
                break
            # If these tags are true, finish
            if did_terminate_normally and have_EET_details:
                break
            counter += 1

        # Third, indicate the stage of completion of this EET job.
        if (did_terminate_normally and have_EET_details):
            return 'C'
        else:
            return 'NC'

    elif software_type == "ORCA":

        raise Exception('add ORCA stuff')




