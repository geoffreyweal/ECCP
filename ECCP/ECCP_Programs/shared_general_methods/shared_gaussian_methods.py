'''
Geoffrey Weal, processing_OPV_Dimer_data.py, 9/3/22

This script contains methods for processing_OPV_Dimer_data.py

'''
import os
from ECCP.ECCP_Programs.shared_general_methods.shared_general_methods import get_lastline, reverse_readline

# -----------------------------------------------------------------

def folder_contains_RE_files(root):
    '''
    This method is designed to determine if a folder contains files associated with a reorganisation energy calculation.

    Parameters
    ----------
    root : str.
        This is the path to the folder you want to analyse.

    Returns
    -------
    True if this folder is used for performing a reorganisation energy calculation. False otherwise
    '''

    # First, go through each file in the root directory and see if anything of them are what we would expect for a reorganisation energy calculations.
    for file in os.listdir(root):
        if os.path.isfile(root+'/'+file):
            # Ground structure files
            if file in ['eGS_gGS_main_opt_preopt.gjf', 'eGS_gGS_main_opt.gjf', 'eGS_gGS_freq.gjf', 'eES_gGS.gjf']:
                return True, 'GS'
            if file in ['eGS_gGS_main_opt_preopt.log', 'eGS_gGS_main_opt.log', 'eGS_gGS_freq.log', 'eES_gGS.log']:
                return True, 'GS'
            # Excited structure files
            if file in ['eES_gES_main_opt_preopt.gjf', 'eES_gES_main_opt.gjf', 'eES_gES_freq.gjf', 'eGS_gES.gjf']:
                return True, 'ES'
            if file in ['eES_gES_main_opt_preopt.log', 'eES_gES_main_opt.log', 'eES_gES_freq.log', 'eGS_gES.log']:
                return True, 'ES'

    # Second, if got to this point, could not file any reorganisation energy calculation files, so return False.
    return False, None

# -----------------------------------------------------------------

def found_a_gaussian_job_that_has_run(root, files):
    '''
    This method is designed to determine if a Gaussian job is running or is currently running.

    Parameters
    ----------
    files : list
        This is a list of all the filenames to look through to see it contains Gaussian input and output files.

    Returns
    -------
    True if both the input .gjf file and the output .log files are found. 
    '''
    does_gjf_exist = False
    does_log_exist = False

    for file in files:
        if file.endswith('.gjf'):
            does_gjf_exist = True
            continue
        elif file.endswith('.log'):
            does_log_exist = True
            continue

    return (does_gjf_exist and does_log_exist)

# ==============================================================================================================================
# ==============================================================================================================================
# ==============================================================================================================================

def did_gaussian_job_complete(log_filepath):
    """
    This method will check to see if the gaussian job has completed successfully. 

    Parameters
    ----------
    log_filepath : str.
        This is the path to the Gaussian job log file.

    Returns
    -------
    True if the log file indicates the Gaussian log file indicates that their has been normal termination, otherwise return False
    """
    did_gaussian_job_terminate_normally = False

    if os.path.exists(log_filepath):
        counter = 0
        for line in reverse_readline(log_filepath):
            # Check if the gaussian file has terminated normally
            if 'Normal termination of Gaussian' in line:
                did_gaussian_job_terminate_normally = True
            # If not found the termination signal after 20 lines from the end of file, 
            # the job probably did not terminate properly
            if counter >= 20 and not did_gaussian_job_terminate_normally:
                break
            # If these tags are true, finish
            if did_gaussian_job_terminate_normally:
                break
            counter += 1

    return did_gaussian_job_terminate_normally

# -----------------------------------------------------------------

multiplier = 1.0
def did_gaussian_opt_job_complete(log_filepath, get_most_converged_image=False, get_total_no_of_images=False):
    """
    This method will check to see if the gaussian optimisation job has completed successfully. 

    Parameters
    ----------
    log_filepath : str.
        This is the path to the Gaussian job log file.
    get_most_converged_image : bool
        This tag indicates if you want to get the most converged image. If True, Yes. If False, just the most recently converged image. 
    get_total_no_of_images : bool
        This tag will indicate if you want to obtain the total number of images that have been created during this optimisation

    Returns
    -------
    did_finish_successfully : bool.
        Did the optimisation finish successfully
    has_fully_converged : bool. 
        This tag indicates if the Gaussian job converged properly (such that either all force and displacement convergence criteria were met all force measurements was 100 time lower than the force convergence thresholds).
    last_converged_image_index/most_converged_image_index
        This is the index of the image from the Gaussian output file to use for further calculations and analysis
    total_no_of_images
    """

    # Preamble, if the log file does not exist, return False.
    if not os.path.exists(log_filepath):
        if get_total_no_of_images:
            return False, None, None, None
        else:
            return False, None, None

    # First, create some settings for determining if a image has force converged
    has_fully_converged                 = False
    has_max_force_converged             = False
    has_RMS_force_converged             = False
    most_recently_converged_image_index = None

    # Second, create variables if you want to find the most converged image in the optimisation.
    if get_most_converged_image:
        lowest_maximum_force = float('inf')
        most_converged_image = None

    # Third, create variables for if you want to obtain the number of images in the optimisation output file. 
    if get_total_no_of_images:
        normal_termination = False
        total_no_of_images = 0

    # Fourth, go through the output file of the gaussian job backward to determine if the job finished or not
    reverse_image_index = -1
    for line in reverse_readline(log_filepath):
        # Check if the gaussian file has terminated normally, or if EET header is found:
        if ('Normal termination of Gaussian' in line) or (('Stationary point found' in line) or ('Optimization completed' in line)):
            most_recently_converged_image_index = reverse_image_index
            most_converged_image = reverse_image_index
            has_fully_converged = True
            if get_total_no_of_images:
                normal_termination = True
            break
        # Indicate if the optimisation have force converged to the degree I (Geoff) is happy with.
        if 'Maximum Force' in line:
            _, _, value, threshold, congerved = line.rstrip().split()
            if value == '********':
                value = float('inf')
            else:
                value = float(value)
            threshold = float(threshold)
            if value*multiplier < threshold:
                has_max_force_converged = True
            # ==============================================================================
            # If both the max force and RMS force has converged, indicate this
            if (has_max_force_converged and has_RMS_force_converged):
                # If most recently converged image index has not been given, give it
                if most_recently_converged_image_index is None:
                    most_recently_converged_image_index = reverse_image_index
                # If you want to get the most converged image index, perform the following
                if get_most_converged_image:
                    if value < lowest_maximum_force:
                        lowest_maximum_force = value
                        most_converged_image = reverse_image_index
                # If you are not recording the get_most_converged_image or get_total_no_of_images, we dont need to do anything more, so break out of loop
                if not (get_most_converged_image or get_total_no_of_images):
                    break
            # Reset the convergence booleans
            has_max_force_converged = False
            has_RMS_force_converged = False
            # move down the image index as we are on two the next (previous) image in the optimisation. 
            reverse_image_index += -1
            if get_total_no_of_images:
                total_no_of_images  +=  1
            # ==============================================================================
        if 'RMS     Force' in line:
            _, _, value, threshold, congerved = line.rstrip().split()
            if value == '********':
                value = float('inf')
            else:
                value = float(value)
            threshold = float(threshold)
            if value*multiplier < threshold:
                has_RMS_force_converged = True

    # Fourth, if you want to get the total no of images but you had a normal termination, use this method to count the number of images you had
    if get_total_no_of_images and normal_termination:
        total_no_of_images = 0
        for line in reverse_readline(log_filepath):
            if 'Maximum Force' in line:
                total_no_of_images += 1

    # Fifth, record data to return about the results about either the most converged result, or the most recently converged result
    if get_most_converged_image:
        # Return most converged result
        # If most_converged_image is None, no converged image was found, job has not completed
        if most_converged_image is None:
            to_return = [False, False, None]
        else:
            to_return = [True, has_fully_converged, most_converged_image]
    else:
        # Return most recently converged result
        # If most_converged_image is None, no converged image was found, job has not completed
        if most_recently_converged_image_index is None:
            to_return = [False, False, None]
        else:
            to_return = [True, has_fully_converged, most_recently_converged_image_index]

    # Sixth, add info for returning about the total number of images
    if get_total_no_of_images:
        to_return.append(total_no_of_images)

    # Seventh, return result as a tuple
    return tuple(to_return)

# -----------------------------------------------------------------

def did_freq_gaussian_job_complete(freq_log_filepath):
    """
    This method will check to see if the gaussian frequency job has completed successfully. 

    Parameters
    ----------
    freq_log_filepath : str.
        This is the path to the Gaussian freq job log file.
    """

    # First, determine if the job finished successfully
    #        This is to check that the force convergence using the full Hessian
    #        rather than the approximate from the optimisation is all good.
    did_gaussian_job_terminate_normally = did_gaussian_job_complete(freq_log_filepath) 
    if not did_gaussian_job_terminate_normally:
        return did_gaussian_job_terminate_normally, None, None

    # Second, check that the force convergence using the full Hessian
    #         rather than the approximate from the optimisation is all good.
    is_freq_force_convergence_good, _, _ = did_gaussian_opt_job_complete(freq_log_filepath) 
    if not is_freq_force_convergence_good:
        return did_gaussian_job_terminate_normally, is_freq_force_convergence_good, None

    # Second, obtain the frequencies from the output file.
    frequency_start_line = 'Harmonic frequencies (cm**-1), IR intensities (KM/Mole), Raman scattering'
    frequency_end_line   = ''
    start_recording = False
    all_frequencies = []
    with open(freq_log_filepath) as logFILE:
        for line in logFILE:
            # 2.1: Check if to start frequency finding or finishing looking for frequencies
            if frequency_start_line in line:
                if start_recording:
                    break
                else:
                    start_recording = True
                    continue
            # 2.2: Look and record freqencies from file.
            if 'Frequencies ---' in line:
                recorded_frequencies = line.replace('Frequencies ---','')
                for recorded_frequency in recorded_frequencies.split():
                    all_frequencies.append(float(recorded_frequency))
    all_frequencies.sort()

    # Third, determine the number of frequencies are negative
    no_of_negative_frequencies = sum([int(frequency < 0.0) for frequency in all_frequencies])

    return did_gaussian_job_terminate_normally, is_freq_force_convergence_good, no_of_negative_frequencies

# ==============================================================================================================================
# ==============================================================================================================================
# ==============================================================================================================================

gaussian_temp_types = ['.d2e','.int','.rwf','.skr']
other_files_to_remove = ['core.']
def gaussian_temp_files_to_remove(root, files, remove_chk_file=True, remove_fort7_file=False, print_to_display=True):
    '''
    This method is designed to remove temporary gaussian files. 

    These are files that will not need to be used again, and many are very large so removing is very advantagous.

    Parameters
    ----------
    root : str
        This is the path to the gaussian temp files
    files : list
        This is a list of all the filenames to look through to see it contains Gaussian temp files.
    remove_fort7_file : bool.
        If True, remove the fort.7 file if found. If False, do not remove this file.
    print_to_display : bool.
        Print what is being remove to screen. Default: True
    '''

    # For determining if to remove wfn file as multiwfn file has already been processed for ATC calculation
    wfn_file = None
    have_chg_file = False

    # Get files to remove
    temp_files_to_remove = []
    for file in files:
        for gaussian_temp_type in gaussian_temp_types:
            if file.endswith(gaussian_temp_type):
                temp_files_to_remove.append(file)
        for other_file_to_remove in other_files_to_remove:
            if other_file_to_remove in file:
                temp_files_to_remove.append(file)
        if remove_chk_file:
            if file.endswith('.chk'):
                temp_files_to_remove.append(file)
        else:
            if file in ['gaussian_freq.chk', 'gaussian_sp.chk']:
                temp_files_to_remove.append(file)
        if file.endswith('.wfn'): 
            wfn_file = file
        if file.endswith('.chg'): 
            have_chg_file = True
        if (file == 'gmon.out'):
            temp_files_to_remove.append(file)
        if (file == 'fort.7') and remove_fort7_file:
            temp_files_to_remove.append(file)

    # Determine if to remove wfn file
    if have_chg_file and (wfn_file is not None):
        temp_files_to_remove.append(wfn_file)

    # Remove files if they exist
    if len(temp_files_to_remove) > 0:
        if print_to_display:
            print('Tidying files in: '+str(root))
            print('Removing Gaussian temp files and other unnecessary files: '+str(temp_files_to_remove))
        for temp_file_to_remove in temp_files_to_remove:
            path_to_temp_file_to_remove = root+'/'+temp_file_to_remove
            if os.path.exists(path_to_temp_file_to_remove):
                os.remove(path_to_temp_file_to_remove)
    else:
        if print_to_display:
            print('No files need to be removed to tidy: '+str(root))

def remove_slurm_output_files(root):
    """
    This method will remove the slurm output files made, including the slurm-XXX.out and slurm-XXX.err files, where XXX is the job number.

    Parameters
    ----------
    root : str
        This is the path to the gaussian temp files.
    """
    for file in os.listdir(root):
        if file.startswith('slurm-') and (file.endswith('.out') or file.endswith('.err')):
            os.remove(root+'/'+file)

# -----------------------------------------------------------------




