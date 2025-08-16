'''
analyse_RE_output.py, Geoffrey Weal, 29/12/22

This method is designed to check if a RE Gaussian job has completed or not.
'''
import os
from ECCP.ECCP_Programs.shared_general_methods.shared_gaussian_methods import did_gaussian_job_complete, did_gaussian_opt_job_complete, did_freq_gaussian_job_complete
from ECCP.ECCP_Programs.shared_general_methods.shared_orca_methods     import did_orca_job_complete,     did_orca_opt_job_complete,     did_freq_orca_job_complete

def analyse_RE_output(software_type, path, calculation_type):
    """
    This method is designed to check if an ATC Gaussian job has completed or not.

    Parameters
    ----------
    path : str.
        This is the path to the RE job files
    calculation_type : str.
        This is the type of reorganisation energy we are performing, either ground_structure or excited_structure.

    Returns
    -------
    The results of the job: str.
        * 'NBY': Not begun yet.
        * 'NC' : Not complete.
        * 'C'  : Complete.
    """

    if software_type == 'Gaussian':
        return analyse_RE_output_Gaussian(path, calculation_type)
    elif software_type == "ORCA":
        return analyse_RE_output_ORCA(path, calculation_type)
    raise Exception('Error, you need to use with Gaussian or ORCA software. software_type = '+str(software_type))

# ----------------------------------------------------------------------------------------

def analyse_RE_output_Gaussian(path, calculation_type):
    """
    This method is designed to check if an ATC Gaussian job has completed or not.

    Parameters
    ----------
    path : str.
        This is the path to the RE job files
    calculation_type : str.
        This is the type of reorganisation energy we are performing, either ground_structure or excited_structure.

    Returns
    -------
    The results of the job: str.
        * 'NBY': Not begun yet.
        * 'NC' : Not complete.
        * 'C'  : Complete.
    """

    # First, if the output.log file does not exist, return this.
    if   calculation_type == 'ground_structure':

        # 1.1.1: This is the ground structure log files to check.
        main_preopt_gjf_name    = 'eGS_gGS_main_opt_preopt.gjf'
        main_preopt_name        = 'eGS_gGS_main_opt_preopt.log'
        main_opt_name           = 'eGS_gGS_main_opt.log'

        # 1.1.2: Get paths to files
        path_to_main_opt_preopt_gjf = path+'/'+main_preopt_gjf_name
        path_to_main_opt_preopt     = path+'/'+main_preopt_name
        path_to_main_opt            = path+'/'+main_opt_name

        # 1.1.3: Check if these files exist
        perform_preopt = os.path.exists(path_to_main_opt_preopt_gjf)
        preopt_begun   = os.path.exists(path_to_main_opt_preopt)
        main_opt_begun = os.path.exists(path_to_main_opt)

        # 1.1.4: If none of the log files exist, this job has not begun.
        if perform_preopt:
            if not preopt_begun:
                return 'NBY', None # Not begun yet.
        else:
            if not main_opt_begun:
                return 'NBY', None # Not begun yet.

    elif calculation_type == 'excited_structure':

        # 1.2.1: This is the excited structure log files to check.
        main_preopt_gjf_name    = 'eES_gES_main_opt_preopt.gjf'
        main_preopt_name        = 'eES_gES_main_opt_preopt.log'
        main_opt_name           = 'eES_gES_main_opt.log'

        # 1.2.2: Get paths to files
        path_to_main_opt_preopt_gjf = path+'/'+main_preopt_gjf_name
        path_to_main_opt_preopt     = path+'/'+main_preopt_name
        path_to_main_opt            = path+'/'+main_opt_name

        # 1.2.3: Check if these files exist
        perform_preopt = os.path.exists(path_to_main_opt_preopt_gjf)
        preopt_begun   = os.path.exists(path_to_main_opt_preopt)
        main_opt_begun = os.path.exists(path_to_main_opt)

        # 1.2.4: If none of the log files exist, this job has not begun.
        if perform_preopt:
            if not preopt_begun:
                return 'NBY', None # Not begun yet.
        else:
            if not main_opt_begun:
                return 'NBY', None # Not begun yet.
            
    else:
        raise Exception('Error, log name must be either "GS_GS.log" or "ES_ES.log". Path to where issue is = '+str(path)+'. Check this out.')

    # ========================================================================
    # Second, get the names of all the log files that are made for a RE job.
    
    # 2.1: convert GS --> ES and ES --> GS
    if   calculation_type == 'ground_structure':
        GS_or_ES_type = 'GS'
    elif calculation_type == 'excited_structure':
        GS_or_ES_type = 'ES'
    else:
        raise Exception('Error, log name must be either "GS_GS.log" or "ES_ES.log". Path to where issue is = '+str(path)+'. Check this out.')

    # 2.2: Get the path to the frequency file
    freq_file         = 'e'+GS_or_ES_type                   +'_g'+GS_or_ES_type+'_freq.log'
    single_point_file = 'e'+convert_GS_and_ES(GS_or_ES_type)+'_g'+GS_or_ES_type+'.log'

    # 2.3: Get the path to the single point file
    freq_log_filepath = path+'/'+freq_file
    sp_log_filepath   = path+'/'+single_point_file
    # ========================================================================

    # Third, check if the pre optimisation finished on the GS_GS or ES_ES files.
    did_main_preopt_finish, has_main_preopt_fully_converged, main_preopt_converged_image_index, total_no_of_images_made_during_main_preopt = did_gaussian_opt_job_complete(path_to_main_opt_preopt, get_most_converged_image=True, get_total_no_of_images=True) if os.path.exists(path_to_main_opt_preopt) else (None, None, None, None)

    # Fourth, check if the main optimisation finished on the GS_GS or ES_ES files.
    did_main_opt_finish,    has_main_opt_fully_converged,    main_opt_converged_image_index,    total_no_of_images_made_during_main_opt    = did_gaussian_opt_job_complete(path_to_main_opt,        get_most_converged_image=True, get_total_no_of_images=True) if os.path.exists(path_to_main_opt)        else (None, None, None, None)

    # Fifth, check if the freq calc finished on the GS_GS_freq.log or ES_ES_freq.log files.
    did_freq_finish, is_freq_force_convergence_good, no_of_negative_frequencies = did_freq_gaussian_job_complete(freq_log_filepath) if os.path.exists(freq_log_filepath) else (None, None, None)

    # Sixth, check if the single point calc finished for the GS_ES.log or ES_GS.log files.
    did_sp_finish = did_gaussian_job_complete(sp_log_filepath) if os.path.exists(sp_log_filepath) else None

    # Seventh, determine if calculations have completed or not.
    re_results_for_determining_completion = (did_main_opt_finish, did_freq_finish, did_sp_finish)
    if all([(re_result_for_determining_completion == True) for re_result_for_determining_completion in re_results_for_determining_completion]):
        if no_of_negative_frequencies == 0:
            has_completed = 'C' 
        else:
            has_completed = 'NC'
    else:
        has_completed = 'NC'

    # Eighth, return results of calculations. 
    if (has_completed == 'NBY') or (has_completed == 'C'):
        re_details = None
    elif has_completed == 'NC':
        re_details = (calculation_type, did_main_preopt_finish, has_main_preopt_fully_converged, main_preopt_converged_image_index, total_no_of_images_made_during_main_preopt, did_main_opt_finish, has_main_opt_fully_converged, main_opt_converged_image_index, total_no_of_images_made_during_main_opt, did_freq_finish, is_freq_force_convergence_good, no_of_negative_frequencies, did_sp_finish)
    else:
        return Exception('Huh?')

    # Ninth, return if the ECCP program completed and reorganisation energy information. 
    return has_completed, re_details

# ----------------------------------------------------------------------------------------

def analyse_RE_output_ORCA(path, calculation_type):
    """
    This method is designed to check if an ATC ORCA job has completed or not.

    Parameters
    ----------
    path : str.
        This is the path to the RE job files
    calculation_type : str.
        This is the type of reorganisation energy we are performing, either ground_structure or excited_structure.

    Returns
    -------
    The results of the job: str.
        * 'NBY': Not begun yet.
        * 'NC' : Not complete.
        * 'C'  : Complete.
    """

    # First, if the output.log file does not exist, return this.
    if   calculation_type == 'ground_structure':
        # 1.1.1: This is the ground structure log files to check.
        main_opt_name      = 'eGS_gGS_main_opt.out'

        # 1.1.2: Get paths to files
        path_to_main_opt   = path+'/'+main_opt_name

        # 1.1.3: If none of the log files exist, this job has not begun.
        if not os.path.exists(path_to_main_opt):
            return 'NBY', None # Not begun yet.

    elif calculation_type == 'excited_structure':
        # 1.2.1: This is the excited structure log files to check.
        main_opt_name       = 'eES_gES_main_opt.out'

        # 1.2.2: Get paths to files
        path_to_main_opt    = path+'/'+main_opt_name

        # 1.2.3: If none of the log files exist, this job has not begun.
        if not os.path.exists(path_to_main_opt):
            return 'NBY', None # Not begun yet.
            
    else:
        raise Exception('Error, output file name must be either "GS_GS.out" or "ES_ES.out". Path to where issue is = '+str(path)+'. Check this out.')

    # ========================================================================
    # Second, get the names of all the output files that are made for a RE job.
    
    # 2.1: convert GS --> ES and ES --> GS
    if   calculation_type == 'ground_structure':
        GS_or_ES_type = 'GS'
    elif calculation_type == 'excited_structure':
        GS_or_ES_type = 'ES'
    else:
        raise Exception('huh?')

    # 2.2: Get the path to the frequency file
    freq_file         = 'e'+GS_or_ES_type                   +'_g'+GS_or_ES_type+'_freq.out'
    single_point_file = 'e'+convert_GS_and_ES(GS_or_ES_type)+'_g'+GS_or_ES_type+'.out'

    # 2.3: Get the path to the single point file
    freq_output_filepath = path+'/'+freq_file
    sp_output_filepath   = path+'/'+single_point_file
    # ========================================================================

    # Third, check if the main optimisation finished on the GS_GS.out or ES_ES.out files.
    did_main_opt_finish, has_main_opt_fully_converged, main_opt_converged_image_index, total_no_of_images_made_during_main_opt = did_orca_opt_job_complete(path_to_main_opt, get_most_converged_image=True, get_total_no_of_images=True) if os.path.exists(path_to_main_opt) else (None, None, None, None)

    # Fourth, check if the freq calc finished on the GS_GS_freq.out or ES_ES_freq.out files.
    did_freq_finish, is_freq_force_convergence_good, no_of_negative_frequencies = did_freq_orca_job_complete(freq_output_filepath) if os.path.exists(freq_output_filepath) else (None, None, None)

    # Fifth, check if the single point calc finished for the GS_ES.out or ES_GS.out files.
    did_sp_finish = did_orca_job_complete(sp_output_filepath) if os.path.exists(sp_output_filepath) else None

    # Sixth, determine if calculations have completed or not.
    if   calculation_type == 'ground_structure': 
        re_results_for_determining_completion = (did_main_opt_finish, did_freq_finish, did_sp_finish)
    elif calculation_type == 'excited_structure':
        re_results_for_determining_completion = (did_main_opt_finish, did_freq_finish, did_sp_finish)
    if all([(re_result_for_determining_completion == True) for re_result_for_determining_completion in re_results_for_determining_completion]):
        if no_of_negative_frequencies == 0:
            has_completed = 'C' 
        else:
            has_completed = 'NC'
    else:
        has_completed = 'NC'

    # Seventh, return results of calculations. 
    if (has_completed == 'NBY') or (has_completed == 'C'):
        re_details = None
    elif has_completed == 'NC':
        if   calculation_type == 'ground_structure':
            re_details = ('ground_structure',  did_main_opt_finish, has_main_opt_fully_converged, main_opt_converged_image_index, total_no_of_images_made_during_main_opt, did_freq_finish, is_freq_force_convergence_good, no_of_negative_frequencies, did_sp_finish)
        elif calculation_type == 'excited_structure':
            re_details = ('excited_structure', did_main_opt_finish, has_main_opt_fully_converged, main_opt_converged_image_index, total_no_of_images_made_during_main_opt, did_freq_finish, is_freq_force_convergence_good, no_of_negative_frequencies, did_sp_finish)
        else:
            raise Exception('Error, log name must be either "eGS_gGS.log" or "eES_gES.log". logfile_name = '+str(logfile_name)+'. Check this out.')
    else:
        return Exception('Huh?')

    # Eighth, return if the ECCP program completed and reorganisation energy information. 
    return has_completed, re_details

# ----------------------------------------------------------------------------------------

def convert_GS_and_ES(input_type):
    if   input_type == 'GS':
        return 'ES'
    elif input_type == 'ES':
        return 'GS'
    else:
        raise Exception('Error in def convert_GS_and_ES, in Did_Complete_Main.py. Can only be ES or GS input_type.')

# ----------------------------------------------------------------------------------------




