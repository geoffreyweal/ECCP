'''
did_finish_calc_on_system.py, Geoffrey Weal, 29/12/22

This method is designed to determine if the ECCP Gaussian job being analysed finished successfully or not.
'''
from ECCP.ECCP_Programs.Did_Complete_Main_methods.determine_what_the_job_is import determine_what_the_job_is
from ECCP.ECCP_Programs.Did_Complete_Main_methods.analyse_ATC_output        import analyse_ATC_output
from ECCP.ECCP_Programs.Did_Complete_Main_methods.analyse_RE_output         import analyse_RE_output
from ECCP.ECCP_Programs.Did_Complete_Main_methods.analyse_EET_output        import analyse_EET_output
from ECCP.ECCP_Programs.Did_Complete_Main_methods.analyse_ICT_output        import analyse_ICT_output

def did_finish_calc_on_system(path, software_type, input_file_name, output_file_name):
    """
    This method will go through the output.log file from a ECCP Gaussian Job and will determine if it completed successfully or not.
 
    Parameters
    ----------
    path : str.
        This is the path to the output.log file.
    software_type : str.
        This is the type of software that was used to perform these calculations.
    input_file_name : str.
        This is the name of the input file. 
    output_file_name : str.
        This is the name of the output file. 

    Returns
    -------
    job_type : str. or None
        This is the type of job that is being analysed
    job_completion_stage : str.
        This indicates what stage of completion the job has got to. 
        * 'NBY': Not begun yet.
        * 'NC' : Not complete.
        * 'C'  : Complete.
    """

    # First, determine what calculations we are looking at.
    job_type = determine_what_the_job_is(path, software_type, output_file_name, input_file_name, break_string='Input orientation:')

    # Second, if job_type is not found, return None.
    re_details = None
    if job_type is None:
        return None, None, re_details
    elif job_type == 'ATC':
        path_to_output_file = path+'/'+output_file_name
        path_to_output_CHG  = path+'/'+'output.chg'
        job_completion_stage = analyse_ATC_output(software_type, path_to_output_file, path_to_output_CHG)
    elif job_type == 'RE':
        if   (input_file_name == 'eGS_gGS_main_opt_preopt.gjf') or (input_file_name == 'eGS_gGS_main_opt.gjf') or (input_file_name == 'eGS_gGS_main_opt_preopt.inp') or (input_file_name == 'eGS_gGS_main_opt.inp'):
            calculation_type = 'ground_structure'
        elif (input_file_name == 'eES_gES_main_opt_preopt.gjf') or (input_file_name == 'eES_gES_main_opt.gjf') or (input_file_name == 'eES_gES_main_opt_preopt.inp') or (input_file_name == 'eES_gES_main_opt.inp') or (input_file_name == 'gaussian_parameters_ES.txt') or (input_file_name == 'orca_parameters_ES.txt'):
            calculation_type = 'excited_structure'
        else:
            raise Exception('Error here.')
        job_completion_stage, re_details = analyse_RE_output(software_type, path, calculation_type)
    elif job_type == 'EET':
        path_to_output_file = path+'/'+output_file_name
        job_completion_stage = analyse_EET_output(software_type, path_to_output_file)
    elif job_type == 'ICT':
        path_to_output_file = path+'/'+output_file_name
        job_completion_stage = analyse_ICT_output(software_type, path_to_output_file)
    else:
        raise Exception('Error. Something happened here?')

    # Third, return the type of ECCP Gaussian job is being analysed, as well as it's stage of completion. 
    return job_type, job_completion_stage, re_details

# -------------------------------------------------------------------------------

