'''
determine_what_the_job_is.py, Geoffrey Weal, 29/12/22

This method is designed to determine what ECCP Gaussian job we are analysing
'''
import os

def determine_what_the_job_is(path, software_type, output_file_name, input_file_name, break_string='Input orientation:'):
    """
    This method is designed to determine the type of ECCP job that we are examining. 

    Parameters
    ----------
    path : str.
        This is the path to the gjf file and the output files.
    software_type : str.
        This is the type of software that was used to perform these calculations.
    output_file_name : str.
        This is the name of the output file.
    input_file_name : str.
        This is the name of the input file.
    break_string : str.
        This is the input line from the output file to break out of the for loop on.

    Returns
    -------
    job_type : str. or None
        This is the type of Gaussian job that is being performed, base on the files in the folder and the Gaussian log file. ATC: Atomic Transition Charge; RE: Reorganisation Energy; EET: Electronic Energy Transfer; ICT: Intermolecular Charge Transfer; None: None of the following job types were recognised in the gaussian log file.
    extra_info : str. or None
        This is any extra info that will be required. 
    """

    # First, initialise the variables.
    job_type = None

    # Second, if the output file doesn't exist, return false to all.
    path_to_output = path+'/'+output_file_name
    path_to_input  = path+'/'+input_file_name
    if   os.path.exists(path_to_output):
        path_to_file = path_to_output
    elif os.path.exists(path_to_input):
        path_to_file = path_to_input
    elif os.path.exists(path+'/'+'gaussian_parameters_ES.txt'):
        return 'RE'
    elif os.path.exists(path+'/'+'orca_parameters_ES.txt'):
        return 'RE'
    else:
        return job_type

    if software_type == 'Gaussian':

        # Second, if there are files called 'GS_GS.gjf' or 'ES_ES.gjf' in this folder, we are dealing with a reorganisation energy file in Gaussian.
        if input_file_name in ['eGS_gGS_main_opt_preopt.gjf', 'eGS_gGS_main_opt.gjf', 'eES_gES_main_opt_preopt.gjf', 'eES_gES_main_opt.gjf']:
            job_type = 'RE'
            looking_at_RE_output = input_file_name.split('_')[0]
            return job_type 

        # Third, look through the output file to determine what type of calculation is being performed. 
        with open(path_to_file,'r') as outputLOG:
            path_including_hash = False
            for line in outputLOG:
                # Get input of line from output file
                if   path_including_hash:
                    if '----------------------' in line:
                        path_including_hash = False
                    else:
                        current_line += line
                        continue
                elif '#' in line:
                    current_line  = line
                    path_including_hash = True
                    continue
                else:
                    current_line = line
                # Perform task with input given
                if 'density=(transition=1)' in current_line:
                    job_type = 'ATC'
                    break
                if 'eet' in current_line:
                    job_type = 'EET'
                    break

                if break_string == 'empty': # huh? this line doesnt make sense? Geoff
                    if line.strip():
                        break
                else:
                    if (break_string in current_line):
                        break

        # Fourth, return the output that determines what 
        return job_type

    if software_type == 'ORCA':

        # Second, if there are files called 'GS_GS.inp' or 'ES_ES.inp' in this folder, we are dealing with a reorganisation energy file in ORCA.
        if input_file_name in ['eGS_gGS_main_preopt.inp', 'eGS_gGS_main_opt.inp', 'eES_gES_main_preopt.inp', 'eES_gES_main_opt.inp']:
            job_type = 'RE'
            looking_at_RE_output = input_file_name.split('_')[0]
            return job_type #, looking_at_RE_output

        # Third, look through the output file to determine what type of calculation is being performed. 
        with open(path_to_file,'r') as outputLOG:
            path_including_hash = False
            for line in outputLOG:
                # Get input of line from output file
                raise Exception('Need to create this once you know what EET and ATC files look like.')

        # Fourth, return the output that determines what 
        return job_type

# ============================================================================================================



