'''
Geoffrey Weal, ECCP_reset_uncompleted_jobs.py, 18/4/22

This program is designed to reset jobs that did not complete.
'''
import os, shutil
from tqdm import tqdm
from ase.io import read, write

from ECCP.ECCP_Programs.shared_general_methods.shared_gaussian_methods import folder_contains_RE_files
from ECCP.ECCP_Programs.shared_general_methods.shared_gaussian_methods import did_gaussian_job_complete, did_gaussian_opt_job_complete
from ECCP.ECCP_Programs.shared_general_methods.shared_gaussian_methods import gaussian_temp_files_to_remove
from ECCP.ECCP_Programs.shared_general_methods.shared_gaussian_methods import remove_slurm_output_files
from ECCP.ECCP_Programs.shared_general_methods.shared_gaussian_methods import found_a_gaussian_job_that_has_run
#from ECCP.Subsidiary_Programs.can_read_data_from_checkpoint_file       import can_read_data_from_checkpoint_file

class CLICommand:
    """Will reset jobs that did not complete. Only run this program if you know all your other jobs have finished, as this program will break and also reset any jobs that are still running.
    """

    @staticmethod
    def add_arguments(parser):
        pass

    @staticmethod
    def run(args):
        Run_method()

def Run_method():
    """
    This method will reset jobs that did not complete. 

    Only run this program if you know all your other jobs have finished, as this program will break and also reset any jobs that are still running.
    """

    # First, ask the user if any jobs are running, as this program will not know if jobs are running still or ended without completing. 
    print('----------------------------------------------')
    print('This program will reset all ECCP jobs that did not complete')
    print()
    print('IMPORTANT: This program can not recognise if a Gaussian/ORCA job is still running')
    print('ONLY USE THIS PROGRAM IF YOU ARE NOT CURRENTLY RUNNING ANY ECCP GAUSSIAN/ORCA JOBS')
    print()
    while True:
        value = input("Are you happy to continue (y/[n]): ")
        if (value is None) or (value == ''):
            value = 'no'
        if   (value.lower() == 'n') or (value.lower() == 'no'):
            exit('Will exit this program without resetting.'+'\n'+'----------------------------------------------')
        elif (value.lower() == 'y') or (value.lower() == 'yes'):
            print('Will proceed to resetting uncomplete ECCP Gaussian/ORCA jobs.')
            break
        else:
            print('Input must be either:')
            print('\tN/n/NO/No/no: Do not perform resetting.')
            print('\tY/y/YES/Yes/yes: Do not perform resetting.')
            print('Try again.')
            print()
    print('----------------------------------------------')

    # Second, setup all the initial variables.
    current_path = os.getcwd()
    print('----------------------------------------------')
    print('Resetting uncompleted jobs from the root path: '+str(current_path))
    print('----------------------------------------------')

    # Third, go through each subdirectory in the parent directory. 
    original_path = os.getcwd()
    pbar = tqdm(os.walk(current_path), bar_format='')
    jobs_that_have_been_reset = []
    for root, dirs, files in pbar:

        pbar.set_description('Reset: '+str(len(jobs_that_have_been_reset))+'; Currently in: '+str(root.replace(original_path+'/','')))

        # 3.1: sort the directory just to make tidying happen in alphabetical order.
        dirs.sort()

        # 3.2: What type of calculation type are we dealing with
        contains_RE_files, RE_type = folder_contains_RE_files(root)
        if contains_RE_files:

            # 3.2.1: We are looking at a reorganisation energy calculation
            if RE_type == 'GS':

                # 3.2.1.1.1: Determine bool statements to determine if to reset files.
                reset_GS_GS_main_preopt_job = os.path.exists(root+'/eGS_gGS_main_opt_preopt.log') and not did_gaussian_opt_job_complete(root+'/eGS_gGS_main_opt_preopt.log')[0]
                reset_GS_GS_main_opt_job    = os.path.exists(root+'/eGS_gGS_main_opt.log')        and not did_gaussian_opt_job_complete(root+'/eGS_gGS_main_opt.log')[0]
                reset_GS_GS_freq_job        = os.path.exists(root+'/eGS_gGS_freq.log')            and not did_gaussian_job_complete    (root+'/eGS_gGS_freq.log')
                reset_GS_ES_job             = os.path.exists(root+'/eES_gGS.log')                 and not did_gaussian_job_complete    (root+'/eES_gGS.log')

                # 3.2.1.1.2: Remove the output files of any non-completed gaussian jobs, and slurm files
                if reset_GS_GS_main_preopt_job:
                    gjf_was_updated = update_gif_file_from_previous_outputLOG(root, 'eGS_gGS_main_opt_preopt.log', 'eGS_gGS_main_opt_preopt.gjf')
                    rename_gaussian_output_file(root, output_name='eGS_gGS_main_opt_preopt.log', gjf_was_updated=gjf_was_updated)
                    remove_slurm_output_files(root)
                    jobs_that_have_been_reset.append((root, 'eGS_gGS_main_opt_preopt'))
                if reset_GS_GS_main_opt_job:
                    gjf_was_updated = update_gif_file_from_previous_outputLOG(root, 'eGS_gGS_main_opt.log', 'eGS_gGS_main_opt.gjf')
                    rename_gaussian_output_file(root, output_name='eGS_gGS_main_opt.log', gjf_was_updated=gjf_was_updated)
                    remove_slurm_output_files(root)
                    jobs_that_have_been_reset.append((root, 'eGS_gGS_main_opt'))
                if reset_GS_GS_freq_job:
                    remove_gaussian_file(root, output_name='eGS_gGS_freq.gjf')
                    remove_gaussian_file(root, output_name='eGS_gGS_freq.log')
                    remove_slurm_output_files(root)
                    jobs_that_have_been_reset.append((root, 'eGS_gGS_freq'))
                if reset_GS_ES_job:
                    remove_gaussian_file(root, output_name='eES_gGS.gjf')
                    remove_gaussian_file(root, output_name='eES_gGS.log')
                    remove_slurm_output_files(root)
                    jobs_that_have_been_reset.append((root, 'eES_gGS'))

            elif RE_type == 'ES':

                # 3.2.1.2.1: Determine bool statements to determine if to reset files.
                reset_ES_ES_main_preopt_job = os.path.exists(root+'/eES_gES_main_opt_preopt.log') and not did_gaussian_opt_job_complete(root+'/eES_gES_main_opt_preopt.log')[0]
                reset_ES_ES_main_opt_job    = os.path.exists(root+'/eES_gES_main_opt.log')        and not did_gaussian_opt_job_complete(root+'/eES_gES_main_opt.log')[0]
                reset_ES_ES_freq_job        = os.path.exists(root+'/eES_gES_freq.log')            and not did_gaussian_job_complete    (root+'/eES_gES_freq.log')
                reset_ES_GS_job             = os.path.exists(root+'/eGS_gES.log')                 and not did_gaussian_job_complete    (root+'/eGS_gES.log')

                # 3.2.1.2.2: Did this excited state reorganisation energy finish?
                if reset_ES_ES_main_preopt_job:
                    gjf_was_updated = update_gif_file_from_previous_outputLOG(root, 'eES_gES_main_opt_preopt.log', 'eES_gES_main_opt_preopt.gjf')
                    rename_gaussian_output_file(root, output_name='eES_gES_main_opt_preopt.log', gjf_was_updated=gjf_was_updated)
                    remove_slurm_output_files(root)
                    jobs_that_have_been_reset.append((root, 'eES_gES_main_opt_preopt'))
                if reset_ES_ES_main_opt_job:
                    gjf_was_updated = update_gif_file_from_previous_outputLOG(root, 'eES_gES_main_opt.log', 'eES_gES_main_opt.gjf')
                    rename_gaussian_output_file(root, output_name='eES_gES_main_opt.log', gjf_was_updated=gjf_was_updated)
                    remove_slurm_output_files(root)
                    jobs_that_have_been_reset.append((root, 'eES_gES_main_opt'))
                if reset_ES_ES_freq_job:
                    remove_gaussian_file(root, output_name='eES_gES_freq.gjf')
                    remove_gaussian_file(root, output_name='eES_gES_freq.log')
                    remove_slurm_output_files(root)
                    jobs_that_have_been_reset.append((root, 'eES_gES_freq'))
                if reset_ES_GS_job:
                    remove_gaussian_file(root, output_name='eGS_gES.gjf')
                    remove_gaussian_file(root, output_name='eGS_gES.log')
                    remove_slurm_output_files(root)
                    jobs_that_have_been_reset.append((root, 'eGS_gES'))

            else:
                raise Exception('huh?')

            # 3.2.2: Clean up the temp files while we are at it for any ECCP reorganisation energy calcs, completed or uncompleted.
            gaussian_temp_files_to_remove(root, files, remove_chk_file=False, remove_fort7_file=True, print_to_display=False) 

            # 3.2.3: Do not need to move further down the subdirectories anymore, remove all dirs and files lists.
            dirs[:] = []
            files[:] = []

        if found_a_gaussian_job_that_has_run(root, files): 

            # 3.3.1: We are looking at a non-reorganisation energy ECCP Gaussian calculation.
            #        If the output.log file shows that the program finished successfully, remove all temp files. 
            if not did_gaussian_job_complete(root+'/output.log'):
                #print('Resetting: '+str(root))
                
                # 3.3.1.1: Remove temp files as well as any results files like output.log files
                remove_gaussian_file(root, output_name='output.log')
                remove_slurm_output_files(root)
                jobs_that_have_been_reset.append((root, 'output'))

            # 3.3.2: Clean up the temp files while we are at it for any ECCP calcs, completed or uncompleted.
            gaussian_temp_files_to_remove(root, files, remove_chk_file=False, remove_fort7_file=False, print_to_display=False)

            # 3.3.3: Do not need to move further down the subdirectories anymore, remove all dirs and files lists.
            dirs[:] = []
            files[:] = []

        pbar.set_description('Reset: '+str(len(jobs_that_have_been_reset))+'; Currently in: '+str(root.replace(original_path+'/','')))

    # Fourth, print out which jobs have finished. 
    print('----------------------------------------------')
    if len(jobs_that_have_been_reset) == 0:
        print('No Jobs were reset')
    else:
        print('The following jobs were reset:')
        print()
        for job, log_name in jobs_that_have_been_reset:
            print(str(job)+': '+str(log_name))
    print('----------------------------------------------')

# --------------------------------------------------------------------------------------------------

old_suffix = 'old'
def update_gif_file_from_previous_outputLOG(root, output_name, previous_gjf_name):
    """
    This method will update the gif file will the previous self-consistancy field completed geometric step.

    Parameters
    ----------
    root : str
        This is the path to the gaussian output file.
    output_name : str.
        This is the name of the output file.
    previous_gjf_name str.
        This is the name of the gjf file that produced the output file given by output_name.

    Returns
    -------
    Return True if this method updated the gjf file, False if not.
    """
    
    # First, obtain the structure of the previously completed gemoetry step.
    try:
        previously_completed_geometric_step = read(root+'/'+output_name,index=-1)
    except:
        print('Log file does not contain a starting configuration.')
        print('This is likely because Gaussian had only just begun for a few secnods before being cancelled.')
        print('The gjf file will not need to be updated.')
        return False

    # Second, obtain the initial lines of the gif file or the previous gaussian job.
    gjf_input_lines_start = []
    gjf_input_lines_end   = []
    old_checkpoint_filename = None
    checkpoint_filename = None
    with open(root+'/'+previous_gjf_name,'r') as previousGJF:
        atomic_positions_number_of_blank_lines = 2
        no_of_blank_lines = 0
        found_second_blank_line_first_line = True
        for line in previousGJF:
            line = line.rstrip()
            # Read if there is already a line to read input from the checkpoint.
            if ('# Geom=Check Guess=Read' in line) or ('# Geom=Check Guess=TCheck' in line): 
                line = '# Geom=Check Guess=TCheck ! Will read in the geometry and electronic details from the checkpoint file'
            # Record lines from gjf about how to run the gaussian job (link0).
            if   no_of_blank_lines < atomic_positions_number_of_blank_lines:
                gjf_input_lines_start.append(line)
            elif no_of_blank_lines > atomic_positions_number_of_blank_lines:
                gjf_input_lines_end.append(line)
            elif found_second_blank_line_first_line:
                gjf_input_lines_start.append(line)
                found_second_blank_line_first_line = False
            # Determine where the breaks between the link0, title, and atom positions lines are.
            if (line == "") or line.isspace():
                no_of_blank_lines += 1

    # Third, rename the previous gjf file as a old file.
    # 3.1: Get all the names of the previous gjf files.
    previous_gjf_names = []
    for file in os.listdir(root):
        if os.path.isfile(root+'/'+file) and file.startswith(previous_gjf_name) and (not file == previous_gjf_name):
            if 'original' in file:
                continue
            previous_gjf_names.append(file)

    # 3.2: Check that there is a consecutive order of numbering in previous_gjf_names
    previous_gjf_old_numbers = sorted([int(x.replace(previous_gjf_name+'.'+old_suffix,'')) for x in previous_gjf_names])
    if not (previous_gjf_old_numbers == list(range(1,len(previous_gjf_old_numbers)+1))):
        raise Exception('huh?')

    # 3.3: Rename gif file to next old file
    os.rename(root+'/'+previous_gjf_name, root+'/'+previous_gjf_name+'.'+old_suffix+str(len(previous_gjf_old_numbers)+1))

    # Fourth, create the new gjf file with updated geometric structure.
    with open(root+'/'+previous_gjf_name,'w') as currentGJF: 
        currentGJF.write('\n'.join(gjf_input_lines_start)+'\n')
        for atom in previously_completed_geometric_step:
            currentGJF.write(str(atom.symbol)+'\t'+str(atom.x)+'\t'+str(atom.y)+'\t'+str(atom.z)+'\n')
        currentGJF.write('\n'+'\n'.join(gjf_input_lines_end)+'\n\n')

    # Fifth, return True, as this method has updated the gjf file
    return True

# --------------------------------------------------------------------------------------------------

def rename_gaussian_output_file(root, output_name='output.log', gjf_was_updated=True):
    """
    This method will remove all files that were created during the Gaussian calculation run by the slurm job.

    Parameters
    ----------
    root : str
        This is the path to the gaussian output file.
    output_name : str.
        This is the name of the output file 
    """

    if not gjf_was_updated:
        os.remove(root+'/'+output_name)
        return 
        '''
        if '.' in output_name:
            output_name = output_name.split('.')
            output_name.insert(-1,'nogeometryupdate')
            output_name = '.'.join(output_name)
        else:
            output_name+'.'+'nogeometryupdate'
        '''

    if output_name in os.listdir(root):

        # 3.1: Get all the names of the previous log files.
        previous_log_names = []
        for file in os.listdir(root):
            if os.path.isfile(root+'/'+file) and file.startswith(output_name) and (not file == output_name):
                previous_log_names.append(file)

        # 3.2: Check that there is a consecutive order of numbering in previous_log_names
        previous_log_old_numbers = sorted([int(x.replace(output_name+'.'+old_suffix,'')) for x in previous_log_names])
        if not (previous_log_old_numbers == list(range(1,len(previous_log_old_numbers)+1))):
            raise Exception('huh?')

        # 3.3: Rename log file to next old file
        os.rename(root+'/'+output_name, root+'/'+output_name+'.'+old_suffix+str(len(previous_log_old_numbers)+1))

def remove_gaussian_file(root, output_name='output.log'):
    """
    This method will remove all files that were created during the Gaussian calculation run by the slurm job.

    Parameters
    ----------
    root : str
        This is the path to the gaussian output file.
    output_name : str.
        This is the name of the output file 
    """
    if output_name in os.listdir(root):
        os.remove(root+'/'+output_name)

# --------------------------------------------------------------------------------------------------



