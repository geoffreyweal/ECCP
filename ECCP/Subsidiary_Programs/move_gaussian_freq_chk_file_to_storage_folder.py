#!/usr/bin/env python3
'''
Geoffrey Weal, move_gaussian_freq_chk_file_to_storage_folder.py, 29/12/22

This program is designed to move the frequency checkpoint file from the RE folder to the RE_Checkpoint folder, once the frequency single point calculation has finished. 
'''
import os, shutil
from ECCP.ECCP_Programs.shared_general_methods.shared_general_methods import reverse_readline

All_RE_Checkpoint_Files_folder    = 'All_RE_Checkpoint_Files'
Unique_RE_Checkpoint_Files_folder = 'Unique_RE_Checkpoint_Files'

def make_folder(path_to_dir):
    # This method allow you to make a new directory if it doesnt already exist.
    if not os.path.exists(path_to_dir):
        os.makedirs(path_to_dir)

def has_freq_file_completed(path_to_outputLOG):
    # This method check if the frequency calculation has completed from the log file.
    did_terminate_normally = False
    for line in reverse_readline(path_to_outputLOG):
        # Check if the gaussian file has terminated normally
        if 'Normal termination of Gaussian' in line:
            did_terminate_normally = True
            break
        if counter >= 20 and not did_terminate_normally:
            break
    return did_terminate_normally

# First, get the path to the current directory and the overall path to the gaussian_freq.chk file.
current_path = os.getcwd()
gaussian_freq_chk_file = 'gaussian_freq.chk'
checkpoint_file_to_move_path = current_path+'/'+gaussian_freq_chk_file

# Second, check if this directory contains a gaussian_freq.chk file, and the GS_GS_freq.log or ES_ES_freq.log file has completed successfully.
if not os.path.exists(gaussian_freq_chk_file):
    raise Exception('Error: No '+str(gaussian_freq_chk_file)+' exists.'+'. Finishing without moving checkpoint file.')
GS_GS_freq_file = 'eGS_gGS_freq.log'
ES_ES_freq_file = 'eES_gES_freq.log'
calculation_type = None
if os.path.exists(GS_GS_freq_file):
    if not has_freq_file_completed(current_path+'/'+GS_GS_freq_file):
        raise Exception('Error: '+str(GS_GS_freq_file)+' indicates that the frequency calculation did not complete successfully.'+'. Finishing without moving checkpoint file.')
    calculation_type = 'eGS_gGS'
elif os.path.exists(ES_ES_freq_file):
    if not has_freq_file_completed(current_path+'/'+ES_ES_freq_file):
        raise Exception('Error: '+str(ES_ES_freq_file)+' indicates that the frequency calculation did not complete successfully.'+'. Finishing without moving checkpoint file.')
    calculation_type = 'eES_gES'
else:
    raise Exception('Error: No frequency calculation log file exists exists (either '+str(GS_GS_freq_file)+' or '+str(ES_ES_freq_file)+').'+'. Finishing without moving checkpoint file.')

# Third, obtain the path to the storage folder to move the gaussian_freq.chk file to.
current_path_segments = current_path.split('/')
overall_RE_Gaussian_Jobs_found = False
for index in range(len(current_path_segments)):
    if current_path_segments[index] == 'All_RE_Gaussian_Jobs':
        current_path_segments[index] = str(All_RE_Checkpoint_Files_folder)
        overall_RE_Gaussian_Jobs_found = True
    if current_path_segments[index] == 'Unique_RE_Gaussian_Jobs':
        current_path_segments[index] = str(Unique_RE_Checkpoint_Files_folder)
        overall_RE_Gaussian_Jobs_found = True
if not overall_RE_Gaussian_Jobs_found:
    raise Exception('Error: Current path not in the "Unique_RE_Gaussian_Jobs" or "All_RE_Gaussian_Jobs" folder. Current Directory: '+str(current_path)+'. Finishing without moving checkpoint file.')
storage_directory = '/'.join(current_path_segments)

# Fourth, remove the 'ground_structure' or 'excited_structure' from our storage folder path and include thr GS or ES name in the checkpoint filename.
if calculation_type == 'eGS_gGS':
    folder_structure_name = 'ground_structure'
    new_checkpoint_filename = 'eGS_gGS_gaussian.chk'
elif calculation_type == 'eES_gES':
    folder_structure_name = 'excited_structure'
    new_checkpoint_filename = 'eES_gES_gaussian.chk'
else:
    raise Exception('Huh?')
if not storage_directory.endswith(folder_structure_name):
    raise Exception('Error: This checkpoint file should be stored in a folder called "'+str(folder_structure_name)+'". Not sure why it is not. Finishing without moving checkpoint file.')
storage_directory = storage_directory[:len(storage_directory)-len(folder_structure_name)]
if storage_directory.endswith('/'):
    storage_directory = storage_directory[:-1]

# Fifth, check that a frequency checkpoint file does not already exist. 
move_the_checkpoint_file_to_path = storage_directory+'/'+new_checkpoint_filename
if os.path.exists(move_the_checkpoint_file_to_path):
    to_string  = 'Error: '+str(move_the_checkpoint_file_to_path)+' already exists. We do not want to overwrite anything, so will finish here without moving anything.'+'\n'
    to_string += 'Checkpoint file you want to move: '+str(checkpoint_file_to_move_path)
    raise Exception(to_string)

# Sixth, create the storage_directory if it doesn't already exist.
make_folder(storage_directory)

# Seventh, move the checkpoint file from current_path to storage_directory
#          and change the name from gaussian_freq_chk_file to new_checkpoint_filename
#          Moving: checkpoint_file_to_move_path   --->   move_the_checkpoint_file_to_path
print('Moving: '+str(checkpoint_file_to_move_path)+'  -->  '+str(move_the_checkpoint_file_to_path))
shutil.move(checkpoint_file_to_move_path, move_the_checkpoint_file_to_path)
print('Successfully Moved: '+str(checkpoint_file_to_move_path)+'  -->  '+str(move_the_checkpoint_file_to_path))

# ----------------------------------------------------------





