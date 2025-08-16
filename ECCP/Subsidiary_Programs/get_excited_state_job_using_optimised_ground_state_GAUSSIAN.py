#!/usr/bin/env python3
"""
get_excited_state_job_using_optimimsed_ground_state_GAUSSIAN.py, Geoffrey Weal, 25/8/24

This program is designed to obtain the optimised ground state structure and prepare it as the starting structure for the excited state calculation. 
"""
import os, sys
from ase.io import read
from ECCP.ECCP.write_molecules_to_disk_methods.write_RE_gaussian_files_methods.make_initial_gaussian_RE_optimisation_gjf_file import make_initial_gaussian_RE_optimisation_gjf_file

# First, indicate to the user what is happening to the terminal
print('Creating input file for eES_gES file using the optimised ground state.')

# Second, obtain the path to the optimised ground state structure, given by the first input
path_to_ground_state_structure = sys.argv[1]

# Third, usually the molecule name given here would be its actual name. However, as this is just for the title bar of the gjf file, we will just called it molecule.
molecule_name 'molecule' #sys.args[2]

# Third obtain the "gaussian_parameters_ES" dictionary from the file in the "excited_structure" folder. 
with open('../excited_structure/gaussian_parameters_ES.txt', 'r') as FILE:
	gaussian_parameters_ES = FILE.read()

# Fourth, take the raw file and evaluated it into a dictinoary containing all the information about the gaussian inputs for performing the excited state calculation. 
gaussian_parameters_ES = eval(gaussian_parameters_ES)

# Fifth, read the locally optimimsed ground state structure from the output.log file. 
molecule = read(path_to_ground_state_structure, index=-1)

# Sixth, set the names for "excited_structure_ES_name" and "excited_structure_ES_DFT_main_opt" names.
excited_structure_ES_name         = 'eES_gES'
excited_structure_ES_DFT_main_opt = 'eES_gES_main_opt.gjf'

# Seventh, obtain the path to the excited state calculation by getting the current path and replacing the final "ground_structure" value to "excited_structure"
excited_structure_calc_folder = os.getcwd()[::-1].replace('ground'[::-1], 'excited'[::-1])[::-1]

# Eighth, Create the gaussian .gjf file for optimising the excited structure.
make_initial_gaussian_RE_optimisation_gjf_file(excited_structure_calc_folder, excited_structure_ES_DFT_main_opt, molecule, True,  molecule_name+' - ES_ES - opt', gaussian_parameters_ES)

# Ninth, we now have the excited state gjf file, so we do not need gaussian_parameters_ES.txt. So we can remove it
os.remove('../excited_structure/gaussian_parameters_ES.txt')

# Tenth, send a final message to the terminal saying that the excited state structure gjf file has been created, and that the excited state optimisation job will proceed. 
print('Created input file for eES_gES file using the optimised ground state.')
print('Will now submit this eES_gES file to slurm.')