"""
write_RE_FC_gaussian_files.py, Geoffrey Weal, 8/5/22

This script is designed to write the gaussian files and submit.sl files required for performing Gaussian jobs for performing reorganisation energy (RE) calculations.
"""
from copy   import deepcopy
from SUMELF import make_folder

from ECCP.ECCP.write_molecules_to_disk_methods.shared_methods import change_folder_name_components

from ECCP.ECCP.write_molecules_to_disk_methods.write_RE_gaussian_files_methods.make_initial_gaussian_RE_optimisation_gjf_file import make_initial_gaussian_RE_optimisation_gjf_file
from ECCP.ECCP.write_molecules_to_disk_methods.write_RE_gaussian_files_methods.make_RE_gaussian_submitSL_preopt               import make_RE_gaussian_submitSL_preopt
from ECCP.ECCP.write_molecules_to_disk_methods.write_RE_gaussian_files_methods.make_RE_gaussian_submitSL                      import make_RE_gaussian_submitSL
from ECCP.ECCP.write_molecules_to_disk_methods.write_RE_gaussian_files_methods.make_RE_freq_gaussian_submitSL                 import make_RE_freq_gaussian_submitSL

def write_RE_gaussian_files(molecule, molecule_name, SolventsList, gaussian_jobs_path, calc_parameters_for_REs, submission_information_for_REs, run_excited_state_from_optimised_ground_structure=False):
	"""
	This method will write information the Gaussian files to disk.

	Parameters
	----------
	molecule : ase.Atoms.
		This is the molecule. 
	molecule_name : str.
		This is the name of the molecule. 
	SolventsList : list of int
		These are the indices of the molecules in the molecules list that have been identified at solvents.
	gaussian_jobs_path : str.
		This is the path to save gaussian jobs to.
	calc_parameters_for_REs : list
		This dictionary contain all the information required for the Gaussian input RE file.
	submission_information_for_REs : list
		This dictionary contain all the information required for the submit.sl script. 
	run_excited_state_from_optimised_ground_structure : bool.
		This boolean indicates if you want to run the excited state calculation from the optimised ground state calculation for reorganisation energy calculations. True if you do, False if you want to run the excited state calculation from the original structure (Default: False).
	"""

	# First, make a copy of the gaussian_parameters and submission_information dictionaries.
	gaussian_parameters    = deepcopy(calc_parameters_for_REs)
	submission_information = deepcopy(submission_information_for_REs)
	del gaussian_parameters['calc_software']

	# Second, determine if some critical tags that are needed are in the submission_information dictionary. 
	got_cpu = 'cpus_per_task' in submission_information
	got_mem = 'mem' in submission_information
	got_time = 'time' in submission_information
	if not (got_cpu and got_time):
		print('Error: You need to specify the following in your submission_information dictionary:')
		if not got_cpu:
			print('\t* cpus_per_task')
		if not got_mem:
			print('\t* mem')
		if not got_time:
			print('\t* time')
		print('See https://github.com/geoffreyweal/ECCP/ for more information about these tags.')
		print('submission_information = '+str(submission_information))
		exit('This program will finish without completing.')

	# Third, copy some tag information that is in the submission_information dictionary to the gaussian_parameters dictionary.
	gaussian_parameters['nprocshared'] = submission_information['cpus_per_task']

	# Fourth, give the name of the folder to place gaussian files to.
	quantum_chemistry_program = 'GAUSSIAN'
	functional                = change_folder_name_components(gaussian_parameters['method'])
	basis_set                 = change_folder_name_components(gaussian_parameters['basis'])
	funct_and_basis_name      = 'F_'+functional+'_B_'+basis_set
	calc_folder               = str(gaussian_jobs_path)+'/'+str(molecule_name)+'/'+quantum_chemistry_program+'_'+str(funct_and_basis_name)

	# =============================================================================
	# Fifth, for those temporary files that I have control over where they get saved to, 
	# indicate to the .gjf file where to save those files.

	# 5.1: setup the folder names.
	ground_state_foldername  = 'ground_structure'
	excited_state_foldername = 'excited_structure'
	opt_calc_type  = 'opt'
	freq_calc_type = 'freq'

	# 5.2: Make a copy of the gaussian_parameters and submission_information for the ground state.
	gaussian_parameters_GS      = deepcopy(gaussian_parameters)
	gaussian_parameters_GS_freq = deepcopy(gaussian_parameters)
	gaussian_parameters_ES      = deepcopy(gaussian_parameters)
	gaussian_parameters_ES_freq = deepcopy(gaussian_parameters)

	# 5.3: Make a copy of the gaussian_parameters and submission_information for the excited state.
	submission_information_GS      = deepcopy(submission_information)
	submission_information_GS_freq = deepcopy(submission_information)
	submission_information_ES      = deepcopy(submission_information)
	submission_information_ES_freq = deepcopy(submission_information)

	# 5.4: Make the datasets to be used in 5.5
	data_sets = []
	data_sets.append((gaussian_parameters_GS,      submission_information_GS,      ground_state_foldername,  opt_calc_type,  False))
	data_sets.append((gaussian_parameters_ES,      submission_information_ES,      excited_state_foldername, opt_calc_type,  False))
	data_sets.append((gaussian_parameters_GS_freq, submission_information_GS_freq, ground_state_foldername,  freq_calc_type, True))
	data_sets.append((gaussian_parameters_ES_freq, submission_information_ES_freq, excited_state_foldername, freq_calc_type, True))
	
	# 5.5: Make the temporary folders and Gaussian files.
	for gaussian_parameters_set, submission_information_set, GorE_state_foldername, calc_type, performing_freq_calc in data_sets:

		# 5.5.1, provide the name and filepath for each of the scratch files.
		scratch_dir_given = ('temp_folder_path' in gaussian_parameters_set) # was gaussian_scratch_name
		if scratch_dir_given:
			temp_folder_path = gaussian_parameters_set['temp_folder_path']+'/'+calc_folder+'/'+GorE_state_foldername+'/'+calc_type
			del gaussian_parameters_set['temp_folder_path']
			submission_information_set['temp_folder_path'] = temp_folder_path

		# 5.5.2: Make the checkpoint file.
		gaussian_parameters_set['chk'] = 'gaussian.chk'

		# 5.5.3: Make path folder and file details for the other gaussian checkpoint files.
		for suffix in ['rwf','int','d2e','skr']:
			if scratch_dir_given: # A scratch path is given.
				gaussian_parameters_set[suffix] = temp_folder_path+'/'+'gaussian.'+str(suffix)
			else: # A scratch path has not been given.
				# Default name given called gaussian.suffix, whether scratch_dir_given is True or False
				gaussian_parameters_set[suffix] = 'gaussian.'+str(suffix)

	# =============================================================================

	# Sixth, the names for the Ground and Excited States.
	ground_structure_foldername       = 'ground_structure'
	excited_structure_foldername      = 'excited_structure'
	ground_structure_calc_folder  = calc_folder +'/'+ ground_structure_foldername
	excited_structure_calc_folder = calc_folder +'/'+ excited_structure_foldername

	# Seventh, write the folder to place gaussian files to.
	make_folder(ground_structure_calc_folder)
	make_folder(excited_structure_calc_folder)

	# Eighth, write the name of the gjf files to make
	ground_structure_GS_name          = 'eGS_gGS'
	ground_structure_GS_DFT_main_opt  = 'eGS_gGS_main_opt.gjf'
	ground_structure_ES               = 'eES_gGS.gjf'

	excited_structure_ES_name         = 'eES_gES'
	excited_structure_ES_DFT_main_opt = 'eES_gES_main_opt.gjf'
	excited_structure_GS              = 'eGS_gES.gjf'

	# =============================================================================
	# Ninth, create the gaussian .gjf file for optimising the ground and excited structures.

	# 9.1: Create the gaussian .gjf file for optimising the ground structure.
	running_pre_optimisation_ground_structure = make_initial_gaussian_RE_optimisation_gjf_file(ground_structure_calc_folder,   ground_structure_GS_DFT_main_opt,  molecule, False, molecule_name+' - GS_GS - opt', gaussian_parameters_GS)

	# 9.2: Create the gaussian .gjf file for optimising the excited structure.
	if run_excited_state_from_optimised_ground_structure:
		with open(excited_structure_calc_folder+'/'+'gaussian_parameters_ES.txt','w') as FILE:
			FILE.write(str(gaussian_parameters_ES))
	else:
		running_pre_optimisation_excited_structure = make_initial_gaussian_RE_optimisation_gjf_file(excited_structure_calc_folder, excited_structure_ES_DFT_main_opt, molecule, True,  molecule_name+' - ES_ES - opt', gaussian_parameters_ES)

	# =============================================================================
	# Tenth, create the submit .sl file for optimising the ground and excited structures.

	# 10.1: Create the submit .sl file for optimising the ground structure.
	# * Make the initial optimimsation files.
	if running_pre_optimisation_ground_structure:
		make_RE_gaussian_submitSL_preopt    (ground_structure_GS_name, ground_structure_GS_DFT_main_opt, ground_structure_ES, ground_structure_calc_folder, False,          functional, basis_set, gaussian_parameters_GS,      **submission_information_GS)
	else:
		make_RE_gaussian_submitSL           (ground_structure_GS_name, ground_structure_GS_DFT_main_opt, ground_structure_ES, ground_structure_calc_folder, False,    True, functional, basis_set, gaussian_parameters_GS,      **submission_information_GS)
	# * Perform the frequency calculation for the optimised ground state structure (set perform_TD = False). 
	make_RE_freq_gaussian_submitSL          (ground_structure_GS_name, ground_structure_calc_folder, False, True,                                                                                                               functional, basis_set, gaussian_parameters_GS_freq, **submission_information_GS_freq)

	# 10.2: Create the submit .sl file for optimising the excited structure.
	# * Make the initial optimimsation files.
	if run_excited_state_from_optimised_ground_structure:
		make_RE_gaussian_submitSL           (excited_structure_ES_name, excited_structure_ES_DFT_main_opt, excited_structure_GS, excited_structure_calc_folder, True, True,  functional, basis_set, gaussian_parameters_ES,      **submission_information_ES)
	else:
		if running_pre_optimisation_excited_structure:
			make_RE_gaussian_submitSL_preopt(excited_structure_ES_name, excited_structure_ES_DFT_main_opt, excited_structure_GS, excited_structure_calc_folder, True,        functional, basis_set, gaussian_parameters_ES,      **submission_information_ES)
		else:
			make_RE_gaussian_submitSL       (excited_structure_ES_name, excited_structure_ES_DFT_main_opt, excited_structure_GS, excited_structure_calc_folder, True, False, functional, basis_set, gaussian_parameters_ES,      **submission_information_ES)
	# * Perform the frequency calculation for the optimised excited state structure (set perform_TD = True). 
	make_RE_freq_gaussian_submitSL          (excited_structure_ES_name, excited_structure_calc_folder, True, False,                                                          functional, basis_set, gaussian_parameters_ES_freq, **submission_information_ES_freq)

	# =============================================================================

# ------------------------------------------------------------------------------------------------------------------------------
