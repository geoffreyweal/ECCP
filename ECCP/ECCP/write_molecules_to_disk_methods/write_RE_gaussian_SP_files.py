"""
write_RE_FC_gaussian_SP_files.py, Geoffrey Weal, 8/5/22

This script is designed to create the single point Gaussian calculation files (.gjf files) for performing reorganisation energy (RE) calculations.
"""
from copy                                                     import deepcopy
from SUMELF                                                   import make_folder
from ECCP.ECCP.write_molecules_to_disk_methods.shared_methods import change_folder_name_components, convert_dict_for_bash_input
from ECCP.ECCP.write_molecules_to_disk_methods.shared_methods import slurmSL_header, load_gaussian_programs, make_gaussian_temp_folder, remove_gaussian_temp_files

def write_RE_gaussian_SP_files(molecule, molecule_name, SolventsList, gaussian_jobs_path, calc_parameters_for_RE_SPs, submission_information_for_RE_SPs):
	"""
	This method will write information the Gaussian files to disk for performing reorganisation energy calculations.

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
	calc_parameters_for_RE_SPs : list
		This dictionary contain all the information required for the Gaussian input single point RE file.
	submission_information_for_RE_SPs : list
		This dictionary contain all the information required for the submit.sl script. 
	"""

	# First, make a copy of the gaussian_parameters and submission_information dictionaries.
	gaussian_parameters    = deepcopy(calc_parameters_for_RE_SPs)
	submission_information = deepcopy(submission_information_for_RE_SPs)
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

	# 5.1: Setup the folder names.
	ground_state_foldername  = 'ground_structure'
	excited_state_foldername = 'excited_structure'
	SP_calc_type             = 'SP'

	# 5.2: Make a copy of the gaussian_parameters and submission_information for gound and excited states
	gaussian_parameters_GS = deepcopy(gaussian_parameters)
	gaussian_parameters_ES = deepcopy(gaussian_parameters)

	submission_information_GS = deepcopy(submission_information)
	submission_information_ES = deepcopy(submission_information)

	# 5.3: Make the datasets to be used in 5.4
	data_sets = []
	data_sets.append((gaussian_parameters_GS, submission_information_GS, ground_state_foldername))
	data_sets.append((gaussian_parameters_ES, submission_information_ES, excited_state_foldername))

	# 5.4: Make the temporary folders and Gaussian files.
	for gaussian_parameters_set, submission_information_set, GorE_state_foldername in data_sets:

		# 5.4.1, provide the name and filepath for each of the scratch files.
		scratch_dir_given = ('temp_folder_path' in gaussian_parameters_set) # was gaussian_scratch_name
		if scratch_dir_given:
			temp_folder_path = gaussian_parameters_set['temp_folder_path']+'/'+calc_folder+'/'+GorE_state_foldername+'/'+SP_calc_type
			del gaussian_parameters_set['temp_folder_path']
			submission_information_set['temp_folder_path'] = temp_folder_path

		# 5.4.2: Make the checkpoint file.
		gaussian_parameters_set['chk'] = 'gaussian.chk'

		# 5.4.3: Make path folder and file details for the other gaussian checkpoint files.
		for suffix in ['rwf','int','d2e','skr']:
			if scratch_dir_given: # A scratch path is given.
				gaussian_parameters_set[suffix] = temp_folder_path+'/'+'gaussian.'+str(suffix)
			else: # A scratch path has not been given.
				# Default name given called gaussian.suffix, whether scratch_dir_given is True or False
				gaussian_parameters_set[suffix] = 'gaussian.'+str(suffix)
	# =============================================================================

	# Sixth, get the names for the Ground and Excited States
	ground_structure_foldername       = 'ground_structure'
	excited_structure_foldername      = 'excited_structure'
	ground_structure_calc_folder  = calc_folder +'/'+ ground_structure_foldername
	excited_structure_calc_folder = calc_folder +'/'+ excited_structure_foldername

	# Seventh, write the folder to place gaussian files to.
	make_folder(ground_structure_calc_folder)
	make_folder(excited_structure_calc_folder)

	# Eighth, write the name of the gjf files to make
	ground_structure_GS  = 'eGS_gGS_main_opt.gjf'
	ground_structure_ES  = 'eES_gGS.gjf'
	excited_structure_ES = 'eES_gES_main_opt.gjf'
	excited_structure_GS = 'eGS_gES.gjf'

	# Ninth, create the submit .sl file for optimised ground structure.
	make_RE_gaussian_submitSL(ground_structure_GS,  ground_structure_ES,  ground_structure_calc_folder,   True, functional, basis_set, gaussian_parameters_GS, **submission_information_GS)

	# Tenth, create the submit .sl file for optimised excited structure.
	make_RE_gaussian_submitSL(excited_structure_ES, excited_structure_GS, excited_structure_calc_folder, False, functional, basis_set, gaussian_parameters_ES, **submission_information_ES)

# ------------------------------------------------------------------------------------------------------------------------------

def make_RE_gaussian_submitSL(optimisation_filename,single_point_filename,local_path,perform_TD,functional,basis_set,gaussian_parameters,cpus_per_task,mem,time,partition='parallel',constraint=None,nodelist=None,exclude=None,email='',python_version='python/3.8.1',gaussian_version='gaussian/g16',temp_folder_path=None):
	"""
	This method will write the individual submit.sl files so all Gaussian jobs can be run individually in slurm (in 'parallel'). 

	Parameters
	----------
	optimisation_filename : str. 
		This is the name of the optimisation file 
	single_point_filename : str. 
		This is the name of the single point calculation file 
	local_path : str. 
		This is the location to save this submit.sl file to
	perform_TD : bool.
		This boolean indicates if you are running TD calculations.
	functional : str. 
		This is the functional you are going to use in your Gaussian calculation.
	basis_set : str. 
		This is the basis set you are going to use in your Gaussian calculation.
	gaussian_parameters : dict.
		This dictionary contains all the input parameters required for creating the gaussian input (.gjf) file.
	cpus_per_task : int
		This is the number of cpus you want to use for Gaussian jobs.
	mem : str.
		This is the amount of memory you want to use for Gaussian jobs.
	time : str.
		This is the amount of time you want to use for Gaussian jobs.
	partition : str.
		This is the partition to run this job on. Default: 'parallel'
	constraint : str.
		This is the slurm constraint. If you dont give this, this wont be set. Default: None
	nodelist : str.
		This is the slurm nodelist. If you dont give this, this wont be set. Default: None
	exclude : str.
		This is the slurm exclude nodes list. If you dont give this, this wont be set. Default: None
	email : str.
		This is the email to email about how this job is going. If you dont give this, this wont be set. Default: ''
	python_version : str.
		This is the version of python you want to load/use in slurm. Default: 'python/3.8.1'
	gaussian_version : str.
		This is the version of Gaussian you want to load/use in slurm. Default: 'gaussian/g16'
	temp_folder_path : str.
		This is the path to the scratch directory to save Gaussian temp files to. If you dont give this, Gaussian temp files will be saves to the default scratch directory. Default: None
	"""

	# Get version of Gaussian to use
	gaussian_version_suffix = str(gaussian_version.split('/')[-1])
	
	# Make changes to gaussian_parameters to change gaussian filenames to include "freq" in their name.
	gaussian_parameters_for_SP = deepcopy(gaussian_parameters)
	for gaussian_file in ['chk', 'd2e', 'int', 'rwf', 'skr']:
		if gaussian_file in gaussian_parameters_for_SP:
			gaussian_filename = gaussian_parameters_for_SP[gaussian_file]
			gaussian_filename = gaussian_filename.split('.')
			gaussian_filename = gaussian_filename[0]+'_sp.'+gaussian_filename[1]
			gaussian_parameters_for_SP[gaussian_file] = gaussian_filename

	# Remove pre_method and pre-basis if they exist on gaussian_parameters
	if ('pre_method' in gaussian_parameters_for_SP):
		del gaussian_parameters_for_SP['pre_method']
	if ('pre_basis'  in gaussian_parameters_for_SP):
		del gaussian_parameters_for_SP['pre_basis']

	# Remove the entry for oldchk as we do not want to use a previous checkpoint file for this calculation.
	if 'oldchk' in gaussian_parameters_for_SP:
		del gaussian_parameters_for_SP['oldchk']

	# If not performing an excited state calculation, remove the 'td_settings' entry so there is not confusion in the file. 
	if (not perform_TD) and ('td_settings' in gaussian_parameters_for_SP):
		del gaussian_parameters_for_SP['td_settings']

	# create name for job
	optimisation_name = optimisation_filename.replace('.gjf','')
	single_point_name = single_point_filename.replace('.gjf','')
	name = '-'.join(local_path.split('/')[-4:-1])+'-ReorgE_SP-'+str(single_point_name)

	# writing the submit.sl script
	with open(local_path+'/'+str(single_point_name)+"_submit.sl", "w") as submitSL:
		slurmSL_header(submitSL, name, mem, partition, constraint, nodelist, time, email, cpus_per_task=cpus_per_task, exclude=exclude)
		make_gaussian_temp_folder(submitSL, temp_folder_path)
		load_gaussian_programs(submitSL, gaussian_version, python_version)
		submitSL.write('# ============================\n')
		submitSL.write('# Prevent the single point job from running if it is already running or has already run.\n')
		submitSL.write('\n')
		submitSL.write('if ! [[ -f '+str(single_point_name)+'.log'+' ]]\n')
		submitSL.write('then\n')
		submitSL.write('\t\n')
		submitSL.write('\t# ============================\n')
		submitSL.write('\t# Settings for creating new gjf files\n')
		submitSL.write('\t\n')
		submitSL.write('\tgaussian_parameters='+str(convert_dict_for_bash_input(gaussian_parameters_for_SP))+'\n')
		submitSL.write('\t\n')
		submitSL.write('\t# ============================\n')
		submitSL.write('\t# Extract optimised structure and place it into single point calculation gaussian input (.gjf) file.\n')
		submitSL.write('\t\n')
		submitSL.write('\tget_single_point_RE_Gaussian_input_file.py '+str(optimisation_name)+'.log '+str(single_point_name)+'.gjf '+str(perform_TD)+' "${gaussian_parameters}"\n')
		submitSL.write('\t\n')
		submitSL.write('\t# ============================\n')
		submitSL.write('\t# Performing Gaussian Calculation\n')
		submitSL.write('\t\n')
		submitSL.write('\tsrun '+gaussian_version_suffix+' < '+str(single_point_filename)+' > '+str(single_point_name)+'.log\n')
		submitSL.write('\t\n')
		submitSL.write('\t# ============================\n')
		remove_gaussian_temp_files(submitSL, gaussian_parameters_for_SP, temp_folder_path, remove_chk_file=True, remove_temp_folder=True, prefix='\t')
		submitSL.write('\t# ============================\n')
		submitSL.write('\techo "End of job"\n')
		submitSL.write('\t# ============================\n')
		submitSL.write('\t\n')
		submitSL.write('fi\n')

# ------------------------------------------------------------------------------------------------------------------------------





