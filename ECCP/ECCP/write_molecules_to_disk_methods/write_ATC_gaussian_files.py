"""
write_ATC_gaussian_files.py, Geoffrey Weal, 8/5/22

This script is designed to write the gaussian files and submit.sl files required for performing Gaussian jobs for performing atomic transition charge (ATC) calculations.
"""
from copy                                                                          import deepcopy
from SUMELF                                                                        import make_folder
from SUMELF                                                                        import check_molecule_against_file
from ECCP.ECCP.write_molecules_to_disk_methods.write_methods.gaussian_modified_ATC import write_gaussian_in as write_gaussian_in_ATC
from ECCP.ECCP.write_molecules_to_disk_methods.shared_methods                      import change_folder_name_components, input_commands_for_multiwfn
from ECCP.ECCP.write_molecules_to_disk_methods.shared_methods                      import slurmSL_header, load_gaussian_programs, make_gaussian_temp_folder, remove_gaussian_temp_files

def write_ATC_gaussian_files(molecule, molecule_name, environment_about_molecule, SolventsList, gaussian_jobs_path, calc_parameters_for_ATCs, submission_information_for_ATCs):
	"""
	This method will write information the Gaussian files to disk.

	Parameters
	----------
	molecule : ase.Atoms.
		This is the molecule. 
	molecule_name : str.
		This is the name of the molecule. 
	environment_about_molecule : ase.Atoms or None
		This is the environment surrounding the molecule. If None, no environment was given.
	SolventsList : list of int
		These are the indices of the molecules in the molecules list that have been identified at solvents.
	gaussian_jobs_path : str.
		This is the path to save gaussian jobs to.
	calc_parameters_for_ATCs : list
		This dictionary contain all the information required for the Gaussian input ATC file.
	submission_information_for_ATCs : list
		This dictionary contain all the information required for the submit.sl script. 
	"""

	# First, make a copy of the gaussian_parameters and submission_information dictionaries.
	gaussian_parameters    = deepcopy(calc_parameters_for_ATCs)
	submission_information = deepcopy(submission_information_for_ATCs)
	del gaussian_parameters['calc_software']

	# Second, determine if some critical tags that are needed are in the submission_information dictionary. 
	got_cpu  = 'cpus_per_task' in submission_information
	got_mem  = 'mem' in submission_information
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

	# Fifth, provide the name and filepath for each of the scratch files.
	scratch_dir_given = ('temp_folder_path' in gaussian_parameters) # was gaussian_scratch_name
	if scratch_dir_given:
		temp_folder_path = gaussian_parameters['temp_folder_path']+'/'+calc_folder
		del gaussian_parameters['temp_folder_path']
		submission_information['temp_folder_path'] = temp_folder_path

	# =============================================================================
	# Sixth, for those temporary files that I have control over where they get saved to, 
	# indicate to the .gjf file where to save those files.

	# 7.1: Make the checkpoint file.
	gaussian_parameters['chk'] = 'gaussian.chk'

	# 7.2: Make path folder and file details for the other gaussian checkpoint files.
	for suffix in ['rwf','int','d2e','skr']:
		if scratch_dir_given: # A scratch path is given.
			gaussian_parameters[suffix] = temp_folder_path+'/'+'gaussian.'+str(suffix)
		else: # A scratch path has not been given.
			# Default name given called gaussian.suffix, whether scratch_dir_given is True or False
			gaussian_parameters[suffix] = 'gaussian.'+str(suffix)

	# =============================================================================

	# Seventh, write the folder to place gaussian files to.
	make_folder(calc_folder)

	# Eighth, set the path to the gaussian input file.
	path_to_gaussian_file = calc_folder+'/'+molecule_name+'.gjf'

	# Ninth, check if there already exists this molecule on file, check if the molecules are the same.
	check_molecule_against_file(molecule, path_to_gaussian_file)

	# Tenth, create the gaussian .gjf file.
	gaussian_parameters['wfn_filename'] = 'output.wfn'
	with open(path_to_gaussian_file, 'w') as fd:
		write_gaussian_in_ATC(fd, molecule, environment=environment_about_molecule, molecule_name=molecule_name, **gaussian_parameters)
	del gaussian_parameters['wfn_filename']

	# Eleventh, create the submit.sl file to submit this gaussian job to slurm.
	make_ATC_gaussian_submitSL(molecule_name+'.gjf',calc_folder,functional,basis_set,gaussian_parameters,**submission_information)

# ------------------------------------------------------------------------------------------------------------------------------

def make_ATC_gaussian_submitSL(gaussian_filename,local_path,functional,basis_set,gaussian_parameters,cpus_per_task,mem,time,partition='parallel',constraint=None,nodelist=None,exclude=None,email='',python_version='python/3.8.1',gaussian_version='gaussian/g16',log_filename='output.log',wfn_filename='output.wfn',temp_folder_path=None,remove_chk_file=True):
	"""
	This method will write the submit.sl file in parallel

	Parameters
	----------
	gaussian_filename : str. 
		This is the name of the gaussian file.
	local_path : str. 
		This is the location to save this submit.sl file to
	functional : str. 
		This is the functional you are going to use in your Gaussian calculation.
	basis_set : str. 
		This is the basis set you are going to use in your Gaussian calculation.
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
	log_filename : str.
		This is the name of the Gaussian output file. Default: 'output.log'
	wfn_filename : str.
		This is the name of the Gaussian Wavefunction file. Default: 'output.wfn'
	temp_folder_path : str. or None
		This is the path to the scratch directory to save Gaussian temp files to. If you dont give this, Gaussian temp files will be saves to the default scratch directory. Default: None
	remove_chk_file : bool.
		This variable indicates if you want to remove the chk file afterwards. Default: False
	"""

	# First, create name for job.
	name = '-'.join(local_path.split('/')[-3:])+'-ATC'

	# Second, write the submit.sl script
	with open(local_path+'/'+"submit.sl", "w") as submitSL:
		slurmSL_header(submitSL, name, mem, partition, constraint, nodelist, time, email, cpus_per_task=cpus_per_task, exclude=exclude)
		make_gaussian_temp_folder(submitSL, temp_folder_path)
		load_gaussian_programs(submitSL, gaussian_version, python_version)
		submitSL.write('# ============================\n')
		submitSL.write('# Prevent the ATC job from running if it is already running or has already run.\n')
		submitSL.write('\n')
		submitSL.write('if ! [[ -f '+str(log_filename)+' ]]\n')
		submitSL.write('then\n')
		submitSL.write('\t\n')
		submitSL.write('\t# ----------------------------\n')
		submitSL.write('\t# Perform Gaussian Calculation\n')
		submitSL.write('\t\n')
		submitSL.write('\tsrun '+str(gaussian_version.split('/')[-1])+' < '+str(gaussian_filename)+' > '+str(log_filename)+'\n')
		submitSL.write('\t\n')
		submitSL.write('\t# ----------------------------\n')
		remove_gaussian_temp_files(submitSL, gaussian_parameters, temp_folder_path, remove_chk_file=remove_chk_file, remove_temp_folder=True, prefix='\t')
		submitSL.write('\t# ----------------------------\n')
		submitSL.write('\t# Performing Multiwfn Calculation\n')
		submitSL.write('\t\n')
		submitSL.write('\tif [[ -f "'+str(wfn_filename)+'" ]]\n')
		submitSL.write('\tthen\n')
		submitSL.write('\t\techo "found '+str(wfn_filename)+', will perform Multiwfn ATC calculation"\n')
		submitSL.write('\t\tsubmit_slurm_job.py multiwfn_submit.sl\n')
		#submitSL.write('\t\t'+input_commands_for_multiwfn(wfn_filename=wfn_filename)+'\n')
		submitSL.write('\telse\n')
		submitSL.write('\t\techo "'+str(wfn_filename)+' file not found. Will not perform Multiwfn ATC calculation"\n')
		submitSL.write('\t\techo "'+str(wfn_filename)+' file not found. Will not perform Multiwfn ATC calculation" 1>&2\n')
		submitSL.write('\tfi\n')
		submitSL.write('\t\n')
		submitSL.write('\t# ----------------------------\n')
		submitSL.write('\techo "End of job"\n')
		submitSL.write('\t# ----------------------------\n')
		submitSL.write('\t\n')
		submitSL.write('fi\n')
		
# ------------------------------------------------------------------------------------------------------------------------------


