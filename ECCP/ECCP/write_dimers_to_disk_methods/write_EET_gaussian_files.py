"""
write_Dimer_gaussian_files.py, Geoffrey Weal, 8/5/22

This script is designed to write the gaussian files and submit.sl files required for performing Gaussian jobs for performing electronic energy transfer (Dimer) calculations.
"""
from copy                                                                       import deepcopy
from SUMELF                                                                     import make_folder
from SUMELF                                                                     import check_molecule_against_file
from ECCP.ECCP.write_dimers_to_disk_methods.write_methods.gaussian_modified_EET import write_gaussian_in as write_gaussian_in_EET
from ECCP.ECCP.write_molecules_to_disk_methods.shared_methods                   import change_folder_name_components
from ECCP.ECCP.write_molecules_to_disk_methods.shared_methods                   import slurmSL_header, load_gaussian_programs, make_gaussian_temp_folder, remove_gaussian_temp_files

def write_EET_gaussian_files(original_molecule_1, original_molecule_2, full_dimer_name, environment_about_dimer, gaussian_jobs_path, fragmentlist, calc_parameters_for_EETs, submission_information_for_EETs):
	"""
	This method will write information the Gaussian files to disk.

	Parameters
	----------
	original_molecule_1 : ase.Atoms
		This is the first molecule in the dimer.
	original_molecule_2 : ase.Atoms
		This is the second molecule in the dimer.

	full_dimer_name : str.
		This is the full name of the dimer. 
	environment_about_dimer : ase.Atoms
		This is the enviroment of neighbouring molecules that surround the dimer. 

	gaussian_jobs_path : str.
		This is the path to save gaussian jobs to

	fragmentlist : list
		This is a list that indicates the assignment of fragments, which associates each atom with the molecule it goes with in a dimer.

	calc_parameters_for_EETs : list
		This list contain all the information required to be included in the EET Gaussian file(s).
	submission_information_for_EETs : list
		This list contain all the information required to be included in the EET submit.sl script.
	"""

	# ******************************************************************************************
	# This is included as a test for debugging to see if the eigendata from EET is the same as for the Eigendata calcs.
	# for index in range(len(calc_parameters_for_EETs)):
	#	calc_parameters_for_EETs[index]['show_eigendata_in_output_file'] = True
	# ******************************************************************************************

	# This for loop will create all the various gaussian settings for the same dimer. 
	# This is important for testing a functionals and basis sets.

	# First, make a copy of the gaussian_parameters and submission_information dictionaries.
	gaussian_parameters    = deepcopy(calc_parameters_for_EETs)
	submission_information = deepcopy(submission_information_for_EETs)
	del gaussian_parameters['calc_software']

	# Second, determine if some critical tags that are needed are in the submission_information dictionary. 
	got_cpu = 'cpus_per_task' in submission_information
	got_mem = 'mem' in submission_information
	got_time = 'time' in submission_information
	if not (got_cpu and got_time):
		print()
		print('Error: You need to specify the following in your submission_information dictionary:')
		if not got_cpu:
			print('\t* cpus_per_task')
		if not got_mem:
			print('\t* mem')
		if not got_time:
			print('\t* time')
		print('See https://github.com/geoffreyweal/ECCP/ for more information about these tags.')
		exit('This program will finish without completing.')

	# Third, copy some tag information that is in the submission_information dictionary to the gaussian_parameters dictionary.
	gaussian_parameters['nprocshared'] = submission_information['cpus_per_task']

	# Fourth, include the fragment list in the gaussian parameters
	gaussian_parameters['fragmentlist'] = fragmentlist

	# Fifth, give the name of the folder to place gaussian files to.
	quantum_chemistry_program = 'GAUSSIAN'
	functional                = change_folder_name_components(gaussian_parameters['method'])
	basis_set                 = change_folder_name_components(gaussian_parameters['basis'])
	funct_and_basis_name      = 'F_'+functional+'_B_'+basis_set
	calc_folder               = str(gaussian_jobs_path)+'/'+str(full_dimer_name)+'/'+quantum_chemistry_program+'_'+str(funct_and_basis_name)

	# Sixth, provide the name and filepath for each of the scratch files.
	scratch_dir_given = ('temp_folder_path' in gaussian_parameters) # was gaussian_scratch_name
	if scratch_dir_given:
		temp_folder_path = gaussian_parameters['temp_folder_path']+'/'+calc_folder
		del gaussian_parameters['temp_folder_path']
		submission_information['temp_folder_path'] = temp_folder_path

	# =============================================================================
	# Seventh, for those temporary files that I have control over where they get saved to, 
	# indicate to the .gjf file where to save those files.
	
	# 7.1, provide the name and filepath for each of the scratch files.
	scratch_dir_given = ('temp_folder_path' in gaussian_parameters) # was gaussian_scratch_name
	if scratch_dir_given:
		temp_folder_path = gaussian_parameters['temp_folder_path']+'/'+calc_folder+'/'+GorE_state_foldername+'/'+calc_type
		del gaussian_parameters['temp_folder_path']
		submission_information['temp_folder_path'] = temp_folder_path

	# 7.2: Make the checkpoint file.
	gaussian_parameters['chk'] = 'gaussian.chk'

	# 7.3: Make path folder and file details for the other gaussian checkpoint files.
	for suffix in ['rwf','int','d2e','skr']:
		if scratch_dir_given: # A scratch path is given.
			gaussian_parameters[suffix] = temp_folder_path+'/'+'gaussian.'+str(suffix)
		else: # A scratch path has not been given.
			# Default name given called gaussian.suffix, whether scratch_dir_given is True or False
			gaussian_parameters[suffix] = 'gaussian.'+str(suffix)

	# =============================================================================

	# Eighth, write the folder to place gaussian files to.
	make_folder(calc_folder)

	# =============================================================================

	# Ninth, make a copy of the first molecule and set the periodic boundary settings (pbc) to False.
	molecule_1 = original_molecule_1.copy()
	molecule_1.set_pbc(False)

	# Tenth, make a copy of the second molecule and set the periodic boundary settings (pbc) to False.
	molecule_2 = original_molecule_2.copy()
	molecule_2.set_pbc(False)

	# =============================================================================

	# Eleventh, create the dimer ase.Atoms object
	dimer = molecule_1 + molecule_2

	# Twelfth, get the path to the dimer to save to
	path_to_dimer = calc_folder+'/'+full_dimer_name+'.gjf'

	# Thirteenth, If this dimer has already been written to disk, check that the atoms in the dimer here is the same as that on file. 
	check_molecule_against_file(dimer, path_to_dimer)

	# Fourteenth, create the gaussian .gjf file.
	with open(path_to_dimer, 'w') as fd:
		write_gaussian_in_EET(fd, molecule_1, molecule_2, environment_about_dimer, run_EET=True, full_dimer_name=full_dimer_name, **gaussian_parameters)

	# Fifteenth, indicate the name of the output file.
	submission_information['log_filename'] = 'output.log'

	# Sixteenth, create the submit.sl file to submit this gaussian job to slurm.
	make_submitSL(full_dimer_name+'.gjf',calc_folder,'EET',functional,basis_set,gaussian_parameters,**submission_information)

# ----------------------------------------------------------------------------------------------------------------------------------

def make_submitSL(gaussian_filename,local_path,suffix,functional,basis_set,gaussian_parameters,cpus_per_task,mem,time,partition='parallel',constraint=None,nodelist=None,exclude=None,email='',python_version='python/3.8.1',gaussian_version='gaussian/g16',log_filename='output.log',temp_folder_path=None,remove_chk_file=True):
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
	gaussian_parameters : dict:
		here
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
	temp_folder_path : str. or None
		This is the path to the scratch directory to save Gaussian temp files to. If you dont give this, Gaussian temp files will be saves to the default scratch directory. Default: None
	remove_chk_file : bool.
		This variable indicates if you want to remove the chk file afterwards. Default: False
	"""

	only_get_FRAG1_STATE1_to_FRAG2_STATE1_EET = True
	
	# create name for job
	name = '-'.join(local_path.split('/')[-3:])+'-'+suffix
	# writing the submit.sl script
	with open(local_path+'/'+"submit.sl", "w") as submitSL:
		slurmSL_header(submitSL, name, mem, partition, constraint, nodelist, time, email, cpus_per_task=cpus_per_task, exclude=exclude)
		make_gaussian_temp_folder(submitSL, temp_folder_path)
		load_gaussian_programs(submitSL, gaussian_version, python_version=None)
		if only_get_FRAG1_STATE1_to_FRAG2_STATE1_EET:
			submitSL.write('# ----------------------------\n')
			submitSL.write('# Execute a program that will monitor when Frag=2,State=1 <=> Frag=1,State=1 has been calculated by Gaussian.\n')
			submitSL.write('# Once Gaussian has obtained this EET value, this program will finish the job as no more combinations of \n')
			submitSL.write('# fragment state EET combinations need to be calculated. \n')
			submitSL.write('\n')
			submitSL.write('monitor_EET_calculation.py '+str(log_filename)+' ${SLURM_JOB_ID} &\n')
			submitSL.write('\n')
		submitSL.write('# ----------------------------\n')
		submitSL.write('# Perform Gaussian Calculation\n')
		submitSL.write('\n')
		submitSL.write('srun '+str(gaussian_version.split('/')[-1])+' < '+str(gaussian_filename)+' > '+str(log_filename)+'\n')
		submitSL.write('\n')
		submitSL.write('# ----------------------------\n')
		# For those temporary files that I have control over where they get saved to, indicate to the .gjf file where to save those files.
		remove_gaussian_temp_files(submitSL, gaussian_parameters, temp_folder_path, remove_chk_file=remove_chk_file, remove_temp_folder=True)
		submitSL.write('# ----------------------------\n')
		submitSL.write('echo "End of job"\n')
		submitSL.write('# ----------------------------\n')

# ----------------------------------------------------------------------------------------------------------------------------------



