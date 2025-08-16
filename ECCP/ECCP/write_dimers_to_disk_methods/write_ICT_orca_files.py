"""
write_Dimer_gaussian_files.py, Geoffrey Weal, 8/5/22

This script is designed to write the gaussian files and submit.sl files required for performing Gaussian jobs for performing electronic energy transfer (Dimer) calculations.
"""
from copy                                                                       import deepcopy
from SUMELF                                                                     import make_folder
from ECCP.ECCP.write_dimers_to_disk_methods.write_methods.gaussian_modified_ICT import write_gaussian_in as write_gaussian_in_ICT
from ECCP.ECCP.write_molecules_to_disk_methods.shared_methods                   import change_folder_name_components
from ECCP.ECCP.write_molecules_to_disk_methods.shared_methods                   import slurmSL_header, load_orca_programs, make_gaussian_temp_folder, remove_gaussian_temp_files

def write_ICT_orca_files(dimer, molecule_1, molecule_2, full_dimer_name, gaussian_jobs_path, all_gaussian_parameters_for_ICTs, all_submission_information_for_ICTs, get_dimer_icts=True):
	"""
	This method will write information the Gaussian files to disk.

	Parameters
	----------
	molecule_1 : ase.Atoms
		This is the first molecule in the dimer.
	molecule_2 : ase.Atoms
		This is the second molecule in the dimer.
	dimer : ase.Atoms
		This is the dimer. 

	full_dimer_name : str.
		This is the full name of the dimer. 
	gaussian_jobs_path : str.
		This is the path to save gaussian jobs to

	fragmentlist : list
		This is a list that indicates the assignment of fragments, which associates each atom with the molecule it goes with in a dimer.

	all_gaussian_parameters_for_ICTs : list
		This list contain all the information required to be included in the ICTs Gaussian file(s).
	all_submission_information_for_ICTs : list
		This list contain all the information required to be included in the ICTs submit.sl script.

	get_dimer_icts : bool.
		This tag indicates if the user wants to obtain Gaussian files for ICTs (such as overlap orbtials and molecular orbital energies and coefficients). Default: True
	"""

	raise Exception('Sort out how environment is dealt with here. Also add methods for checking molecules and dimers on file')
	raise Exception('Add the molecule graph back to molecule here')

	# This for loop will create all the various gaussian settings for the same dimer. 
	# This is important for testing a functionals and basis sets.
	for original_gaussian_parameters, original_submission_information in zip(all_gaussian_parameters_for_ICTs, all_submission_information_for_ICTs):

		# First, make a copy of the gaussian_parameters and submission_information dictionaries.
		gaussian_parameters    = deepcopy(original_gaussian_parameters)
		submission_information = deepcopy(original_submission_information)
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

		# Fourth, give the name of the folder to place gaussian files to.
		quantum_chemistry_program = 'ORCA'
		functional                = change_folder_name_components(gaussian_parameters['method'])
		basis_set                 = change_folder_name_components(gaussian_parameters['basis'])
		funct_and_basis_name      = 'F_'+functional+'_B_'+basis_set
		calc_folder               = str(gaussian_jobs_path)+'/'+str(full_dimer_name)+'/'+quantum_chemistry_program+'_'+str(funct_and_basis_name)

		# Fifth, provide the name and filepath for each of the scratch files.
		scratch_dir_given = ('temp_folder_path' in gaussian_parameters) # was gaussian_scratch_name
		if scratch_dir_given:
			temp_folder_path = gaussian_parameters['temp_folder_path']+'/'+calc_folder
			del gaussian_parameters['temp_folder_path']
			submission_information['temp_folder_path'] = temp_folder_path

		# =============================================================================
		# Sixth, for those temporary files that I have control over where they get saved to, 
		# indicate to the .gjf file where to save those files.

		# 6.1, provide the name and filepath for each of the scratch files.
		scratch_dir_given = ('temp_folder_path' in gaussian_parameters) # was gaussian_scratch_name
		if scratch_dir_given:
			temp_folder_path = gaussian_parameters['temp_folder_path']+'/'+calc_folder+'/'+GorE_state_foldername+'/'+calc_type
			del gaussian_parameters['temp_folder_path']
			submission_information['temp_folder_path'] = temp_folder_path

		# 6.2: Make the checkpoint file.
		gaussian_parameters['chk'] = 'gaussian.chk'

		# 6.3: Make path folder and file details for the other gaussian checkpoint files.
		for suffix in ['rwf','int','d2e','skr']:
			if scratch_dir_given: # A scratch path is given.
				gaussian_parameters[suffix] = temp_folder_path+'/'+'gaussian.'+str(suffix)
			else: # A scratch path has not been given.
				# Default name given called gaussian.suffix, whether scratch_dir_given is True or False
				gaussian_parameters[suffix] = 'gaussian.'+str(suffix)

		# =============================================================================

		# Seventh, write the folder to place gaussian files to.
		make_folder(calc_folder)

		# ----------------------------------------------------------------
		# Eighth, create the gaussian .gjf file

		# 8.1: Get the names of the molecules
		dimer_name, molecule1_name, molecule2_name = full_dimer_name.split('_')

		# 8.2: Write the folder to save ICTs of the dimer and both molecules to.
		dimer_path = calc_folder+'/Dimer'
		make_folder(dimer_path)

		# 8.3: Save the .gjf file of the dimer.
		with open(dimer_path+'/'+full_dimer_name+'.gjf','w') as fd:
			dimer_for_input = dimer.copy()
			dimer_for_input.set_pbc(False)
			write_gaussian_in_ICTs(fd, dimer_for_input, get_icts=get_dimer_icts, full_dimer_name=full_dimer_name, **gaussian_parameters)

		# 8.4: Write the .gjf file for the first molecule.
		molecule1_path = calc_folder+'/Monomer_1'
		make_folder(molecule1_path)
		with open(molecule1_path+'/'+molecule1_name+'.gjf','w') as fd:
			molecule_1_for_input = molecule_1.copy()
			molecule_1_for_input.set_pbc(False)
			write_gaussian_in_ICTs(fd, molecule_1_for_input, get_icts=get_dimer_icts, full_dimer_name='Monomer 1 - '+molecule1_name, **gaussian_parameters)

		# 8.5: Write the .gjf file for the second molecule.
		molecule2_path = calc_folder+'/Monomer_2'
		make_folder(molecule2_path)
		with open(molecule2_path+'/'+molecule2_name+'.gjf','w') as fd:
			molecule_2_for_input = molecule_2.copy()
			molecule_2_for_input.set_pbc(False)
			write_gaussian_in_ICTs(fd, molecule_2_for_input, get_icts=get_dimer_icts, full_dimer_name='Monomer 2 - '+molecule2_name, **gaussian_parameters)

		# ----------------------------------------------------------------

		# Ninth, indicate the name of the output file
		submission_information['log_filename'] = 'output.log'

		# Tenth, create the submit.sl file to submit this gaussian job to slurm
		make_submitSL(full_dimer_name+'.gjf',dimer_path,    'Eigen',functional,basis_set,gaussian_parameters,**submission_information)
		make_submitSL(molecule1_name +'.gjf',molecule1_path,'MO',   functional,basis_set,gaussian_parameters,**submission_information)
		make_submitSL(molecule2_name +'.gjf',molecule2_path,'MO',   functional,basis_set,gaussian_parameters,**submission_information)

# ----------------------------------------------------------------------------------------------------------------------------------

def make_submitSL(gaussian_filename,local_path,suffix,functional,basis_set,gaussian_parameters,cpus_per_task,mem,time,partition='parallel',constraint=None,nodelist=None,email='',python_version='python/3.8.1',orca_version='ORCA/5.0.3',gcc_version='GCC/11.2.0',openmpi_version='OpenMPI/4.1.1',log_filename='output.log',temp_folder_path=None,remove_chk_file=False):
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
	
	# create name for job
	name = '-'.join(local_path.split('/')[-3:])+'-'+suffix
	# writing the submit.sl script
	with open(local_path+'/'+"submit.sl", "w") as submitSL:
		slurmSL_header(submitSL, name, cpus_per_task, mem, partition, constraint, nodelist, time, email)
		make_gaussian_temp_folder(submitSL, temp_folder_path)
		load_orca_programs(submitSL, orca_version, gcc_version, openmpi_version, python_version=None)
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


