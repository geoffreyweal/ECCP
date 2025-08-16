"""
write_RE_FC_gaussian_files.py, Geoffrey Weal, 8/5/22

This script is designed to write the gaussian files and submit.sl files required for performing Gaussian jobs for performing reorganisation energy (RE) calculations.
"""
from ase                                                                          import Atoms
from copy                                                                         import deepcopy
from SUMELF                                                                       import make_folder
from SUMELF                                                                       import check_molecule_against_file
from ECCP.ECCP.write_molecules_to_disk_methods.write_methods.gaussian_modified_FC import write_gaussian_in as write_gaussian_in_FC
from ECCP.ECCP.write_molecules_to_disk_methods.shared_methods                     import change_folder_name_components, convert_dict_for_bash_input
from ECCP.ECCP.write_molecules_to_disk_methods.shared_methods                     import slurmSL_header, load_gaussian_programs, make_gaussian_temp_folder, remove_gaussian_temp_files

def write_FC_gaussian_files(molecule, molecule_name, SolventsList, gaussian_jobs_path, calc_parameters_for_FCs, submission_information_for_FCs):
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
	calc_parameters_for_FCs : list
		This dictionary contain all the information required for the Gaussian input FC file.
	submission_information_for_FCs : list
		This dictionary contain all the information required for the submit.sl script. 
	"""

	raise Exception('Need to make sure this emthod works and is complient')

	# First, check if there already exists this molecule on file, check if the molecules are the same.
	check_molecule_against_file(molecule, calc_folder+'/'+molecule_name+'.gjf') # Do we need this method

	# First, make a copy of the gaussian_parameters and submission_information dictionaries.
	gaussian_parameters    = deepcopy(calc_parameters_for_FCs)
	submission_information = deepcopy(submission_information_for_FCs)
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

	# Fifth, provide the name and filepath for each of the scratch files.
	scratch_dir_given = ('temp_folder_path' in gaussian_parameters) # was gaussian_scratch_name
	if scratch_dir_given:
		temp_folder_path = gaussian_parameters['temp_folder_path']+'/'+calc_folder
		del gaussian_parameters['temp_folder_path']
		submission_information['temp_folder_path'] = temp_folder_path

	# =============================================================================
	# Sixth, for those temporary files that I have control over where they get saved to, 
	# indicate to the .gjf file where to save those files.

	# 6.1: check if 'remove_chk_file' is in gaussian_parameters, otherwise give this.
	if 'remove_chk_file' in gaussian_parameters:
		del gaussian_parameters['remove_chk_file']
	if 'remove_chk_file' in submission_information:
		del submission_information['remove_chk_file']
	
	# 6.2: Make path folder and file details for the other gaussian checkpoint files.
	for suffix in ['chk','rwf','int','d2e','skr']:
		if scratch_dir_given: # A scratch path is given.
			gaussian_parameters[suffix] = temp_folder_path+'/'+'gaussian.'+str(suffix)
		else: # A scratch path has not been given.
			# Default name given called gaussian.suffix, whether scratch_dir_given is True or False
			gaussian_parameters[suffix] = 'gaussian.'+str(suffix)

	# =============================================================================
	# Seventh, set some settings to the gaussian_parameters so it is ready for performing the franck-condon calculation.

	# 7.1: Set the oldchk name to 'eGS_gGS_gaussian.chk'
	gaussian_parameters['oldchk'] = 'eGS_gGS_gaussian.chk'

	# 7.2: Add the end lines to the addsec
	addsec = "Final=Source=Chk\nPrint=(Spectra=All,Matrix=JK,HuangRhys)\n\neES_gES_gaussian.chk"
	gaussian_parameters['addsec'] = addsec

	# =============================================================================

	# Eighth, the names for the Franck-Condon Folder
	#frank_condon_foldername  = 'frank_condon'
	frank_condon_calc_folder  = calc_folder # +'/'+ frank_condon_foldername

	# Ninth, write the folder to place gaussian files to.
	make_folder(frank_condon_calc_folder)

	# Tenth, write the name of the gjf files to make
	frank_condon_Gaussian_filename = 'FC.gjf'

	# Eleventh, create the gaussian .gjf file for optimising the ground structure.
	with open(frank_condon_calc_folder+'/'+frank_condon_Gaussian_filename,'w') as fd:
		write_gaussian_in_FC (fd, Atoms(), molecule_name=molecule_name+' - FC', **gaussian_parameters)

	# Twelfth, create the submit .sl file for optimising the ground structure, and to perform the frequency calculation for the optimised ground state structure. 
	make_FC_gaussian_submitSL(frank_condon_Gaussian_filename, frank_condon_calc_folder, functional, basis_set, gaussian_parameters, **submission_information)

# ------------------------------------------------------------------------------------------------------------------------------

def make_FC_gaussian_submitSL(gaussian_input_filename,local_path,functional,basis_set,gaussian_parameters,cpus_per_task,mem,time,partition='parallel',constraint=None,nodelist=None,exclude=None,email='',python_version='python/3.8.1',gaussian_version='gaussian/g16',temp_folder_path=None):
	"""
	This method will write the submit.sl file in parallel

	Parameters
	----------
	gaussian_input_filename : str. 
		This is the name of the optimisation file 
	local_path : str. 
		This is the location to save this submit.sl file to
	perform_TD_and_freq : bool.
		This boolean indicates if you are running TD and Freq calculations.
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
	temp_folder_path : str. or None
		This is the path to the scratch directory to save Gaussian temp files to. If you dont give this, Gaussian temp files will be saves to the default scratch directory. Default: None
	"""

	# create name for job.
	gaussian_input_name = gaussian_input_filename.replace('.gjf','')
	name = '-'.join(local_path.split('/')[-4:-1])+'-ReorgE-'+str(gaussian_input_name)
	
	# writing the submit.sl script
	with open(local_path+'/'+str(gaussian_input_name)+"_submit.sl", "w") as submitSL:
		slurmSL_header(submitSL, name, mem, partition, constraint, nodelist, time, email, cpus_per_task=cpus_per_task, exclude=exclude)
		make_gaussian_temp_folder(submitSL, temp_folder_path)
		load_gaussian_programs(submitSL, gaussian_version, python_version)
		submitSL.write('# ----------------------------\n')
		submitSL.write('# Move checkpoint files over\n')
		submitSL.write('\n')
		submitSL.write('copy_checkpoint_files_for_franck_condon_calculation.py\n')
		submitSL.write('\n')
		submitSL.write('# ----------------------------\n')
		submitSL.write('# Perform geometry optimisation calculation\n')
		submitSL.write('\n')
		submitSL.write('srun '+str(gaussian_version.split('/')[-1])+' < '+str(gaussian_input_filename)+' > '+str(gaussian_input_name)+'.log\n')
		submitSL.write('\n')
		submitSL.write('# ----------------------------\n')
		remove_gaussian_temp_files(submitSL, gaussian_parameters, temp_folder_path, remove_chk_file=True, remove_temp_folder=True)
		submitSL.write('# ----------------------------\n')
		submitSL.write('echo "End of job"\n')
		submitSL.write('# ----------------------------\n')

# ------------------------------------------------------------------------------------------------------------------------------




