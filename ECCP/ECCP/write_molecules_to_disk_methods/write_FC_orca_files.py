"""
write_FC_orca_files.py, Geoffrey Weal, 8/5/22

This script is designed to write the ORCA files and submit.sl files required for performing ORCA jobs for performing reorganisation energy (RE) calculations.
"""
from ase                                                                      import Atoms
from copy                                                                     import deepcopy
from SUMELF                                                                   import make_folder
from SUMELF                                                                   import check_molecule_against_file
from ECCP.ECCP.write_molecules_to_disk_methods.write_methods.orca_modified_FC import write_orca_in_FC
from ECCP.ECCP.write_molecules_to_disk_methods.shared_methods                 import change_folder_name_components, convert_dict_for_bash_input
from ECCP.ECCP.write_molecules_to_disk_methods.shared_methods                 import slurmSL_header, load_orca_programs, make_orca_temp_folder, remove_orca_temp_files

def write_FC_orca_files(molecule, molecule_name, SolventsList, orca_jobs_path, calc_parameters_for_FCs, submission_information_for_FCs):
	"""
	This method will write information the ORCA files to disk.

	Parameters
	----------
	molecule : ase.Atoms.
		This is the molecule. 
	molecule_name : str.
		This is the name of the molecule. 
	SolventsList : list of int
		These are the indices of the molecules in the molecules list that have been identified at solvents.
	orca_jobs_path : str.
		This is the path to save ORCA jobs to.
	calc_parameters_for_FCs : list
		This dictionary contain all the information required for the ORCA input FC file.
	submission_information_for_FCs : list
		This dictionary contain all the information required for the submit.sl script. 
	"""

	raise Exception('Need to make sure this emthod works and is complient')
	raise Exception('Make sure that info from molecule graph is added to molecule before writing the orca file.')
	
	# First, check if there already exists this molecule on file, check if the molecules are the same.
	check_molecule_against_file(molecule, calc_folder+'/'+molecule_name+'.gjf') # Do we need this method

	# First, make a copy of the orca_parameters and submission_information dictionaries.
	orca_parameters    = deepcopy(calc_parameters_for_FCs)
	submission_information = deepcopy(submission_information_for_FCs)
	del orca_parameters['calc_software']

	# Second, determine if some critical tags that are needed are in the submission_information dictionary. 
	got_cpu = 'ntasks' in submission_information
	got_mem = 'mem' in submission_information
	got_time = 'time' in submission_information
	if not (got_cpu and got_time):
		print('Error: You need to specify the following in your submission_information dictionary:')
		if not got_cpu:
			print('\t* ntasks')
		if not got_mem:
			print('\t* mem')
		if not got_time:
			print('\t* time')
		print('See https://github.com/geoffreyweal/ECCP/ for more information about these tags.')
		print('submission_information = '+str(submission_information))
		exit('This program will finish without completing.')

	# Third, copy some tag information that is in the submission_information dictionary to the orca_parameters dictionary.
	orca_parameters['NPROCS'] = submission_information['ntasks']

	# Fourth, give the name of the folder to place ORCA files to.
	quantum_chemistry_program = 'ORCA'
	functional                = change_folder_name_components(gaussian_parameters['method'])
	basis_set                 = change_folder_name_components(gaussian_parameters['basis'])
	funct_and_basis_name      = 'F_'+functional+'_B_'+basis_set
	calc_folder               = str(gaussian_jobs_path)+'/'+str(molecule_name)+'/'+quantum_chemistry_program+'_'+str(funct_and_basis_name)

	# Fifth, provide the name and filepath for each of the scratch files.
	scratch_dir_given = ('temp_folder_path' in orca_parameters) # was orca_scratch_name
	if scratch_dir_given:
		temp_folder_path = orca_parameters['temp_folder_path']+'/'+orca_folder
		del orca_parameters['temp_folder_path']
		submission_information['temp_folder_path'] = temp_folder_path

	# =============================================================================
	# Sixth, for those temporary files that I have control over where they get saved to, 
	# indicate to the .inp file where to save those files.

	# 6.1: check if 'remove_chk_file' is in orca_parameters, otherwise give this.
	if 'remove_chk_file' in orca_parameters:
		del orca_parameters['remove_chk_file']
	if 'remove_chk_file' in submission_information:
		del submission_information['remove_chk_file']
	
	# 6.2: Make path folder and file details for the other ORCA checkpoint files.
	for suffix in ['chk','rwf','int','d2e','skr']:
		if scratch_dir_given: # A scratch path is given.
			orca_parameters[suffix] = temp_folder_path+'/'+'orca.'+str(suffix)
		else: # A scratch path has not been given.
			# Default name given called orca.suffix, whether scratch_dir_given is True or False
			orca_parameters[suffix] = 'orca.'+str(suffix)

	# =============================================================================
	# Seventh, set some settings to the orca_parameters so it is ready for performing the franck-condon calculation.

	# 7.1: Set the oldchk name to 'eGS_gGS_orca.chk'
	orca_parameters['oldchk'] = 'eGS_gGS_orca.chk'

	# 7.2: Add the end lines to the addsec
	addsec = "Final=Source=Chk\nPrint=(Spectra=All,Matrix=JK,HuangRhys)\n\neES_gES_orca.chk"
	orca_parameters['addsec'] = addsec

	# =============================================================================

	# Eighth, the names for the Franck-Condon Folder
	#frank_condon_foldername  = 'frank_condon'
	frank_condon_orca_folder  = orca_folder # +'/'+ frank_condon_foldername

	# Ninth, write the folder to place orca files to.
	make_folder(frank_condon_ORCA_folder)

	# Tenth, write the name of the inp files to make
	frank_condon_orca_filename = 'FC.inp'

	# Eleventh, create the ORCA .inp file for optimising the ground structure.
	with open(frank_condon_orca_folder+'/'+frank_condon_ORCA_filename,'w') as fd:
		write_orca_in_FC (fd, Atoms(), molecule_name=molecule_name+' - FC', **orca_parameters)

	# Twelfth, create the submit .sl file for optimising the ground structure, and to perform the frequency calculation for the optimised ground state structure. 
	make_FC_orca_submitSL(frank_condon_ORCA_filename, frank_condon_ORCA_folder, functional, basis_set, orca_parameters, **submission_information)

# ------------------------------------------------------------------------------------------------------------------------------

def make_FC_orca_submitSL(orca_input_filename,local_path,functional,basis_set,orca_parameters,ntasks,mem,time,partition='parallel',constraint=None,nodelist=None,email='',python_version='python/3.8.1',orca_version='ORCA/5.0.3',gcc_version='GCC/11.2.0',openmpi_version='OpenMPI/4.1.1',temp_folder_path=None):
	"""
	This method will write the submit.sl file in parallel

	Parameters
	----------
	orca_input_filename : str. 
		This is the name of the optimisation file 
	local_path : str. 
		This is the location to save this submit.sl file to
	perform_TD_and_freq : bool.
		This boolean indicates if you are running TD and Freq calculations.
	functional : str. 
		This is the functional you are going to use in your ORCA calculation.
	basis_set : str. 
		This is the basis set you are going to use in your ORCA calculation.
	orca_parameters : dict.
		This dictionary contains all the input parameters required for creating the ORCA input file.
	ntasks : int
		This is the number of cpus you want to use for ORCA jobs.
	mem : str.
		This is the amount of memory you want to use for ORCA jobs.
	time : str.
		This is the amount of time you want to use for ORCA jobs.
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
	orca_version : str.
		This is the version of ORCA you want to load/use in slurm. Default: 'ORCA/5.0.3'
	gcc_version : str.
		This is the version of GCC you want to load/use in slurm. Default: 'GCC/11.2.0'
	openmpi_version : str.
		This is the version of OpenMPI you want to load/use in slurm. Default: 'OpenMPI/4.1.1'
	temp_folder_path : str. or None
		This is the path to the scratch directory to save ORCA temp files to. If you dont give this, ORCA temp files will be saves to the default scratch directory. Default: None
	"""

	# create name for job.
	orca_input_name = orca_input_filename.replace('.inp','')
	name = '-'.join(local_path.split('/')[-4:-1])+'-ReorgE-'+str(orca_input_name)
	
	# writing the submit.sl script
	with open(local_path+'/'+str(orca_input_name)+"_submit.sl", "w") as submitSL:
		slurmSL_header(submitSL, name, mem, partition, constraint, nodelist, time, email, ntasks=ntasks)
		make_orca_temp_folder(submitSL, temp_folder_path)
		load_orca_programs(submitSL, orca_version, gcc_version, openmpi_version, python_version)
		submitSL.write('# ----------------------------\n')
		submitSL.write('# ORCA under MPI requires that it be called via its full absolute path\n')
		submitSL.write('\n')
		submitSL.write('orca_exe=$(which orca)\n')
		submitSL.write('\n')
		submitSL.write('# ----------------------------\n')
		submitSL.write('# Move checkpoint files over\n')
		submitSL.write('\n')
		submitSL.write('copy_checkpoint_files_for_franck_condon_calculation.py\n')
		submitSL.write('\n')
		submitSL.write('# ----------------------------\n')
		submitSL.write('# Perform geometry optimisation calculation\n')
		submitSL.write('\n')
		submitSL.write('${orca_exe} '+str(orca_input_filename)+' > '+str(orca_input_name)+'.out\n')
		submitSL.write('\n')
		submitSL.write('# ----------------------------\n')
		remove_orca_temp_files(submitSL, orca_parameters, temp_folder_path, remove_chk_file=True, remove_temp_folder=True)
		submitSL.write('# ----------------------------\n')
		submitSL.write('echo "End of job"\n')
		submitSL.write('# ----------------------------\n')

# ------------------------------------------------------------------------------------------------------------------------------




