"""
write_RE_orca_SP_files.py, Geoffrey Weal, 8/5/22

This script is designed to create the single point ORCA calculation files (.inp files) for performing reorganisation energy (RE) calculations.
"""
from copy                                                     import deepcopy
from SUMELF                                                   import make_folder
from ECCP.ECCP.write_molecules_to_disk_methods.shared_methods import change_folder_name_components, convert_dict_for_bash_input
from ECCP.ECCP.write_molecules_to_disk_methods.shared_methods import slurmSL_header, load_orca_programs, make_orca_temp_folder, remove_orca_temp_files

def write_RE_orca_SP_files(molecule, molecule_name, SolventsList, orca_jobs_path, calc_parameters_for_RE_SPs, submission_information_for_RE_SPs):
	"""
	This method will write information the ORCA files to disk for performing reorganisation energy calculations.

	Parameters
	----------
	molecule : ase.Atoms.
		This is the molecule. 
	molecule_name : str.
		This is the name of the molecule. 
	SolventsList : list of int
		These are the indices of the molecules in the molecules list that have been identified at solvents.
	orca_jobs_path : str.
		This is the path to save orca jobs to.
	calc_parameters_for_RE_SPs : list
		This dictionary contain all the information required for the ORCA input single point RE file.
	submission_information_for_RE_SPs : list
		This dictionary contain all the information required for the submit.sl script. 
	"""

	#raise Exception('Check this works')
	#raise Exception('Make sure that info from molecule graph is added to molecule before writing the orca file.')

	# First, make a copy of the orca_parameters and submission_information dictionaries.
	orca_parameters        = deepcopy(calc_parameters_for_RE_SPs)
	submission_information = deepcopy(submission_information_for_RE_SPs)
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
	functional                = change_folder_name_components(orca_parameters['method'])
	basis_set                 = change_folder_name_components(orca_parameters['basis'])
	funct_and_basis_name      = 'F_'+functional+'_B_'+basis_set
	calc_folder               = str(orca_jobs_path)+'/'+str(molecule_name)+'/'+quantum_chemistry_program+'_'+str(funct_and_basis_name)

	# =============================================================================
	# Fifth, for those temporary files that I have control over where they get saved to, 
	# indicate to the .inp file where to save those files.

	# 5.1: Setup the folder names.
	ground_state_foldername  = 'ground_structure'
	excited_state_foldername = 'excited_structure'
	SP_calc_type             = 'SP'

	# 5.2: Make a copy of the orca_parameters and submission_information for gound and excited states
	orca_parameters_GS = deepcopy(orca_parameters)
	orca_parameters_ES = deepcopy(orca_parameters)

	submission_information_GS = deepcopy(submission_information)
	submission_information_ES = deepcopy(submission_information)

	# =============================================================================

	# Sixth, get the names for the Ground and Excited States
	ground_structure_foldername   = 'ground_structure'
	excited_structure_foldername  = 'excited_structure'
	ground_structure_orca_folder  = calc_folder +'/'+ ground_structure_foldername
	excited_structure_orca_folder = calc_folder +'/'+ excited_structure_foldername

	# Seventh, write the folder to place ORCA files to.
	make_folder(ground_structure_orca_folder)
	make_folder(excited_structure_orca_folder)

	# Eighth, write the name of the inp files to make
	ground_structure_GS  = 'eGS_gGS_main_opt.inp'
	ground_structure_ES  = 'eES_gGS.inp'
	excited_structure_ES = 'eES_gES_main_opt.inp'
	excited_structure_GS = 'eGS_gES.inp'

	# Ninth, create the submit .sl file for optimised ground structure.
	make_RE_orca_submitSL(ground_structure_GS,  ground_structure_ES,  ground_structure_orca_folder,   True, functional, basis_set, orca_parameters_GS, **submission_information_GS)

	# Tenth, create the submit .sl file for optimised excited structure.
	make_RE_orca_submitSL(excited_structure_ES, excited_structure_GS, excited_structure_orca_folder, False, functional, basis_set, orca_parameters_ES, **submission_information_ES)

# ------------------------------------------------------------------------------------------------------------------------------

def make_RE_orca_submitSL(optimisation_filename,single_point_filename,local_path,perform_TD,functional,basis_set,orca_parameters,ntasks,mem,time,partition='parallel',constraint=None,nodelist=None,email='',python_version='python/3.8.1',orca_version='ORCA/5.0.3',gcc_version='GCC/11.2.0',openmpi_version='OpenMPI/4.1.1',temp_folder_path=None):
	"""
	This method will write the individual submit.sl files so all ORCA jobs can be run individually in slurm (in 'parallel'). 

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
	temp_folder_path : str.
		This is the path to the scratch directory to save ORCA temp files to. If you dont give this, ORCA temp files will be saves to the default scratch directory. Default: None
	"""

	# Get version of ORCA to use
	orca_version_suffix = str(orca_version.split('/')[-1])
	
	# Make changes to orca_parameters to change orca filenames to include "freq" in their name.
	orca_parameters_for_SP = deepcopy(orca_parameters)
	
	# If not performing an excited state calculation, remove the 'td_settings' entry so there is not confusion in the file. 
	if (not perform_TD) and ('td_settings' in orca_parameters_for_SP):
		del orca_parameters_for_SP['td_settings']

	# create name for job
	optimisation_name = optimisation_filename.replace('.inp','')
	single_point_name = single_point_filename.replace('.inp','')
	name = '-'.join(local_path.split('/')[-4:-1])+'-ReorgE_SP-'+str(single_point_name)

	# writing the submit.sl script
	with open(local_path+'/'+str(single_point_name)+"_submit.sl", "w") as submitSL:
		slurmSL_header(submitSL, name, mem, partition, constraint, nodelist, time, email, ntasks=ntasks)
		make_orca_temp_folder(submitSL, temp_folder_path)
		load_orca_programs(submitSL, orca_version, gcc_version, openmpi_version, python_version)
		submitSL.write('# ============================\n')
		submitSL.write('# Prevent the single point job from running if it is already running or has already run.\n')
		submitSL.write('\n')
		submitSL.write('if ! [[ -f '+str(single_point_name)+'.out'+' ]]\n')
		submitSL.write('then\n')
		submitSL.write('\t# ============================\n')
		submitSL.write('\t# ORCA under MPI requires that it be called via its full absolute path\n')
		submitSL.write('\t\n')
		submitSL.write('\torca_exe=$(which orca)\n')
		submitSL.write('\t\n')
		submitSL.write('\t# ============================\n')
		submitSL.write('\t# Settings for creating new inp files\n')
		submitSL.write('\t\n')
		submitSL.write('\torca_parameters='+str(convert_dict_for_bash_input(orca_parameters_for_SP))+'\n')
		submitSL.write('\t\n')
		submitSL.write('\t# ============================\n')
		submitSL.write('\t# Extract optimised structure and place it into single point calculation ORCA input file.\n')
		submitSL.write('\t\n')
		submitSL.write('\tget_single_point_RE_ORCA_input_file.py '+str(optimisation_name)+'.out '+str(single_point_name)+'.inp '+str(perform_TD)+' "${orca_parameters}"\n')
		submitSL.write('\t\n')
		submitSL.write('\t# ============================\n')
		submitSL.write('\t# Performing ORCA Calculation\n')
		submitSL.write('\t\n')
		submitSL.write('\t${orca_exe} '+str(single_point_filename)+' > '+str(single_point_name)+'.out\n')
		submitSL.write('\t\n')
		#remove_orca_temp_files(submitSL, orca_parameters_for_SP, temp_folder_path, remove_chk_file=True, remove_temp_folder=True)
		submitSL.write('\t# ============================\n')
		submitSL.write('\techo "End of job"\n')
		submitSL.write('\t# ============================\n')
		submitSL.write('fi\n')

# ------------------------------------------------------------------------------------------------------------------------------





