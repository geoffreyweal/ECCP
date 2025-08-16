"""
write_RE_orca_files.py, Geoffrey Weal, 8/5/22

This script is designed to write the ORCA files and submit.sl files required for performing ORCA jobs for performing reorganisation energy (RE) calculations.
"""
from copy                                                                     import deepcopy
from ase.io                                                                   import write
from SUMELF                                                                   import make_folder
from SUMELF                                                                   import check_molecule_against_file
from SUMELF                                                                   import obtain_graph
from ECCP.ECCP.write_molecules_to_disk_methods.write_methods.orca_modified_RE import write_orca_in_RE
from ECCP.ECCP.write_molecules_to_disk_methods.shared_methods                 import change_folder_name_components, convert_dict_for_bash_input
from ECCP.ECCP.write_molecules_to_disk_methods.shared_methods                 import slurmSL_header, load_orca_programs, make_orca_temp_folder, remove_orca_temp_files

def write_RE_orca_files(molecule, molecule_name, SolventsList, orca_jobs_path, calc_parameters_for_REs, submission_information_for_REs):
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
	calc_parameters_for_REs : list
		This dictionary contain all the information required for the ORCA input RE file.
	submission_information_for_REs : list
		This dictionary contain all the information required for the submit.sl script. 
	"""

	#raise Exception('Check this works')
	#raise Exception('Make sure that info from molecule graph is added to molecule before writing the orca file.')

	# First, make a copy of the orca_parameters and submission_information dictionaries.
	orca_parameters        = deepcopy(calc_parameters_for_REs)
	submission_information = deepcopy(submission_information_for_REs)
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

	# 5.1: Make a copy of the orca_parameters and submission_information for ground and excited states
	orca_parameters_GS_opt_main        = deepcopy(orca_parameters)
	orca_parameters_ES_opt_main        = deepcopy(orca_parameters)
	submission_information_GS_opt_main = deepcopy(submission_information)
	submission_information_ES_opt_main = deepcopy(submission_information)

	# 5.2: Make a copy of the orca_parameters and submission_information for ground and excited states for freq calcs
	orca_parameters_GS_freq        = deepcopy(orca_parameters)
	orca_parameters_ES_freq        = deepcopy(orca_parameters)
	submission_information_GS_freq = deepcopy(submission_information)
	submission_information_ES_freq = deepcopy(submission_information)

	# =============================================================================

	# Sixth, the names for the Ground and Excited States.
	ground_structure_foldername   = 'ground_structure'
	excited_structure_foldername  = 'excited_structure'
	ground_structure_orca_folder  = calc_folder +'/'+ ground_structure_foldername
	excited_structure_orca_folder = calc_folder +'/'+ excited_structure_foldername

	# Seventh, write the folder to place ORCA files to.
	make_folder(ground_structure_orca_folder)
	make_folder(excited_structure_orca_folder)

	# Eighth, write the name of the inp files to make
	ground_structure_GS_name              = 'eGS_gGS'
	ground_structure_GS_DFT_main_opt      = 'eGS_gGS_main_opt.inp'
	ground_structure_GS_DFT_main_opt_xyz  = 'eGS_gGS_main_opt_input.xyz'
	ground_structure_GS_SP                = 'eES_gGS.inp'
	#ground_structure_GS_SP_xyz            = 'eES_gGS_input.xyz'

	excited_structure_ES_name             = 'eES_gES'
	excited_structure_ES_DFT_main_opt     = 'eES_gES_main_opt.inp'
	excited_structure_ES_DFT_main_opt_xyz = 'eES_gES_main_opt_input.xyz'
	excited_structure_ES_SP               = 'eGS_gES.inp'
	#excited_structure_ES_SP_xyz           = 'eGS_gES.xyz'

	# =============================================================================
	# Tenth, create the ORCA .inp file for optimising the ground and excited structures.

	# 10.1: Get the path to the ground state ORCA file
	path_to_ground_state_ORCA_file = ground_structure_orca_folder+'/'+ground_structure_GS_DFT_main_opt
	path_to_ground_state_XYZ_file  = ground_structure_orca_folder+'/'+ground_structure_GS_DFT_main_opt_xyz

	# 10.2: Check if there already exists this molecule on file, check if the molecules are the same.
	check_molecule_against_file(molecule, path_to_ground_state_XYZ_file)

	# 10.3: Create the ORCA input file for optimising the ground structure.
	with open(path_to_ground_state_ORCA_file, 'w') as fd:
		molecule_for_input = molecule.copy()
		molecule_for_input.set_pbc(False)
		write_orca_in_RE(fd, molecule_for_input, ground_structure_GS_DFT_main_opt_xyz, perform_opt=True, perform_CalcAll=False, perform_TD=False, perform_freq=False, perform_pop=False, **orca_parameters_GS_opt_main)

	# 10.4: Create the xyz file for optimising the ground structure. 
	#magnetic_moments = molecule.get_initial_magnetic_moments()
	molecule_copy, molecule_copy_graph = obtain_graph(molecule.copy(), name='molecule_copy', no_of_cpus=1)
	molecule_copy.set_initial_charges(None); molecule_copy.set_cell(None)
	#molecule_copy.set_initial_magnetic_moments(magnetic_moments)
	write(path_to_ground_state_XYZ_file, molecule_copy)	

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	# 10.5: Get the path to the excited state ORCA file
	path_to_excited_state_ORCA_file = excited_structure_orca_folder+'/'+excited_structure_ES_DFT_main_opt
	path_to_excited_state_XYZ_file  = excited_structure_orca_folder+'/'+excited_structure_ES_DFT_main_opt_xyz

	# 10.6: Check if there already exists this molecule on file, check if the molecules are the same.
	check_molecule_against_file(molecule, path_to_excited_state_XYZ_file)

	# 10.7: Create the ORCA input file for optimising the excited structure.
	with open(path_to_excited_state_ORCA_file, 'w') as fd:
		molecule_for_input = molecule.copy()
		molecule_for_input.set_pbc(False)
		write_orca_in_RE(fd, molecule_for_input, excited_structure_ES_DFT_main_opt_xyz, perform_opt=True, perform_CalcAll=False, perform_TD=True,  perform_freq=False, perform_pop=False, **orca_parameters_ES_opt_main)

	# 10.8: Create the xyz file for optimising the excited structure. 
	#magnetic_moments = molecule.get_initial_magnetic_moments()
	molecule_copy, molecule_copy_graph = obtain_graph(molecule.copy(), name='molecule_copy', no_of_cpus=1)
	molecule_copy.set_initial_charges(None); molecule_copy.set_cell(None)
	#molecule_copy.set_initial_magnetic_moments(magnetic_moments)
	write(path_to_excited_state_XYZ_file, molecule_copy)

	# =============================================================================
	# Eleventh, create the submit .sl file for optimising the ground and excited structures.

	# 11.1: Create the submit .sl file for optimising the ground structure.
	# Make the initial optimimsation files.
	main_opt_chk_filepath = make_RE_orca_submitSL(ground_structure_GS_name,  ground_structure_GS_DFT_main_opt,  ground_structure_GS_SP,  ground_structure_orca_folder, False,        functional, basis_set,                          **submission_information_GS_opt_main)
	# Perform the frequency calculation for the optimised ground state structure (set perform_TD = False). 
	make_RE_freq_orca_submitSL                   (ground_structure_GS_name,                                                              ground_structure_orca_folder, False, True,  functional, basis_set, orca_parameters_GS_freq, **submission_information_GS_freq)

	# 11.2: Create the submit .sl file for optimising the excited structure.
	# Make the initial optimimsation files.
	main_opt_chk_filepath = make_RE_orca_submitSL(excited_structure_ES_name, excited_structure_ES_DFT_main_opt, excited_structure_ES_SP, excited_structure_orca_folder, True,        functional, basis_set,                          **submission_information_ES_opt_main)
	# Perform the frequency calculation for the optimised excited state structure (set perform_TD = True). 
	make_RE_freq_orca_submitSL                   (excited_structure_ES_name,                                                            excited_structure_orca_folder, True, False, functional, basis_set, orca_parameters_ES_freq, **submission_information_ES_freq)

	# =============================================================================

# ------------------------------------------------------------------------------------------------------------------------------

def make_RE_orca_submitSL(main_calculation_type_name, optimisation_filename_DFT_main_opt, single_point_filename, local_path, perform_excited_state_calc, functional, basis_set, ntasks, mem, time, partition='parallel', constraint=None, nodelist=None, email='', python_version='python/3.8.1', orca_version='ORCA/5.0.3',gcc_version='GCC/11.2.0',openmpi_version='OpenMPI/4.1.1', temp_folder_path=None):
	"""
	This method will write the submit.sl file in parallel

	Parameters
	----------
	main_calculation_type_name : str.
		This is the name of the calculation type (either GS_GS or ES_ES).
	optimisation_filename_PM6 : str. 
		This is the name of the optimisation file for optimising with PM6 force field function
	optimisation_filename_DFT_main_opt : str. 
		This is the name of the optimisation file for performing the main optimisation process with DFT
	single_point_filename : str. 
		This is the name of the single point calculation file.
	local_path : str. 
		This is the location to save this submit.sl file to
	perform_excited_state_calc : bool.
		This boolean indicates if you are performing a excited state optimisation.
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

	# Get version of ORCA to use
	orca_version_suffix = str(orca_version.split('/')[-1])

	# create names for job.
	optimisation_name_DFT_main_opt   = optimisation_filename_DFT_main_opt.replace('.inp','')
	single_point_name                = single_point_filename.replace('.inp','')
	name = '-'.join(local_path.split('/')[-4:-1])+'-ReorgE_main-'+str(main_calculation_type_name)

	# writing the submit.sl script
	with open(local_path+'/'+str(main_calculation_type_name)+"_main_opt_submit.sl", "w") as submitSL:
		slurmSL_header(submitSL, name, mem, partition, constraint, nodelist, time, email, ntasks=ntasks)
		make_orca_temp_folder(submitSL, temp_folder_path)
		load_orca_programs(submitSL, orca_version, gcc_version, openmpi_version, python_version)
		
		# Perform the main DFT optimisation

		submitSL.write('# ============================\n')
		submitSL.write('# Prevent the optimisation job from running if it is already running or has already run.\n')
		submitSL.write('\n')
		submitSL.write('if ! [[ -f '+str(optimisation_name_DFT_main_opt)+'.out'+' ]]\n')
		submitSL.write('then\n')
		submitSL.write('\t# ============================\n')
		#submitSL.write('\t# Make temporary folder\n')
		#submitSL.write('\t\n')
		#submitSL.write('\ttempdir=$(pwd)/orca_temp_folder\n')
		#submitSL.write('\tmkdir -p $tempdir\n')
		#submitSL.write('\tmy $scratchdir = "$tempdir"')
		#submitSL.write('\t\n')
		#submitSL.write('\t# ============================\n')
		submitSL.write('\t# ORCA under MPI requires that it be called via its full absolute path\n')
		submitSL.write('\t\n')
		submitSL.write('\torca_exe=$(which orca)\n')
		submitSL.write('\t\n')
		submitSL.write('\t# ============================\n')
		submitSL.write('\t# Perform the main geometry optimisation calculation with the desired functional and basis set.\n')
		submitSL.write('\t\n')
		submitSL.write('\techo "Performing main optimisation calculation"\n')
		submitSL.write('\t${orca_exe} '+str(optimisation_filename_DFT_main_opt)+' > '+str(optimisation_name_DFT_main_opt)+'.out\n')
		submitSL.write('\techo "Finished main optimisation calculation"\n')
		submitSL.write('\t\n')
		#remove_orca_temp_files(submitSL, orca_parameters_main_opt, temp_folder_path, remove_temp_folder=True)
		submitSL.write('\t# ----------------------------\n')
		submitSL.write('\t# Submit the Frequency and single point calculations to slurm.\n')
		submitSL.write('\t\n')
		submitSL.write('\tsubmit_slurm_job.py '+str(main_calculation_type_name)+'_freq_submit.sl\n')
		submitSL.write('\tsubmit_slurm_job.py '+str(single_point_name)+'_submit.sl\n')
		submitSL.write('\t\n')
		submitSL.write('\t# ----------------------------\n')
		submitSL.write('\techo "End of job"\n')
		submitSL.write('\t# ----------------------------\n')
		submitSL.write('fi\n')

# ------------------------------------------------------------------------------------------------------------------------------

def make_RE_freq_orca_submitSL(freq_calc_filename, local_path, perform_TD, perform_raman, functional, basis_set, orca_parameters, ntasks, mem, time, partition='parallel', constraint=None, nodelist=None, email='', python_version='python/3.8.1', orca_version='ORCA/5.0.3', gcc_version='GCC/11.2.0', openmpi_version='OpenMPI/4.1.1', temp_folder_path=None, remove_chk_file=False):
	"""
	This method will write the submit.sl file in parallel for performing frequency calculations. 

	Parameters
	----------
	freq_calc_filename : str. 
		This is the name of the frequency calculation file.
	local_path : str. 
		This is the location to save this submit.sl file to
	perform_TD : bool.
		This tag indicates if the frequency calculation will be performed with TD or not.
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
	remove_chk_file : bool.
		This variable indicates if you want to remove the chk file afterwards. Default: False
	"""

	# Get version of ORCA to use
	orca_version_suffix = str(orca_version.split('/')[-1])

	# If not performing an excited state calculation, remove the 'td_settings' entry so there is not confusion in the file. 
	if (not perform_TD) and ('td_settings' in orca_parameters):
		del orca_parameters['td_settings']

	# create name for job
	freq_calc_name = freq_calc_filename.replace('.inp','')
	name = '-'.join(local_path.split('/')[-4:-1])+'-ReorgE_freq-'+str(freq_calc_name)

	# writing the submit.sl script
	with open(local_path+'/'+str(freq_calc_name)+"_freq_submit.sl", "w") as submitSL:
		slurmSL_header(submitSL, name, mem, partition, constraint, nodelist, time, email, ntasks=ntasks)
		make_orca_temp_folder(submitSL, temp_folder_path)
		load_orca_programs(submitSL, orca_version, gcc_version, openmpi_version, python_version)
		submitSL.write('# ============================\n')
		submitSL.write('# Prevent the optimisation job from running if it is already running or has already run.\n')
		submitSL.write('\n')
		submitSL.write('if ! [[ -f '+str(freq_calc_name)+'_freq.out'+' ]]\n')
		submitSL.write('then\n')
		submitSL.write('\t# ============================\n')
		submitSL.write('\t# ORCA under MPI requires that it be called via its full absolute path\n')
		submitSL.write('\t\n')
		submitSL.write('\torca_exe=$(which orca)\n')
		submitSL.write('\t\n')
		submitSL.write('\t# ============================\n')
		submitSL.write('\t# Settings for creating new inp files\n')
		submitSL.write('\t\n')
		submitSL.write('\torca_parameters='+str(convert_dict_for_bash_input(orca_parameters))+'\n')
		submitSL.write('\t\n')
		submitSL.write('\t# ============================\n')
		submitSL.write('\t# Extract optimised structure and place it into single point calculation ORCA input file.\n')
		submitSL.write('\t\n')
		submitSL.write('\tget_freq_RE_ORCA_input_file.py '+str(freq_calc_name)+'_main_opt.out '+str(freq_calc_name)+'_freq.inp '+str(perform_TD)+' '+str(perform_raman)+' "${orca_parameters}"\n')
		submitSL.write('\t\n')
		submitSL.write('\t# ============================\n')
		submitSL.write('\t# Perform geometry optimisation calculation\n')
		submitSL.write('\t\n')
		submitSL.write('\t${orca_exe} '+str(freq_calc_name)+'_freq.inp > '+str(freq_calc_name)+'_freq.out\n')
		submitSL.write('\t\n')
		if not remove_chk_file:
			# This method will move the frequency calculation checkpoint file to its own mass folder point.
			submitSL.write('\t# ============================\n')
			submitSL.write('\t# Move the frequency calculation checkpoint file (orca_freq.chk) file to its own mass folder\n')
			submitSL.write('\t\n')
			submitSL.write('\tmove_orca_freq_chk_file_to_storage_folder.py\n')
			submitSL.write('\t\n')
		#remove_orca_temp_files(submitSL, orca_parameters, temp_folder_path, remove_chk_file=remove_chk_file, remove_temp_folder=True)
		submitSL.write('\t# ============================\n')
		submitSL.write('\techo "End of job"\n')
		submitSL.write('\t# ============================\n')
		submitSL.write('fi\n')

# ------------------------------------------------------------------------------------------------------------------------------



