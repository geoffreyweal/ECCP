"""
make_RE_gaussian_submitSL.py, Geoffrey Weal, 8/5/22

This method will write the submit.sl file in parallel
"""
from copy                                                     import deepcopy
from ECCP.ECCP.write_molecules_to_disk_methods.shared_methods import slurmSL_header, load_gaussian_programs, make_gaussian_temp_folder, remove_gaussian_temp_files

def make_RE_gaussian_submitSL(main_calculation_type_name, optimisation_filename_DFT_main_opt, single_point_filename, local_path, perform_excited_state_calc, run_excited_state_from_optimised_ground_structure, functional, basis_set, gaussian_parameters, cpus_per_task, mem, time, partition='parallel', constraint=None, nodelist=None, exclude=None, email='', python_version='python/3.8.1', gaussian_version='gaussian/g16', temp_folder_path=None):
	"""
	This method will write the submit.sl file in parallel

	Parameters
	----------
	main_calculation_type_name : str.
		This is the name of the calculation type (either eGS_gGS or eES_gES).
	optimisation_filename_DFT_main_opt : str. 
		This is the name of the optimisation file for performing the main optimisation process with DFT
	single_point_filename : str. 
		This is the name of the single point calculation file.
	local_path : str. 
		This is the location to save this submit.sl file to
	perform_excited_state_calc : bool.
		This boolean indicates if you are performing a excited state optimisation.
	run_excited_state_from_optimised_ground_structure : bool.
		This tag indicates you would like to run the excited state calculation immediately (True, .sl file), or only run them after optimising the ground state calculation and use the optimised ground state structure as the initial input (not active, .sl.inactive file)
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

	# Get version of Gaussian to use
	gaussian_version_suffix = str(gaussian_version.split('/')[-1])

	# Make changes to gaussian_parameters to change gaussian filenames to include "main_opt" in their name.
	gaussian_parameters_main_opt = deepcopy(gaussian_parameters)
	for gaussian_file in ['chk', 'd2e', 'int', 'rwf', 'skr']:
		if gaussian_file in gaussian_parameters_main_opt:
			gaussian_filename = gaussian_parameters_main_opt[gaussian_file]
			gaussian_filename = gaussian_filename.split('.')
			gaussian_filename = gaussian_filename[0]+'_main_opt.'+gaussian_filename[1]
			gaussian_parameters_main_opt[gaussian_file] = gaussian_filename

	# If performing ground state optimisation (not excited state), remove component to indicate td setting from gaussian_parameters. 
	if (not perform_excited_state_calc) and ('td_settings' in gaussian_parameters_main_opt):
		del gaussian_parameters_main_opt['td_settings']

	# create names for job.
	optimisation_name_DFT_main_opt   = optimisation_filename_DFT_main_opt.replace('.gjf','')
	single_point_name                = single_point_filename.replace('.gjf','')
	name = '-'.join(local_path.split('/')[-4:-1])+'-ReorgE_main-'+str(main_calculation_type_name)

	#import pdb; pdb.set_trace()
	# Either active or do not active the submit file by changing the filename
	if run_excited_state_from_optimised_ground_structure and (main_calculation_type_name == 'eES_gES'):
		path_to_submit_file = local_path+'/eES_gES_main_opt_submit.sl.inactive'
	else:
		path_to_submit_file = local_path+'/'+str(main_calculation_type_name)+"_main_opt_submit.sl"

	# writing the submit.sl script
	with open(path_to_submit_file, "w") as submitSL:
		slurmSL_header(submitSL, name, mem, partition, constraint, nodelist, time, email, cpus_per_task=cpus_per_task, exclude=exclude)
		make_gaussian_temp_folder(submitSL, temp_folder_path)
		load_gaussian_programs(submitSL, gaussian_version, python_version)
		
		submitSL.write('# ============================\n')
		submitSL.write('# Prevent the optimisation job from running if it is already running or has already run.\n')
		submitSL.write('\n')
		submitSL.write('if ! [[ -f '+str(optimisation_name_DFT_main_opt)+'.log'+' ]]\n')
		submitSL.write('then\n')

		# Perform the proper DFT optimisation

		submitSL.write('\t# ============================\n')
		submitSL.write('\t# Perform the main geometry optimisation calculation with the desired functional and basis set.\n')
		submitSL.write('\t\n')
		submitSL.write('\techo "Performing main optimisation calculation"\n')
		submitSL.write('\tsrun '+gaussian_version_suffix+' < '+str(optimisation_filename_DFT_main_opt)+' > '+str(optimisation_name_DFT_main_opt)+'.log\n')
		submitSL.write('\techo "Finished main optimisation calculation"\n')
		submitSL.write('\t\n')
		submitSL.write('\t# ============================\n')

		# Perform final tasks and submit for frequency and single point analysis

		remove_gaussian_temp_files(submitSL, gaussian_parameters_main_opt, temp_folder_path, remove_chk_file=True, remove_temp_folder=True, prefix='\t')
		submitSL.write('\t# ----------------------------\n')
		submitSL.write('\t# Submit the Frequency and single point calculations to slurm.\n')
		submitSL.write('\t\n')
		submitSL.write('\tsubmit_slurm_job.py '+str(main_calculation_type_name)+'_freq_submit.sl\n')
		submitSL.write('\tsubmit_slurm_job.py '+str(single_point_name)+'_submit.sl\n')
		submitSL.write('\t\n')

		# Create and submit the excited state calculation using the ground state calculation if desired

		if run_excited_state_from_optimised_ground_structure and (main_calculation_type_name == 'eGS_gGS'):
			submitSL.write('\t# ----------------------------\n')
			submitSL.write('\t# Submit the excited state calculation using the optimised ground state to slurm.\n')
			submitSL.write('\t\n')
			submitSL.write('\tget_excited_state_job_using_optimised_ground_state_GAUSSIAN.py '+str(optimisation_name_DFT_main_opt)+f'.log\n')
			submitSL.write('\t\n')
			submitSL.write('\tcd ../excited_structure\n')
			submitSL.write('\tmv eES_gES_main_opt_submit.sl.inactive eES_gES_main_opt_submit.sl\n')
			submitSL.write('\tsubmit_slurm_job.py eES_gES_main_opt_submit.sl\n')
			submitSL.write('\t\n')
			submitSL.write('\t# ----------------------------\n')

		# End of sl file: print to terminal the job has finished. 

		submitSL.write('\t# ----------------------------\n')
		submitSL.write('\techo "End of job"\n')
		submitSL.write('\t# ----------------------------\n')

		submitSL.write('fi\n')

# ------------------------------------------------------------------------------------------------------------------------------
