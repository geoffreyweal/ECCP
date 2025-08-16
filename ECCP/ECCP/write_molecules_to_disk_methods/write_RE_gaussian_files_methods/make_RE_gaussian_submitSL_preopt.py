"""
make_RE_gaussian_submitSL_preopt.py, Geoffrey Weal, 8/5/22

This method will write the submit.sl file in parallel
"""
from copy                                                     import deepcopy
from ECCP.ECCP.write_molecules_to_disk_methods.shared_methods import convert_dict_for_bash_input
from ECCP.ECCP.write_molecules_to_disk_methods.shared_methods import slurmSL_header, load_gaussian_programs, make_gaussian_temp_folder, remove_gaussian_temp_files


def make_RE_gaussian_submitSL_preopt(main_calculation_type_name, optimisation_filename_DFT_main_opt, single_point_filename, local_path, perform_excited_state_calc, functional, basis_set, gaussian_parameters, cpus_per_task, mem, time, partition='parallel', constraint=None, nodelist=None, exclude=None, email='', python_version='python/3.8.1', gaussian_version='gaussian/g16', temp_folder_path=None):
	"""
	This method will write the submit.sl file in parallel

	Parameters
	----------
	main_calculation_type_name : str.
		This is the name of the calculation type (either GS_GS or ES_ES).
	optimisation_filename_DFT_main_opt : str. 
		This is the name of the optimisation file for performing the main optimisation process with DFT
	single_point_filename : str. 
		This is the name of the single point calculation file.
	local_path : str. 
		This is the location to save this submit.sl file to
	perform_excited_state_calc : bool.
		This boolean indicates if you are performing a excited state optimisation.
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
	gaussian_parameters_main_preopt = deepcopy(gaussian_parameters)
	for gaussian_file in ['chk', 'd2e', 'int', 'rwf', 'skr']:
		if gaussian_file in gaussian_parameters_main_preopt:
			gaussian_filename = gaussian_parameters_main_preopt[gaussian_file]
			gaussian_filename = gaussian_filename.split('.')
			gaussian_filename = gaussian_filename[0]+'_main_preopt.'+gaussian_filename[1]
			gaussian_parameters_main_preopt[gaussian_file] = gaussian_filename

	# If performing ground state optimisation (not excited state), remove component to indicate td setting from gaussian_parameters. 
	if (not perform_excited_state_calc) and ('td_settings' in gaussian_parameters_main_preopt):
		del gaussian_parameters_main_preopt['td_settings']

	# Make changes to gaussian_parameters to change gaussian filenames to include "main_opt" in their name.
	gaussian_parameters_main_opt = deepcopy(gaussian_parameters)
	for gaussian_file in ['chk', 'd2e', 'int', 'rwf', 'skr']:
		if gaussian_file in gaussian_parameters_main_opt:
			gaussian_filename = gaussian_parameters_main_opt[gaussian_file]
			gaussian_filename = gaussian_filename.split('.')
			gaussian_filename = gaussian_filename[0]+'_main_opt.'+gaussian_filename[1]
			gaussian_parameters_main_opt[gaussian_file] = gaussian_filename

	if ('pre_method' in gaussian_parameters_main_opt):
		del gaussian_parameters_main_opt['pre_method']
	if ('pre_basis'  in gaussian_parameters_main_opt):
		del gaussian_parameters_main_opt['pre_basis']

	# If performing ground state optimisation (not excited state), remove component to indicate td setting from gaussian_parameters. 
	if (not perform_excited_state_calc) and ('td_settings' in gaussian_parameters_main_opt):
		del gaussian_parameters_main_opt['td_settings']

	# create names for job.
	optimisation_name_DFT_main_opt   = optimisation_filename_DFT_main_opt.replace('.gjf','')
	single_point_name                = single_point_filename.replace('.gjf','')
	name = '-'.join(local_path.split('/')[-4:-1])+'-ReorgE_main-'+str(main_calculation_type_name)

	# writing the submit.sl script
	with open(local_path+'/'+str(main_calculation_type_name)+"_main_opt_submit.sl", "w") as submitSL:
		slurmSL_header(submitSL, name, mem, partition, constraint, nodelist, time, email, cpus_per_task=cpus_per_task, exclude=exclude)
		make_gaussian_temp_folder(submitSL, temp_folder_path)
		load_gaussian_programs(submitSL, gaussian_version, python_version)
		
		# Perform the preoptimisation optimisation.

		submitSL.write('# ----------------------------\n')
		submitSL.write('\n')
		submitSL.write('if [ ! -f '+str(optimisation_name_DFT_main_opt)+'_preopt.log'+' ]; then\n')

		# Perform the proper DFT optimisation

		submitSL.write('\t# ============================\n')
		submitSL.write('\t# Perform the pre geometry optimisation calculation with the desired functional and basis set.\n')
		submitSL.write('\t\n')
		submitSL.write('\techo "Performing pre optimisation calculation"\n')
		submitSL.write('\tsrun '+gaussian_version_suffix+' < '+str(optimisation_filename_DFT_main_opt.replace('.gjf','_preopt.gjf'))+' > '+str(optimisation_name_DFT_main_opt)+'_preopt.log\n')
		submitSL.write('\techo "Finished pre optimisation calculation"\n')
		submitSL.write('\t\n')
		submitSL.write('\t# ============================\n')
		submitSL.write('\t\n')

		# Perform final tasks and submit for frequency and single point analysis

		remove_gaussian_temp_files(submitSL, gaussian_parameters_main_preopt, temp_folder_path, remove_chk_file=True, remove_temp_folder=True, prefix='\t')
		submitSL.write('\t# ----------------------------\n')
		submitSL.write('\t# Submit the Frequency and single point calculations to slurm.\n')
		submitSL.write('\t\n')
		submitSL.write('\tsubmit_slurm_job.py '+str(main_calculation_type_name)+'_main_opt_submit.sl\n')
		submitSL.write('\t\n')
		submitSL.write('\t# ----------------------------\n')
		submitSL.write('\techo "End of job"\n')
		submitSL.write('\t# ----------------------------\n')
		submitSL.write('\t\n')

		# Run the process for the main optimisation.

		submitSL.write('elif [ -f '+str(optimisation_name_DFT_main_opt)+'_preopt.log'+' ] && [ ! -f '+str(optimisation_name_DFT_main_opt)+'.log'+' ]; then\n')
		submitSL.write('\t# ============================\n')
		submitSL.write('\t# Perform the geometry optimisation calculation with the desired functional and basis set.\n')
		submitSL.write('\t\n')
		submitSL.write('\tdid_converged=$(has_gaussian_optimisation_converged.py '+str(optimisation_name_DFT_main_opt)+'_preopt.log)\n')
		#submitSL.write('\thas_gaussian_optimisation_converged.py '+str(optimisation_name_DFT_main_opt)+'_preopt.log\n')
		#submitSL.write('\tdid_converged=$?\n')
		submitSL.write('\tif [ $did_converged -eq 0 ]; then\n')
		submitSL.write('\t\t\n')

		# Perform the main DFT optimisation

		submitSL.write('\t\t# ============================\n')
		submitSL.write('\t\t# Settings for creating new gjf files\n')
		submitSL.write('\t\t\n')
		submitSL.write('\t\tgaussian_parameters='+str(convert_dict_for_bash_input(gaussian_parameters_main_opt))+'\n')
		submitSL.write('\t\t\n')
		submitSL.write('\t\t# ============================\n')
		submitSL.write('\t\t# Extract optimised structure and place it into frequency calculation gaussian input (.gjf) file.\n')
		submitSL.write('\t\t\n')
		submitSL.write('\t\tget_optimisation_RE_Gaussian_input_file.py '+str(optimisation_name_DFT_main_opt)+'_preopt.log '+str(optimisation_filename_DFT_main_opt)+' "${gaussian_parameters}"\n')

		submitSL.write('\t\t# ============================\n')
		submitSL.write('\t\t# Perform the main geometry optimisation calculation with the desired functional and basis set.\n')
		submitSL.write('\t\t\n')
		submitSL.write('\t\techo "Performing main optimisation calculation"\n')
		submitSL.write('\t\tsrun '+gaussian_version_suffix+' < '+str(optimisation_filename_DFT_main_opt)+' > '+str(optimisation_name_DFT_main_opt)+'.log\n')
		submitSL.write('\t\techo "Finished main optimisation calculation"\n')
		submitSL.write('\t\t\n')
		submitSL.write('\t\t# ============================\n')

		# Perform final tasks and submit for frequency and single point analysis

		remove_gaussian_temp_files(submitSL, gaussian_parameters_main_opt, temp_folder_path, remove_chk_file=True, remove_temp_folder=True, prefix='\t\t')
		submitSL.write('\t\t# ----------------------------\n')
		submitSL.write('\t\t# Submit the Frequency and single point calculations to slurm.\n')
		submitSL.write('\t\t\n')
		submitSL.write('\t\tsubmit_slurm_job.py '+str(main_calculation_type_name)+'_freq_submit.sl\n')
		submitSL.write('\t\tsubmit_slurm_job.py '+str(single_point_name)+'_submit.sl\n')
		submitSL.write('\t\t\n')
		submitSL.write('\t\t# ----------------------------\n')
		submitSL.write('\t\techo "End of job"\n')
		submitSL.write('\t\t# ----------------------------\n')
		submitSL.write('\t\n')

		submitSL.write('\telse\n')
		submitSL.write('\t\t\n')
		submitSL.write('\t\t# ----------------------------\n')
		submitSL.write('\t\t# Indicate that the preoptimisation did not converge.\n')
		submitSL.write('\t\t\n')
		submitSL.write('\t\techo "The preoptimisation did not converge."\n')
		submitSL.write('\t\t\n')
		submitSL.write('\tfi\n')
		submitSL.write('\n')

		submitSL.write('fi\n')

# ------------------------------------------------------------------------------------------------------------------------------
