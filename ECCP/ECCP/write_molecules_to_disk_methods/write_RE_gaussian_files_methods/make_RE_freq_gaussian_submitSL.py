"""
make_RE_freq_gaussian_submitSL.py, Geoffrey Weal, 8/5/22

This method will write the submit.sl file in parallel for performing frequency calculations. 
"""
from copy                                                     import deepcopy
from ECCP.ECCP.write_molecules_to_disk_methods.shared_methods import convert_dict_for_bash_input
from ECCP.ECCP.write_molecules_to_disk_methods.shared_methods import slurmSL_header, load_gaussian_programs, make_gaussian_temp_folder, remove_gaussian_temp_files

def make_RE_freq_gaussian_submitSL(freq_calc_filename, local_path, perform_TD, perform_raman, functional, basis_set, gaussian_parameters, cpus_per_task, mem, time, partition='parallel', constraint=None, nodelist=None, exclude=None, email='', python_version='python/3.8.1', gaussian_version='gaussian/g16', temp_folder_path=None, remove_chk_file=False):
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
	remove_chk_file : bool.
		This variable indicates if you want to remove the chk file afterwards. Default: False
	"""

	# Get version of Gaussian to use
	gaussian_version_suffix = str(gaussian_version.split('/')[-1])

	# Make changes to gaussian_parameters to change gaussian filenames to include "freq" in their name.
	gaussian_parameters_for_freq = deepcopy(gaussian_parameters)
	for gaussian_file in ['chk', 'd2e', 'int', 'rwf', 'skr']:
		if gaussian_file in gaussian_parameters_for_freq:
			gaussian_filename = gaussian_parameters_for_freq[gaussian_file]
			gaussian_filename = gaussian_filename.split('.')
			gaussian_filename = gaussian_filename[0]+'_freq.'+gaussian_filename[1]
			gaussian_parameters_for_freq[gaussian_file] = gaussian_filename

	# Remove pre_method and pre-basis if they exist on gaussian_parameters
	if ('pre_method' in gaussian_parameters_for_freq):
		del gaussian_parameters_for_freq['pre_method']
	if ('pre_basis'  in gaussian_parameters_for_freq):
		del gaussian_parameters_for_freq['pre_basis']

	# Remove oldchk in gaussian_parameters_for_freq if it exists, as we will want to recalculate the full Hessian rather than the one from the optimisation that may not use the full hessian.
	if 'oldchk' in gaussian_parameters_for_freq:
		del gaussian_parameters_for_freq['oldchk']

	# If not performing an excited state calculation, remove the 'td_settings' entry so there is not confusion in the file. 
	if (not perform_TD) and ('td_settings' in gaussian_parameters_for_freq):
		del gaussian_parameters_for_freq['td_settings']

	# create name for job
	freq_calc_name = freq_calc_filename.replace('.gjf','')
	name = '-'.join(local_path.split('/')[-4:-1])+'-ReorgE_freq-'+str(freq_calc_name)

	# writing the submit.sl script
	with open(local_path+'/'+str(freq_calc_name)+"_freq_submit.sl", "w") as submitSL:
		slurmSL_header(submitSL, name, mem, partition, constraint, nodelist, time, email, cpus_per_task=cpus_per_task, exclude=exclude)
		make_gaussian_temp_folder(submitSL, temp_folder_path)
		load_gaussian_programs(submitSL, gaussian_version, python_version)
		submitSL.write('# ============================\n')
		submitSL.write('# Prevent the frequency job from running if it is already running or has already run.\n')
		submitSL.write('\n')
		submitSL.write('if ! [[ -f '+str(freq_calc_name)+'_freq.log'+' ]]\n')
		submitSL.write('then\n')
		submitSL.write('\t\n')
		submitSL.write('\t# ============================\n')
		submitSL.write('\t# Settings for creating new gjf files\n')
		submitSL.write('\t\n')
		submitSL.write('\tgaussian_parameters='+str(convert_dict_for_bash_input(gaussian_parameters_for_freq))+'\n')
		submitSL.write('\t\n')
		submitSL.write('\t# ============================\n')
		submitSL.write('\t# Extract optimised structure and place it into frequency calculation gaussian input (.gjf) file.\n')
		submitSL.write('\t\n')
		submitSL.write('\tget_freq_RE_Gaussian_input_file.py '+str(freq_calc_name)+'_main_opt.log '+str(freq_calc_name)+'_freq.gjf '+str(perform_TD)+' '+str(perform_raman)+' "${gaussian_parameters}"\n')
		submitSL.write('\t\n')
		submitSL.write('\t# ============================\n')
		submitSL.write('\t# Perform frequency calculation\n')
		submitSL.write('\t\n')
		submitSL.write('\tsrun '+gaussian_version_suffix+' < '+str(freq_calc_name)+'_freq.gjf > '+str(freq_calc_name)+'_freq.log\n')
		submitSL.write('\t\n')
		submitSL.write('\t# ============================\n')
		if not remove_chk_file:
			# This method will move the frequency calculation checkpoint file to its own mass folder point.
			submitSL.write('\t# Move the frequency calculation checkpoint file (gaussian_freq.chk) file to its own mass folder\n')
			submitSL.write('\t\n')
			submitSL.write('\tmove_gaussian_freq_chk_file_to_storage_folder.py\n')
			submitSL.write('\t\n')
			submitSL.write('\t# ============================\n')
		remove_gaussian_temp_files(submitSL, gaussian_parameters_for_freq, temp_folder_path, remove_chk_file=remove_chk_file, remove_temp_folder=True, prefix='\t')
		submitSL.write('\t# ============================\n')
		submitSL.write('\techo "End of job"\n')
		submitSL.write('\t# ============================\n')
		submitSL.write('\t\n')
		submitSL.write('fi\n')

# ------------------------------------------------------------------------------------------------------------------------------

