'''
Geoffrey Weal, ECCP_submit_gaussian_jobs_to_slurm.py, 15/6/2023

This program contains methods for submitting Gaussian jobs to slurm if appropriate to do so.
'''
from ECCP.ECCP_Programs.shared_general_methods.shared_gaussian_methods import did_gaussian_opt_job_complete
#from ECCP.ECCP_Programs.shared_general_methods.shared_gaussian_methods import did_gaussian_job_complete

def general_gaussian_submission(filenames, are_RE_jobs_running_currently):
	"""
	This method is design to determine if you can submit Gaussian ATC or EET calculations.

	Parameters
	----------
	filenames : list of str.
		These are all the filenames of files in the directory you want to look at.
	are_RE_jobs_running_currently : bool
		Has the user already indicates these jobs may already be running on slurm.

	Returns
	-------
	submission_filenames : all the filenames for submission scripts to submit to slurm.
	"""
	
	# First, create the submission_filenames list
	submission_filenames = []

	# Second, if output.log is already in the directory, this job has already run so dont submit it. 
	#         --> Only submit this job if output.log does not already exist.
	if ('output.log' not in filenames): 
		submission_filenames.append('submit.sl')

	# Third, if performing ATC calculations and Gaussian has run but not MultiWFN, submit multiwfn_submit.sl to slurm
	if (not are_RE_jobs_running_currently) and ('output.log' in filenames) and ('output.wfn' in filenames) and ('output.chg' not in filenames): 
		submission_filenames.append('multiwfn_submit.sl')

	# Fourth, return submission_filenames
	return submission_filenames

def RE_GStructure_gaussian_submission(filenames, dirpath, are_RE_jobs_running_currently):
	"""
	This method is design to determine if you can submit Gaussian ground structure reorganisation energy calculations.

	Parameters
	----------
	filenames : list of str.
		These are all the filenames of files in the directory you want to look at.
	dirpath : str.
		The path to the directory the Gaussian and submission files are in.
	are_RE_jobs_running_currently : bool
		Has the user already indicates these jobs may already be running on slurm.

	Returns
	-------
	submission_filenames : all the filenames for submission scripts to submit to slurm.
	"""

	# First, create the submission_filenames list
	submission_filenames = []

	# Second, determine if to submit the GS_GS_submit.sl
	will_perform_preopt = 'eGS_gGS_main_opt_preopt.gjf' in filenames
	begun_main_preopt   = 'eGS_gGS_main_opt_preopt.log' in filenames
	begun_main_opt      = 'eGS_gGS_main_opt.log'        in filenames
	begun_freq          = 'eGS_gGS_freq.log'            in filenames
	begun_SP            = 'eES_gGS.log'                 in filenames

	# Third, determine times to allow GS_GS_submit.sl submission
	has_main_preopt_step_finished = (begun_main_preopt and has_RE_opt_finished('eGS_gGS_main_opt_preopt', dirpath))
	has_main_opt_step_finished    = (begun_main_opt    and has_RE_opt_finished('eGS_gGS_main_opt', dirpath))

	# Fourth, determine what submission files to submit to slurm.
	if will_perform_preopt and not has_main_preopt_step_finished:

		if not begun_main_preopt:

			# 4.1: The pre-optimisation has not begun, so submit the main optimisation submit script (which will run the preoptimisation). 
			if begun_freq or begun_SP:
				raise Exception('Not sure what is happening, the frequency or GS_ES file has run/is running, but the optimisation has not run? optimisation file may have been deleted.') 
			submission_filenames.append('eGS_gGS_main_opt_submit.sl')

	elif not begun_main_opt:

		# 4.2: The main optimisation has not begun, so submit the main optimisation submit script. 
		if begun_freq or begun_SP:
			raise Exception('Not sure what is happening, the frequency or GS_ES file has run/is running, but the optimisation has not run?') 
		submission_filenames.append('eGS_gGS_main_opt_submit.sl')

	elif (not are_RE_jobs_running_currently) and has_main_opt_step_finished:

		# 4.3: The optimisation has completed successfully, so look at submitting frequency and single point Gaussian file.
		if not begun_freq: # or (begun_freq and not did_gaussian_job_complete(dirpath+'/'+'GS_GS_freq.log')):
			submission_filenames.append('eGS_gGS_freq_submit.sl')
		if not begun_SP: # or (begun_SP and not did_gaussian_job_complete(dirpath+'/'+'GS_ES.log')):
			submission_filenames.append('eES_gGS_submit.sl')

	# Fifth, return submission_filenames
	return submission_filenames

def RE_EStructure_gaussian_submission(filenames, dirpath, are_RE_jobs_running_currently):
	"""
	This method is design to determine if you can submit Gaussian excited structure reorganisation energy calculations.

	Parameters
	----------
	filenames : list of str.
		These are all the filenames of files in the directory you want to look at.
	dirpath : str.
		The path to the directory the Gaussian and submission files are in.
	are_RE_jobs_running_currently : bool
		Has the user already indicates these jobs may already be running on slurm.

	Returns
	-------
	submission_filenames : all the filenames for submission scripts to submit to slurm.
	"""

	# First, create the submission_filenames list
	submission_filenames = []

	# Second, if the excited structure folder contains 'eES_gES_main_opt_submit.sl.inactive', ignore submitting this as the file is unactive.
	#         It will be activated once the ground state optimisation has completed. 
	if 'eES_gES_main_opt_submit.sl.inactive' in filenames:
		return submission_filenames

	# Third, determine if to submit the ES_ES_submit.sl
	will_perform_preopt = 'eES_gES_main_opt_preopt.gjf' in filenames
	begun_main_preopt   = 'eES_gES_main_opt_preopt.log' in filenames
	begun_main_opt      = 'eES_gES_main_opt.log'        in filenames
	begun_freq          = 'eES_gES_freq.log'            in filenames
	begun_SP            = 'eGS_gES.log'                 in filenames

	# Fourth, determine times to allow ES_ES_submit.sl submission
	has_main_preopt_step_finished = (begun_main_preopt and has_RE_opt_finished('eES_gES_main_opt_preopt', dirpath))
	has_main_opt_step_finished    = (begun_main_opt    and has_RE_opt_finished('eES_gES_main_opt', dirpath))

	# Fifth, determine what submission files to submit to slurm.
	if will_perform_preopt and not has_main_preopt_step_finished:

		if not begun_main_preopt:

			# 4.1: The preopt or main optimisation has not begun, so submit the main optimisation submit script (this script will also run the preopt calculation if desired). 
			if begun_freq or begun_SP:
				raise Exception('Not sure what is happening, the frequency or GS_ES file has run/is running, but the optimisation has not run? optimisation file may have been deleted.') 
			submission_filenames.append('eES_gES_main_opt_submit.sl')

	elif not begun_main_opt:

		# 4.2: The preopt or main optimisation has not begun, so submit the main optimisation submit script (this script will also run the preopt calculation if desired). 
		if begun_freq or begun_SP:
			raise Exception('Not sure what is happening, the frequency or GS_ES file has run/is running, but the optimisation has not run? optimisation file may have been deleted.') 
		submission_filenames.append('eES_gES_main_opt_submit.sl')

	elif (not are_RE_jobs_running_currently) and has_main_opt_step_finished:
		
		# 4.3: The optimisation has completed successfully, so look at submitting frequency and single point Gaussian file.
		if not begun_freq: # or (begun_freq and not did_gaussian_job_complete(dirpath+'/'+'ES_ES_freq.log')):
			submission_filenames.append('eES_gES_freq_submit.sl')
		if not begun_SP: # or (begun_SP and not did_gaussian_job_complete(dirpath+'/'+'ES_GS.log')):
			submission_filenames.append('eGS_gES_submit.sl')

	# Sixth, return submission_filenames
	return submission_filenames

# ==============================================================================================================

def has_RE_opt_finished(reorg_energy_name, root):
	'''
	This method will check if reorganisation energy calculation have finished successfully or not. 

	Parameters
	----------
	reorg_energy_name : str.
		This is the name of the reorganisation energy process
	root : str.
		This is the path to the log file for this process

	Returns
	-------
	True if the Gaussian job finished successfully and found a stationary point/local minimum.
	'''
	path_to_outputLOG = root + '/' + reorg_energy_name + '.log'
	did_main_opt_finish, has_fully_converged, converged_image_index =  did_gaussian_opt_job_complete(path_to_outputLOG)
	return did_main_opt_finish

# ==============================================================================================================



