'''
Geoffrey Weal, determine_quantum_computing_software_type.py, 16/06/2021

This program is designed
'''
from ECCP.ECCP_Programs.shared_general_methods.shared_orca_methods import is_orca_output_file

def determine_quantum_computing_software_type(dirpath, filenames):
	"""
	This method is designed to determine if the job(s) you are trying are Gaussian or ORCA jobs.

	Parameters
	----------
	dirpath : str.
		The string to the path that a potential Gaussian or ORCA job is in.
	filenames : list of str.
		These are the files in the dirpath directory.

	Returns
	-------
	If the job in this directory contains a Gaussian or ORCA job. 
	"""

	# First, define booleans for if we are looking at a Gaussian or ORCA job
	found_gaussian_file = False
	found_orca_file     = False

	# Second, look through all files in filenames for keywrods in the title that indicate it contains Gaussian or ORCA files.
	for filename in filenames:
		if filename.endswith('.gjf') or filename.endswith('.log') or (filename == 'gaussian_parameters_ES.txt'):
			found_gaussian_file = True
		if filename.endswith('.inp') or is_orca_output_file(dirpath, filename) or (filename == 'orca_parameters_ES.txt'):
			found_orca_file = True

	# Third,  return whether we are wanting to submit a Gaussian or ORCA job(s). 
	if found_gaussian_file and found_orca_file:
		raise Exception('Error: Both Gaussian and ORCA files were found when trying to submit jobs to slurm. Check this. Path: '+str(dirpath))
	elif found_gaussian_file:
		return 'Gaussian'
	elif found_orca_file:
		return 'ORCA'
	raise Exception('Error: Could not find either Gaussian nor ORCA files when trying to submit jobs to slurm. Check this. Path: '+str(dirpath))




