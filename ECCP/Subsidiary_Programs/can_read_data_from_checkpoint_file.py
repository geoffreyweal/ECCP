'''
Geoffrey Weal, read_data_from_checkpoint_file.py, 12/1/23

This program is designed to 
'''
import os
from subprocess import run, PIPE

max_counter = 500
def can_read_data_from_checkpoint_file(checkpoint_filepath):
	"""
	This method is designed to determine if information can be read from the checkpoint file

	Parameters
	----------
	checkpoint_filepath : str.
		This is the path to the checkpoint file

	Returns
	-------
	True if the calculation geometry and SCF data can be obtained from the checkpoint file. False otherwise.
	"""

	raise Exception('This method may be defuncted now.')

	# First, determine if the checkpoint file exists
	if not os.path.exists(checkpoint_filepath):
		return False

	# Second, use the chkchk utility provided by Gaussian to determine if the checkpoint file can be read from.
	# 2.1: Obtain the data from the chkchk utility.
	execute_chkchk = run(['chkchk', checkpoint_filepath], stdout=PIPE, stderr=PIPE)
	stdout_result = execute_chkchk.stdout.decode('utf-8')
	stderr_result = execute_chkchk.stderr.decode('utf-8')

	# 2.2: If there was an issue with chkchk, return False and don't worry about the checkpoint file.
	if not execute_chkchk.returncode == 0:
		return False

	# 2.3: Read lines from chkchk output and see if this checkpoint file can be used.
	atomic_coordinates_present = False
	SCF_restart_data_present   = False
	job_successfully_finished = False
	for line in stdout_result.split('\n'):
		if 'Atomic coordinates present.' in line:
			atomic_coordinates_present = True
		if 'SCF restart data present.' in line:
			SCF_restart_data_present = True
		if 'This checkpoint is from a job which completed successfully.' in line:
			job_successfully_finished = True

	# Third, return the result from 2.
	return atomic_coordinates_present and (SCF_restart_data_present or job_successfully_finished)


