'''
analyse_ATC_output.py, Geoffrey Weal, 29/12/22

This method is designed to check if an ATC Gaussian job has completed or not.
'''
import os
from ECCP.ECCP_Programs.shared_general_methods.shared_general_methods import reverse_readline

def analyse_ATC_output(software_type, path_to_outputLOG, path_to_outputCHG):
	"""
	This method is designed to check if an ATC Gaussian job has completed or not.

	Parameters
	----------
    software_type : str.
        This is the type of software that was used to perform these calculations.
	path_to_outputLOG : str.
		This is the path to the output.log/.out file for this ATC job.
	path_to_outputCHG : str.
		This is the path to the output.chg file for this ATC job.

	Returns
	-------
	The results of the job: str.
		* 'NBY' : Not begun yet.
		* 'NC'  : Not complete.
		* 'NCh' : Quantum calculation has complete, but MultiWFN has not obtained the ATC .chg file.
		* 'C'   : Complete.
	"""

	if software_type == "Gaussian":

		# First, if the output.log file does not exist, return this
		if not os.path.exists(path_to_outputLOG):
			return 'NBY' # Not begun yet.

		# Second, check to see if this ATC calculation has finished successfully.
		did_terminate_normally = False
		counter = 0
		for line in reverse_readline(path_to_outputLOG):
			# Check if the gaussian file has terminated normally
			if 'Normal termination of Gaussian' in line:
				did_terminate_normally = True
				break
			if counter >= 20 and not did_terminate_normally:
				break

		# Third, determine if the job completed
		if not did_terminate_normally:
			return 'NC' # Not complete

		# Fourth, determine if the chg file completed successfully. 
		if not os.path.exists(path_to_outputCHG):
			return 'NCh' # Not complete

		# Fifth, if at this point, this indicates the program finished successfully. 
		return 'C' # Completed

	elif software_type == "ORCA":

		raise Exception('add ORCA stuff')




