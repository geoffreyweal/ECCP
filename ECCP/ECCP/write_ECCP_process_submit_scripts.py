"""
write_ECCP_process_submit_scripts.py, Geoffrey Weal, 8/5/22

This script is designed to create submit.sl scripts for submitting ECCP jobs to slurm for processing ATC, RE, EET, Eigendata, and ICT data.

These jobs can take a long time to run if one is wanting to obtain overlap orbtials and MO coefficients and energies due to the sizes of these matrices, so these submit scripts can be useful for getting slurm to do this for us.
"""

from copy import deepcopy

def write_ECCP_process_ATC_submit_script(Unique_ATC_Gaussian_Jobs_folder, all_submission_information_for_ATCs):
	"""
	This method is designed to create a submit.sl script for performing ECCP process_ATC program in slurm.

	Parameters
	----------
	Unique_ATC_Gaussian_Jobs_folder : str.
		This is the path to the "Unique_ATC_Gaussian_Jobs" folder
	all_submission_information_for_ATCs : list of dict.
		This is the information about how to make a submit script for ATC jobs. This is used to make a basis submit.sl file. Only all_submission_information_for_ATCs[0] is used in this method. 
	"""
	submission_information_for_ATCs_for_ECCP_prep = deepcopy(all_submission_information_for_ATCs[0])
	if 'ntasks' in submission_information_for_ATCs_for_ECCP_prep:
		submission_information_for_ATCs_for_ECCP_prep['cpus_per_task'] = submission_information_for_ATCs_for_ECCP_prep['ntasks']
		del submission_information_for_ATCs_for_ECCP_prep['ntasks']
	submission_information_for_ATCs_for_ECCP_prep['cpus_per_task'] = 1
	submission_information_for_ATCs_for_ECCP_prep['mem'] = '10GB'
	if 'remove_chk_file' in submission_information_for_ATCs_for_ECCP_prep:
		del submission_information_for_ATCs_for_ECCP_prep['remove_chk_file']
	make_submitSL(Unique_ATC_Gaussian_Jobs_folder, 'ECCP_process_ATC', 'ECCP -T process_ATC', **submission_information_for_ATCs_for_ECCP_prep)

def write_ECCP_process_RE_submit_script(Unique_RE_FC_Gaussian_Jobs_folder,  all_submission_information_for_REs):
	"""
	This method is designed to create a submit.sl script for performing ECCP process_RE program in slurm.

	Parameters
	----------
	Unique_RE_FC_Gaussian_Jobs_folder : str.
		This is the path to the "Unique_RE_FC_Gaussian_Jobs" folder
	all_submission_information_for_REs_and_FCs : list of dict.
		This is the information about how to make a submit script for RE jobs. This is used to make a basis submit.sl file. Only all_submission_information_for_REs_and_FCs[0] is used in this method. 
	"""
	submission_information_for_REs_for_ECCP_prep = deepcopy(all_submission_information_for_REs[0])
	if 'ntasks' in submission_information_for_REs_for_ECCP_prep:
		submission_information_for_REs_for_ECCP_prep['cpus_per_task'] = submission_information_for_REs_for_ECCP_prep['ntasks']
		del submission_information_for_REs_for_ECCP_prep['ntasks']
	submission_information_for_REs_for_ECCP_prep['cpus_per_task'] = 1
	submission_information_for_REs_for_ECCP_prep['mem'] = '10GB'
	if 'remove_chk_file' in submission_information_for_REs_for_ECCP_prep:
		del submission_information_for_REs_for_ECCP_prep['remove_chk_file']
	make_submitSL(Unique_RE_FC_Gaussian_Jobs_folder,  'ECCP_process_RE',  'ECCP -T process_RE',  **submission_information_for_REs_for_ECCP_prep)

def write_ECCP_process_FC_submit_script(Unique_RE_FC_Gaussian_Jobs_folder,  all_submission_information_for_FCs):
	"""
	This method is designed to create a submit.sl script for performing ECCP process_FC program in slurm.

	Parameters
	----------
	Unique_RE_FC_Gaussian_Jobs_folder : str.
		This is the path to the "Unique_RE_FC_Gaussian_Jobs" folder
	all_submission_information_for_REs_and_FCs : list of dict.
		This is the information about how to make a submit script for RE jobs. This is used to make a basis submit.sl file. Only all_submission_information_for_REs_and_FCs[0] is used in this method. 
	"""
	submission_information_for_FCs_for_ECCP_prep = deepcopy(all_submission_information_for_FCs[0])
	if 'ntasks' in submission_information_for_FCs_for_ECCP_prep:
		submission_information_for_FCs_for_ECCP_prep['cpus_per_task'] = submission_information_for_FCs_for_ECCP_prep['ntasks']
		del submission_information_for_FCs_for_ECCP_prep['ntasks']
	submission_information_for_FCs_for_ECCP_prep['cpus_per_task'] = 1
	submission_information_for_FCs_for_ECCP_prep['mem'] = '10GB'
	if 'remove_chk_file' in submission_information_for_FCs_for_ECCP_prep:
		del submission_information_for_FCs_for_ECCP_prep['remove_chk_file']
	make_submitSL(Unique_RE_FC_Gaussian_Jobs_folder,  'ECCP_process_FC',  'ECCP -T process_FC',  **submission_information_for_FCs_for_ECCP_prep)

def write_ECCP_process_EET_submit_script(Unique_EET_Gaussian_Jobs_folder, all_submission_information_for_EETs):
	"""
	This method is designed to create a submit.sl script for performing ECCP process_EET program in slurm.

	Parameters
	----------
	Unique_EET_Gaussian_Jobs_folder : str.
		This is the path to the "Unique_EET_Gaussian_Jobs_folder" folder
	all_submission_information_for_EETs : list of dict.
		This is the information about how to make a submit script for EET jobs. This is used to make a basis submit.sl file. Only all_submission_information_for_EETs[0] is used in this method. 
	"""
	submission_information_for_EETs_for_ECCP_prep = deepcopy(all_submission_information_for_EETs[0])
	if 'ntasks' in submission_information_for_EETs_for_ECCP_prep:
		submission_information_for_EETs_for_ECCP_prep['cpus_per_task'] = submission_information_for_EETs_for_ECCP_prep['ntasks']
		del submission_information_for_EETs_for_ECCP_prep['ntasks']
	submission_information_for_EETs_for_ECCP_prep['cpus_per_task'] = 1
	submission_information_for_EETs_for_ECCP_prep['mem'] = '10GB'
	if 'remove_chk_file' in submission_information_for_EETs_for_ECCP_prep:
		del submission_information_for_EETs_for_ECCP_prep['remove_chk_file']
	make_submitSL(Unique_EET_Gaussian_Jobs_folder, 'ECCP_process_EET', 'ECCP -T process_EET', **submission_information_for_EETs_for_ECCP_prep)

def write_ECCP_process_Eigendata_submit_script(Unique_Eigendata_Gaussian_Jobs_folder, all_submission_information_for_Eigendata):
	"""
	This method is designed to create a submit.sl script for performing ECCP process_Eigendata program in slurm.

	Parameters
	----------
	Unique_Eigendata_Gaussian_Jobs_folder : str.
		This is the path to the "Unique_Eigendata_Gaussian_Jobs_folder" folder
	all_submission_information_for_Eigendata : list of dict.
		This is the information about how to make a submit script for obtaining eigendata. This is used to make a basis submit.sl file. Only all_submission_information_for_Eigendata[0] is used in this method. 
	"""
	submission_information_for_Eigendata_for_ECCP_prep = deepcopy(all_submission_information_for_Eigendata[0])
	if 'ntasks' in submission_information_for_Eigendata_for_ECCP_prep:
		submission_information_for_Eigendata_for_ECCP_prep['cpus_per_task'] = submission_information_for_Eigendata_for_ECCP_prep['ntasks']
		del submission_information_for_Eigendata_for_ECCP_prep['ntasks']
	submission_information_for_Eigendata_for_ECCP_prep['cpus_per_task'] = 1
	submission_information_for_Eigendata_for_ECCP_prep['mem'] = '10GB'
	if 'remove_chk_file' in submission_information_for_Eigendata_for_ECCP_prep:
		del submission_information_for_Eigendata_for_ECCP_prep['remove_chk_file']
	make_submitSL(Unique_Eigendata_Gaussian_Jobs_folder, 'ECCP_process_Eigendata', 'ECCP -T process_Eigendata', **submission_information_for_Eigendata_for_ECCP_prep)

def write_ECCP_process_ICT_submit_script(Unique_ICT_Gaussian_Jobs_folder, all_submission_information_for_ICT):
	"""
	This method is designed to create a submit.sl script for performing ECCP process_ICT program in slurm.

	Parameters
	----------
	Unique_ICT_Gaussian_Jobs_folder : str.
		This is the path to the "Unique_ICT_Gaussian_Jobs_folder" folder
	all_submission_information_for_ICT : list of dict.
		This is the information about how to make a submit script for obtaining ICT. This is used to make a basis submit.sl file. Only all_submission_information_for_ICTs[0] is used in this method. 
	"""
	submission_information_for_ICT_for_ECCP_prep = deepcopy(all_submission_information_for_ICT[0])
	if 'ntasks' in submission_information_for_ICT_for_ECCP_prep:
		submission_information_for_ICT_for_ECCP_prep['cpus_per_task'] = submission_information_for_ICT_for_ECCP_prep['ntasks']
		del submission_information_for_ICT_for_ECCP_prep['ntasks']
	submission_information_for_ICT_for_ECCP_prep['cpus_per_task'] = 1
	submission_information_for_ICT_for_ECCP_prep['mem'] = '10GB'
	if 'remove_chk_file' in submission_information_for_ICT_for_ECCP_prep:
		del submission_information_for_ICT_for_ECCP_prep['remove_chk_file']
	make_submitSL(Unique_ICT_Gaussian_Jobs_folder, 'ECCP_process_ICT', 'ECCP -T process_ICT', **submission_information_for_ICT_for_ECCP_prep)

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def make_submitSL(local_path,submit_name,input_string,cpus_per_task,mem,time,partition='parallel',constraint=None,nodelist=None,exclude=None,email='',python_version='Python/3.9.5',gaussian_version='gaussian/g16',orca_version='ORCA/5.0.4',log_filename='gaussian.log',gaussian_scratch_name=None):
	"""
	This method will create the submit.sl file for submitting ECCP jobs to slurm.

	"""
	with open(local_path+'/'+submit_name+"_submit.sl", "w") as submitSL:
		submitSL.write('#!/bin/bash -e\n')
		submitSL.write('#SBATCH --job-name=' + str(submit_name) + '\n')
		submitSL.write('#SBATCH --cpus-per-task=' + str(cpus_per_task) + '\n')
		submitSL.write('#SBATCH --mem=' + str(mem) + '\n')
		submitSL.write('#SBATCH --partition=' + str(partition) + '\n')
		if (constraint is not None):
			submitSL.write('#SBATCH --constraint=' + str(constraint) + '\n')
		if (nodelist is not None):
			submitSL.write('#SBATCH --nodelist=' + str(nodelist) + '\n')
		if (exclude is not None):
			submitSL.write('#SBATCH --exclude=' + str(exclude) + '\n')
		submitSL.write('#SBATCH --time=' + str(time) + '     # Walltime\n')
		submitSL.write('#SBATCH --output=slurm-%j.out      # %x and %j are replaced by job name and ID'+'\n')
		submitSL.write('#SBATCH --error=slurm-%j.err'+'\n')
		if not email == '':
		    submitSL.write('#SBATCH --mail-user=' + str(email) + '\n')
		    submitSL.write('#SBATCH --mail-type=ALL\n')
		submitSL.write('\n')
		submitSL.write('module load '+str(python_version)+'\n')
		submitSL.write('\n')
		submitSL.write(input_string+'\n')
		submitSL.write('\n')
		submitSL.write('echo "FINISHED RUNNING THIS SCRIPT"\n')

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------






