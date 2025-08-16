"""
shared_methods.py, Geoffrey Weal, 8/5/22

This script contains methods that are used by the write_molecules_to_disk.py script and other scripts in the write_molecules_to_disk_methods folder.
"""
# -----------------------------------------------------------------

def change_folder_name_components(input_string):
	"""
	This method will change certain characters that should be used in file and folder names into ones that can be.

	Parameters
	----------
	input_string : str.
		This is the string to modify

	Returns
	-------
	Returns a string that is safe to use in a file or folder name. 
	"""
	return str(input_string).replace('-','_').replace('+','plus').replace('(','').replace(')','').replace(',','_')

# -----------------------------------------------------------------
# 'I': 1.98, 
no_preset_vdW_radius = {'Br': 1.85, 'Se': 1.90, 'Te': 2.06, 'As': 1.85}
def input_commands_for_multiwfn(wfn_filename='output.wfn',molecule_symbols=[]):
	"""
	This method is designed to provide the commands for performing the ATC calculations in Multiwfn.

	Parameters
	----------
	wfn_filename : str.
		This is the name of the wfn file
	molecule_symbols : list of str.
		This list contains all the elements in the molecule.

	Returns
	-------
	Returns a string of the commands needed to run Multiwfn in the terminal. This is required to run Multiwfn in slurm to perform atomic transition charge (ATC) calculations using the output.wfn file created with Gaussian. 
	"""
	
	# First, obtain the commands needed to run ATC calculations in Multiwfn in the terminal.
	multiwfn_echo_commands = [7,12,5,3,1]
	for element in sorted(set(molecule_symbols)):
		if element in no_preset_vdW_radius.keys():
			multiwfn_echo_commands.append(no_preset_vdW_radius[element])
			#multiwfn_echo_commands.append('')
	multiwfn_echo_commands += ['y',0,0,'q']
	multiwfn_echo_commands = [str(command) for command in multiwfn_echo_commands]
	multiwfn_echo_commands = '\\n'.join(multiwfn_echo_commands)

	# Second, get the string to run in the terminal
	multiwfn_commands = 'echo -e "'+str(multiwfn_echo_commands)+'" | Multiwfn '+str(wfn_filename)

	# Third, return the command string to run ATC calculations in Multiwfn in the terminal.
	return multiwfn_commands

# -----------------------------------------------------------------

def convert_dict_for_bash_input(gaussian_parameters):
	"""
	This method will convert the gaussian_parameters into a string that can be placed in a .sh file or submit.sl file for slurm.

	Parameters
	----------
	gaussian_parameters : dict.
		This is the dictionary that contains the inputs required to run the Gaussian calculation.

	Returns
	-------
	toString : str.
		This is the gaussian_parameters dictionary is a string form that can be placed in a .sh file or submit.sl file for slurm.
	"""
	toString = '"{'
	counter = 0
	for key, value in gaussian_parameters.items():
		toString += "'"+str(key)+"': "
		if isinstance(value,str):
			toString += "'"+str(value)+"'"
		else:
			toString += str(value)
		if counter < len(gaussian_parameters)-1:
			toString += ", "
		counter += 1
	toString += '}"'
	return toString

# -----------------------------------------------------------------

def slurmSL_header(submitSL, name, mem, partition, constraint, nodelist, time, email, cpus_per_task=None, ntasks=None, nodes=None, exclude=None):
	"""
	This method provides the header information for the slurm file. 

	Parameters
	----------
	submitSL : open()
		This is the slurm submit file to add commands to.
	name : str. 
		This is the name of the job.
	cpus_per_task : int. 
		This is the number of cpus assigned to this/these Gaussian jobs.
	mem : str. 
		This is the amount of memory assigned to this job.
	partition : str. 
		This is the partition to assign this job to.
	constraint : str. 
		Not sure 
	nodelist : str. 
		These are the nodes to assign this job to.
	time : str. 
		This is the amount of time to assign this job to.
	email : str. 
		This is the email to send job information to.
	"""
	submitSL.write('#!/bin/bash -e\n')
	submitSL.write('#SBATCH --job-name=' + str(name) + '\n')
	if (cpus_per_task is not None) and (ntasks is not None):
		raise Exception('Error: Your slurm submission dictionary can include only either cpus_per_task or ntasks. Remove one of these and try again: cpus_per_task='+str(cpus_per_task)+', ntasks='+str(ntasks))
	elif (cpus_per_task is not None):
		submitSL.write('#SBATCH --cpus-per-task=' + str(cpus_per_task) + '\n')
	elif (ntasks is not None):
		submitSL.write('#SBATCH --ntasks=' + str(ntasks) + '\n')
	else:
		raise Exception('Error: Your slurm submission dictionary needs to include either either cpus_per_task or ntasks.')

	submitSL.write('#SBATCH --mem=' + str(mem) + '\n')
	submitSL.write('#SBATCH --partition=' + str(partition) + '\n')
	if (constraint is not None):
		submitSL.write('#SBATCH --constraint=' + str(constraint) + '\n')
	if (nodelist is not None):
		submitSL.write('#SBATCH --nodelist=' + str(nodelist) + '\n')
	if (not exclude == None):
		submitSL.write('#SBATCH --exclude=' + str(exclude) + '\n')
	if (nodes is not None):
		submitSL.write('#SBATCH --nodes=' + str(nodes) + '\n')
	submitSL.write('#SBATCH --time=' + str(time) + '     # Walltime\n')
	submitSL.write('#SBATCH --output=slurm-%j.out      # %x and %j are replaced by job name and ID'+'\n')
	submitSL.write('#SBATCH --error=slurm-%j.err'+'\n')
	if not email == '':
		submitSL.write('#SBATCH --mail-user=' + str(email) + '\n')
		submitSL.write('#SBATCH --mail-type=ALL\n')
	submitSL.write('\n')

def load_gaussian_programs(submitSL, gaussian_version=None, python_version=None):
	"""
	This method will allow you to load Gaussian and python in slurm.

	Parameters
	----------
	submitSL : open()
		This is the slurm submit file to add commands to.
	gaussian_version : str.
		This is the version of Gaussian you want to use. Default: None. 
	python_version : str.
		This is the version of Python you want to use. Default: None. 
	"""

	# Before running, do not write anything if you are loading Gaussian or Python
	if not gaussian_version and not python_version:
		return

	# First, check what programs you are going to load.
	running_gaussian = not gaussian_version is None
	running_python   = not python_version is None

	# Second, add bar.
	if running_gaussian or running_python:
		submitSL.write('# ----------------------------\n')
	
	# Third, add lines to indicate what programs are being loaded.
	if running_gaussian:
		submitSL.write('# Load Gaussian\n')
	elif running_python:
		submitSL.write('# Load Python\n')
	elif running_gaussian and running_python:
		submitSL.write('# Load Gaussian and Python\n')
	submitSL.write('\n')

	# Fourth, write lines in slurm file to load progams
	if running_gaussian:
		submitSL.write('module load '+str(gaussian_version)+'\n')
	if running_python:
		submitSL.write('module load '+str(python_version)+'\n')
	submitSL.write('\n')

def load_orca_programs(submitSL, orca_version=None, gcc_version=None, openmpi_version=None, python_version=None):
	"""
	This method will allow you to load Gaussian and python in slurm.

	Parameters
	----------
	submitSL : open()
		This is the slurm submit file to add commands to.
	orca_version : str.
		This is the version of ORCA you want to use. Default: None. 
	gcc_version : str.
		This is the version of GCC you want to use. Default: None. 
	openmpi_version : str.
		This is the version of OpenMPI you want to use. Default: None. 
	python_version : str.
		This is the version of Python you want to use. Default: None. 
	"""

	# Before running, do not write anything if you are loading Gaussian or Python
	if not orca_version and not python_version:
		return

	# First, check what programs you are going to load.
	running_orca    = not orca_version is None
	running_python  = not python_version is None

	# Second, add bar.
	if running_orca or running_python:
		submitSL.write('# ----------------------------\n')
	
	# Third, add lines to indicate what programs are being loaded.
	if running_orca:
		submitSL.write('# Load ORCA\n')
	elif running_python:
		submitSL.write('# Load Python\n')
	elif running_gaussian and running_python:
		submitSL.write('# Load ORCA and Python\n')
	submitSL.write('\n')

	# Fourth, write lines in slurm file to load progams
	if running_orca:
		submitSL.write('module load '+str(gcc_version)+'\n')
		submitSL.write('module load '+str(openmpi_version)+'\n')
		submitSL.write('module load '+str(orca_version)+'\n')
	if running_python:
		submitSL.write('module load '+str(python_version)+'\n')
	submitSL.write('\n')

def make_gaussian_temp_folder(submitSL, temp_folder_path):
	"""
	This method will write commands for create the temp folder to write Gaussian Temp files to.

	Parameters
	----------
	submitSL : open()
		This is the slurm submit file to add commands to.
	temp_folder_path : str.
		This is the path to the folder to store temp files in.
	"""
	if temp_folder_path is not None:
		submitSL.write('# ----------------------------\n')
		submitSL.write('# Create temp folder for Gaussian files.\n')
		submitSL.write('\n')
		submitSL.write('mkdir -p '+str(temp_folder_path)+'\n')
		submitSL.write('\n')

def remove_gaussian_temp_files(submitSL, gaussian_parameters, temp_folder_path, remove_chk_file=False, remove_temp_folder=False, remove_chk_filepath=None, prefix='',spaces_in_code=True):
	"""
	This method will write removal commands for removing Gaussian Temp files.

	Parameters
	----------
	submitSL : open()
		This is the slurm submit file to add commands to.
	gaussian_parameters : dict.
		This is the dictionary that provides the information for the Gaussian gjf file.
	temp_folder_path : str.
		This is the path to the folder to store temp files in.
	remove_chk_file : bool.
		This boolean indicates if you want to add commands for removing the checkpoint file. Default: False.
	remove_temp_folder : bool.
		This boolean indicates if you want to add commands for removing the temp folder. Default: False.
	prefix : str.
		if you want to add a prefix to each line of your code, do it here. Default: ''.
	spaces_in_code : bool.
		This indicates if you want to include spaces between lines in your slurm file. Default: True.
	"""

	# First, write what you are doing to the slurm file.
	submitSL.write(prefix+'# Remove temp files\n')
	if spaces_in_code:
		submitSL.write(prefix+'\n')

	# Second, remove the chk file if you would like to.
	if remove_chk_file:
		submitSL.write(prefix+'rm -fvr '+str(gaussian_parameters['chk'])+'\n')

	# Third, remove the other temp gaussian files.
	for suffix in ['rwf','int','d2e','skr']:
		submitSL.write(prefix+'rm -fvr '+str(gaussian_parameters[suffix]+'\n'))
	if spaces_in_code:
		submitSL.write(prefix+'\n')

	# Fourth, remove the old chkpoint file if desired
	if remove_chk_filepath is not None:
		submitSL.write(prefix+'# Remove the old checkpoint file.\n')
		submitSL.write(prefix+'rm -fvr '+str(remove_chk_filepath)+'\n')
		if spaces_in_code:
			submitSL.write(prefix+'\n')

	# Fifth, remove the temp folder if desired.
	if remove_temp_folder and (temp_folder_path is not None):
		submitSL.write(prefix+'# ----------------------------\n')
		submitSL.write(prefix+'# Remove the temp folder.\n')
		if spaces_in_code:
			submitSL.write(prefix+'\n')
		submitSL.write(prefix+'rm -fvr '+str(temp_folder_path)+'\n')
		if spaces_in_code:
			submitSL.write(prefix+'\n')

# ------------------------------------------------------------------------------------------------------------------------------

def make_orca_temp_folder(submitSL, temp_folder_path):
	"""
	This method will write commands for create the temp folder to write ORCA Temp files to.

	Parameters
	----------
	submitSL : open()
		This is the slurm submit file to add commands to.
	temp_folder_path : str.
		This is the path to the folder to store temp files in.
	"""
	if temp_folder_path is not None:
		submitSL.write('# ----------------------------\n')
		submitSL.write('# Create temp folder for ORCA files.\n')
		submitSL.write('\n')
		submitSL.write('mkdir -p '+str(temp_folder_path)+'\n')
		submitSL.write('\n')

def remove_orca_temp_files(submitSL, gaussian_parameters, temp_folder_path, remove_temp_folder=False, remove_gbw_file=False):
	"""
	"""
	pass

# ------------------------------------------------------------------------------------------------------------------------------


