"""
write_ECCP_prepare_submit_scripts.py, Geoffrey Weal, 8/5/22

This script is designed to write the submit.sl files required for running Multiwfn calculations.
"""
from copy                                                                          import deepcopy
from SUMELF                                                                        import make_folder
from SUMELF                                                                        import add_graph_to_ASE_Atoms_object
from ECCP.ECCP.write_molecules_to_disk_methods.write_methods.gaussian_modified_ATC import write_gaussian_in as write_gaussian_in_ATC
from ECCP.ECCP.write_molecules_to_disk_methods.shared_methods                      import change_folder_name_components, input_commands_for_multiwfn
from ECCP.ECCP.write_molecules_to_disk_methods.shared_methods                      import slurmSL_header

def write_ATC_multiwfn_files(molecule, molecule_name, environment_about_molecule, SolventsList, gaussian_jobs_path, all_calc_parameters, all_submission_information_for_multiwfn, submit_ATCs_in_parallel=True):
	"""
	This method will write information the Gaussian files to disk.

	Parameters
	----------
	molecule : ase.Atoms.
		This is the molecule. 
	molecule_name : str.
		This is the name of the molecule. 
	environment_about_molecule : ase.Atoms or None
		This is the environment surrounding the molecule. If None, no environment was given.
	SolventsList : list of int
		These are the indices of the molecules in the molecules list that have been identified at solvents.
	gaussian_jobs_path : str.
		This is the path to save gaussian jobs to.
	all_submission_information_for_multiwfn : list
		This list contain all the information required to be included in the multiwfn_submit.sl script. 
	submit_ATCs_in_parallel : bool.
		If True, write submit.sl files for submitting individual jobs in slurm ('parallel'). If False, write a single submit.sl file for submitting this job to run ATC jobs one after another. This should be set to false if you know your ATC jobs will take less than 5 minutes to run each job. Default: True. 
	"""

	# This for loop will create all the various gaussian settings for the same molecule. 
	# This is important for testing a functionals and basis sets.
	for calc_parameters, submission_information in zip(all_calc_parameters, all_submission_information_for_multiwfn):

		# First, determine if some critical tags that are needed are in the submission_information dictionary. 
		if 'ntasks' in submission_information:
			submission_information['cpus_per_task'] = submission_information['ntasks']
			del submission_information['ntasks']
		got_cpu = 'cpus_per_task' in submission_information
		got_mem = 'mem' in submission_information
		got_time = 'time' in submission_information
		if not (got_cpu and got_mem and got_time):
			print('Error: You need to specify the following in your submission_information dictionary:')
			if not got_cpu:
				print('\t* cpus_per_task')
			if not got_mem:
				print('\t* mem')
			if not got_time:
				print('\t* time')
			print('See https://github.com/geoffreyweal/ECCP/ for more information about these tags.')
			print('submission_information = '+str(submission_information))
			exit('This program will finish without completing.')

		# Second, make sure that calc_parameters['calc_software'] is either 'gaussian' or 'orca'
		if   calc_parameters['calc_software'].lower() == 'gaussian':
			calc_software = 'gaussian'
		elif calc_parameters['calc_software'].lower() == 'orca':
			calc_software = 'orca'
		else:
			to_string  = 'Error: calc_software in the "calc_parameters" dictionary needs to be either "gaussian" or "orca"\n'
			to_string += "calc_parameters['calc_software'] = "+str(calc_parameters['calc_software'])+'\n'
			to_string += 'Check this'
			raise Exception(to_string)

		# Third, give the name of the folder to place gaussian files to.
		quantum_chemistry_program = change_folder_name_components(calc_parameters['calc_software']).upper()
		functional                = change_folder_name_components(calc_parameters['method'])
		basis_set                 = change_folder_name_components(calc_parameters['basis'])
		funct_and_basis_name      = 'F_'+functional+'_B_'+basis_set
		calc_folder               = str(gaussian_jobs_path)+'/'+str(molecule_name)+'/'+quantum_chemistry_program+'_'+str(funct_and_basis_name)

		# Fourth, write the folder to place gaussian files to.
		make_folder(calc_folder)

		# Fifth, create the submit.sl file to submit this gaussian job to slurm.
		if submit_ATCs_in_parallel:
			if calc_software == 'gaussian':
				wavefunction_filename = 'output.wfn'
			elif calc_software == 'orca':
				wavefunction_filename = 'orca.molden.input'	
			else:
				to_string  = 'Error: calc_software in the "calc_parameters" dictionary needs to be either "gaussian" or "orca"\n'
				to_string += "calc_parameters['calc_software'] = "+str(calc_parameters['calc_software'])+'\n'
				to_string += 'Check this'
				raise Exception(to_string)
			make_ATC_multiwfn_submitSL(molecule_name+'.gjf',molecule.get_chemical_symbols(),calc_folder,functional,basis_set,wavefunction_filename=wavefunction_filename,**submission_information)

	# Sixth, if desired, write the submit.sl file for all the functionals and basis sets in series.
	if not submit_ATCs_in_parallel:
		calc_folder = str(gaussian_jobs_path)+'/'+str(molecule_name)
		make_ATC_multiwfn_submitSL_in_series(molecule_name+'.gjf',calc_folder,all_calc_parameters,**submission_information)

# ------------------------------------------------------------------------------------------------------------------------------

def make_ATC_multiwfn_submitSL(gaussian_filename,molecule_symbols,local_path,functional,basis_set,cpus_per_task,mem,time,partition='parallel',constraint=None,nodelist=None,exclude=None,nodes=1,email='',log_filename='output.log',wavefunction_filename='output.wfn'):
	"""
	This method will write the submit.sl file in parallel

	Parameters
	----------
	gaussian_filename : str. 
		This is the name of the gaussian file.
	molecule_symbols : list of str.
		This list contains all the elements in the molecule.
	local_path : str. 
		This is the location to save this submit.sl file to.
	functional : str. 
		This is the functional you are going to use in your Gaussian calculation.
	basis_set : str. 
		This is the basis set you are going to use in your Gaussian calculation.
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
	log_filename : str.
		This is the name of the log file obtained by Gaussian. Default: 'output.log'
	wavefunction_filename : str.
		This is the name of the Gaussian Wavefunction file. Default: 'output.wfn'
	"""
	# create name for job
	name = '-'.join(local_path.split('/')[-3:])+'-ATC-Multiwfn'
	# writing the submit.sl script
	with open(local_path+'/'+"multiwfn_submit.sl", "w") as submitSL:
		slurmSL_header(submitSL, name, mem, partition, constraint, nodelist, time, email, cpus_per_task=cpus_per_task, nodes=nodes, exclude=exclude)
		submitSL.write('# ============================\n')
		submitSL.write('# Perform Multiwfn on '+str(wavefunction_filename)+'\n')
		submitSL.write('\n')
		submitSL.write(input_commands_for_multiwfn(wfn_filename=wavefunction_filename,molecule_symbols=molecule_symbols))
		submitSL.write('\n')
		submitSL.write('# Remove '+str(wavefunction_filename)+'\n')
		submitSL.write('rm -fvr '+str(wavefunction_filename)+'\n')
		submitSL.write('\n')
		submitSL.write('echo "Done performing Multiwfn"\n')

# ------------------------------------------------------------------------------------------------------------------------------



