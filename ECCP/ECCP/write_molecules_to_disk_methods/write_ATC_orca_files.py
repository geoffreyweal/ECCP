"""
write_ATC_orca_files.py, Geoffrey Weal, 8/5/22

This script is designed to write the ORCA files and submit.sl files required for performing ORCA jobs for performing atomic transition charge (ATC) calculations.
"""
from copy                                                                      import deepcopy
from ase.io                                                                    import write
from SUMELF                                                                    import make_folder
from SUMELF                                                                    import check_molecule_against_file
from SUMELF                                                                    import obtain_graph
from ECCP.ECCP.remove_unwanted_entries                                         import remove_unwanted_entries
from ECCP.ECCP.write_molecules_to_disk_methods.write_methods.orca_modified_ATC import write_orca_in_ATC
from ECCP.ECCP.write_molecules_to_disk_methods.shared_methods                  import change_folder_name_components, input_commands_for_multiwfn
from ECCP.ECCP.write_molecules_to_disk_methods.shared_methods                  import slurmSL_header, load_orca_programs, make_orca_temp_folder, remove_orca_temp_files

def write_ATC_orca_files(molecule, molecule_name, environment_about_molecule, SolventsList, orca_jobs_path, calc_parameters_for_ATCs, submission_information_for_ATCs):
	"""
	This method will write information the ORCA files to disk.

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
	orca_jobs_path : str.
		This is the path to save ORCA jobs to.
	calc_parameters_for_ATCs : list
		This dictionary contain all the information required for the ORCA input ATC file.
	submission_information_for_ATCs : list
		This dictionary contain all the information required for the submit.sl script. 
	"""

	#raise Exception('Check this is working properly')
	#raise Exception('Make sure that info from molecule graph is added to molecule before writing the orca file.')
	
	# First, check if there already exists this molecule on file, check if the molecules are the same.
	#check_molecule_against_file(molecule, calc_folder+'/'+molecule_name+'.inp')

	# First, make a copy of the orca_parameters and submission_information dictionaries.
	orca_parameters        = deepcopy(calc_parameters_for_ATCs)
	submission_information = deepcopy(submission_information_for_ATCs)
	del orca_parameters['calc_software']

	# Second, determine if some critical tags that are needed are in the submission_information dictionary. 
	got_cpu  = 'ntasks' in submission_information
	got_mem  = 'mem' in submission_information
	got_time = 'time' in submission_information
	if not (got_cpu and got_time):
		print('Error: You need to specify the following in your submission_information dictionary:')
		if not got_cpu:
			print('\t* ntasks')
		if not got_mem:
			print('\t* mem')
		if not got_time:
			print('\t* time')
		print('See https://github.com/geoffreyweal/ECCP/ for more information about these tags.')
		print('submission_information = '+str(submission_information))
		exit('This program will finish without completing.')

	# Third, copy some tag information that is in the submission_information dictionary to the orca_parameters dictionary.
	orca_parameters['NPROCS'] = submission_information['ntasks']
	
	# Fourth, give the name of the folder to place orca files to.
	quantum_chemistry_program = 'ORCA'
	functional                = change_folder_name_components(calc_parameters_for_ATCs['method'])
	basis_set                 = change_folder_name_components(calc_parameters_for_ATCs['basis'])
	funct_and_basis_name      = 'F_'+functional+'_B_'+basis_set
	calc_folder               = str(orca_jobs_path)+'/'+str(molecule_name)+'/'+quantum_chemistry_program+'_'+str(funct_and_basis_name)

	# Fifth, write the folder to place ORCA files to.
	make_folder(calc_folder)

	# Sixth, create the orca input file. 
	with open(calc_folder+'/'+molecule_name+'.inp','w') as fd:
		write_orca_in_ATC(fd, molecule, environment_about_molecule, molecule_name, **orca_parameters)

	# Seventh, create the xyz file for this ATC structure. 
	#magnetic_moments = molecule.get_initial_magnetic_moments()
	molecule_copy, molecule_copy_graph = obtain_graph(molecule.copy(), name='molecule_copy', no_of_cpus=1)
	molecule_copy.set_initial_charges(None); molecule_copy.set_cell(None)
	#molecule_copy.set_initial_magnetic_moments(magnetic_moments)
	write(calc_folder+'/'+molecule_name+'.xyz', molecule_copy)

	# Eigth, create the submit.sl file to submit this orca job to slurm.
	make_ATC_orca_submitSL(molecule_name+'.inp',calc_folder,functional,basis_set,orca_parameters,**submission_information)

# ------------------------------------------------------------------------------------------------------------------------------

def make_ATC_orca_submitSL(orca_filename,local_path,functional,basis_set,orca_parameters,ntasks,mem,time,partition='parallel',constraint=None,nodelist=None,email='',python_version='python/3.8.1',orca_version='ORCA/5.0.3',gcc_version='GCC/11.2.0',openmpi_version='OpenMPI/4.1.1',out_filename='output.out',gbw_filename='output.gbw',temp_folder_path=None,remove_chk_file=True):
	"""
	This method will write the submit.sl file in parallel

	Parameters
	----------
	orca_filename : str. 
		This is the name of the orca file.
	local_path : str. 
		This is the location to save this submit.sl file to
	functional : str. 
		This is the functional you are going to use in your ORCA calculation.
	basis_set : str. 
		This is the basis set you are going to use in your ORCA calculation.
	ntasks : int
		This is the number of cpus you want to use for ORCA jobs.
	mem : str.
		This is the amount of memory you want to use for ORCA jobs.
	time : str.
		This is the amount of time you want to use for ORCA jobs.
	partition : str.
		This is the partition to run this job on. Default: 'parallel'
	constraint : str.
		This is the slurm constraint. If you dont give this, this wont be set. Default: None
	nodelist : str.
		This is the slurm nodelist. If you dont give this, this wont be set. Default: None
	email : str.
		This is the email to email about how this job is going. If you dont give this, this wont be set. Default: ''
	python_version : str.
		This is the version of python you want to load/use in slurm. Default: 'python/3.8.1'
	orca_version : str.
		This is the version of ORCA you want to load/use in slurm. Default: 'ORCA/5.0.3'
	gcc_version : str.
		This is the version of GCC you want to load/use in slurm. Default: 'GCC/11.2.0'
	openmpi_version : str.
		This is the version of OpenMPI you want to load/use in slurm. Default: 'OpenMPI/4.1.1'
	out_filename : str.
		This is the name of the ORCA output file. Default: 'output.out'
	gbw_filename : str.
		This is the name of the ORCA Wavefunction file. Default: 'output.gbw'
	temp_folder_path : str. or None
		This is the path to the scratch directory to save ORCA temp files to. If you dont give this, ORCA temp files will be saves to the default scratch directory. Default: None
	remove_chk_file : bool.
		This variable indicates if you want to remove the chk file afterwards. Default: False
	"""
	# create name for job
	name = '-'.join(local_path.split('/')[-3:])+'-ATC'
	# writing the submit.sl script
	with open(local_path+'/'+"submit.sl", "w") as submitSL:
		slurmSL_header(submitSL, name, mem, partition, constraint, nodelist, time, email, ntasks=ntasks)
		make_orca_temp_folder(submitSL, temp_folder_path)
		load_orca_programs(submitSL, orca_version, gcc_version, openmpi_version, python_version)
		submitSL.write('# ----------------------------\n')
		submitSL.write('# ORCA under MPI requires that it be called via its full absolute path\n')
		submitSL.write('\n')
		submitSL.write('orca_exe=$(which orca)\n')
		submitSL.write('\n')
		submitSL.write('# ----------------------------\n')
		submitSL.write('# Perform ORCA Calculation\n')
		submitSL.write('\n')
		submitSL.write('${orca_exe} '+str(orca_filename)+' > '+str(out_filename)+'\n')
		submitSL.write('\n')
		submitSL.write('# ----------------------------\n')
		submitSL.write('# Perform orca_2mkl to export natural orbitals corresponding to transition density matrix (TDM) be converting the information in the gbw file into a molden.input file.\n')
		submitSL.write('\n')
		submitSL.write(f'orca_2mkl {gbw_filename} -molden\n')
		submitSL.write('\n')
		submitSL.write('# ----------------------------\n')
		remove_orca_temp_files(submitSL, orca_parameters, temp_folder_path, remove_temp_folder=True, remove_gbw_file=True)
		submitSL.write('# ----------------------------\n')
		submitSL.write('# Performing Multiwfn Calculation\n')
		submitSL.write('\n')
		submitSL.write('if [[ -f "'+str(gbw_filename)+'" ]]\n')
		submitSL.write('then\n')
		submitSL.write('\techo "found '+str(gbw_filename)+', will perform Multiwfn ATC calculation"\n')
		submitSL.write('\tsubmit_slurm_job.py multiwfn_submit.sl\n')
		submitSL.write('else\n')
		submitSL.write('\techo "'+str(gbw_filename)+' file not found. Will not perform Multiwfn ATC calculation"\n')
		submitSL.write('\techo "'+str(gbw_filename)+' file not found. Will not perform Multiwfn ATC calculation" 1>&2\n')
		submitSL.write('fi\n')
		submitSL.write('\n')
		submitSL.write('# ----------------------------\n')
		submitSL.write('echo "End of job"\n')
		submitSL.write('# ----------------------------\n')

# ------------------------------------------------------------------------------------------------------------------------------


