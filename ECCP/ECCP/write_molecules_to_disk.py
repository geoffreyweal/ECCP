"""
make_molecule.py, Geoffrey Weal, 17/2/22

This script will save the individual molecules that are found within a crystal structure to disk.
"""
import numpy as np
from ase import Atoms
from ase.io import write
from copy import deepcopy

from SUMELF                                                               import make_folder, add_graph_to_ASE_Atoms_object
from SUMELF                                                               import check_molecule_against_file

from ECCP.ECCP.invariance_methods.are_environments_equivalent             import get_environment

from ECCP.ECCP.write_molecules_to_disk_methods.write_ATC_gaussian_files   import write_ATC_gaussian_files
from ECCP.ECCP.write_molecules_to_disk_methods.write_ATC_orca_files       import write_ATC_orca_files
from ECCP.ECCP.write_molecules_to_disk_methods.write_ATC_multiwfn_files   import write_ATC_multiwfn_files

from ECCP.ECCP.write_molecules_to_disk_methods.write_RE_gaussian_files    import write_RE_gaussian_files
from ECCP.ECCP.write_molecules_to_disk_methods.write_RE_orca_files        import write_RE_orca_files
from ECCP.ECCP.write_molecules_to_disk_methods.write_RE_gaussian_SP_files import write_RE_gaussian_SP_files
from ECCP.ECCP.write_molecules_to_disk_methods.write_RE_orca_SP_files     import write_RE_orca_SP_files

from ECCP.ECCP.write_molecules_to_disk_methods.write_FC_gaussian_files    import write_FC_gaussian_files
from ECCP.ECCP.write_molecules_to_disk_methods.write_FC_orca_files        import write_FC_orca_files

def write_molecules_to_disk(structurally_unique_molecules_indices, conformationally_unique_molecules_indices, molecules, molecule_graphs, neighbouring_molecules_about_molecules, SolventsList, molecules_folderpath, molecules_with_environment_folderpath, atc_calc_jobs_path, re_calc_jobs_path, fc_calc_jobs_path, all_calc_parameters_for_ATCs, all_calc_parameters_for_multiwfn, all_calc_parameters_for_REs, all_calc_parameters_for_FCs, all_submission_information_for_ATCs, all_submission_information_for_multiwfn, all_submission_information_for_REs, all_submission_information_for_FCs, get_molecule_atcs=True, get_molecule_res=True, get_molecule_fcs=True, run_excited_state_from_optimised_ground_structure=False):
	"""
	This method will save the individual molecules that are found within a crystal structure to disk.

	This method will also save the ATC files of these molecules if desired. 

	Parameters
	----------
	structurally_unique_molecules_indices : list of ints
		This is the indices of structurally unique molecules to save files ATC files of.
	conformationally_unique_molecules_indices : list of ints
		This is the indices of structurally unique molecules to save files RE and FC files of.
	molecules : list of ase.Atoms
		This is a list of all the molecules in the crystal
	molecule_graphs : list of networkx.graph
		These are all the associated graphs that describe the bonding network for each molecule in the crystal.
	neighbouring_molecules_about_molecules : list
		This is a list of all the details about the molecules that neighbour each of the molecules in the crystal. 
	SolventsList : list of int
		These are the indices of the molecules in the molecules list that have been identified at solvents.

	molecules_folderpath : str.  
		This is the folder to save xyz files of the molecule to.
	molecules_with_environment_folderpath : str.  
		This is the folder to save xyz files of the molecule with its environment surround it.
	atc_calc_jobs_path : str.  
		This is the folder to save calc jobs of the molecules to in order to perform ATC calculations. 
	re_calc_jobs_path : str.  
		This is the folder to save calc jobs of the molecules to in order to perform RE calculations. 
	fc_calc_jobs_path : str.  
		This is the folder to save calc jobs of the molecules to in order to perform FC (and HR) calculations. 

	all_calc_parameters_for_ATCs : list of dict.
		These are dictionaries that contain all the information to create the calc input files for ATC calculations.
	all_calc_parameters_for_multiwfn : list of dict.
		These are dictionaries that contain all the information to create the multiwfn input files for ATC calculations.
	all_calc_parameters_for_REs : list of dict.
		These are dictionaries that contain all the information to create the calc input files for performing reorganisation energy calculations.
	all_calc_parameters_for_FCs : list of dict.
		These are dictionaries that contain all the information to create the calc input files for performing franck-condon (and huang-rhys) factor calculations.

	all_submission_information_for_ATCs : list of dict.
		These are dictionaries that contain all the information to create submit.sl file to submit calc jobs to slurm.
	all_submission_information_for_multiwfn : list of dict.
		These are dictionaries that contain all the information to create multiwfn_submit.sl file to submit multiwfn jobs to slurm.
	all_submission_information_for_REs : list of dict.
		These are dictionaries that contain all the information to create submit.sl file to submit calc jobs to slurm for obtaining reorganisation energies.
	all_submission_information_for_FCs : list of dict.
		These are dictionaries that contain all the information to create submit.sl file to submit calc jobs to slurm for obtaining franck-condon (and huang-rhys) factor calculations.

	get_molecule_atcs : bool.
		This tag indicates if you want to perform ATC calculations on your molecules. 
	get_molecule_res  : bool.
		This tag indicates if you want to perform RE  calculations on your molecules. 
	get_molecule_fcs  : bool.
		This tag indicates if you want to perform FC  calculations on your molecules. 
	
	run_excited_state_from_optimised_ground_structure : bool.
		This boolean indicates if you want to run the excited state calculation from the optimised ground state calculation for reorganisation energy calculations. True if you do, False if you want to run the excited state calculation from the original structure (Default: False).
	"""

	# First, make the folder to place molecule files in.
	make_folder(molecules_folderpath)
	if len(neighbouring_molecules_about_molecules) > 0:
		make_folder(molecules_with_environment_folderpath)

	# Second, for all the structurally unique molecules in the crystal, performing the following:
	for mol_name in structurally_unique_molecules_indices:

		# ==================================================================================================
		# Third, obtain the following variables
		molecule       = molecules[mol_name]
		molecule_graph = molecule_graphs[mol_name]

		# Fourth, Add details from the molecules graph back to the molecule.
		molecule_copy = molecule.copy()
		add_graph_to_ASE_Atoms_object(molecule_copy, deepcopy(molecule_graph))

		# Fifth, obtain the name of the molecule filename. 
		molecule_name = 'molecule_'+str(mol_name)
		if mol_name in SolventsList:
			molecule_name += 'S'

		# Sixth, if there already exists this molecule on file, check if the molecules are the same.
		check_molecule_against_file(molecule_copy, molecules_folderpath+'/'+molecule_name+'.xyz')

		# Seventh, write the molecules to molecules_folderpath as xyz files.
		write(molecules_folderpath+'/'+molecule_name+'.xyz', molecule_copy)

		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  

		# Eighth, write the molecules with their environments to file
		if len(neighbouring_molecules_about_molecules) > 0:

			# 8.1: Write the molecule file including it's environment.
			environment_about_molecule = get_environment(mol_name, neighbouring_molecules_about_molecules, molecules)

			# 8.2: Include no information about the neighbourlist for environmental molecules.
			environment_about_molecule.set_array('NeighboursList', np.array(['-']*len(environment_about_molecule)))

			# 8.3: Add the molecules that are apart of the environment to the molecule of interest. 
			molecule_with_environment = molecule.copy() + environment_about_molecule

			# 8.4: Indicate which atoms are apart of the molecule and which atoms are apart of the environment. 
			#molecule_with_environment.set_tags([1]*len(molecule) + [2]*len(environment_about_molecule))

			raise Exception('Need to check how this works, and include method for checking Atoms object against that on file if it exists.')

			# 8.5: If there already exists this molecule on file, check if the molecules are the same.
			check_molecule_against_file(molecule_copy, molecules_folderpath+'/'+molecule_name+'.xyz')

			# 8.6: Write the molecule with its environment to disk. 
			write(molecules_with_environment_folderpath+'/'+molecule_name+'.xyz', molecule_with_environment)

		else:

			# 8.7: Do not include the environment around the molecule. 
			environment_about_molecule = None

		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

		# Ninth, write the DFT calculation files related to the ATC if desired.
		if get_molecule_atcs and (all_calc_parameters_for_ATCs is not None):

			# 9.1: Write the ATC input files.
			for calc_parameters_for_ATCs, submission_information_for_ATCs in zip(all_calc_parameters_for_ATCs, all_submission_information_for_ATCs):

				# 9.1.1: Make sure that the 'calc_software' tag has been added, as this tells ECCP what DFT program you want to use. 
				if not 'calc_software' in calc_parameters_for_ATCs:
					raise Exception("Error: You need to specify a value for 'calc_software' in calc_parameters_for_ATCs.")

				# 9.1.2: Create the files for running the program in Gaussian or ORCA. 
				if   calc_parameters_for_ATCs['calc_software'].lower() == 'gaussian':
					write_ATC_gaussian_files(molecule_copy, molecule_name, environment_about_molecule, SolventsList, atc_calc_jobs_path, calc_parameters_for_ATCs, submission_information_for_ATCs)
				elif calc_parameters_for_ATCs['calc_software'].lower() == 'orca':
					write_ATC_orca_files    (molecule_copy, molecule_name, environment_about_molecule, SolventsList, atc_calc_jobs_path, calc_parameters_for_ATCs, submission_information_for_ATCs)
				else:
					raise Exception("Error: calc_parameters_for_ATCs['calc_software'] needs to be either Gaussian or ORCA. calc_parameters_for_ATCs['calc_software'] = "+str(calc_parameters_for_ATCs['calc_software']))

			# 9.2: Write the multiwfn input files for the molecules.
			for index in range(len(all_submission_information_for_ATCs)):
				if 'log_filename' in all_submission_information_for_ATCs[index]:
					all_submission_information_for_multiwfn[index]['log_filename'] = all_submission_information_for_ATCs[index]['log_filename']
			write_ATC_multiwfn_files(molecule_copy, molecule_name, environment_about_molecule, SolventsList, atc_calc_jobs_path, all_calc_parameters_for_multiwfn, all_submission_information_for_multiwfn)

	# Tenth, For all the conformationally unique molecules in the crystal, performing the following:
	#        Note: If a molecule is structurally unique, it will also be conformationally unique.
	#              A conformationally unique molecule might not be necessarily structurally unique.
	for mol_name in conformationally_unique_molecules_indices:

		# Eleventh, obtain the following variables
		molecule       = molecules[mol_name]
		molecule_graph = molecule_graphs[mol_name]

		# Twelfth, Add details from the molecules graph back to the molecule.
		molecule_copy = molecule.copy()
		add_graph_to_ASE_Atoms_object(molecule_copy, deepcopy(molecule_graph))

		# Thirteenth, write the molecules to molecules_folderpath as xyz files.
		molecule_name = 'molecule_'+str(mol_name)
		if mol_name in SolventsList:
			molecule_name += 'S'

		# Fourteenth, write the reorganisation energy (RE) input file.
		if get_molecule_res and (all_calc_parameters_for_REs is not None):

			# 14.1: Write the RE input files.
			for calc_parameters_for_REs, submission_information_for_REs in zip(all_calc_parameters_for_REs, all_submission_information_for_REs):

				# 14.1.1: Make sure that the 'calc_software' tag has been added, as this tells ECCP what DFT program you want to use. 
				if not 'calc_software' in calc_parameters_for_REs:
					raise Exception("Error: You need to specify a value for 'calc_software' in calc_parameters_for_REs.")

				# 14.1.2: Create the files for running the program in Gaussian or ORCA. 
				if   calc_parameters_for_REs['calc_software'].lower() == 'gaussian':
					write_RE_gaussian_files   (molecule_copy, molecule_name, SolventsList, re_calc_jobs_path, calc_parameters_for_REs, submission_information_for_REs, run_excited_state_from_optimised_ground_structure)
					write_RE_gaussian_SP_files(molecule_copy, molecule_name, SolventsList, re_calc_jobs_path, calc_parameters_for_REs, submission_information_for_REs)
				elif calc_parameters_for_REs['calc_software'].lower() == 'orca':
					write_RE_orca_files       (molecule_copy, molecule_name, SolventsList, re_calc_jobs_path, calc_parameters_for_REs, submission_information_for_REs, run_excited_state_from_optimised_ground_structure)
					write_RE_orca_SP_files    (molecule_copy, molecule_name, SolventsList, re_calc_jobs_path, calc_parameters_for_REs, submission_information_for_REs)
				else:
					raise Exception("Error: calc_parameters_for_REs['calc_software'] needs to be either Gaussian or ORCA. calc_parameters_for_REs['calc_software'] = "+str(calc_parameters_for_REs['calc_software']))

		# Fifteenth, write the Franck-Condon (FC) input files.
		if get_molecule_fcs and (all_calc_parameters_for_FCs is not None):

			# 15.1: Write the FC input files.
			for calc_parameters_for_FCs, submission_information_for_FCs in zip(all_calc_parameters_for_FCs, all_submission_information_for_FCs):

				# 15.1.1: Make sure that the 'calc_software' tag has been added, as this tells ECCP what DFT program you want to use. 
				if not 'calc_software' in calc_parameters_for_FCs:
					raise Exception("Error: You need to specify a value for 'calc_software' in calc_parameters_for_FCs.")

				# 15.1.2: Create the files for running the program in Gaussian or ORCA. 
				if   calc_parameters_for_FCs['calc_software'].lower() == 'gaussian':
					write_FC_gaussian_files(molecule_copy, molecule_name, SolventsList, fc_calc_jobs_path, calc_parameters_for_FCs, submission_information_for_FCs)
				elif calc_parameters_for_FCs['calc_software'].lower() == 'orca':
					raise Exception('Write this method for ORCA')
					pass	# while figuring this out, ignore this part
					#write_FC_orca_files    (molecule_copy, molecule_name, SolventsList, fc_calc_jobs_path, calc_parameters_for_FCs, submission_information_for_FCs)
				else:
					raise Exception("Error: calc_parameters_for_FCs['calc_software'] needs to be either Gaussian or ORCA. calc_parameters_for_FCs['calc_software'] = "+str(calc_parameters_for_FCs['calc_software']))

