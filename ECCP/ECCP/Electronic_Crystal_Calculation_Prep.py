"""
Electronic_Crystal_Calculation_Prep.py, Geoffrey Weal, 17/2/22

Electronic_Crystal_Calculation_Prep is designed to provide information about the molecules in a OPV crystal structure, as well as dimer between individual molecules.
"""
#from SUMELF import molecules_inbetween_dimers

import os
from copy import deepcopy

from SUMELF                                                                                  import read_crystal
from SUMELF                                                                                  import obtain_graph

from ECCP.ECCP.introduction_message                                                          import introduction_message

from ECCP.ECCP.prepare_calc_params_and_submit_information                                    import prepare_ATC_calc_params_and_submit_information, prepare_RE_calc_params_and_submit_information, prepare_FC_calc_params_and_submit_information, prepare_dimer_calc_params_and_submit_information

from ECCP.ECCP.convert_ECCP_Information_data.read_data_from_ECCP_Information                 import read_data_from_ECCP_Information

from ECCP.ECCP.process_crystal_methods.convert_bonds_to_ignore_file_to_list                  import convert_bonds_to_ignore_file_to_list
from SUMELF                                                                                  import process_crystal, get_spacegroup_molecules
from ECCP.ECCP.remove_unwanted_entries                                                       import remove_unwanted_entries
from SUMELF                                                                                  import is_solvent
from ECCP.ECCP.remove_solvents_from_molecules_dict                                           import remove_solvents_from_molecules_dict

from ECCP.ECCP.get_neighbouring_molecules                                                    import get_neighbouring_molecules
from ECCP.ECCP.convert_ECCP_Information_data.convert_dimer_details                           import convert_dimer_details_to_neighbourhood_molecules_for_dimer_method
from ECCP.ECCP.check_dimer_duplication                                                       import check_dimer_duplication

from ECCP.ECCP.centre_molecules                                                              import centre_molecules

from ECCP.ECCP.get_neighbouring_molecules_about_system                                       import get_neighbouring_molecules_about_molecules, get_neighbouring_molecules_about_dimers

from ECCP.ECCP.write_crystal_structures_to_disk                                              import write_crystal_structures_to_disk

from ECCP.ECCP.get_simple_molecule_graphs                                                    import get_simple_molecule_graphs

from SUMELF                                                                                  import obtain_unique_molecules
from ECCP.ECCP.get_unique_molecules                                                          import get_unique_molecules
from ECCP.ECCP.convert_ECCP_Information_data.convert_unique_molecule_data                    import convert_existing_unique_molecule_data_from_ECCP_Information

from ECCP.ECCP.get_dimers                                                                    import get_dimers
from ECCP.ECCP.convert_ECCP_Information_data.convert_dimer_details_to_all_dimers_info_method import convert_dimer_details_to_all_dimers_info_method
from ECCP.ECCP.utilities.add_dimer_name_to_neighbourhood_molecules_for_dimer_method          import add_dimer_name_to_neighbourhood_molecules_for_dimer_method
from SUMELF                                                                                  import obtain_unique_dimers
from ECCP.ECCP.get_unique_dimers                                                             import get_unique_dimers
from ECCP.ECCP.convert_ECCP_Information_data.convert_unique_dimer_data                       import convert_existing_unique_dimer_data_from_ECCP_Information

from SUMELF                                                                                  import make_folder, remove_folder
from ECCP.ECCP.write_molecules_to_disk                                                       import write_molecules_to_disk
from ECCP.ECCP.write_dimers_to_disk                                                          import write_dimers_to_disk

from ECCP.ECCP.write_ECCP_process_submit_scripts                                             import write_ECCP_process_ATC_submit_script, write_ECCP_process_RE_submit_script, write_ECCP_process_FC_submit_script, write_ECCP_process_EET_submit_script, write_ECCP_process_Eigendata_submit_script, write_ECCP_process_ICT_submit_script

from ECCP.ECCP.write_results_document                                                        import write_results_document

no_of_char_in_divides = 57
divide_string = '.'+'-'*no_of_char_in_divides+'.'

def ECCP(filepath, bonds_to_ignore=None, make_molecule_method='component_assembly_approach', molecule_equivalence_method={'method': 'invariance_method', 'type': 'combination'}, make_dimer_method={'method': 'nearest_atoms_method', 'max_dimer_distance': 8.0}, dimer_equivalence_method={'method': 'invariance_method', 'type': 'combination'}, environment_settings={'include_environment_where_possible': False, 'environment_radius': 8.0}, include_hydrogens_in_neighbour_analysis=False, include_hydrogens_in_uniqueness_analysis=False, remove_solvents=False, atc_file_creation_information=None, re_file_creation_information=None, fc_file_creation_information=None, eet_file_creation_information=None, ict_file_creation_information=None, overall_folder_suffix_name='', run_excited_state_from_optimised_ground_structure=False, no_of_cpus=1):
	"""
	The Electronic Crystal Calculation Prep (ECCP) Program is designed to:

	Part 1: Gather molecules from crystal structure.
		* Load the crystal and get a graph of the crystal, where nodes are atoms and edges are bonds between atoms.
		* Separate the crystal into separate individual molecules (species).
		* If you are supplying custom molecules, make sure that the custom molecules only contain atoms that have been removed from the molecule, and then import the custom molecules.
		* Centre the connected molecule as much as possible so that most of its atoms lay inside the unit cell without loosing its general position in the unit cell.

	Part 2: Determine which molecules neighbour each other in the crystal.
		* Molecules that are within the neighbouring distance of each other will be recorded, based on the settings given by the make_dimer_method and environment_settings dictionaries. 

	Part 3: Writing input files that involve individual molecules.
		* Write all the xyz, atomic transition charge (ATC) calc input files, reorganisation energy (RE), and Franck-Condon/Huang-Rhys (FC/HR) calc input files for all molecules in the crystal.
		* Determine which molecules are structurally unique and which are structurally equivalent to other molecules in the crystal (including same bends and twists).
		* Write all the xyz, atomic transition charge (ATC) calc input files, and reorganisation energy (RE), and Franck-Condon/Huang-Rhys (FC/HR) calc input files for only the unique molecules in the crystal.
	
	Part 4: Writing input files that involve dimers between pairs of molecules in the crystal.
		* Obtain the dimers of the molecules in the crystal across all cells, where the molecules are within some distance of each other.
		* Write all the xyz, electronic energy transfer (EET), and intermolecular charge transfer (ICT) calc input files for all dimers in the crystal.
		* Determine which dimers are unique and which are equivalent.
		* Write all the xyz, electronic energy transfer (EET), and intermolecular charge transfer (ICT) calc input files for only unique dimers in the crystal.

	Part 5: Final information about molecules and dimers.
		* Report all the relavant structural information about molecules and dimers that have been obtained.
		* Report information about the unique and equivalent molecules and dimers.

	While this program is designed for OPV crystals, this program is designed to work for any crystal system. 

	Parameters
	----------
	filepath : list of ase.Atoms objects.
		These are all the individual molecules identified in the crystal that you want to determine dimers for.

	bonds_to_ignore : list or str.
		This is a list of the atom pairs (given as indices) to ignore bonding between. This is given as a list of index pairs or a string of the name of the file that contains this list.

	make_molecule_method : str.
		This is the name of the method you want to use to create the molecule. See https://github.com/geoffreyweal/ECCP for more information. Default: 'component_assembly_approach'. 
	molecule_equivalence_method : dict. 
		This contains the information about the method used to determine which molecule are equivalent. See https://github.com/geoffreyweal/ECCP for more information. Default: {'method': 'invariance_method'}. 
	
	make_dimer_method : dict.
		This dictionary contains information about the method used to obtain dimers. See https://github.com/geoffreyweal/ECCP for more information. Default: {'method': 'nearest_atoms_method', 'rCut': 5.0}. 
	dimer_equivalence_method : dict.
		This contains the information about the method used to determine which dimers are equivalent. See https://github.com/geoffreyweal/ECCP for more information. Default: {'method': 'invariance_method'}.
	
	environment_settings : dict.
		This dictionary contains information about how to treat the local environment about a dimer (where applicable). 

	include_hydrogens_in_neighbour_analysis: bool. 
		This tag indicates if you want to include hydrogens when accessing neighbours between molecules. Default: False
	include_hydrogens_in_uniqueness_analysis: bool. 
		This tag indicates if you want to include hydrogens when accessing uniquess between molecules and dimers. Default: False
	remove_solvents: bool.
		This tag indicates if the user wants to include solvents in calculations. Default: False

	atc_file_creation_information : list of dict. or None
		This variable contains the calc Parameters and Submission information required for performing Atomic Transition Charge (ATC) calculations in calc and Multiwfn. Set to None if you do not want to perform ATC calculations. Default: None
	re_file_creation_information : list of dict. or None
		This variable contains the calc Parameters and Submission information required for performing reorganisation energy (RE) calculations in calc. Set to None if you do not want to perform RE calculations. Default: None
	fc_file_creation_information : list of dict. or None
		This variable contains the calc Parameters and Submission information required for performing Franck-Condon (FC) and Huang-Rhys (HR) calculations in calc. Set to None if you do not want to perform FC (and HR) calculations. Default: None

	eet_file_creation_information : list of dict. or None
		This variable contains the calc Parameters and Submission information required for performing Electronic Energy Transfer (EET) calculations in calc. Set to None if you do not want to perform EET calculations. Default: None
	ict_file_creation_information : list of dict. or None
		This variable contains the calc Parameters and Submission information required for performing intermolecular charge transfer calculations. This includes obtaining the eigendata (such as overlap orbtials and molecular orbital energies and coefficients). Set to None if you do not want to obtain eigendata. Default: None

	overall_folder_suffix_name : str.
		This is the suffix to add to the ECCP_Data name if you want to distinguish it in any way. If you set this to something, the overall folder name will be given as 'ECCP_Data_'+str(overall_folder_suffix_name). Default: ''
	run_excited_state_from_optimised_ground_structure : bool.
		This boolean indicates if you want to run the excited state calculation from the optimised ground state calculation for reorganisation energy calculations. True if you do, False if you want to run the excited state calculation from the original structure (Default: False).
	no_of_cpus : int.
		This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.

	Attributes
	----------
	calc_parameters_for_atomic_transition_charges : dict./list of dict.
		This contains all the information required by calc to perform atomic transition charge (ATC) calculations, including the functional and basis set. 
	submission_information_for_atomic_transition_charges : dict./list of dict.
		This contains all the information required for the submit.sl file to perform atomic transition charge (ATC) calculations, required for submitting calc jobs to slurm. 
	submission_information_for_multiwfn : dict./list of dict.
		This contains all the information required for the multiwfn_submit.sl file, required for submitting multiwfn jobs to slurm to get atc chg files. 

	calc_parameters_for_reorganisation_energy : dict./list of dict.
		This contains all the information required by calc to perform reorganisation energy (RE) calculations, including the functional and basis set. 
	submission_information_for_reorganisation_energy : dict./list of dict.
		This contains all the information required for the submit.sl file to perform reorganisation energy (RE) calculations, required for submitting calc jobs to slurm. 

	calc_parameters_for_electronic_energy_transfer_calcs : dict./list of dict.
		This contains all the information required by calc to perform electronic energy transfer (EET) calculations, including the functional and basis set. 
	submission_information_for_electronic_energy_transfer_calcs : dict./list of dict.
		This contains all the information required for the submit.sl file to perform electronic energy transfer (EET) calculations, required for submitting calc jobs to slurm. 

	calc_parameters_for_eigendata_calcs : dict./list of dict.
		This contains all the information required by calc to obtain eigendata (such as overlap orbtials and molecular orbital energies and coefficients), including the functional and basis set. 
	submission_information_for_eigendata_calcs : dict./list of dict.
		This contains all the information required for the submit.sl file to obtain eigendata (such as overlap orbtials and molecular orbital energies and coefficients), required for submitting calc jobs to slurm. 
	"""

	# ----------------------------------------------------------------------------------------- #
	# ----------------------------------------------------------------------------------------- #
	# ----------------------------------------------------------------------------------------- #
	#                                 INITIATION OF ECCP PROGRAM                                #
	# ----------------------------------------------------------------------------------------- #
	# ----------------------------------------------------------------------------------------- #
	# ----------------------------------------------------------------------------------------- #

	# First, we will perform some tasks to intialise the ECCP program for this crystal before beginning.

	# 1.1: Perform introduction message for this program.
	introduction_message()

	# 1.2: Make sure all Calculation and Submission dictionaries are all good to go.
	all_calc_parameters_for_ATCs, all_calc_parameters_for_multiwfn,   all_submission_information_for_ATCs, all_submission_information_for_multiwfn, get_molecule_atcs = prepare_ATC_calc_params_and_submit_information(deepcopy(atc_file_creation_information))
	all_calc_parameters_for_REs,  all_submission_information_for_REs,                                                                               get_molecule_res  = prepare_RE_calc_params_and_submit_information (deepcopy(re_file_creation_information))
	all_calc_parameters_for_FCs,  all_submission_information_for_FCs,                                                                               get_molecule_fcs  = prepare_FC_calc_params_and_submit_information (deepcopy(fc_file_creation_information))
	if get_molecule_fcs and (not get_molecule_res):
		raise Exception('Error: You have provided settings for obtain Franck-Condon (and Huang-Rhys) factors, but not to get reorganisation energies. The calculations needed to obtain reorganisation energies are also required to obtain Franck-Condon (and Huang-Rhys) factors. Include settings for obtaining reorganisation energies and try again.')
	all_calc_parameters_for_EETs, all_submission_information_for_EETs, get_dimer_eets = prepare_dimer_calc_params_and_submit_information(deepcopy(eet_file_creation_information), job_type='EET')
	all_calc_parameters_for_ICTs, all_submission_information_for_ICTs, get_dimer_icts = prepare_dimer_calc_params_and_submit_information(deepcopy(ict_file_creation_information), job_type='ICT')

	# 1.3: Write what crystal we will be focusing on.
	filename = os.path.basename(filepath)
	print('Performing ECCP upon: '+str(filename))
	print('Number of CPUs utilised: '+str(no_of_cpus))

	# 1.4: Names of folders to save data into.
	ECCP_Data_path                    = 'ECCP_Data'+('' if (overall_folder_suffix_name == '') else ('_'+str(overall_folder_suffix_name)))
	ECCP_Information_path             = ECCP_Data_path+'/'+'ECCP_Information'
	All_ATC_Calc_Jobs_folder          = ECCP_Data_path+'/'+'All_ATC_Calc_Jobs'
	Unique_ATC_Calc_Jobs_folder       = ECCP_Data_path+'/'+'Unique_ATC_Calc_Jobs'
	All_RE_Calc_Jobs_folder           = ECCP_Data_path+'/'+'All_RE_Calc_Jobs'
	Unique_RE_Calc_Jobs_folder        = ECCP_Data_path+'/'+'Unique_RE_Calc_Jobs'
	All_RE_Checkpoint_Files_folder    = ECCP_Data_path+'/'+'All_RE_Checkpoint_Files'
	Unique_RE_Checkpoint_Files_folder = ECCP_Data_path+'/'+'Unique_RE_Checkpoint_Files'
	All_FC_Calc_Jobs_folder           = ECCP_Data_path+'/'+'All_FC_Calc_Jobs'
	Unique_FC_Calc_Jobs_folder        = ECCP_Data_path+'/'+'Unique_FC_Calc_Jobs'
	All_EET_Calc_Jobs_folder          = ECCP_Data_path+'/'+'All_EET_Calc_Jobs'
	Unique_EET_Calc_Jobs_folder       = ECCP_Data_path+'/'+'Unique_EET_Calc_Jobs'
	All_Eigendata_Calc_Jobs_folder    = ECCP_Data_path+'/'+'All_Eigendata_Calc_Jobs'
	Unique_Eigendata_Calc_Jobs_folder = ECCP_Data_path+'/'+'Unique_Eigendata_Calc_Jobs'

	# 1.5: Paths to save files to.
	subfolder_name                        = filename.replace('.','_')
	path_to_eccp_folder                   = ECCP_Information_path+'/'+subfolder_name
	all_atc_calc_jobs_path                = str(All_ATC_Calc_Jobs_folder)+'/'+str(subfolder_name)
	unique_atc_calc_jobs_path             = str(Unique_ATC_Calc_Jobs_folder)+'/'+str(subfolder_name)
	all_re_calc_jobs_path                 = str(All_RE_Calc_Jobs_folder)+'/'+str(subfolder_name)
	unique_re_calc_jobs_path              = str(Unique_RE_Calc_Jobs_folder)+'/'+str(subfolder_name)
	all_re_checkpoint_files_path          = str(All_RE_Checkpoint_Files_folder)+'/'+str(subfolder_name)
	unique_re_checkpoint_files_path       = str(Unique_RE_Checkpoint_Files_folder)+'/'+str(subfolder_name)
	all_fc_calc_jobs_path                 = str(All_FC_Calc_Jobs_folder)+'/'+str(subfolder_name)
	unique_fc_calc_jobs_path              = str(Unique_FC_Calc_Jobs_folder)+'/'+str(subfolder_name)
	all_eet_calc_jobs_path                = str(All_EET_Calc_Jobs_folder)+'/'+str(subfolder_name)
	unique_eet_calc_jobs_path             = str(Unique_EET_Calc_Jobs_folder)+'/'+str(subfolder_name)
	all_eigendata_calc_jobs_path          = str(All_Eigendata_Calc_Jobs_folder)+'/'+str(subfolder_name)
	unique_eigendata_calc_jobs_path       = str(Unique_Eigendata_Calc_Jobs_folder)+'/'+str(subfolder_name)
	
	# ----------------------------------------------------------------------------------------- #
	# Second, obtain the data from the ``ECCP_Information`` folder.

	# 2.1: Obtain the data from previous ECCP runs from the ``ECCP_Information`` folder.
	have_ECCP_Information_crystal_file, has_neighbouring_molecules, has_unique_molecules, has_unique_dimers, ECCP_Information_data = read_data_from_ECCP_Information(path_to_eccp_folder, make_dimer_method, environment_settings)

	# 2.2: Record the information about the crystal.xyz file from ECCP_Information if this exists.
	if have_ECCP_Information_crystal_file:
		crystal = ECCP_Information_data['crystal']

	# 2.3: Record information about the dimer details from the ECCP_Information folder if this exists.
	if has_neighbouring_molecules:
		dimer_details = ECCP_Information_data['dimer_details']

	# 2.4: Record information about the structural and conformational equivalent molecule groups from the ECCP_Information folder if this exists.
	if has_unique_molecules:
		structurally_equivalent_molecule_groups     = ECCP_Information_data['structurally_equivalent_molecule_groups']
		conformationally_equivalent_molecule_groups = ECCP_Information_data['conformationally_equivalent_molecule_groups']
		structurally_equivalent_molecule_pairs      = ECCP_Information_data['structurally_equivalent_molecule_pairs']

	# 2.5: Record information about the structural equivalent dimer groups from the ECCP_Information folder if this exists.
	if has_unique_dimers:
		structurally_equivalent_dimer_groups = ECCP_Information_data['structurally_equivalent_dimer_groups']
		structurally_equivalent_dimer_pairs  = ECCP_Information_data['structurally_equivalent_dimer_pairs']

	# ----------------------------------------------------------------------------------------- #

	# Third, remove associated folders for saving data and calc jobs to. 
	make_folder(path_to_eccp_folder)
	if (all_calc_parameters_for_ATCs is not None):
		for paths_to_calc_job in [all_atc_calc_jobs_path, unique_atc_calc_jobs_path]: 
			make_folder(paths_to_calc_job)
	if (all_calc_parameters_for_REs is not None):
		for paths_to_calc_job in [all_re_calc_jobs_path, unique_re_calc_jobs_path]: 	
			make_folder(paths_to_calc_job)
		make_folder(All_RE_Checkpoint_Files_folder)
		make_folder(Unique_RE_Checkpoint_Files_folder)
		#for paths_to_calc_job in [all_re_checkpoint_files_path, unique_re_checkpoint_files_path]: 	
			#make_folder(paths_to_calc_job)
	if (all_calc_parameters_for_FCs is not None):
		for paths_to_calc_job in [all_fc_calc_jobs_path, unique_fc_calc_jobs_path]: 	
			make_folder(paths_to_calc_job)
	if (all_calc_parameters_for_EETs is not None):
		for paths_to_calc_job in [all_eet_calc_jobs_path, unique_eet_calc_jobs_path]: 
			make_folder(paths_to_calc_job)
	if (all_calc_parameters_for_ICTs is not None):
		for paths_to_calc_job in [all_eigendata_calc_jobs_path, unique_eigendata_calc_jobs_path]: 
			make_folder(paths_to_calc_job)

	# Fourth, create the logger for recording messages about about programs run to the logfile.
	import logging
	Log_Format = "%(levelname)s %(asctime)s - %(message)s"
	logging.basicConfig(filename=ECCP_Data_path+'/'+"ECCP_logfile.log", filemode="w", format=Log_Format, level=logging.NOTSET)
	logger = logging.getLogger()

	# ----------------------------------------------------------------------------------------- #
	# ----------------------------------------------------------------------------------------- #
	# ----------------------------------------------------------------------------------------- #
	#                      PART 1: Gather molecules from crystal structure                      #
	# ----------------------------------------------------------------------------------------- #
	# ----------------------------------------------------------------------------------------- #
	# ----------------------------------------------------------------------------------------- #

	print(divide_string)

	# Fifth, read in the crystal file in the ASE.
	if not have_ECCP_Information_crystal_file:
		crystal = read_crystal(filepath)

	# Sixth, make sure that the periodic boundary condition of the crystal is set to true.
	crystal.set_pbc(True)

	# Seventh, remove some of the unnecessary attributes from the crystal
	crystal = remove_unwanted_entries(crystal)

	# Eighth, get the graph of the crystal.
	#crystal, crystal_graph = obtain_graph(crystal, name='crystal', no_of_cpus=no_of_cpus)
	crystal_graph = obtain_graph(crystal, name='crystal', no_of_cpus=no_of_cpus)

	# Ninth, obtain the list of bonds in the crystal to ignore if given.
	if isinstance(bonds_to_ignore,str):
		bonds_to_ignore = convert_bonds_to_ignore_file_to_list(bonds_to_ignore)

	# Tenth, we will now process the crystal and obtain the molecules in the crystal and other molecule and crystal data. 
	print('Processing Crystal: Obtaining the molecules in the crystal')
	molecules, molecule_graphs, SolventsList, symmetry_operations, unitcelllatticevectors = process_crystal(crystal, crystal_graph=crystal_graph, take_shortest_distance=True, no_of_cpus=no_of_cpus, logger=logger, bonds_to_ignore=bonds_to_ignore)

	# Eleventh, remove unwanted entries from molecules
	for molecule_name in sorted(molecules.keys()):
		molecules[molecule_name] = remove_unwanted_entries(molecules[molecule_name])

	# Twelfth, if there are no entries in SolventsList, this may mean that the solvents have not been recorded in the crystal file
	#        This method will determine which molecules are solvents in the crystal. 
	if len(SolventsList) == 0:
		solvent_components = [molecule_name for molecule_name in molecule_graphs.keys() if is_solvent(molecule_graphs[molecule_name])]

	# Thirteenth, remove solvents if you do not want to include solvents in your input files for further calculations in Gaussian/ORCA
	if remove_solvents:
		print('Removing Solvents from the Crystal')
		molecules, molecule_graphs = remove_solvents_from_molecules_dict(molecules, molecule_graphs, SolventsList)
		SolventsList = []

	# Fourteenth, get the extra molecules in the crystal that exist due to symmetry operations in the spacegroup.
	print('Obtaining molecules from the crystal that exist due to spacegroup symmetries.')
	molecules, molecule_graphs, crystal, crystal_graph, SolventsList = get_spacegroup_molecules(molecules, molecule_graphs, SolventsList, symmetry_operations, unitcelllatticevectors)

	# Fifteenth, remove unwanted entries from the crystal and it's associated molecules.
	crystal = remove_unwanted_entries(crystal)
	for molecule_name in sorted(molecules.keys()):
		molecules[molecule_name] = remove_unwanted_entries(molecules[molecule_name])

	# Sixteenth, centre the molecules in the middle of the origin unit cell
	centre_molecules(molecules)
	print('All molecules in the crystal have been identified.')

	# ----------------------------------------------------------------------------------------- #
	# ----------------------------------------------------------------------------------------- #
	# ----------------------------------------------------------------------------------------- #
	#                    PART 2: Determine neighbouring molecules in crystal                    #
	# ----------------------------------------------------------------------------------------- #
	# ----------------------------------------------------------------------------------------- #
	# ----------------------------------------------------------------------------------------- #

	print(divide_string)

	# Seventeenth, get the neighbourhood_molecules information needed for running the dimer method. 
	print('Determining neighbours between molecules in the crystal.')
	if not has_neighbouring_molecules:
		neighbourhood_molecules_for_dimer_method, neighbourhood_molecules_for_environment_method = get_neighbouring_molecules(molecules, molecule_graphs, make_dimer_method=make_dimer_method, environment_settings=environment_settings, include_hydrogens_in_neighbour_analysis=include_hydrogens_in_neighbour_analysis, no_of_cpus=no_of_cpus)
	else:
		print('Obtaining "neighbours between molecules" data from ECCP_Information folder.')
		neighbourhood_molecules_for_dimer_method = convert_dimer_details_to_neighbourhood_molecules_for_dimer_method(dimer_details, molecules, unitcelllatticevectors, make_dimer_method)
		neighbourhood_molecules_for_environment_method = [] # To do when I know how to use this 

	# Eighteenth, perform a check to make sure that a dimer hasn't been entered in twice into neighbourhood_molecules_for_dimer_method
	check_dimer_duplication(neighbourhood_molecules_for_dimer_method)

	# Nineteenth, sort the dimers in neighbourhood_molecules_for_dimer_method. 
	#             Sorting format: (distance between molecules, mol1_index, mol2_index, unit_cell_i, unit_cell_j, unit_cell_k)
	neighbourhood_molecules_for_dimer_method.sort(key=lambda x: (x[4], x[0], x[1], x[2][0], x[2][1], x[2][2]))

	# NOTE: Will probably need order the neighbourhood_molecules_for_environment_method when created.

	print(divide_string)

	# ----------------------------------------------------------------------------------------- #
	# ----------------------------------------------------------------------------------------- #
	# ----------------------------------------------------------------------------------------- #
	#               PART 3: Writing input files that involve individual molecules               #
	# ----------------------------------------------------------------------------------------- #
	# ----------------------------------------------------------------------------------------- #
	# ----------------------------------------------------------------------------------------- #

	# Twentieth, get all the molecules that neighbour each molecule in the unit cell to be recorded into Gaussian/ORCA/Quantum Chemistry program.
	neighbouring_molecules_about_molecules = get_neighbouring_molecules_about_molecules(environment_settings, neighbourhood_molecules_for_environment_method)

	# Twenty-first, remove unwanted entries from molecules
	for mol_name in molecules.keys():
		molecules[mol_name] = remove_unwanted_entries(molecules[mol_name])

	# Twenty-second, write all molecules to disk. This will change the file in the All_Molecules folder to custom molecules if you are using molecules from the molecules
	print('Found '+str(len(molecules))+' molecules in the unit cell.')
	print('Writing molecules to '+str(path_to_eccp_folder))
	all_molecules_folderpath                  = path_to_eccp_folder+'/'+'All_Molecules'
	all_molecules_with_environment_folderpath = path_to_eccp_folder+'/'+'All_Molecules_with_environment'
	all_molecules_names                     = list(molecules.keys())
	write_molecules_to_disk(all_molecules_names, all_molecules_names, molecules, molecule_graphs, neighbouring_molecules_about_molecules, SolventsList, all_molecules_folderpath, all_molecules_with_environment_folderpath, all_atc_calc_jobs_path, all_re_calc_jobs_path, all_fc_calc_jobs_path, all_calc_parameters_for_ATCs=all_calc_parameters_for_ATCs, all_calc_parameters_for_multiwfn=all_calc_parameters_for_multiwfn, all_calc_parameters_for_REs=all_calc_parameters_for_REs, all_calc_parameters_for_FCs=all_calc_parameters_for_FCs, all_submission_information_for_ATCs=all_submission_information_for_ATCs, all_submission_information_for_multiwfn=all_submission_information_for_multiwfn, all_submission_information_for_REs=all_submission_information_for_REs, all_submission_information_for_FCs=all_submission_information_for_FCs, get_molecule_atcs=get_molecule_atcs, get_molecule_res=get_molecule_res, get_molecule_fcs=get_molecule_fcs, run_excited_state_from_optimised_ground_structure=run_excited_state_from_optimised_ground_structure)

	# -----------------------------------------------------------------------------------------

	# Twenty-third, determine if the user wants to record unique molecules to disk.
	obtain_unique_molecules_bool = obtain_unique_molecules(molecule_equivalence_method)

	# Twenty-fourth, determine if the user wants to record unique dimers to disk.
	obtain_unique_dimers_bool = obtain_unique_dimers(dimer_equivalence_method)

	# -----------------------------------------------------------------------------------------
	# Twenty-fifth, create a molecule_graphs list where the graph only records the element of each atom. 
	# * We are wanting to determine unique molecules and/or dimers, the ECCP program will use the list of molecule_graphs to help determine which molecules and dimers are the same or different. 
	# * However, currently the graphs in molecule_graphs potentially contain lots of information (features) about the atoms and bonds in the molecules. 
	#   --> For example, examples of the atom and bond features the graph may contain are:
	#                    Atom Features: E (element), mass, initial_magmom, momentum, involved_in_no_of_rings, is_spiro_atom, is_H_acceptor, is_H_donor, no_of_neighbouring_non_cord_H, hybridisation
	#                    Bond Features: involved_in_no_of_rings, is_conjugated, bond_type, bond_type_from_sybyl_type, 'is_cyclic
	# * For comparing unique molecules and/or dimers, we do not want the ECCP to compare all these atoms and bond features.
	#   --> We only want it to compare the elements of the atoms between molecules.
	# 
	# * So that we only compare the elements of each atom between molecules/dimers, we will create a new list of molecule graphs that only contains the element of each atom.
	#   --> Note: we also include "no_of_neighbouring_non_cord_H" as an atom feature in the simple graph as this helps for comparing atoms that do have hydrogen attahced to it, but CCDC does not give them spatial components.
	#
	# * "simple_molecule_graphs" is only used for the get_unique_molecules and get_unique_dimers methods. 

	if obtain_unique_molecules_bool or obtain_unique_dimers_bool:
		simple_molecule_graphs = get_simple_molecule_graphs(molecule_graphs)

	# -----------------------------------------------------------------------------------------

	# Twenty-sixth, if there are unique molecules: 
	if obtain_unique_molecules_bool:
		print(divide_string)

		# Twenty-seventh, get the unique molecules if desired.
		print('Getting Unique Molecules')
		if not has_unique_molecules:
			structurally_unique_molecules_names, structurally_equivalent_molecule_groups, conformationally_unique_molecules_names, conformationally_equivalent_molecule_groups, structurally_equivalent_molecule_pairs = get_unique_molecules(molecules, simple_molecule_graphs, crystal, molecule_equivalence_method=molecule_equivalence_method, neighbouring_molecules_about_molecules=neighbouring_molecules_about_molecules, include_hydrogens_in_uniqueness_analysis=include_hydrogens_in_uniqueness_analysis, no_of_cpus=no_of_cpus)
		else:
			print('Reading Unique Molecules from Equivalence_Group_Information/Structurally_Equivalent_Molecule_Groups.txt and Equivalence_Group_Information/Conformationally_Equivalent_Molecule_Groups.txt')
			structurally_unique_molecules_names, structurally_equivalent_molecule_groups, conformationally_unique_molecules_names, conformationally_equivalent_molecule_groups = convert_existing_unique_molecule_data_from_ECCP_Information(structurally_equivalent_molecule_groups, conformationally_equivalent_molecule_groups, len(molecules))
			# Note: We have already got structurally_equivalent_molecule_pairs. 
		print('No of structurally unique individual molecules: '+str(len(structurally_unique_molecules_names)))
		print('No of conformationally unique individual molecules: '+str(len(conformationally_unique_molecules_names)))

		# Twenty-eighth, Remove unwanted entries from molecules
		for mol_name in molecules.keys():
			molecules[mol_name] = remove_unwanted_entries(molecules[mol_name])

		# Twenty-ninth, save the ATC files of unique molecules if desired.
		print('Writing unique molecules to '+str(path_to_eccp_folder))
		unique_molecules_folderpath                  = path_to_eccp_folder+'/'+'Unique_Molecules'
		unique_molecules_with_environment_folderpath = path_to_eccp_folder+'/'+'Unique_Molecules_with_environment'
		write_molecules_to_disk(structurally_unique_molecules_names, conformationally_unique_molecules_names, molecules, molecule_graphs, neighbouring_molecules_about_molecules, SolventsList, unique_molecules_folderpath, unique_molecules_with_environment_folderpath, unique_atc_calc_jobs_path, unique_re_calc_jobs_path, unique_fc_calc_jobs_path, all_calc_parameters_for_ATCs=all_calc_parameters_for_ATCs, all_calc_parameters_for_multiwfn=all_calc_parameters_for_multiwfn, all_calc_parameters_for_REs=all_calc_parameters_for_REs, all_calc_parameters_for_FCs=all_calc_parameters_for_FCs, all_submission_information_for_ATCs=all_submission_information_for_ATCs, all_submission_information_for_multiwfn=all_submission_information_for_multiwfn, all_submission_information_for_REs=all_submission_information_for_REs, all_submission_information_for_FCs=all_submission_information_for_FCs, get_molecule_atcs=get_molecule_atcs, get_molecule_res=get_molecule_res, get_molecule_fcs=get_molecule_fcs, run_excited_state_from_optimised_ground_structure=run_excited_state_from_optimised_ground_structure)
		if get_molecule_res:
			write_ECCP_process_RE_submit_script(Unique_RE_Calc_Jobs_folder,  all_submission_information_for_REs)
		if get_molecule_fcs:
			write_ECCP_process_FC_submit_script(Unique_FC_Calc_Jobs_folder,  all_submission_information_for_FCs)
	else:
		structurally_unique_molecules_names, structurally_equivalent_molecule_groups, conformationally_unique_molecules_names, conformationally_equivalent_molecule_groups = None, None, None, None

	# Thirtieth, remove any unwanted enteries in the crystal object
	crystal = remove_unwanted_entries(crystal)

	# Thirty-first, save the molecule files to disk.
	print('Writing crystal structure to disk.')
	if not have_ECCP_Information_crystal_file:
		write_crystal_structures_to_disk(crystal, molecules, path_to_eccp_folder, crystal_graph, molecule_graphs, structurally_equivalent_molecule_groups, conformationally_equivalent_molecule_groups)

	# ----------------------------------------------------------------------------------------- #
	# ----------------------------------------------------------------------------------------- #
	# ----------------------------------------------------------------------------------------- #
	# PART 4: Writing input files that involve dimers between pairs of molecules in the crystal #
	# ----------------------------------------------------------------------------------------- #
	# ----------------------------------------------------------------------------------------- #
	# ----------------------------------------------------------------------------------------- #

	print(divide_string)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	# Thirty-first, obtain all dimers of molecules in the crystal. 
	print('Getting All Dimers')
	
	# 32.1: Make sure that the neighbourhood_molecules_for_dimer_method list is sorted by shortest_distance
	neighbourhood_molecules_for_dimer_method.sort(key=lambda x: (x[4], x[0], x[1], x[2][0], x[2][1], x[2][2]))
	
	# 32.2: Obtain the information about all the dimers recorded by ECCP (either in the current ECCP run or from the ECCP_Information file).
	if not has_neighbouring_molecules:
		all_dimers_info = get_dimers(molecules, molecule_graphs, neighbourhood_molecules_for_dimer_method=neighbourhood_molecules_for_dimer_method, no_of_cpus=no_of_cpus)
	else:
		print('Obtaining dimer details from ECCP_Information folder.')
		all_dimers_info = convert_dimer_details_to_all_dimers_info_method(dimer_details, neighbourhood_molecules_for_dimer_method)
	
	# 32.3: Include the name of the dimer in the dimer enteries in neighbourhood_molecules_for_dimer_method
	neighbourhood_molecules_for_dimer_method = add_dimer_name_to_neighbourhood_molecules_for_dimer_method(neighbourhood_molecules_for_dimer_method)
	
	# 32.4: Record the total number of dimers in the list.
	print('Total no. of dimers: '+str(len(all_dimers_info)))

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	# Thirty-third, get all the molecules that neighbour each molecule in the unit cell to be recorded into Gaussian/ORCA/Quantum Chemistry program.
	neighbouring_molecules_about_dimers = get_neighbouring_molecules_about_dimers(environment_settings, neighbourhood_molecules_for_dimer_method, neighbourhood_molecules_for_environment_method)

	# Thirty-fourth, delete neighbourhood_molecules_for_dimer_method and neighbourhood_molecules_for_environment_method. Not needed anymore.
	del neighbourhood_molecules_for_dimer_method
	del neighbourhood_molecules_for_environment_method

	# Thirty-fifth, remove unwanted entries from molecules
	for mol_name in molecules.keys():
		molecules[mol_name] = remove_unwanted_entries(molecules[mol_name])

	# Thirty-sixth, write all dimers and their calc jobs to disk
	print('Writing all dimers to disk.')
	all_dimers_folderpath                  = path_to_eccp_folder+'/'+'All_Dimers'
	all_dimers_with_environment_folderpath = path_to_eccp_folder+'/'+'All_Dimers_with_environment'
	all_dimers_info_names                  = list(all_dimers_info.keys())
	write_dimers_to_disk(all_dimers_info_names, all_dimers_info, molecules, molecule_graphs, neighbouring_molecules_about_dimers, SolventsList, all_dimers_folderpath, all_dimers_with_environment_folderpath, all_eet_calc_jobs_path, all_eigendata_calc_jobs_path, all_calc_parameters_for_EETs=all_calc_parameters_for_EETs, all_submission_information_for_EETs=all_submission_information_for_EETs, all_calc_parameters_for_ICTs=all_calc_parameters_for_ICTs, all_submission_information_for_ICTs=all_submission_information_for_ICTs, get_dimer_eets=get_dimer_eets, get_dimer_icts=get_dimer_icts)

	# Thirty-seventh, if there are unique dimers: 
	if obtain_unique_dimers_bool:
		print(divide_string)

		# Thirty-eighth, obtain the unique, non-symmetric dimers of molecules in the crystal.  
		print('Getting Unique Dimers')
		if not has_unique_dimers:
			unique_dimers_names, structurally_equivalent_dimer_groups, structurally_equivalent_dimer_pairs = get_unique_dimers(all_dimers_info, molecules, simple_molecule_graphs, dimer_equivalence_method=dimer_equivalence_method, neighbouring_molecules_about_dimers=neighbouring_molecules_about_dimers, include_hydrogens_in_uniqueness_analysis=include_hydrogens_in_uniqueness_analysis, no_of_cpus=no_of_cpus)
		else:
			print('Reading Unique Dimers from Structurally_Equivalent_Dimer_Groups.txt')
			unique_dimers_names, structurally_equivalent_dimer_groups = convert_existing_unique_dimer_data_from_ECCP_Information(structurally_equivalent_dimer_groups, len(all_dimers_info))
			# Note: We have already got structurally_equivalent_dimer_pairs. 
		print('No of unique dimers: '+str(len(unique_dimers_names)))

		# Thirty-ninth, remove unwanted entries from molecules
		for mol_name in molecules.keys():
			molecules[mol_name] = remove_unwanted_entries(molecules[mol_name])

		# Fortieth, save the unique on-symmetric dimers to disk.
		print('Writing unique dimers to '+str(path_to_eccp_folder))
		unique_dimers_folderpath                  = path_to_eccp_folder+'/'+'Unique_Dimers'
		unique_dimers_with_environment_folderpath = path_to_eccp_folder+'/'+'Unique_Dimers_with_environment'
		write_dimers_to_disk(unique_dimers_names, all_dimers_info, molecules, molecule_graphs, neighbouring_molecules_about_dimers, SolventsList, unique_dimers_folderpath, unique_dimers_with_environment_folderpath, unique_eet_calc_jobs_path, unique_eigendata_calc_jobs_path, all_calc_parameters_for_EETs=all_calc_parameters_for_EETs, all_submission_information_for_EETs=all_submission_information_for_EETs, all_calc_parameters_for_ICTs=all_calc_parameters_for_ICTs, all_submission_information_for_ICTs=all_submission_information_for_ICTs, get_dimer_eets=get_dimer_eets, get_dimer_icts=get_dimer_icts)
		if get_dimer_eets:
			write_ECCP_process_EET_submit_script(Unique_EET_Calc_Jobs_folder, all_submission_information_for_EETs)
		if get_dimer_icts:
			write_ECCP_process_Eigendata_submit_script(Unique_Eigendata_Calc_Jobs_folder, all_submission_information_for_ICTs)
			write_ECCP_process_ICT_submit_script      (Unique_Eigendata_Calc_Jobs_folder, all_submission_information_for_ICTs)

	# ----------------------------------------------------------------------------------------- #
	# ----------------------------------------------------------------------------------------- #
	# ----------------------------------------------------------------------------------------- #
	#                    PART 5: Final information about molecules and dimers                   #
	# ----------------------------------------------------------------------------------------- #
	# ----------------------------------------------------------------------------------------- #
	# ----------------------------------------------------------------------------------------- #

	# Forty-first, write the result to txt files.

	# 41.1: Write the required infromation for creating the results document (ECCP_Information.txt) to the document_info dict.
	document_info = {'molecules': molecules, 'SolventsList': SolventsList, 'obtain_unique_molecules_bool': obtain_unique_molecules_bool, 'all_dimers_info': all_dimers_info, 'obtain_unique_dimers_bool': obtain_unique_dimers_bool, 'path_to_eccp_folder': path_to_eccp_folder, 'filename': filename}
	document_info['make_dimer_method']           = make_dimer_method
	document_info['environment_settings']        = environment_settings

	# 41.2: Write the information for recording structurally and conformationally equivalent molecules groups. 
	document_info['structurally_equivalent_molecule_groups']     = None if not obtain_unique_molecules_bool else structurally_equivalent_molecule_groups
	document_info['conformationally_equivalent_molecule_groups'] = None if not obtain_unique_molecules_bool else conformationally_equivalent_molecule_groups
	document_info['structurally_equivalent_molecule_pairs']      = None if not obtain_unique_molecules_bool else structurally_equivalent_molecule_pairs

	# 41.3: Write the information for recording structurally equivalent dimers groups. 
	document_info['structurally_equivalent_dimer_groups']       = None if not obtain_unique_dimers_bool else structurally_equivalent_dimer_groups
	document_info['structurally_equivalent_dimer_pairs']        = None if not obtain_unique_dimers_bool else structurally_equivalent_dimer_pairs
	
	# 41.4: Write the results document (ECCP_Information.txt). 
	write_results_document(**document_info)
	
	# Forty-second, conclude with finishing remarks.
	print('Finishing running ECCP upon: '+str(filename))
	print(divide_string)
	print(divide_string)
	print(divide_string)

	# ----------------------------------------------------------------------------------------- #
	# ----------------------------------------------------------------------------------------- #
	# ----------------------------------------------------------------------------------------- #
	#                                    END OF ECCP PROGRAM                                    #
	# ----------------------------------------------------------------------------------------- #
	# ----------------------------------------------------------------------------------------- #
	# ----------------------------------------------------------------------------------------- #


