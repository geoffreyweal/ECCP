"""
read_data_from_ECCP_Information.py, Geoffrey Weal, 28/3/24

This method is designed to record all the information from the ECCP Information folder.
"""

from SUMELF import what_ECCP_Information_files_do_we_have
from SUMELF import run_which_ECCP_operations
from SUMELF import get_ECCP_Information_data
from SUMELF import check_ECCP_Information_details
from SUMELF import get_crystal_file_from_ECCP_Information
from SUMELF import get_dimer_details_data
from SUMELF import get_equivalent_molecule_group_data
from SUMELF import get_equivalent_dimer_group_data
from SUMELF import check_consistancy_between_files

def read_data_from_ECCP_Information(path_to_eccp_folder, make_dimer_method, environment_settings):
	"""
	This method is designed to record all the information from the ECCP Information folder.

	Parameters
	----------
	path_to_eccp_folder : str.
		This is the path to the ECCP_Information data. 
	make_dimer_method : dict.
		This dictionary contains information about the method used to obtain dimers. See https://github.com/geoffreyweal/ECCP for more information. Default: {'method': 'nearest_atoms_method', 'rCut': 5.0}. 
	environment_settings : dict.
		This dictionary contains information about how to treat the local environment about a dimer (where applicable). 

	Returns
	-------
	have_ECCP_Information_crystal_file : bool.
		This boolean indicates if you want to load the crystal.xyz file from the ECCP_Information folder. 
	has_neighbouring_molecules : bool.
		This boolean indicates if you want to load dimer information from All_Dimer_Information.txt. This folder indicates how atoms neighbour each other in the crystal. 
	has_unique_molecules : bool.
		This boolean indicates if you have information about the struccturally and conformationally equivalent molecule groups to load into memory from file. 
	has_unique_dimers : bool.
		This boolean indicates if you have information about the struccturally equivalent dimer groups to load into memory from file. 
	ECCP_Information_data : dict.
		This dictionary contains all the data to load from the ECCP_Information folder. 
	"""

	# First, check if we can load the data from the ECCP_Information folder for this crystal. 
	contains_what_ECCP_Information_files = what_ECCP_Information_files_do_we_have(path_to_eccp_folder)

	# Second, determine which ECCP operations to run, and which can be read from file. 
	has_neighbouring_molecules, has_unique_molecules, has_unique_dimers = run_which_ECCP_operations(**contains_what_ECCP_Information_files)

	# Third, if any of booleans above are true, then we expect to have the ECCP_Information.txt and crystal.xyz file as well.
	have_ECCP_Information_file = have_ECCP_Information_crystal_file = any([has_neighbouring_molecules, has_unique_molecules, has_unique_dimers])

	# Fourth, create a dictionary for holding data from the ECCP_Information folder. 
	#         * This is used to check for consistency between ECCP_Information files. 
	ECCP_Information_data = {}

	# Fifth, if ECCP_Information.txt exists, obtain the information from it and perform checks. 
	if have_ECCP_Information_file:

		# 5.1: Get the information from the ECCP_Information.txt file.
		eccp_information = get_ECCP_Information_data(path_to_eccp_folder)

		# 5.2: Check to make sure that the ECCP_Information.txt file is consistent with the Run_ECCP.py file.
		check_ECCP_Information_details(eccp_information, make_dimer_method, environment_settings)

		# 5.3: Add eccp_information to ECCP_Information_data.
		ECCP_Information_data['eccp_information'] = eccp_information

	# Sixth, if crystal.xyz exists in the ECCP_Information folder, get it and performed checks.
	if have_ECCP_Information_crystal_file:
		
		# 6.1: Obtain the crystal file from the crystal.xyz file from the ECCP_Information folder. 
		crystal = get_crystal_file_from_ECCP_Information(path_to_eccp_folder)

		# 6.2: Add crystal to ECCP_Information_data.
		ECCP_Information_data['crystal'] = crystal

	# Seventh, if you have information about the neighbours around each molecule in the crystal, get it from file. 
	if has_neighbouring_molecules:

		# 7.1: Get dimer_details from file.
		dimer_details = get_dimer_details_data(path_to_eccp_folder, eccp_information)

		# 7.2: Add eccp_information to ECCP_Information_data.
		ECCP_Information_data['dimer_details'] = dimer_details

	# Eighth, if the structural and conformational equivalent molecule groups have been recorded, read them from file. 
	if has_unique_molecules:

		# 8.1: Get structural and conformational equivalent molecule groups from file.
		structurally_equivalent_molecule_groups, conformationally_equivalent_molecule_groups, structurally_equivalent_molecule_pairs = get_equivalent_molecule_group_data(path_to_eccp_folder, eccp_information)

		# 8.2: Add eccp_information to ECCP_Information_data.
		ECCP_Information_data['structurally_equivalent_molecule_groups']     = structurally_equivalent_molecule_groups
		ECCP_Information_data['conformationally_equivalent_molecule_groups'] = conformationally_equivalent_molecule_groups
		ECCP_Information_data['structurally_equivalent_molecule_pairs']      = structurally_equivalent_molecule_pairs

	# Ninth, if the structural equivalent dimer groups have been recorded, read them from file. 
	if has_unique_dimers:

		# 9.1: Get structural equivalent dimer groups from file. 
		structurally_equivalent_dimer_groups, structurally_equivalent_dimer_pairs = get_equivalent_dimer_group_data(path_to_eccp_folder, eccp_information)

		# 9.2: Add eccp_information to ECCP_Information_data.
		ECCP_Information_data['structurally_equivalent_dimer_groups'] = structurally_equivalent_dimer_groups
		ECCP_Information_data['structurally_equivalent_dimer_pairs']  = structurally_equivalent_dimer_pairs

	# Tenth, check that all the files are consistent with each other. 
	check_consistancy_between_files(have_ECCP_Information_file, have_ECCP_Information_crystal_file, has_neighbouring_molecules, has_unique_molecules, has_unique_dimers, **ECCP_Information_data)

	# Eleventh, return booleans for determining what data is contained in the ECCP_Information folder, and all the information from it.
	return have_ECCP_Information_crystal_file, has_neighbouring_molecules, has_unique_molecules, has_unique_dimers, ECCP_Information_data




