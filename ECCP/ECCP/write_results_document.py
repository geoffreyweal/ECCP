"""
write_results_document.py, Geoffrey Weal, 18/2/22

This script includes methods for writing the results from the ECCP method.
"""
import os
from SUMELF import make_folder
from ECCP.ECCP.write_results_document_method.write_all_dimer_information import write_all_dimer_information

def write_results_document(molecules : list, SolventsList : list, obtain_unique_molecules_bool : bool, all_dimers_info : list, obtain_unique_dimers_bool : bool, make_dimer_method : dict, environment_settings : dict, structurally_equivalent_molecule_groups : dict, conformationally_equivalent_molecule_groups : dict, structurally_equivalent_molecule_pairs : dict, structurally_equivalent_dimer_groups : dict, structurally_equivalent_dimer_pairs : list, path_to_eccp_folder : str, filename : str):
	"""
	This method will write the results from the strutural analysis of the molecules and dimers in the crystal as performed by this ECCP program. 

	Parameters
	----------
	molecules : list of ase.Atoms objects.
		These are all the individual molecules identified in the crystal that you want to determine dimers for.
	SolventsList : list of int
		These are the indices of the molecules in the molecules list that have been identified at solvents.
	obtain_unique_molecules_bool : bool.
		This tag indicates if the user wanted to obtain the unique molecules.

	all_dimers_info : list 
		This is all the information about all the dimers that have been idenified in the crystal using the settings as given by the user.
	obtain_unique_dimers_bool : bool.
		This tag indicates if the user wanted to obtain the unique dimers.

	make_dimer_method : dict.
		This contains all the information about how the dimers were created
	environment_settings : dict.
		This contains all the information about how molecules surrounding the molecules and dimers are included in calculations. 

	structurally_equivalent_molecule_groups : dict.
		This dictionary contains the structurally equivalent molecule groups, given as --> representative structurally unique molecule: list of structurally equivalent molecules.
	conformationally_equivalent_molecule_groups : dict.
		This dictionary contains the conformationally equivalent molecule groups, given as --> representative conformationally unique molecule: list of conformationally equivalent molecules.
	structurally_equivalent_molecule_pairs : list
		This list contains all the pairs of structurally equivalent molecules. 

	structurally_equivalent_dimer_groups : dict.
		This dictionary contains the structurally equivalent dimer groups, given as --> representative structurally unique dimer: list of structurally equivalent dimers.
	structurally_equivalent_dimer_pairs : dict
		This list contains all the pairs of structurally equivalent dimers. 

	path_to_eccp_folder : str.
		This is the path to this crystal in the ECCP folder
	filename : str.
		This is the name of the crystal, given by its folder name.
	"""

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# First, obtain the lists of unique and equivalent molecules. 

	if obtain_unique_molecules_bool:

		# 1.1: Obtain the structurally unique molecule group indices.
		structurally_unique_molecule_indices         = tuple(sorted(structurally_equivalent_molecule_groups.keys()))

		# 1.2: Obtain the structurally equivalent molecule group indices.
		structurally_equivalent_molecule_indices     = tuple(sorted([j for sub in structurally_equivalent_molecule_groups.values() for j in sub]))

		# 1.3: Obtain the conformationally unique molecule group indices.
		conformationally_unique_molecule_indices     = tuple(sorted(conformationally_equivalent_molecule_groups.keys()))

		# 1.4: Obtain the conformationally equivalent molecule group indices.
		conformationally_equivalent_molecule_indices = tuple(sorted([j for sub in conformationally_equivalent_molecule_groups.values() for j in sub]))

	else:

		structurally_unique_molecule_indices = None
		structurally_equivalent_molecule_indices = None
		conformationally_unique_molecule_indices = None
		conformationally_equivalent_molecule_indices = None

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# Second, obtain the lists of structurally unique and equivalent dimers. 

	if obtain_unique_dimers_bool:

		# 2.1: Obtain the structurally unique dimer indices.
		structurally_unique_dimer_indices     = tuple(sorted(structurally_equivalent_dimer_groups.keys()))

		# 2.2: Obtain the structurally equivalent dimer indices.
		structurally_equivalent_dimer_indices = tuple(sorted([j for sub in structurally_equivalent_dimer_groups.values() for j in sub]))

	else:

		structurally_unique_dimer_indices = None
		structurally_equivalent_dimer_indices = None

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Third, get a list of all the names of the molecules that are solvents.
	#        * Given as a list of strings
	molecule_names_that_are_solvents = [str(mol_no) for mol_no in sorted(SolventsList)]

	# Fourth, get a list of all the names of the dimers that contain solvents.
	#        * Given as a list of strings
	dimers_that_contain_solvents = [] 
	for dimer_no, (molecule1_name, molecule2_name, _, _, _, _) in sorted(all_dimers_info.items(), key=lambda x: x[0]):
		if (molecule1_name in SolventsList) or (molecule2_name in SolventsList):
			dimers_that_contain_solvents.append(dimer_no)
	dimers_that_contain_solvents = [str(dimer_no) for dimer_no in sorted(dimers_that_contain_solvents)]

	# Fifth, make the ECCP_Information.txt file.
	make_ECCP_Information_file(path_to_eccp_folder, filename, make_dimer_method, environment_settings, molecules, molecule_names_that_are_solvents, obtain_unique_molecules_bool, structurally_unique_molecule_indices, structurally_equivalent_molecule_indices, conformationally_unique_molecule_indices, conformationally_equivalent_molecule_indices, all_dimers_info, dimers_that_contain_solvents, obtain_unique_dimers_bool, structurally_unique_dimer_indices, structurally_equivalent_dimer_indices)

	# Sixth, create a folder for placing equivalency groups into.
	equivalency_group_folder_name = 'Equivalence_Group_Information'
	make_folder(path_to_eccp_folder+'/'+equivalency_group_folder_name)

	# Seventh, record all the information about the structurally and conformationally unique and equivalent molecule groups in the crystal.
	if obtain_unique_molecules_bool:
		write_molecule_equivalence_groups(path_to_eccp_folder+'/'+equivalency_group_folder_name, structurally_equivalent_molecule_groups, conformationally_equivalent_molecule_groups)
		write_molecule_equivalence_pairs (path_to_eccp_folder+'/'+equivalency_group_folder_name, structurally_equivalent_molecule_pairs)

	# Eighth, record all the information about the dimers created during this ECCP program
	write_all_dimer_information(path_to_eccp_folder, all_dimers_info)

	# Ninth, record all the information about the structurally unique and equivalent dimer groups created during this ECCP program.
	if obtain_unique_dimers_bool:
		write_dimer_equivalence_groups(path_to_eccp_folder+'/'+equivalency_group_folder_name, structurally_equivalent_dimer_groups)
		write_dimer_equivalence_pairs (path_to_eccp_folder+'/'+equivalency_group_folder_name, structurally_equivalent_dimer_pairs)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def make_ECCP_Information_file(path_to_eccp_folder, filename, make_dimer_method, environment_settings, molecules, molecule_names_that_are_solvents, obtain_unique_molecules_bool, structurally_unique_molecule_indices, structurally_equivalent_molecule_indices, conformationally_unique_molecule_indices, conformationally_equivalent_molecule_indices, all_dimers_info, dimers_that_contain_solvents, obtain_unique_dimers_bool, structurally_unique_dimer_indices, structurally_equivalent_dimer_indices):
	"""
	This method is designed to create the ECCP_Information.txt file.

	Parameters
	----------
	path_to_eccp_folder : str.
		This is the path to this crystal in the ECCP folder
	filename : str.
		This is the name of the crystal, given by its folder name.

	make_dimer_method : dict.
		This contains all the information about how the dimers were created
	environment_settings : dict.
		This contains all the information about how molecules surrounding the molecules and dimers are included in calculations. 

	molecules : list of ase.Atoms objects.
		These are all the individual molecules identified in the crystal that you want to determine dimers for.
	molecule_names_that_are_solvents : list of ints
		This is a list of the names of all the molecules that are solvents.

	obtain_unique_molecules_bool : bool.
		This tag indicates if the user wanted to obtain the unique molecules.

	structurally_unique_molecule_indices : tuple
		This list contains all the structurally unique molecules in the crystal. 
	structurally_equivalent_molecule_indices : tuple
		This list contains all the structurally equivalent molecules in the crystal. 
	conformationally_unique_molecule_indices : tuple
		This list contains all the conformationally unique molecules in the crystal. 
	conformationally_equivalent_molecule_indices : tuple
		This list contains all the conformationally equivalent molecules in the crystal. 

	all_dimers_info : list 
		This is all the information about all the dimers that have been idenified in the crystal using the settings as given by the user.
	dimers_that_contain_solvents : list of ints
		This is a list of the names of all the molecules that are solvents.

	obtain_unique_dimers_bool : bool.
		This tag indicates if the user wanted to obtain the unique dimers.

	structurally_unique_dimer_indices : tuple
		This list contains all the structurally unique dimers recorded.
	structurally_equivalent_dimer_indices : tuple
		This list contains all the structurally equivalent dimers recorded.
	"""
	data_filename = 'ECCP_Information.txt'
	with open(path_to_eccp_folder+'/'+data_filename,'w') as ECCP_Information_TXT:
		ECCP_Information_TXT.write('ECCP Information about: '+str(filename)+'\n')
		ECCP_Information_TXT.write('\n')
		ECCP_Information_TXT.write('make_dimer_method: '+str(make_dimer_method)+'\n')
		ECCP_Information_TXT.write('environment_settings: '+str(environment_settings)+'\n')
		ECCP_Information_TXT.write('\n')
		ECCP_Information_TXT.write('---------------------------------------------------------------\n')
		ECCP_Information_TXT.write('---------------------------------------------------------------\n')
		ECCP_Information_TXT.write('\n')
		ECCP_Information_TXT.write('Number of molecules obtained from crystal: '+str(len(molecules))+'\n')
		ECCP_Information_TXT.write('\n')
		ECCP_Information_TXT.write('Molecules that are solvents: '+(', '.join(molecule_names_that_are_solvents) if (len(molecule_names_that_are_solvents) > 0) else 'None') +'\n')
		ECCP_Information_TXT.write('\n')
		ECCP_Information_TXT.write('---------------------------------------------------------------\n')
		if obtain_unique_molecules_bool:
			ECCP_Information_TXT.write('\n')
			ECCP_Information_TXT.write('Number of structurally unique molecules obtained from crystal: '+(str(len(structurally_unique_molecule_indices)) if structurally_unique_molecule_indices is not None else '-')+'\n')
			ECCP_Information_TXT.write('Number of structurally equivalent molecules obtained from crystal: '+(str(len(structurally_equivalent_molecule_indices)) if structurally_equivalent_molecule_indices is not None else '-')+'\n')
			ECCP_Information_TXT.write('\n')
			ECCP_Information_TXT.write('Number of conformationally unique molecules obtained from crystal: '+(str(len(conformationally_unique_molecule_indices)) if conformationally_unique_molecule_indices is not None else '-')+'\n')
			ECCP_Information_TXT.write('Number of conformationally equivalent molecules obtained from crystal: '+(str(len(conformationally_equivalent_molecule_indices)) if conformationally_equivalent_molecule_indices is not None else '-')+'\n')
			ECCP_Information_TXT.write('\n')
			ECCP_Information_TXT.write('---------------------------------------------------------------\n')
		ECCP_Information_TXT.write('---------------------------------------------------------------\n')
		ECCP_Information_TXT.write('\n')
		ECCP_Information_TXT.write('Number of dimers obtained from crystal: '+str(len(all_dimers_info))+'\n')
		ECCP_Information_TXT.write('\n')
		ECCP_Information_TXT.write('Dimers that contain solvents: '+(', '.join(dimers_that_contain_solvents) if (len(dimers_that_contain_solvents) > 0) else 'None') +'\n')
		ECCP_Information_TXT.write('\n')
		ECCP_Information_TXT.write('---------------------------------------------------------------\n')
		if obtain_unique_dimers_bool:
			ECCP_Information_TXT.write('\n')
			ECCP_Information_TXT.write('Number of structurally unique dimers obtained from crystal: '+(str(len(structurally_unique_dimer_indices)) if structurally_unique_dimer_indices is not None else '-')+'\n')
			ECCP_Information_TXT.write('Number of structurally equivalent dimers obtained from crystal: '+(str(len(structurally_equivalent_dimer_indices)) if structurally_equivalent_dimer_indices is not None else '-')+'\n')
			ECCP_Information_TXT.write('\n')
			ECCP_Information_TXT.write('---------------------------------------------------------------\n')
		ECCP_Information_TXT.write('---------------------------------------------------------------\n')

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def write_molecule_equivalence_groups(path_to_equivalency_group_folder, structurally_equivalent_molecule_groups, conformationally_equivalent_molecule_groups):
	"""
	This method will create a file of all the structurally and conformationally equivalent molecule groups in the crystal.

	Parameters
	----------
	path_to_equivalency_group_folder : str.
		This is the path to this crystal in the ECCP folder, specifically to the folder containing equivalency group information. 
	structurally_equivalent_molecule_groups : dict.
		This dictionary contains the structurally equivalent molecule groups, where the key is the representative structurally unique molecule, and the value is the list of structurally equivalent molecules.
	conformationally_equivalent_molecule_groups : dict.
		This dictionary contains the conformationally equivalent molecule groups, where the key is the representative conformationally unique molecule, and the value is the list of conformationally equivalent molecules.
	"""

	# Write all the information about the structurally equivalent molecule groups to file.
	data_filename = 'Structurally_Equivalent_Molecule_Groups.txt'
	with open(path_to_equivalency_group_folder+'/'+data_filename,'w') as Molecule_Information_TXT:
		for unique_mol_name, symmetric_mol_names in sorted(structurally_equivalent_molecule_groups.items()):
			structurally_equivalent_molecule_group = sorted([unique_mol_name] + [symmetric_mol_name for symmetric_mol_name in symmetric_mol_names])
			Molecule_Information_TXT.write(str(structurally_equivalent_molecule_group)+'\n')

	# Write all the information about the conformationally equivalent molecule groups to file.
	data_filename = 'Conformationally_Equivalent_Molecule_Groups.txt'
	with open(path_to_equivalency_group_folder+'/'+data_filename,'w') as Molecule_Information_TXT:
		for unique_mol_name, symmetric_mol_names in sorted(conformationally_equivalent_molecule_groups.items()):
			conformationally_equivalent_molecule_group = sorted([unique_mol_name] + [symmetric_mol_name for symmetric_mol_name in symmetric_mol_names])
			Molecule_Information_TXT.write(str(conformationally_equivalent_molecule_group)+'\n')

def write_molecule_equivalence_pairs(path_to_equivalency_group_folder, structurally_equivalent_molecule_pairs):
	"""
	This method will create a file of all the structurally equivalent molecule pairs in the crystal.

	Parameters
	----------
	path_to_equivalency_group_folder : str.
		This is the path to this crystal in the ECCP folder, specifically to the folder containing equivalency group information. 
	structurally_equivalent_molecule_pairs : list
		This list hold all the pairs of structrually equivalent molecules. 
	"""
	data_filename = 'Structurally_Equivalent_Molecule_Pairs.txt'
	with open(path_to_equivalency_group_folder+'/'+data_filename,'w') as Molecule_Information_TXT:
		Molecule_Information_TXT.write(str([(mol1_name, mol2_name) for mol1_name, mol2_name in sorted(structurally_equivalent_molecule_pairs)])+'\n')

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def write_dimer_equivalence_groups(path_to_equivalency_group_folder, structurally_equivalent_dimer_groups):
	"""
	This method will create a file of all the structurally equivalent dimer groups obtained in this ECCP run.

	Parameters
	----------
	path_to_equivalency_group_folder : str.
		This is the path to this crystal in the ECCP folder, specifically to the folder containing equivalency group information. 
	structurally_equivalent_molecule_groups : dict.
		This dictionary contains the structurally equivalent dimer groups, where the key is the representative structurally unique dimer, and the value is the list of structurally equivalent dimers.
	"""

	# Write all the information about the structurally equivalent dimer groups to file.
	data_filename = 'Structurally_Equivalent_Dimer_Groups.txt'
	with open(path_to_equivalency_group_folder+'/'+data_filename,'w') as Dimer_Information_TXT:
		for unique_dimer_name, symmetric_dimer_names in sorted(structurally_equivalent_dimer_groups.items()):
			conformationally_equivalent_dimer_group = sorted([unique_dimer_name] + [symmetric_dimer_name for symmetric_dimer_name in symmetric_dimer_names])
			Dimer_Information_TXT.write(str(conformationally_equivalent_dimer_group)+'\n')

def write_dimer_equivalence_pairs(path_to_equivalency_group_folder, structurally_equivalent_dimer_pairs):
	"""
	This method will create a file of all the structurally equivalent dimer pairs obtained in this ECCP run.

	Parameters
	----------
	path_to_equivalency_group_folder : str.
		This is the path to this crystal in the ECCP folder, specifically to the folder containing equivalency group information. 
	structurally_equivalent_dimer_pairs : list
		This list hold all the pairs of structrually equivalent dimers. 
	"""
	data_filename = 'Structurally_Equivalent_Dimer_Pairs.txt'
	with open(path_to_equivalency_group_folder+'/'+data_filename,'w') as Dimer_Information_TXT:
		Dimer_Information_TXT.write(str([(dimer1_name, dimer2_name) for dimer1_name, dimer2_name in sorted(structurally_equivalent_dimer_pairs)])+'\n')

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------

