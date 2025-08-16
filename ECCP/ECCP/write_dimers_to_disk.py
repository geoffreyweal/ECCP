"""
write_dimers_to_disk.py, Geoffrey Weal, 22/2/22

write_dimers_to_disk is designed to save data about the dimer to disk.
"""
import numpy as np
import networkx as nx
from ase.io import write
from copy import deepcopy

from SUMELF                                                          import make_folder, make_dimer, add_graph_to_ASE_Atoms_object

from ECCP.ECCP.invariance_methods.are_environments_equivalent        import get_environment

from ECCP.ECCP.write_dimers_to_disk_methods.write_EET_gaussian_files import write_EET_gaussian_files
from ECCP.ECCP.write_dimers_to_disk_methods.write_ICT_gaussian_files import write_ICT_gaussian_files

from ECCP.ECCP.write_dimers_to_disk_methods.write_EET_orca_files     import write_EET_orca_files
from ECCP.ECCP.write_dimers_to_disk_methods.write_ICT_orca_files     import write_ICT_orca_files

def write_dimers_to_disk(all_dimers_info_names, all_dimers_info, molecules, molecule_graphs, neighbouring_molecules_about_dimers, SolventsList, dimers_folderpath, dimers_with_environment_folderpath, eet_calc_jobs_path, ict_calc_jobs_path, all_calc_parameters_for_EETs=None, all_submission_information_for_EETs=None, all_calc_parameters_for_ICTs=None, all_submission_information_for_ICTs=None, get_dimer_eets=True, get_dimer_icts=True):
	"""
	This method will save dimer files to disk.

	Parameters
	----------
	all_dimers_info_names : list
		A list of the names of the dimers you want to write to disk.
	all_dimers_info : dict
		This is the dictionary of dimers that have been obtained by this program. This method will only save information about those dimers given in all_dimers_info_names.
	molecules : list of ase.Atoms
		These are the ase.Atoms objects of all the molecules that are used to make up the dimers in your crystal.
	molecule_graphs : list of networkx.graph
		These are all the associated graphs that describe the bonding network for each molecule in the crystal.
	neighbouring_molecules_about_dimers : list
		This is a list of all the details about the molecules that neighbour each of the dimer in the crystal. 
	SolventsList : list of int
		These are the names of the molecules in the molecules list that have been identified at solvents.
		
	dimers_folderpath : str.
		This is the path to save the dimer to.
	dimers_with_environment_folderpath : str.  
		This is the folder to save xyz files of the molecule with its environment surround it.
	eet_calc_jobs_path : str.
		This is the path to save eet calc jobs to.
	ict_calc_jobs_path : str.
		This is the path to save ict calc jobs to.

	all_calc_parameters_for_EETs : list
		This list contain all the information required to be included in the EET calc file(s).
	all_submission_information_for_EETs : list
		This list contain all the information required to be included in the EET submit.sl script. 
	all_calc_parameters_for_ICTs : list
		This list contain all the information required to be included in the ICT calc file(s).
	all_submission_information_for_ICTs : list
		This list contain all the information required to be included in the ICT submit.sl script. 

	get_dimer_eets : bool.
		This tag indicates if the user wants to obtain calc files for running EET jobs on the dimers. Default: True
	get_dimer_icts : bool.
		This tag indicates if the user wants to obtain calc files for running ICT jobs on the dimers, This includes obtaining eigendata (such as overlap orbtials and molecular orbital energies and coefficients). Default: True

	submit_EETs_in_parallel : bool.
		This tag indicates if the user wants to submit EET jobs of different calc parameters from all_calc_parameters in series or parallel. Default: True.
	submit_ICTs_in_parallel : bool.
		This tag indicates if the user wants to submit ICT jobs of different calc parameters from all_calc_parameters in series or parallel. Default: True.

	"""

	# Prestep, make a note if EET has been requested at least once.
	have_requested_EET = False

	# First, make the dimers_folderpath folder if it doesn't already exist
	make_folder(dimers_folderpath)
	if len(neighbouring_molecules_about_dimers) > 0:
		make_folder(dimers_with_environment_folderpath)

	# Second, for each dimer information in the list of dimers
	for dimer_name in sorted(all_dimers_info_names):

		# Second, get information from the index-th position on the all_dimers_info list
		mol_name1, mol_name2, unit_cell_displacement, displacement, move_centre_of_mass_by, shortest_distance = all_dimers_info[dimer_name]

		# Third, indicate which monomers are solvents
		Solvent_name_1 = 'S' if (mol_name1 in SolventsList) else ''
		Solvent_name_2 = 'S' if (mol_name2 in SolventsList) else ''

		# Fourth, make all names needed
		dimer_name_prefix = 'Dimer'+str(dimer_name)
		full_dimer_name   = dimer_name_prefix+'_M'+str(mol_name1)+Solvent_name_1+'_M'+str(mol_name2)+Solvent_name_2
		molecule1_name    = dimer_name_prefix+'_Monomer1_M'+str(mol_name1)+Solvent_name_1
		molecule2_name    = dimer_name_prefix+'_Monomer2_M'+str(mol_name2)+Solvent_name_2

		# Fifth, copy the molecules just to prevent anything being overwritten in the original_file
		dimer, molecule1, molecule2 = make_dimer(molecules, mol_name1, mol_name2, displacement, move_centre_of_mass_by)

		# Sixth, mark out which indicies for with which molecule in the dimer. 
		molecule_1_indices = list(range(len(molecule1)))
		molecule_2_indices = list(range(len(molecule1),len(molecule1)+len(molecule2)))

		# Seventh, record the fragments that each atom in the dimer belongs to for calc. 
		fragment1_list = [1]*len(molecule1)
		fragment2_list = [2]*len(molecule2)

		# Eighth, make the Dimer
		dimer         = molecule1 + molecule2
		fragmentlist = fragment1_list + fragment2_list
		dimer.set_tags(fragmentlist)

		# Ninth, add the graphs (bonding system) of the molecules in the dimer together.
		molecule1_graph   = deepcopy(molecule_graphs[mol_name1])
		molecule2_graph   = deepcopy(molecule_graphs[mol_name2])
		molecule2_mapping = {m2_atom_index: m2_atom_index+len(molecule1_graph) for m2_atom_index in tuple(molecule2_graph.nodes.keys())}
		molecule2_graph   = nx.relabel_nodes(molecule2_graph, molecule2_mapping)
		dimer_graph       = nx.compose(molecule1_graph,molecule2_graph)

		# Tenth, add the graph of the dimer back to the dimer before it is saved as an xyz file.
		add_graph_to_ASE_Atoms_object(dimer, dimer_graph)
		
		# Eleventh, assign the bonding system to it and write it as a xyz file. 
		write(dimers_folderpath+'/'+full_dimer_name+'.xyz', dimer)

		# Twelfth, write the molecules with their environments to file
		if len(neighbouring_molecules_about_dimers) > 0:
			raise Exception('Check if this is working when you first run this.')
			details_of_dimer_to_get_environ_for = (mol_name1, mol_name2, unit_cell_displacement)
			environment_about_dimer = get_environment(details_of_dimer_to_get_environ_for, neighbouring_molecules_about_dimers, molecules)
			environment_about_dimer.set_array('NeighboursList', np.array(['-']*len(environment_about_dimer)))
			dimer_with_environment = dimer.copy() + environment_about_dimer
			fragmentlist += [3]*len(environment_about_dimer)
			dimer_with_environment.set_tags(fragmentlist)
			write(dimers_with_environment_folderpath+'/'+full_dimer_name+'.xyz', dimer_with_environment)
		else:
			environment_about_dimer = None

		# Thirteenth, write the EET files for the dimer
		if get_dimer_eets and (all_calc_parameters_for_EETs is not None):
			for calc_parameters_for_EETs, submission_information_for_EETs in zip(all_calc_parameters_for_EETs, all_submission_information_for_EETs):
				if not 'calc_software' in calc_parameters_for_EETs:
					raise Exception("Error: You need to specify a value for 'calc_software' in calc_parameters_for_EETs.")
				if   calc_parameters_for_EETs['calc_software'].lower() == 'gaussian':
					write_EET_gaussian_files(molecule1, molecule2, full_dimer_name, environment_about_dimer, eet_calc_jobs_path, fragmentlist, calc_parameters_for_EETs, submission_information_for_EETs)
				elif calc_parameters_for_EETs['calc_software'].lower() == 'orca':
					if not have_requested_EET:
						print('Note: There is no EET function in ORCA')
					#raise Exception('Need to write this part.')
					#write_EET_orca_files    (molecule1, molecule2, full_dimer_name, environment_about_dimer, eet_calc_jobs_path, fragmentlist, calc_parameters_for_EETs, submission_information_for_EETs)
				else:
					raise Exception("Error: calc_parameters_for_EETs['calc_software'] needs to be either Gaussian or ORCA. calc_parameters_for_EETs['calc_software'] = "+str(calc_parameters_for_EETs['calc_software']))
			have_requested_EET = True
		
		# Fourteenth, write the ICT files for the dimer
		if get_dimer_icts  and (all_calc_parameters_for_ICTs is not None):
			raise Exception('Need to write this part.')
			'''
			for original_calc_parameters, original_submission_information in zip(all_calc_parameters_for_EETs, all_submission_information_for_EETs):
				if   all_calc_parameters_for_ICTs['calc_software'] == 'gaussian':
					write_ICTs_gaussian_files(dimer, molecule1, molecule2, environment_about_dimer, full_dimer_name, icts_calc_jobs_path, all_calc_parameters_for_ICTs, all_submission_information_for_ICTs)
				elif all_calc_parameters_for_ICTs['calc_software'] == 'orca':
					write_ICTs_orca_files    (dimer, molecule1, molecule2, environment_about_dimer, full_dimer_name, icts_calc_jobs_path, all_calc_parameters_for_ICTs, all_submission_information_for_ICTs)
				else:
					raise Exception('Error here.')
			'''
			#write_ICTs_gaussian_files(dimer, molecule1, molecule2, full_dimer_name, environment_about_dimer, icts_calc_jobs_path, all_calc_parameters_for_ICTs, all_submission_information_for_ICTs)

		# Fiftheenth, add the graph of each molecule in the dimer back to the molecule itself before it is saved.
		add_graph_to_ASE_Atoms_object(molecule1, deepcopy(molecule_graphs[mol_name1]))
		add_graph_to_ASE_Atoms_object(molecule2, deepcopy(molecule_graphs[mol_name2]))

		# Sixteenth, save the individual xyz files of each molecule in the dimer
		make_folder(dimers_folderpath+'/'+full_dimer_name)
		write(dimers_folderpath+'/'+full_dimer_name+'/'+molecule1_name+'.xyz', molecule1)
		write(dimers_folderpath+'/'+full_dimer_name+'/'+molecule2_name+'.xyz', molecule2)

		# Seventeenth, save information about the specific dimer. 
		which_indices_go_with_which_molecule(molecule1_name, molecule2_name, molecule_1_indices, molecule_2_indices, dimers_folderpath+'/'+full_dimer_name, full_dimer_name)

# for testing xyz integration grid
from math import pi
import numpy as np
def get_rotation(roll_angle, pitch_angle, yaw_angle):

	roll_angle *= pi/180.0
	pitch_angle *= pi/180.0
	yaw_angle *= pi/180.0

	roll_matrix  = np.array([[1, 0, 0], [0, np.cos(roll_angle), -np.sin(roll_angle)], [0, np.sin(roll_angle), np.cos(roll_angle)]])
	pitch_matrix = np.array([[np.cos(pitch_angle), 0, np.sin(pitch_angle)], [0, 1, 0], [-np.sin(pitch_angle), 0, np.cos(pitch_angle)]])
	yaw_matrix   = np.array([[np.cos(yaw_angle), -np.sin(yaw_angle), 0], [np.sin(yaw_angle), np.cos(yaw_angle), 0], [0, 0, 1]])

	return yaw_matrix @ pitch_matrix @ roll_matrix

# ----------------------------------------------------------------------------------------------------------------------------------

def which_indices_go_with_which_molecule(molecule1_name, molecule2_name, molecule_1_indices, molecule_2_indices, path_to, full_dimer_name):
	"""
	This method will write information about molecules in the dimer to disk

	Parameters
	----------
	molecule1_name : str.
		This is the name of the first molecule in the dimer.
	molecule2_name : str.
		This is the name of the second molecule in the dimer.

	molecule_1_indices : list
		This list contains the indicies of atoms that are associated with the first molecule in the dimer. 
	molecule_2_indices : list
		This list contains the indicies of atoms that are associated with the second molecule in the dimer. 

	path_to : str.
		This is the path to place this information file in.
	full_dimer_name : str.
		This is the full name of the dimer. 
	"""
	data_filename = full_dimer_name+'_Information.txt'
	with open(path_to+'/'+data_filename,'w') as Dimer_Information_TXT:
		Dimer_Information_TXT.write('Dimer Information about: '+str(full_dimer_name)+'\n')
		Dimer_Information_TXT.write('\n')
		Dimer_Information_TXT.write('Atom indices in '+str(molecule1_name)+': '+str(molecule_1_indices)+'\n')
		Dimer_Information_TXT.write('Atom indices in '+str(molecule2_name)+': '+str(molecule_2_indices)+'\n')

# ----------------------------------------------------------------------------------------------------------------------------------

