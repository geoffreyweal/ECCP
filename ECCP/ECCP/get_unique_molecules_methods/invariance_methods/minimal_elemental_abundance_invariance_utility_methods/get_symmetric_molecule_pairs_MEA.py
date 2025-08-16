"""
get_symmetric_molecule_pairs_MEA.py, Geoffrey Weal, 5/4/22

This script is designed to use the procrustes analysis to determine if dimers are rotationally, translationally, and reflectively invariant.
"""
import sys, warnings
from tqdm import tqdm

import multiprocessing as mp

from ECCP.ECCP.get_unique_molecules_methods.invariance_methods.minimal_elemental_abundance_invariance_utility_methods.are_molecules_variant_MEA                import are_molecules_variant_MEA
from ECCP.ECCP.get_unique_molecules_methods.invariance_methods.minimal_elemental_abundance_invariance_utility_methods.are_molecules_variant_from_comprehensive import are_molecules_variant_from_comprehensive

def get_symmetric_molecule_pairs_MEA(unique_molecules_names, non_hydrogen_molecules_elements, non_hydrogen_molecules_positions, all_no_of_H_on_atoms_in_molecule, max_distance_disparity, neighbouring_molecules_about_molecules, non_hydrogen_molecules, non_hydrogen_graphs, no_of_cpus=1):
	"""
	This method looks for the elements in lowest abundance and aligns them between molecules to determine if they are the same or not.

	This method works best if your molecule contains one or a few elements that are in low abundance in the molecule.

	Parameters
	----------
	unique_molecules_names : list of ints
		This is a list of unique molecules based on the symmetry of the crystal.
	non_hydrogen_molecules_elements : dict of str. 
		These are the list of the elements of the atoms for each molecule in non_hydrogen_molecules
	non_hydrogen_molecules_positions : dict of 2D np.array
		These are the list of the positions of the atoms for each molecule in non_hydrogen_molecules
	all_no_of_H_on_atoms_in_molecule : dict of int
		This is a list of all the hydrogens that are bound to each "heavy" atom in the molecule
	max_distance_disparity : float
		This is the maximum disparity between two molecules to be considered invariant. If max_distance_disparity is given as None, the default value will be given. Default: 0.01 Å.
	neighbouring_molecules_about_molecules : dict.
		This is the information about the molecules that surround (in the vicinity of) each moleule in the crystal.
	non_hydrogen_molecules : dict of ase.Atoms
		These are the original ase.Atoms object with their hydrogens removed. This is for debugging. 
	non_hydrogen_graphs : dict. of networkx.Graph
		These are the non-hydrogen graphs for the corresponding molecules in non_hydrogen_molecules. Here, hydrogens has been removed as nodes and instead been append to each "heavy" atom node as a node variable. 
	no_of_cpus : int.
		This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.

	Returns
	-------
	symmetric_molecule_pairs : dict
		This dictionary stores which molecules are symmetric to each other, as well as which indices in molecule 1 map onto which indices in molecule 2. 
	"""

	# Prestep: Set the number of cores for each comparison of two molecules
	no_of_cpus_for_comparing_two_molecules = no_of_cpus

	# First, determine the number of molecule pairwise comparisons that will be made. 
	nn = int((len(unique_molecules_names)*(len(unique_molecules_names)-1))/2)

	# Second, obtain the pairs of symmetric molecules in the crystal. 
	if True: #no_of_cpus == 1: # If the user only wants to use 1 cpu, perform tasks without using multiprocessing

		# 2.1: Initialise the symmetric_molecule_pairs dictionary. This will hold all the information about symmetric dimers, as well as which names in molecule 1 map onto which names in molecule 2. 
		symmetric_molecule_pairs = {}

		# 2.2: Create a progress bar for running this task
		pbar = tqdm(get_inputs(unique_molecules_names, non_hydrogen_molecules_elements, non_hydrogen_molecules_positions, all_no_of_H_on_atoms_in_molecule, max_distance_disparity, neighbouring_molecules_about_molecules, non_hydrogen_molecules, non_hydrogen_graphs, symmetric_molecule_pairs, no_of_cpus_for_comparing_two_molecules), total=nn, unit='molecule pair')

		# 2.3: For each comparison of molecules.
		for input_data in pbar:

			# 2.3.1: Obtain the names of the two molecules being compared. 
			mol_name1 = input_data[0][0]
			mol_name2 = input_data[1][0]

			# 2.3.2: Write a message to tell the user what two molecules are being comapred. 
			pbar.set_description('Comparing molecules '+str(mol_name1)+' and '+str(mol_name2))

			# 2.3.3: Compare the two molecules. 
			compare_if_two_molecules_are_symmetric_single_process(input_data)

		# 2.4: Close the progress bar. 
		pbar.close()

	else:

		# 2.5: Create the manager to save lists to
		with mp.Manager() as manager:

			# 2.6: Initialise the manager dictionary to save symmetric molecule information to during the compare_if_two_molecules_are_symmetric_single_process method
			symmetric_molecule_pairs = manager.dict()

			# 2.7: Write a message for the user
			print('Comparing dimers (Please wait until after 100%, as the process will still be running.)', file=sys.stderr)

			# 2.8: Perform compare_two_dimers on each way that dimer 2 can be mapped onto dimer 1 using multiprocessing. 
			#process_map(compare_if_two_molecules_are_symmetric_single_process, get_inputs(unique_molecules_names, non_hydrogen_molecules_elements, non_hydrogen_molecules_positions, all_no_of_H_on_atoms_in_molecule, max_distance_disparity, neighbouring_molecules_about_molecules, non_hydrogen_molecules, non_hydrogen_graphs, symmetric_molecule_pairs), total=nn, unit='molecule pair', desc='Comparing molecules', max_workers=no_of_cpus)
			pool = mp.Pool(processes=no_of_cpus)
			pool.map_async(compare_if_two_molecules_are_symmetric_single_process, tqdm(get_inputs(unique_molecules_names, non_hydrogen_molecules_elements, non_hydrogen_molecules_positions, all_no_of_H_on_atoms_in_molecule, max_distance_disparity, neighbouring_molecules_about_molecules, non_hydrogen_molecules, non_hydrogen_graphs, symmetric_molecule_pairs, no_of_cpus_for_comparing_two_molecules), total=nn, unit='molecule pair', desc='Comparing molecules'))
			pool.close()
			pool.join()

			# 2.9: Convert symmetric_molecule_pairs from a manager.dict object to a dict object.
			symmetric_molecule_pairs = dict(symmetric_molecule_pairs)

	# Third, return symmetric_molecule_pairs
	return symmetric_molecule_pairs

# ---------------------------------------------------------------------------------------------------------------------------------------------------------

def compare_if_two_molecules_are_symmetric_single_process(input_data):
	"""
	This method is designed to determine if two molecules are symmetric to each other. 

	Parameters
	----------
	mol_name1 : int
		This is the name of the first molecule we want to compare.
	molecule1_elements : list of str.
		This is the list of elements that make up molecule 1.
	molecule1_positions : 2D np.array
		These are the positions of the atoms in molecule 1
	no_of_H_on_atoms_in_molecule1 : list of int.
		These are the number of hydrogen atoms bound to each "heavy" atom in molecule 1. 

	mol_name2 : int
		This is the name of the second molecule we want to compare.
	molecule2_elements : list of str.
		This is the list of elements that make up molecule 2.
	molecule2_positions : 2D np.array
		These are the positions of the atoms in molecule 2.
	no_of_H_on_atoms_in_molecule2 : list of int.
		These are the number of hydrogen atoms bound to each "heavy" atom in molecule 2.

	em_indices_m2_to_m1 : dict
		This dictionary indicates how atoms in molecule 2 could be mapped onto molecule 1. 
	max_distance_disparity: float
		This is the maximum that any two "could be equivalent" atoms can be between molecule 1 and molecule 2 for molecules 1 and 2 to be considered variant.
	neighbouring_molecules_about_molecules : dict.
		This is the information about the molecules that surround (in the vicinity of) each moleule in the crystal.
	non_hydrogen_molecules : list of ase.Atoms
		This is the list of molecules in the crystal, not including hydrogens.

	symmetric_molecule_pairs : dict
		This dictionary stores which molecules are symmetric to each other, as well as which indices in molecule 1 map onto which indices in molecule 2. This list is populated during this method. 

	no_of_cpus_for_comparing_two_molecules : int
		This is the number of cpus reserved for comparing two molecules together.
	"""

	# First, extract the input variables from input_data
	molecule_1_information, molecule_2_information, max_distance_disparity, neighbouring_molecules_about_molecules, non_hydrogen_molecules, m1_non_hydrogen_graph, m2_non_hydrogen_graph, symmetric_molecule_pairs, no_of_cpus_for_comparing_two_molecules = input_data
	
	# Second, get the individual variables from molecule_1_information. 
	m1_name, m1_original_elements, m1_original_positions, no_of_H_on_atoms_in_molecule1 = molecule_1_information

	# Third, get the individual variables from molecule_1_information. 
	m2_name, m2_original_elements, m2_original_positions, no_of_H_on_atoms_in_molecule2 = molecule_2_information

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Fourth, make a tuple of the names of the molecules being compared against each other.
	molecules_being_compared = (m1_name, m2_name)

	# Fifth, determine if the two molecules are structurally variant of each other.
	#        * We determine if we can use the Minimal Elemental Abundance Method based on if there are enough atoms in the molecules to use it. 
	#        * If we do not have enough atoms to use this method, default to using the comprehensive method. 
	if (len(m1_original_elements) >= 4) and (len(m2_original_elements) >= 4):

		# 5.1: Determine if the molecules are variants using the Minimal Elemental Abundance Method. 
		positions_are_variant_MEA, mol1_to_mol2_conversion_MEA, align_vectors_gave_a_warning = are_molecules_variant_MEA(m1_original_elements, m1_original_positions, no_of_H_on_atoms_in_molecule1, m2_original_elements, m2_original_positions, no_of_H_on_atoms_in_molecule2, max_distance_disparity, molecules_being_compared, neighbouring_molecules_about_molecules, non_hydrogen_molecules, no_of_cpus=no_of_cpus_for_comparing_two_molecules)

		# 5.2: If the positions between the two molecule are variant:
		if positions_are_variant_MEA:
			
			# 5.3: Found that these molecules are variants using the Minimal Elemental Abundance Method. Append this molecule pair to the symmetric_dimer_pairs list, along with the indices of each molecule that map to each other. 
			symmetric_molecule_pairs[(m1_name, m2_name)] = [(mn_i, mn_j) for mn_i, mn_j in mol1_to_mol2_conversion_MEA.items()]
			return

		# 5.3: If:
		#          1. The Minimal Elemental Abundance (MEA) method could not find a way that these molecules were variant and, 
		#          2. There was a warning issue given when running the Rotation.align_vectors method in scipy, 
		#      We will perform the comprehensive method just to make sure these two molecules are structurally equivalent or not. 
		if len(align_vectors_gave_a_warning) > 0:

			# 5.5: Print to the user that you obtain a warning from Rotation.align_vectors method in scipy, and did not find a way to show that these two molecules are structurally equivalent,
			#      so we will make sure these two molecules are structurally equivalent or not using the comprehensive method.
			to_string = f'We received a warning from from scipy.spatial.transform.Rotation when using Rotation.align_vectors in the Minimal Elemental Abundance (MEA) method. Because we did not find structural equivalency between molecules {m1_name+1} and {m2_name+1} using the MEA method, we will double check this using the comprehensive method.'
			#warnings.warn(to_string)

			# 5.6: Determine if the molecules are variants using the Comprehensive Method. 
			#      * There are not enough different elements to perform the Minimal Elemental Abundance Method, so fall back to the more 
			#        computationally demanding comprehensive method. 
			positions_are_variant_Comp, mol1_to_mol2_conversion_Comp = are_molecules_variant_from_comprehensive(m1_original_elements, m1_original_positions, no_of_H_on_atoms_in_molecule1, m2_original_elements, m2_original_positions, no_of_H_on_atoms_in_molecule2, max_distance_disparity, molecules_being_compared, neighbouring_molecules_about_molecules, non_hydrogen_molecules, m1_non_hydrogen_graph, m2_non_hydrogen_graph, no_of_cpus=no_of_cpus_for_comparing_two_molecules)

			raise Exception('Error: Need to check the next stepo.')
			import pdb; pdb.set_trace()

			# 5.7: If the positions between the two molecule are variant:
			if positions_are_variant_Comp:

				# 5.8: Found that these molecules are variants using the Comprehensive Method. Append this molecule pair to the symmetric_dimer_pairs list, along with the indices of each molecule that map to each other. 
				symmetric_molecule_pairs[(m1_name, m2_name)] = [(m1_i, m2_j) for m1_i, m2_j in mol1_to_mol2_conversion_Comp_Only.items()] 

	else:

		# 5.9: Determine if the molecules are variants using the Comprehensive Method. 
		#      * There are not enough different elements to perform the Minimal Elemental Abundance Method, so fall back to the more 
		#        computationally demanding comprehensive method. 
		positions_are_variant_Comp_Only, mol1_to_mol2_conversion_Comp_Only = are_molecules_variant_from_comprehensive(m1_original_elements, m1_original_positions, no_of_H_on_atoms_in_molecule1, m2_original_elements, m2_original_positions, no_of_H_on_atoms_in_molecule2, max_distance_disparity, molecules_being_compared, neighbouring_molecules_about_molecules, non_hydrogen_molecules, m1_non_hydrogen_graph, m2_non_hydrogen_graph, no_of_cpus=no_of_cpus_for_comparing_two_molecules)

		# 5.10: If the positions between the two molecule are variant:
		if positions_are_variant_Comp_Only:
			
			# 5.11: Found that these molecules are variants using the Comprehensive Method. Append this molecule pair to the symmetric_dimer_pairs list, along with the indices of each molecule that map to each other. 
			symmetric_molecule_pairs[(m1_name, m2_name)] = [(m1_i, m2_j) for m1_i, m2_j in mol1_to_mol2_conversion_Comp_Only.items()] 

def get_inputs(unique_molecules_names, non_hydrogen_molecules_elements, non_hydrogen_molecules_positions, all_no_of_H_on_atoms_in_molecule, max_distance_disparity, neighbouring_molecules_about_molecules, non_hydrogen_molecules, non_hydrogen_graphs, symmetric_molecule_pairs, no_of_cpus_for_comparing_two_molecules):
	"""
	This generator is designed to return all the input methods required for the compare_if_two_molecules_are_symmetric_single_process method. 

	Parameters
	----------
	unique_molecules_names : list of ints
		This is a list of unique molecules based on the symmetry of the crystal.
	non_hydrogen_molecules_elements : list
		These are the element of all the molecules that can make up the dimers.
	non_hydrogen_molecules_positions : list
		These are the positions of all the molecules that can make up the dimers.
	all_no_of_H_on_atoms_in_molecule : list
		These are the number of hydrogens bound to each "Heavy" atom in each molecule that can make up the dimers. 
	max_distance_disparity : float
		This is the maximum difference in relative positions of atoms between dimer 1 and dimer 2 to be considered equivalent/identical
	neighbouring_molecules_about_molecules : dict.
		This is the information about the molecules that surround (in the vicinity of) each moleule in the crystal.
	non_hydrogen_molecules : list of ase.Atoms
		These are Atoms of objects of the molecules in the crystal. These do not contain hydrogen atoms. This is for debugging use only. 
	symmetric_molecule_pairs : dict
		This dictionary stores which molecules are symmetric to each other, as well as which names in molecule 1 map onto which names in molecule 2. 
	no_of_cpus_for_comparing_two_molecules : int
		This is the number of cpus reserved for comparing two molecules together.

	Returns
	-------
	m1_name : int
		This is the name of the first molecule we want to compare.
	m1_original_elements : list of str.
		This is the list of elements that make up molecule 1.
	m1_original_positions : 2D np.array
		These are the positions of the atoms in molecule 1
	no_of_H_on_atoms_in_molecule1 : list of int.
		These are the number of hydrogen atoms bound to each "heavy" atom in molecule 1. 

	m2_name : int
		This is the name of the second molecule we want to compare.
	m2_original_elements : list of str.
		This is the list of elements that make up molecule 2.
	m2_original_positions : 2D np.array
		These are the positions of the atoms in molecule 2.
	no_of_H_on_atoms_in_molecule2 : list of int.
		These are the number of hydrogen atoms bound to each "heavy" atom in molecule 2.

	max_distance_disparity: float
		This is the maximum that any two "could be equivalent" atoms can be between molecule 1 and molecule 2 for molecules 1 and 2 to be considered variant.
	neighbouring_molecules_about_molecules : dict.
		This is the information about the molecules that surround (in the vicinity of) each moleule in the crystal.
	non_hydrogen_molecules : list of ase.Atoms
		This is the list of molecules in the crystal, not including hydrogens.

	symmetric_molecule_pairs : dict
		This dictionary stores which molecules are symmetric to each other, as well as which names in molecule 1 map onto which names in molecule 2. 

	no_of_cpus_for_comparing_two_molecules : int
		This is the number of cpus reserved for comparing two molecules together.
	"""

	# First, for each molecule in unique_molecules_names
	for index1 in range(len(unique_molecules_names)):

		# Second, get the name of first molecule
		m1_name                           = unique_molecules_names[index1]

		# Third, obtain the elements, positions, and number of hydrogens bound to "heavy" atoms in molecule 1.
		m1_original_elements              = non_hydrogen_molecules_elements[m1_name]
		m1_original_positions             = non_hydrogen_molecules_positions[m1_name]
		no_of_H_on_atoms_in_molecule1     = all_no_of_H_on_atoms_in_molecule[m1_name]

		# Fourth, obtain the non-hydrogen graph for molecule 1.
		m1_non_hydrogen_graph             = non_hydrogen_graphs[m1_name]

		# Fifth, collect the above details into a single tuple to keep the data together.
		molecule1_details = (m1_name, m1_original_elements, m1_original_positions, no_of_H_on_atoms_in_molecule1)

		# Sixth, for each other molecule in unique_molecules_names
		for index2 in range(index1+1,len(unique_molecules_names)):

			# Seventh, get name of second molecule
			m2_name                       = unique_molecules_names[index2]

			# Eighth, obtain the elements, positions, and number of hydrogens bound to "heavy" atoms in molecule 2.
			m2_original_elements          = non_hydrogen_molecules_elements[m2_name]
			m2_original_positions         = non_hydrogen_molecules_positions[m2_name]
			no_of_H_on_atoms_in_molecule2 = all_no_of_H_on_atoms_in_molecule[m2_name]

			# Ninth, obtain the non-hydrogen graph for molecule 1.
			m2_non_hydrogen_graph         = non_hydrogen_graphs[m2_name]

			# Tenth, collect the above details into a single tuple to keep the data together.
			molecule2_details = (m2_name, m2_original_elements, m2_original_positions, no_of_H_on_atoms_in_molecule2)

			# Eleventh, yield the input variables for the XXX method.
			yield (molecule1_details, molecule2_details, max_distance_disparity, neighbouring_molecules_about_molecules, non_hydrogen_molecules, m1_non_hydrogen_graph, m2_non_hydrogen_graph, symmetric_molecule_pairs, no_of_cpus_for_comparing_two_molecules)

# ---------------------------------------------------------------------------------------------------------------------------------------------------------

