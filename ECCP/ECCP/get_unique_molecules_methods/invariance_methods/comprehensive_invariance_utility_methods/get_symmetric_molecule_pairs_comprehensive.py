"""
get_symmetric_molecule_pairs_comprehensive.py, Geoffrey Weal, 11/2/24

This script is designed to determine all the spatially symmetric molecules in the unique_molecules_names list. 
"""
from tqdm import tqdm

import multiprocessing as mp

from ECCP.ECCP.get_unique_molecules_methods.invariance_methods.comprehensive_invariance_utility_methods.are_molecules_variant_comprehensive import are_molecules_variant_comprehensive

def get_symmetric_molecule_pairs_comprehensive(unique_molecules_names, non_hydrogen_molecules_elements, non_hydrogen_molecules_positions, all_no_of_H_on_atoms_in_molecule, equivalent_molecule_atom_index_comparisons, max_distance_disparity, neighbouring_molecules_about_molecules, non_hydrogen_molecules, no_of_cpus=1):
	"""
	This method is designed to determine all the spatially symmetric molecules in the unique_molecules_names list. 

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
	no_of_cpus : int.
		This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.

	Returns
	-------
	symmetric_molecule_pairs : dict
		This dictionary stores which molecules are symmetric to each other, as well as which atom indices in molecule 1 map onto which atom indices in molecule 2. 
	"""

	# Prestep: Set the number of cores for each comparison of two molecules
	no_of_cpus_for_comparing_two_molecules = no_of_cpus

	# First, determine the number of molecule pairwise comparisons that will be made. 
	nn = int((len(unique_molecules_names)*(len(unique_molecules_names)-1))/2)

	# Second, obtain the pairs of symmetric molecules in the crystal. 
	if True: #no_of_cpus == 1: # If the user only wants to use 1 cpu, perform tasks without using multiprocessing

		# 2.1: Initialise the symmetric_molecule_pairs dictionary. This will hold all the information about symmetric dimers, as well as which atom indices in molecule 1 map onto which atom indices in molecule 2. 
		symmetric_molecule_pairs = {}

		# 2.2: Create a progress bar for running this task
		pbar = tqdm(get_inputs(unique_molecules_names, non_hydrogen_molecules_elements, non_hydrogen_molecules_positions, all_no_of_H_on_atoms_in_molecule, equivalent_molecule_atom_index_comparisons, max_distance_disparity, neighbouring_molecules_about_molecules, non_hydrogen_molecules, symmetric_molecule_pairs, no_of_cpus=no_of_cpus_for_comparing_two_molecules), total=nn, unit='molecule pair')

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
			print('Comparing molecules (Please wait until after 100%, as the process will still be running.)', file=sys.stderr)

			# 2.8: Perform compare_two_dimers on each way that dimer 2 can be mapped onto dimer 1 using multiprocessing. 
			#process_map(compare_if_two_molecules_are_symmetric_single_process, get_inputs(unique_molecules_names, non_hydrogen_molecules_elements, non_hydrogen_molecules_positions, all_no_of_H_on_atoms_in_molecule, equivalent_molecule_atom_index_comparisons, max_distance_disparity, neighbouring_molecules_about_molecules, non_hydrogen_molecules, symmetric_molecule_pairs), total=nn, unit='molecule pair', desc='Comparing molecules', max_workers=no_of_cpus)
			pool = mp.Pool(processes=no_of_cpus)
			pool.map_async(compare_if_two_molecules_are_symmetric_single_process, tqdm(get_inputs(unique_molecules_names, non_hydrogen_molecules_elements, non_hydrogen_molecules_positions, all_no_of_H_on_atoms_in_molecule, equivalent_molecule_atom_index_comparisons, max_distance_disparity, neighbouring_molecules_about_molecules, non_hydrogen_molecules, symmetric_molecule_pairs, no_of_cpus=1), total=nn, unit='molecule pair', desc='Comparing molecules'))
			pool.close()
			pool.join()

			# 2.9: Convert symmetric_molecule_pairs from a manager.dict object to a dict object.
			symmetric_molecule_pairs = dict(symmetric_molecule_pairs)

	# Third, return information about the symmetric pairs of molecules
	return symmetric_molecule_pairs

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

def compare_if_two_molecules_are_symmetric_single_process(input_data):
	"""
	This method is designed to determine if two molecules are symmetric to each other. 

	Parameters
	----------
	molecule_1_information : list
		This list contains the following information about molecules 1.
			mol_name1 : int
				This is the name of the first molecule we want to compare.
			molecule1_elements : list of str.
				This is the list of elements that make up molecule 1.
			molecule1_positions : 2D np.array
				These are the positions of the atoms in molecule 1
			no_of_H_on_atoms_in_molecule1 : list of int.
				These are the number of hydrogen atoms bound to each "heavy" atom in molecule 1. 

	molecule_2_information : list
		This list contains the following information about molecules 2.
			mol_name2 : int
				This is the name of the second molecule we want to compare.
			molecule2_elements : list of str.
				This is the list of elements that make up molecule 2.
			molecule2_positions : 2D np.array
				These are the positions of the atoms in molecule 2.
			no_of_H_on_atoms_in_molecule2 : list of int.
				These are the number of hydrogen atoms bound to each "heavy" atom in molecule 2.

	em_indices_m2_to_m1 : dict
		This dictionary indicates how atom indices in molecule 2 could be mapped onto molecule 1. 
	max_distance_disparity: float
		This is the maximum that any two "could be equivalent" atoms can be between molecule 1 and molecule 2 for molecules 1 and 2 to be considered variant.
	neighbouring_molecules_about_molecules : dict.
		This is the information about the molecules that surround (in the vicinity of) each moleule in the crystal.
	non_hydrogen_molecules : list of ase.Atoms
		This is the list of molecules in the crystal, not including hydrogens.

	symmetric_molecule_pairs : dict
		This dictionary stores which molecules are symmetric to each other, as well as which atom indices in molecule 1 map onto which atom indices in molecule 2. 

	no_of_cpus : int.
		This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.
	"""

	# First, extract the input variables from input_data
	molecule_1_information, molecule_2_information, em_indices_m2_to_m1, max_distance_disparity, neighbouring_molecules_about_molecules, non_hydrogen_molecules, symmetric_molecule_pairs, no_of_cpus = input_data

	# Second, determine if the two molecules the variants of each other
	#         * If the two molecules are variant, give a list that indicates how the atom indices in molecule 1 relate to the atom indices in molecule 2. 
	is_variant, mol1_to_mol2_conversion_Comp = are_molecules_variant_comprehensive(molecule_1_information, molecule_2_information, em_indices_m2_to_m1, max_distance_disparity, neighbouring_molecules_about_molecules, non_hydrogen_molecules, symmetric_molecule_pairs, no_of_cpus=no_of_cpus)

	# Third, if the molecule is variant, add mol1_to_mol2_conversion_Comp to the symmetric_molecule_pairs dictionary
	if is_variant:

		# 3.1: Get the names of the molecules
		mol_name1 = molecule_1_information[0]
		mol_name2 = molecule_2_information[0]

		# 3.2: Record mol1_to_mol2_conversion_Comp in symmetric_molecule_pairs dictionary with key (mol_name1, mol_name2)
		#import pdb; pdb.set_trace()
		symmetric_molecule_pairs[(mol_name1, mol_name2)] = list(mol1_to_mol2_conversion_Comp)

def get_inputs(unique_molecules_names, non_hydrogen_molecules_elements, non_hydrogen_molecules_positions, all_no_of_H_on_atoms_in_molecule, equivalent_molecule_atom_index_comparisons, max_distance_disparity, neighbouring_molecules_about_molecules, non_hydrogen_molecules, symmetric_molecule_pairs, no_of_cpus):
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
	equivalent_molecule_atom_index_comparisons : dict.
		This is all the ways that two molecules in non_hydrogen_graphs can map onto each other. If you put in (name1, name2), the output will be all the ways that the atom indices in name1 can be mapped (changed) to obtain an equivalent molecule to name2. 
	max_distance_disparity : float
		This is the maximum difference in relative positions of atoms between dimer 1 and dimer 2 to be considered equivalent/identical
	neighbouring_molecules_about_molecules : dict.
		This is the information about the molecules that surround (in the vicinity of) each moleule in the crystal.
	non_hydrogen_molecules : list of ase.Atoms
		These are Atoms of objects of the molecules in the crystal. These do not contain hydrogen atoms. This is for debugging use only. 
	symmetric_molecule_pairs : dict
		This dictionary stores which molecules are symmetric to each other, as well as which indices in molecule 1 map onto which indices in molecule 2. 
	no_of_cpus : int.
		This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.

	Returns
	-------
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
		This dictionary stores which molecules are symmetric to each other, as well as which indices in molecule 1 map onto which indices in molecule 2. 

	no_of_cpus : int.
		This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.
	"""

	# First, for each molecule in unique_molecules_names:
	for index1 in range(len(unique_molecules_names)):

		# Second, get the name of first molecule.
		mol1_name                    = unique_molecules_names[index1]

		# Third, obtain the elements, positions, and number of hydrogens attached to each non-hydrogen atom in the first molecule. 
		molecule1_elements            = non_hydrogen_molecules_elements [mol1_name]
		molecule1_positions           = non_hydrogen_molecules_positions[mol1_name]
		no_of_H_on_atoms_in_molecule1 = all_no_of_H_on_atoms_in_molecule[mol1_name]

		# Fourth, for each other molecule in unique_molecules_names:
		for index2 in range(index1+1, len(unique_molecules_names)):

			# Fifth, get the name of second molecule.
			mol2_name                     = unique_molecules_names[index2]

			# Sixth, obtain the elements, positions, and number of hydrogens attached to each non-hydrogen atom in the second molecule. 
			molecule2_elements            = non_hydrogen_molecules_elements [mol2_name]
			molecule2_positions           = non_hydrogen_molecules_positions[mol2_name]
			no_of_H_on_atoms_in_molecule2 = all_no_of_H_on_atoms_in_molecule[mol2_name]

			# Eighth, get all the indices that are equivalent to eachother in each molecule in each dimer. 
			#         * Important: Here we want to convert mol2 --> mol1 (i.e. map mol2 onto mol1).
			em_indices_m2_to_m1 = equivalent_molecule_atom_index_comparisons[ (mol2_name, mol1_name) ]

			# Ninth, return input variables. 
			yield ((mol1_name, molecule1_elements, molecule1_positions, no_of_H_on_atoms_in_molecule1), (mol2_name, molecule2_elements, molecule2_positions, no_of_H_on_atoms_in_molecule2), em_indices_m2_to_m1, max_distance_disparity, neighbouring_molecules_about_molecules, non_hydrogen_molecules, symmetric_molecule_pairs, no_of_cpus)

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

