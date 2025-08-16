"""
combination_invariance_method.py, Geoffrey Weal, 5/4/22

This script is designed to use the procrustes analysis to determine if dimers are rotationally, translationally, and reflectively invariant.
"""

from collections import Counter

from ECCP.ECCP.get_unique_molecules_methods.invariance_methods.comprehensive_invariance_method               import remove_equivalent_molecules_comprehensive_invariance_method
from ECCP.ECCP.get_unique_molecules_methods.invariance_methods.minimal_elemental_abundance_invariance_method import remove_equivalent_molecules_minimal_elemental_abundance_invariance_method

def remove_equivalent_molecules_combination_invariance_method(unique_molecule_names, molecules, molecule_graphs, neighbouring_molecules_about_molecules={}, max_distance_disparity=0.01, include_hydrogens_in_uniqueness_analysis=False, no_of_cpus=1):
	"""
	This method uses the procrustes analysis to determine if molecules are rotationally, translationally, and reflectively invariant.

	From using the procrustes analysis, this method will determine which molecules are equivalent.

	This analysis will ignore hydrogen atoms. 

	The procrustes analysis use is from scipy. See the website below for more information: 
		* https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.procrustes.html

	Parameters
	----------
	unique_molecule_names : list of ints
		This is a list of unique molecules based on the symmetry of the crystal
	molecules : list of ase.Atoms
		This is the list of molecules in the crystal
	molecule_graphs : list of networkx.Graph
		This is the list of the undirected graph representations of each molecule
	neighbouring_molecules_about_molecules : dict.
		This is the information about the molecules that surround (in the vicinity of) each moleule in the crystal.
	max_distance_disparity : float
		This is the maximum disparity between two molecules to be considered invariant. Default: 0.01 A
	include_hydrogens_in_uniqueness_analysis: bool. 
		This tag indicates if you want to include hydrogens when accessing uniquess between molecules and dimers. Default: False
	no_of_cpus : int.
		This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.

	Returns
	-------
	symmetric_molecule_pairs : dictionary of (int, int): [mapping of atoms]
		This dictionary contains which molecules are the same as each other (tuples of ints), along with a list of the atom indices that match each other between the molecules. 
	"""

	# -------------------------------
	# First, we want to determine what the total number of the lowest abundant elements are.
	# This is because if these are less than about 4 to 6, then we can take advantage that the dimers have 
	# a few elements in low abundance that we can allign dimers onto each other.
	# This is a very fast method in this case.
	#
	# If this is the case, for example a large dimer that contains just carbons and hydrogens, this method
	# will likely take a long time to do and it might be better to perform the comprehensive invariance method
	# Which has been optimised as best as possible to deal with general cases.

	# ----------------------------------------------------------------------------------------------------------------------------
	# First, determine the total number of lowest abundance elements in each molecule of the dimer
	total_number_of_lowest_abundance_elements_in_molecules = []

	for mol_name, molecule in molecules.items():

		# 1.1: count the number of (abundance of) elements, except for hydrogens which we will exclude.
		number_of_elements = Counter(molecule.get_chemical_symbols())
		if 'H' in number_of_elements:
			del number_of_elements['H']

		# 1.2: Determine the total number of lowest abundance elements in each molecule of the dimer
		total_number_of_lowest_abundance_elements = 0
		for element, no_of_that_element in sorted(number_of_elements.items(), key=lambda x: x[1]):
			total_number_of_lowest_abundance_elements += no_of_that_element
			if total_number_of_lowest_abundance_elements >= 2:
				break
		total_number_of_lowest_abundance_elements_in_molecules.append(total_number_of_lowest_abundance_elements)
	
	# ----------------------------------------------------------------------------------------------------------------------------
	# Second, take the worst case scenario, which is the molecule with the greatest total_number_of_lowest_abundance_elements
	greatest_total_number_of_lowest_abundance_elements_in_molecules = max(total_number_of_lowest_abundance_elements_in_molecules)

	# ----------------------------------------------------------------------------------------------------------------------------
	# Third, based on the worst case scenario, select the invarience method you will use

	if greatest_total_number_of_lowest_abundance_elements_in_molecules <= 4: # Continue with quick invariance method
		print('Determining unique molecules using the minimal elemental abundance invarience method')
		symmetric_molecule_pairs = remove_equivalent_molecules_minimal_elemental_abundance_invariance_method(unique_molecule_names, molecules, molecule_graphs, neighbouring_molecules_about_molecules=neighbouring_molecules_about_molecules, max_distance_disparity=max_distance_disparity, include_hydrogens_in_uniqueness_analysis=include_hydrogens_in_uniqueness_analysis, no_of_cpus=no_of_cpus)
	else: # revert back to original comprehensive invariance method
		print('Determining unique molecules using the comprehensive invarience method')
		symmetric_molecule_pairs = remove_equivalent_molecules_comprehensive_invariance_method              (unique_molecule_names, molecules, molecule_graphs, neighbouring_molecules_about_molecules=neighbouring_molecules_about_molecules, max_distance_disparity=max_distance_disparity, include_hydrogens_in_uniqueness_analysis=include_hydrogens_in_uniqueness_analysis, no_of_cpus=no_of_cpus)
	# ----------------------------------------------------------------------------------------------------------------------------

	# Fourth, all symmetric molecule have been obtained and stored in symmetric_molecule_pairs
	return symmetric_molecule_pairs







