"""
invariance_method.py, Geoffrey Weal, 23/2/22

This script is designed to determine if two dimers are varients of each other or are invariant.
"""

from ECCP.ECCP.get_unique_dimers_methods.invariance_methods.comprehensive_invariance_method               import remove_equivalent_dimers_comprehensive_invariance_method
from ECCP.ECCP.get_unique_dimers_methods.invariance_methods.minimal_elemental_abundance_invariance_method import remove_equivalent_dimers_minimal_elemental_abundance_invariance_method
from ECCP.ECCP.get_unique_dimers_methods.invariance_methods.combination_invariance_method                 import remove_equivalent_dimers_combination_invariance_method

def remove_equivalent_dimers_invariance_method(invariance_method_type, dimers, molecules, molecule_graphs, neighbouring_molecules_about_dimers={}, max_distance_disparity=None, include_hydrogens_in_uniqueness_analysis=False, no_of_cpus=1):
	"""
	This method uses the procrustes analysis to determine if dimers are rotationally, translationally, and reflectively invariant.

		* https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.procrustes.html

	Parameters
	----------
	invariance_method_type : str.
		This is the type of invariance method you would like to use. Options given in the instruction manual on the Github page.
	dimers : list
		A list of dimers, consisting of the two molecules involved in the dimer.
	molecules : list of ase.Atoms
		This is the list of molecules that can be used to make dimers
	molecule_graphs : list of networkx.Graph
		This is the list of the undirected graph representations of each molecule
	neighbouring_molecules_about_dimers : dict.
		This is the information about the molecules that surround (in the vicinity of) each dimer in the crystal.
	max_distance_disparity : float
		This is the maximum disparity between two dimers to be considered invariant. 
	include_hydrogens_in_uniqueness_analysis: bool. 
		This tag indicates if you want to include hydrogens when accessing uniquess between molecules and dimers. Default: False
	no_of_cpus : int.
		This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.

	Returns
	-------
	symmetric_dimer_pairs : list
		A list of the indices of equivalent dimers in the dimers list. 
	"""

	# First, choose the invarience method you would like to use.
	if invariance_method_type == 'comprehensive':
		perform_invariance_method = remove_equivalent_dimers_comprehensive_invariance_method
	elif invariance_method_type == 'minimal_elemental_abundance':
		perform_invariance_method = remove_equivalent_dimers_minimal_elemental_abundance_invariance_method
	elif invariance_method_type == 'combination':
		perform_invariance_method = remove_equivalent_dimers_combination_invariance_method
	else:
		print('Error: The input method for the invariance method can be either (given as invariance_method_type):')
		print('\t* comprehensive: The Comprehensive Invariance Method')
		print('\t* minimal_elemental_abundance: The Minimal Elemental Abundance Invariance Method')
		print('\t* combination: The Combination Invariance Method')
		print('See https://github.com/geoffreyweal/ECCP for more information')
		exit('This program will finish without completing')

	# Second, perform the invarance method.
	symmetric_dimer_pairs = perform_invariance_method(dimers, molecules, molecule_graphs, neighbouring_molecules_about_dimers=neighbouring_molecules_about_dimers, max_distance_disparity=max_distance_disparity, include_hydrogens_in_uniqueness_analysis=include_hydrogens_in_uniqueness_analysis, no_of_cpus=no_of_cpus)

	# Third, return symmetric dimers.
	return symmetric_dimer_pairs

