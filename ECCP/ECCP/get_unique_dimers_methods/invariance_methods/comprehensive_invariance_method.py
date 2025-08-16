"""
comprehensive_invariance_method.py, Geoffrey Weal, 23/2/22

This script is designed to use the procrustes analysis to determine if dimers are rotationally, translationally, and reflectively invariant.
"""
from ECCP.ECCP.invariance_methods.extract_non_hydrogen_lists_from_molecule_graphs                                                               import extract_non_hydrogen_lists_from_molecule_graphs
from ECCP.ECCP.invariance_methods.common_comprehensive_invariance_utility_methods.get_equivalent_molecule_names                                 import get_equivalent_molecule_names
from ECCP.ECCP.get_unique_dimers_methods.invariance_methods.comprehensive_invariance_utility_methods.determine_number_of_permutations_of_dimers import determine_number_of_permutations_of_dimers
from ECCP.ECCP.get_unique_dimers_methods.invariance_methods.comprehensive_invariance_utility_methods.get_symmetric_dimer_pairs_comprehensive    import get_symmetric_dimer_pairs_comprehensive

def remove_equivalent_dimers_comprehensive_invariance_method(dimers, molecules, molecule_graphs, neighbouring_molecules_about_dimers={}, max_distance_disparity=None, include_hydrogens_in_uniqueness_analysis=False, no_of_cpus=1):
	"""
	This method uses the procrustes analysis to determine if dimers are rotationally, translationally, and reflectively invariant.

	From using the procrustes analysis, this method will determine which dimers are equivalent

	This analysis will ignore hydrogen atoms. 

	The procrustes analysis use is from scipy. See the website below for more information: 
		* https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.procrustes.html

	Parameters
	----------
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
	symmetric_dimers : list
		A list of the indices of equivalent dimers in the dimers list. 
	"""

	# First, before beginning, if not value for max_distance_disparity is given, set to a default value.
	if max_distance_disparity is None:
		max_distance_disparity = 0.01

	# Second, obtain all molecules without hydrogen.
	molecules_names = list(molecules.keys())
	non_hydrogen_molecules, non_hydrogen_molecules_elements, non_hydrogen_molecules_positions, all_no_of_H_on_atoms_in_molecule, non_hydrogen_graphs = extract_non_hydrogen_lists_from_molecule_graphs(molecules_names, molecules, molecule_graphs, include_hydrogens_in_uniqueness_analysis=include_hydrogens_in_uniqueness_analysis)

	# Third, determine how the names of molecules can be interchanged to give the same molecule.
	print('Comprehensive Dimer Method: Examining equivalent atoms between the '+str(len(non_hydrogen_graphs))+' molecules identified. This can take a while with large and complex molecules and molecules with several branches.')
	equivalent_molecule_names = get_equivalent_molecule_names(non_hydrogen_graphs, include_comparisons_with_itself=True, no_of_cpus=no_of_cpus)

	# Fourth, this will indicate if this method will probably take too long to perform. (May need to check to see if this is giving good predictions).
	determine_number_of_permutations_of_dimers(equivalent_molecule_names)

	# Fifth, determine which pairs of dimers are symmetric.
	symmetric_dimer_pairs = get_symmetric_dimer_pairs_comprehensive(dimers, non_hydrogen_molecules_elements, non_hydrogen_molecules_positions, all_no_of_H_on_atoms_in_molecule, equivalent_molecule_names, max_distance_disparity, neighbouring_molecules_about_dimers, non_hydrogen_molecules, no_of_cpus=no_of_cpus)

	# Sixth, all symmetric dimers have been obtained and stored in symmetric_dimers.
	return symmetric_dimer_pairs

# ------------------------------------------------------------------------------------------------------------------------------------------
