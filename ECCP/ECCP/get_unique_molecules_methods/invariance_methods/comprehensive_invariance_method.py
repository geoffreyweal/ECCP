"""
invariance_method.py, Geoffrey Weal, 23/2/22

This script is designed to use the procrustes analysis to determine if molecules are rotationally, translationally, and reflectively invariant.
"""
from ECCP.ECCP.invariance_methods.extract_non_hydrogen_lists_from_molecule_graphs                                                                  import extract_non_hydrogen_lists_from_molecule_graphs
from ECCP.ECCP.invariance_methods.common_comprehensive_invariance_utility_methods.get_equivalent_molecule_names                                    import get_equivalent_molecule_names
from ECCP.ECCP.get_unique_molecules_methods.invariance_methods.comprehensive_invariance_utility_methods.get_symmetric_molecule_pairs_comprehensive import get_symmetric_molecule_pairs_comprehensive

def remove_equivalent_molecules_comprehensive_invariance_method(unique_molecule_names, molecules, molecule_graphs, neighbouring_molecules_about_molecules={}, max_distance_disparity=0.01, include_hydrogens_in_uniqueness_analysis=False, no_of_cpus=1):
	"""
	This method uses the procrustes analysis to determine if molecules are rotationally, translationally, and reflectively invariant.

	From using the procrustes analysis, this method will determine which molecules are equivalent.

	This analysis will ignore hydrogen atoms. 

	The procrustes analysis use is from scipy. See the website below for more information: 
		* https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.procrustes.html

	Parameters
	----------
	unique_molecule_names : list of ints
		This is a list of unique molecules based on the symmetry of the crystal.

	molecules : list of ase.Atoms
		This is the list of molecules in the crystal.
	molecule_graphs : list of networkx.Graph
		This is the list of the undirected graph representations of each molecule.

	neighbouring_molecules_about_molecules : dict.
		This is the information about the molecules that surround (in the vicinity of) each moleule in the crystal.

	max_distance_disparity : float
		This is the maximum disparity between two molecules to be considered invariant. If max_distance_disparity is given as None, the default value will be given. Default: 0.01 Å.
	include_hydrogens_in_uniqueness_analysis: bool. 
		This tag indicates if you want to include hydrogens when accessing uniquess between molecules and dimers. Default: False
	
	no_of_cpus : int.
		This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.

	Returns
	-------
	symmetric_molecules : list
		A list of the indices of symmetric molecules in the molecules list. 
	"""

	# First, if not value for max_distance_disparity is given, set to a default value.
	if max_distance_disparity is None:
		max_distance_disparity = 0.01

	# Second, obtain all molecules without hydrogen.
	non_hydrogen_molecules, non_hydrogen_molecules_elements, non_hydrogen_molecules_positions, all_no_of_H_on_atoms_in_molecule, non_hydrogen_graphs = extract_non_hydrogen_lists_from_molecule_graphs(unique_molecule_names, molecules, molecule_graphs, include_hydrogens_in_uniqueness_analysis=False)

	# Third, determine how the indices of molecules can be interchanged to give the same molecule.
	print('Examining equivalent atoms between the '+str(len(non_hydrogen_graphs))+' molecules identified. This can take a while with large and complex molecules.')
	equivalent_molecule_indices = get_equivalent_molecule_names(non_hydrogen_graphs, include_comparisons_with_itself=False, no_of_cpus=no_of_cpus)

	# Fourth, determine which pairs of molecules are symmetric.
	print('Comparing translational, rotational, and reflective invarience between molecules. '+str(len(unique_molecule_names))+' molecules to be examined. This can take a while with large and complex molecules.')
	symmetric_molecule_pairs = get_symmetric_molecule_pairs_comprehensive(unique_molecule_names, non_hydrogen_molecules_elements, non_hydrogen_molecules_positions, all_no_of_H_on_atoms_in_molecule, equivalent_molecule_indices, max_distance_disparity, neighbouring_molecules_about_molecules, non_hydrogen_molecules, no_of_cpus=no_of_cpus)

	# Fifth, return the indices of all the equivalent molecules in the crystal.
	return symmetric_molecule_pairs
