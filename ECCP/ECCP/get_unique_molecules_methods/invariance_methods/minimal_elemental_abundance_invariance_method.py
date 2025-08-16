"""
minimal_elemental_abundance_invariance_method.py, Geoffrey Weal, 5/4/22

This script is designed to use the procrustes analysis to determine if dimers are rotationally, translationally, and reflectively invariant.
"""
from ECCP.ECCP.invariance_methods.extract_non_hydrogen_lists_from_molecule_graphs                                                                      import extract_non_hydrogen_lists_from_molecule_graphs
from ECCP.ECCP.get_unique_molecules_methods.invariance_methods.minimal_elemental_abundance_invariance_utility_methods.get_symmetric_molecule_pairs_MEA import get_symmetric_molecule_pairs_MEA

def remove_equivalent_molecules_minimal_elemental_abundance_invariance_method(unique_molecules_indices, molecules, molecule_graphs, neighbouring_molecules_about_molecules={}, max_distance_disparity=0.01, include_hydrogens_in_uniqueness_analysis=False, no_of_cpus=1):
	"""
	This method looks for the elements in lowest abundance and aligns them between molecules to determine if they are the same or not.

	This method works best if your molecule contains one or a few elements that are in low abundance in the molecule.

	Parameters
	----------
	unique_molecules_indices : list of ints
		This is a list of unique molecules based on the symmetry of the crystal.
	molecules : list of ase.Atoms
		This is the list of molecules in the crystal.
	molecule_graphs : list of networkx.Graph
		This is the list of the undirected graph representations of each molecule.
	neighbouring_molecules_about_molecules : dict.
		This is the information about the molecules that surround (in the vicinity of) each moleule in the crystal.
	max_distance_disparity : float
		This is the maximum disparity between two molecules to be considered invariant. 
	include_hydrogens_in_uniqueness_analysis: bool. 
		This tag indicates if you want to include hydrogens when accessing uniquess between molecules and dimers. Default: False
	no_of_cpus : int.
		This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.

	Returns
	-------
	symmetric_molecules_indices : list
		A list of the indices of equivalent molecules in the molecules list. 
	"""

	# First, before beginning, if not value for max_distance_disparity is given, set to a default value.
	if max_distance_disparity is None:
		max_distance_disparity = 0.01

	# Second, remove hydrogens from molecules
	non_hydrogen_molecules, non_hydrogen_molecules_elements, non_hydrogen_molecules_positions, all_no_of_H_on_atoms_in_molecule, non_hydrogen_graphs = extract_non_hydrogen_lists_from_molecule_graphs(unique_molecules_indices, molecules, molecule_graphs, include_hydrogens_in_uniqueness_analysis=False)

	# Third, determine which molecules are invariant.
	print('Comparing translational, rotational, and reflective invarience between molecules. '+str(len(unique_molecules_indices))+' molecules to be examined. This can take a while with large and complex molecules.')
	symmetric_molecule_pairs = get_symmetric_molecule_pairs_MEA(unique_molecules_indices, non_hydrogen_molecules_elements, non_hydrogen_molecules_positions, all_no_of_H_on_atoms_in_molecule, max_distance_disparity, neighbouring_molecules_about_molecules, non_hydrogen_molecules, non_hydrogen_graphs, no_of_cpus=no_of_cpus)

	# Fourth, return the indices of all the equivalent molecules in the crystal.
	return symmetric_molecule_pairs


