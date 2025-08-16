"""
minimal_elemental_abundance_invariance_method.py, Geoffrey Weal, 5/4/22

This method uses the procrustes analysis to determine if dimers are rotationally, translationally, and reflectively invariant.

This method will focus on those elements that are in least abundances and use them to allign the dimers on top of each other.
"""
from ECCP.ECCP.invariance_methods.extract_non_hydrogen_lists_from_molecule_graphs                                                                import extract_non_hydrogen_lists_from_molecule_graphs
from ECCP.ECCP.get_unique_dimers_methods.invariance_methods.minimal_elemental_abundance_invariance_utility_methods.get_symmetric_dimer_pairs_MEA import get_symmetric_dimer_pairs_MEA

def remove_equivalent_dimers_minimal_elemental_abundance_invariance_method(dimers, molecules, molecule_graphs, neighbouring_molecules_about_dimers={}, max_distance_disparity=None, include_hydrogens_in_uniqueness_analysis=False, no_of_cpus=1):
	"""
	This method uses the procrustes analysis to determine if dimers are rotationally, translationally, and reflectively invariant.

	This method will focus on those elements that are in least abundances and use them to allign the dimers on top of each other.

	Parameters
	----------
	dimers : list
		A list of dimers, consisting of the two molecules involved in the dimer.
	molecules : list of ase.Atoms
		This is the list of molecules that can be used to make dimers.
	molecule_graphs : list of networkx.Graph
		This is the list of the undirected graph representations of each molecule.
	neighbouring_molecules_about_dimers : dict.
		This is the information about the molecules that surround (in the vicinity of) each dimer in the crystal.
	max_distance_disparity : float
		This is the maximum disparity between two dimers to be considered invariant. If None given, set max_distance_disparity = 0.01. Default: None 
	include_hydrogens_in_uniqueness_analysis: bool. 
		This tag indicates if you want to include hydrogens when accessing uniquess between molecules and dimers. Default: False
	no_of_cpus : int.
		This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.

	Returns
	-------
	symmetric_dimer_pairs : list
		A list of the indices of equivalent dimers in the dimers list. 
	"""

	# First, before beginning, if not value for max_distance_disparity is given, set to a default value.
	if max_distance_disparity is None:
		max_distance_disparity = 0.01

	# Second, obtain all molecules without hydrogen
	molecules_names = list(molecules.keys())
	non_hydrogen_molecules, non_hydrogen_molecules_elements, non_hydrogen_molecules_positions, all_no_of_H_on_atoms_in_molecule, non_hydrogen_graphs = extract_non_hydrogen_lists_from_molecule_graphs(molecules_names, molecules, molecule_graphs, include_hydrogens_in_uniqueness_analysis=include_hydrogens_in_uniqueness_analysis)

	# Third, determine which pairs of dimers are symmetric.
	symmetric_dimer_pairs = get_symmetric_dimer_pairs_MEA(dimers, non_hydrogen_molecules_elements, non_hydrogen_molecules_positions, all_no_of_H_on_atoms_in_molecule, max_distance_disparity, neighbouring_molecules_about_dimers, non_hydrogen_molecules, non_hydrogen_graphs, no_of_cpus=no_of_cpus)

	# Fourth, return information about symmetric dimer.
	return symmetric_dimer_pairs
