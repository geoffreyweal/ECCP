"""
extract_non_hydrogen_lists_from_molecule_graphs.py, Geoffrey Weal, 11/2/14

This script is designed to remove the hydrogens from the molecules and their respective graphs, as well as to extract the graphs into lists of useful information. 
"""
from copy        import deepcopy
from collections import Counter
from SUMELF      import remove_hydrogens

def extract_non_hydrogen_lists_from_molecule_graphs(unique_molecule_names, molecules, molecule_graphs, include_hydrogens_in_uniqueness_analysis=False):
	"""
	This method is designed to remove the hydrogens from the molecules and their respective graphs, as well as to extract the graphs into lists of useful information. 

	Parameters
	----------
	unique_molecule_names : list of ints
		This is a list of unique molecules based on the symmetry of the crystal.
	molecules : list of ase.Atoms
		This is the list of molecules in the crystal.
	molecule_graphs : list of networkx.Graph
		This is the list of the undirected graph representations of each molecule.
	include_hydrogens_in_uniqueness_analysis: bool. 
		This tag indicates if you want to include hydrogens when accessing uniquess between molecules and dimers. Default: False

	Returns
	-------
	non_hydrogen_molecules : dict of ase.Atoms
		These are the original ase.Atoms object with their hydrogens removed. This is for debugging. 
	non_hydrogen_molecules_elements : dict of str. 
		These are the list of the elements of the atoms for each molecule in non_hydrogen_molecules
	non_hydrogen_molecules_positions : dict of 2D np.array
		These are the list of the positions of the atoms for each molecule in non_hydrogen_molecules
	all_no_of_H_on_atoms_in_molecule : dict of int
		This is a list of all the hydrogens that are bound to each "heavy" atom in the molecule
	non_hydrogen_graphs : dict of networkx.Graphs
		This is a list of all the graph of all moelcules with hydrogens removed as nodes and attached to its respective "heavy" atom as a node variable. 
	"""

	# First, check that there are no duplicate names in unique_molecule_names
	if not (len(unique_molecule_names) == len(set(unique_molecule_names))):
		to_string  = 'Error: There are duplicate molecules in the unique_molecule_names list.\n'
		duplicates = [mol_name for mol_name, count in Counter(unique_molecule_names).items() if (count >= 2)]
		to_string += f'Duplicate entries in unique_molecule_names: {duplicates}\n'
		to_string += f'unique_molecule_names = {unique_molecule_names}\n'
		to_string += 'Check this'
		raise Exception(to_string)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Second, initialise list for storing information about molecules
	non_hydrogen_molecules           = {} # for debugging
	non_hydrogen_molecules_elements  = {}
	non_hydrogen_molecules_positions = {}
	all_no_of_H_on_atoms_in_molecule = {}
	non_hydrogen_graphs              = {}

	# Third, extract information from molecules and their graphs.
	for mol_name in unique_molecule_names:

		# 3.1: Remove hydrogens from molecules if desired
		if include_hydrogens_in_uniqueness_analysis:
			non_hydrogen_molecule, non_hydrogen_graph, no_of_H_on_atoms_in_molecule = molecules[mol_name].copy(), deepcopy(molecule_graphs[mol_name]), []
		else:
			non_hydrogen_molecule, non_hydrogen_graph, no_of_H_on_atoms_in_molecule = remove_hydrogens(molecules[mol_name], graph=molecule_graphs[mol_name], give_no_H_atoms_on_each_atom=True)

		# 3.2: Check that mol_name is not already in non_hydrogen_molecules, non_hydrogen_molecules_elements, non_hydrogen_molecules_positions, all_no_of_H_on_atoms_in_molecule and non_hydrogen_graphs:
		if mol_name in non_hydrogen_molecules:
			raise Exception('Error: '+str(mol_name)+' already in non_hydrogen_molecules. This indicates a programming error. Check this.')
		if mol_name in non_hydrogen_molecules_elements:
			raise Exception('Error: '+str(mol_name)+' already in non_hydrogen_molecules_elements. This indicates a programming error. Check this.')
		if mol_name in non_hydrogen_molecules_positions:
			raise Exception('Error: '+str(mol_name)+' already in non_hydrogen_molecules_positions. This indicates a programming error. Check this.')
		if mol_name in all_no_of_H_on_atoms_in_molecule:
			raise Exception('Error: '+str(mol_name)+' already in all_no_of_H_on_atoms_in_molecule. This indicates a programming error. Check this.')
		if mol_name in non_hydrogen_graphs:
			raise Exception('Error: '+str(mol_name)+' already in non_hydrogen_graphs. This indicates a programming error. Check this.')

		# 3.3: Add data from non_hydrogen_molecule, non_hydrogen_graph, no_of_H_on_atoms_in_molecule to dicts
		non_hydrogen_molecules[mol_name]           = non_hydrogen_molecule
		non_hydrogen_molecules_elements[mol_name]  = non_hydrogen_molecule.get_chemical_symbols()
		non_hydrogen_molecules_positions[mol_name] = non_hydrogen_molecule.get_positions()
		all_no_of_H_on_atoms_in_molecule[mol_name] = list(no_of_H_on_atoms_in_molecule)
		non_hydrogen_graphs[mol_name]              = non_hydrogen_graph

	# Fourth, return extracted data.
	return non_hydrogen_molecules, non_hydrogen_molecules_elements, non_hydrogen_molecules_positions, all_no_of_H_on_atoms_in_molecule, non_hydrogen_graphs
