"""
make_human_friendly_crystal.py, Geoffrey Weal, 20/2/22

This script is designed to create a super crystal that makes it easier for a human to view.
"""
import numpy    as np
import networkx as nx
from copy     import deepcopy
from ase      import Atoms
from SUMELF   import get_cell_corner_points
from networkx import Graph, relabel_nodes, compose

def make_human_friendly_crystal(molecules, super_cell_reach=1, molecule_graphs=None):
	"""
	This method will create a super crystal that makes it easier for a human to view. 

	Parameters
	----------
	molecules : list of ase.Atoms objects.
		These are all the individual molecules identified in the crystal.
	super_cell_reach : int
		This is the size of the super crystal. Default: 1. 
	molecule_graphs : list of networkx.Graph
		This is a list of all the graph of the molecules in the molecules list.

	Returns
	-------
	human_friendly_crystal : ase.Atoms
		This is the human-friendly super crystal. 
	human_friendly_crystal_graph : networkx.Graph
		This is the graph of the human-friendly super crystal. 
	"""

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# First, get one of the names in the molecules dictionary. 
	first_molecule_name = sorted(molecules)[0]

	# Second, get the cell points of the crystal or super crystal.
	crystal_cell_lattice = molecules[first_molecule_name].get_cell()
	cell_points = get_cell_corner_points(crystal_cell_lattice, super_cell_reach=super_cell_reach, bottom_of_range=0)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# Third, perform some checks

	# 3.1: Make sure that the (0., 0., 0.) displacement is first in cell_points
	if not all(cell_points[0] == np.array([0.0,0.0,0.0])):
		to_string  = 'Error: The origin displacement has not been given as the first displacement vector.\n'
		to_string += 'cell_points: '+str(cell_points)+'\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Fourth, obtain the maximum mol name for the molecules.
	max_mol_name = max(molecules.keys())

	# Fifth, initialise the dict. for holding the molecules in the desired order.
	molecules_in_order_for_human_friendly_crystal = {}

	# Sixth, initialise the dict. for holding the molecule graphs in the desired order.
	if molecule_graphs is not None:
		molecule_graphs_in_order_for_human_friendly_crystal = {}

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# Seventh, obtain the molecules that have also been displaced by the displacements in cell_points.

	# 7.1: For each displacement given in the cell_points list:
	for cp_counter, displacement in enumerate(cell_points):

		# 7.2: For each molecule in the molecules list:
		for mol_name, molecule in sorted(molecules.items(), key=lambda x: x[0]):

			# 7.3: Make sure that mol_name is greater than 1, otherwise there will be a problem
			if not mol_name >= 1:
				to_string = f'Error: mol_name is less than 1 (It is {mol_name}).\nCheck this'
				raise Exception(to_string)

			# 7.4: Obtain the name to give this molecule in the crystal
			mol_name_in_crystal = mol_name + (cp_counter*max_mol_name)

			# 7.5: Make a copy of the molecule
			molecule_disp = molecule.copy()

			# 7.6: Apply the displacement vector for each atom in the molecule.
			molecule_disp.set_positions(molecule_disp.get_positions() + displacement)

			# 7.7: Make sure that mol_name_in_crystal is not already in molecules_in_order_for_human_friendly_crystal
			if mol_name_in_crystal in molecules_in_order_for_human_friendly_crystal:
				raise Exception('Error: mol_name_in_crystal already in molecules_in_order_for_human_friendly_crystal. Probably programming issue.')

			# 7.8: Append molecule_disp to molecules_in_order_for_human_friendly_crystal
			molecules_in_order_for_human_friendly_crystal[mol_name_in_crystal] = molecule_disp

			# 7.9: If the graph for the human_friendly_crystal is being recorded
			if molecule_graphs is not None:

				# 7.10: Make a copy of the graph
				molecule_graph_to_append = deepcopy(molecule_graphs[mol_name])

				# 7.11: Make sure that mol_name_in_crystal is not already in molecules_in_order_for_human_friendly_crystal
				if mol_name_in_crystal in molecule_graphs_in_order_for_human_friendly_crystal:
					raise Exception('Error: mol_name_in_crystal already in molecule_graphs_in_order_for_human_friendly_crystal. Probably programming issue.')

				# 7.12: Add the graph to human_friendly_crystal_graph
				molecule_graphs_in_order_for_human_friendly_crystal[mol_name_in_crystal] = molecule_graph_to_append

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Eighth, initiate counter for recording molecule index and atom count.
	atom_number_counter = 0

	# Ninth, initiate the ase.Atoms object for the human-friendly crystal.
	human_friendly_crystal = Atoms(cell=molecules[first_molecule_name].get_cell()); 

	# Tenth, initiate the networkx.Graph object for the human-friendly crystal (if molecule_graphs is given).
	if molecule_graphs is not None:
		human_friendly_crystal_graph = Graph()

	# Eleventh, for each molecule in the molecules_in_order_for_human_friendly_crystal dict
	for mol_name in sorted(molecules_in_order_for_human_friendly_crystal.keys()):

		# 11.1: Obtain the molecule to append to the human_friendly_crystal
		molecule_to_append = molecules_in_order_for_human_friendly_crystal[mol_name]

		# 11.2: Append molecule_to_append to human_friendly_crystal.
		human_friendly_crystal += molecule_to_append.copy()

		# 11.4: If the graph for the human_friendly_crystal is being recorded
		if molecule_graphs is not None:

			# 11.3: Obtain the graph for this molecule. 
			molecule_graph_to_append = deepcopy(molecule_graphs_in_order_for_human_friendly_crystal[mol_name])

			# 11.4: Add the MoleculeList for this individual molecule into molecule_graph_to_append
			for node in molecule_graph_to_append.nodes.values():
				node['MoleculeList'] = mol_name

			# 11.5: Create a relabelling map to relabel the atoms in this molecule for the overall crystal structure.
			mapping = {index: index+atom_number_counter for index in molecule_graph_to_append.nodes}

			# 11.6: Update the names of the nodes in the molecule graph so the indices of atoms in the molecule are those of the same atoms indices in the crystal
			molecule_graph_to_append = relabel_nodes(molecule_graph_to_append, mapping)

			# 11.7: Add the molecule_graph_to_append to the currently being built human_friendly_crystal_graph
			human_friendly_crystal_graph = compose(human_friendly_crystal_graph, molecule_graph_to_append)

		# 11.8: Increase atom_number_counter by the number of atoms in the molecule. 
		atom_number_counter += len(molecule_to_append)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Twelfth, return the human friendly version of the crystal.
	if molecule_graphs is not None:
		return human_friendly_crystal, human_friendly_crystal_graph
	else:
		return human_friendly_crystal



