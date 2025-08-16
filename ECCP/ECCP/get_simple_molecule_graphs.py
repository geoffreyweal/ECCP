"""
get_simple_molecule_graphs.py, Geoffrey Weal, 5/2/24

This script is designed to return the original list of molecule graph, but where atoms only include E (element) and no_of_neighbouring_non_cord_H as node features.
"""

from networkx import Graph

def get_simple_molecule_graphs(molecule_graphs):
	"""
	This method is designed to return the original list of molecule graph, but where atoms only include E (element) and no_of_neighbouring_non_cord_H as node features.

	Parameters
	----------
	molecule_graphs : list of networkx.Graph
		These are the graphs that correspond to the molecules in the molecules list

	Returns
	-------
	simple_molecule_graphs : list of networkx.Graph
		This corresponds to molecule_graphs, but only E and no_of_neighbouring_non_cord_H are included as atom (node) features. These graphs contain no bond (edge) features.
	"""

	# First, initialise the simple_molecule_graphs list
	simple_molecule_graphs = {}

	# Second, set a boolean to record if no_of_neighbouring_non_cord_H is given in the molecule graphs. 
	all_does_contain_no_of_neighbouring_non_cord_H_details     = []
	all_does_not_contain_no_of_neighbouring_non_cord_H_details = []

	# Third, for each graph in molecule_graphs:
	for mol_name, molecule_graph in molecule_graphs.items():

		# 3.1: Initialise a new graph
		simple_molecule_graph = Graph()

		# 3.2: Copy the node names from molecule_graph to simple_molecule_graph
		simple_molecule_graph.add_nodes_from(molecule_graph.nodes.keys())

		# 3.3: Add the 'E' and 'no_of_neighbouring_non_cord_H' features from simple_molecule_graph to simple_molecule_graph for each atom (node).
		for node_index, node_features in molecule_graph.nodes.items():

			# 3.3.1: Check that element 'E' is given in molecule_graph node. If it is not, there is a programming issue. 
			if not 'E' in node_features:
				raise Exception('Error: The element has not been included in your molecule graph. This should not happen in ECCP. This is a programming issue. ')

			# 3.3.2: Copy the element from molecule_graph to simple_molecule_graph
			simple_molecule_graph.nodes[node_index]['E']                             = node_features['E']

			# 3.3.3: If 'no_of_neighbouring_non_cord_H' is given, copy it.
			if 'no_of_neighbouring_non_cord_H' in node_features:
				simple_molecule_graph.nodes[node_index]['no_of_neighbouring_non_cord_H'] = node_features['no_of_neighbouring_non_cord_H']

			# 3.3.4: Record if 'no_of_neighbouring_non_cord_H' was in the node_features for this node in molecule_graph
			if 'no_of_neighbouring_non_cord_H' in node_features:
				all_does_contain_no_of_neighbouring_non_cord_H_details.append((molecule_graph.name, node_index))
			else:
				all_does_not_contain_no_of_neighbouring_non_cord_H_details.append((molecule_graph.name, node_index))

		# 3.4: Copy the edge names from molecule_graph to simple_molecule_graph
		simple_molecule_graph.add_edges_from(molecule_graph.edges.keys())

		# 3.5: Check that mol_name is not already in simple_molecule_graphs.keys().
		if mol_name in simple_molecule_graphs.keys():
			raise Exception('Error: mol_name already in simple_molecule_graphs.keys(). This likely indicates a programming.')

		# 3.6: Append simple_molecule_graph to simple_molecule_graphs
		simple_molecule_graphs[mol_name] = simple_molecule_graph

	# Fourth, make sure that there is consistancy across all atoms for all molecules for the 'no_of_neighbouring_non_cord_H' feature.
	#         * We want either all atom to have the 'no_of_neighbouring_non_cord_H' feature, or none of them have this feature. 
	if (len(all_does_contain_no_of_neighbouring_non_cord_H_details) != 0) and (len(all_does_not_contain_no_of_neighbouring_non_cord_H_details) != 0):
		toString  = 'Error: There is inconsistancy for the "no_of_neighbouring_non_cord_H" feature. Some atoms have them, some do not.'+'\n'
		tostring += 'Nodes containing "no_of_neighbouring_non_cord_H":     '+str(all_does_contain_no_of_neighbouring_non_cord_H_details)+'\n'
		tostring += 'Nodes not containing "no_of_neighbouring_non_cord_H": '+str(all_does_not_contain_no_of_neighbouring_non_cord_H_details)+'\n'
		tostring += 'Check this, as this indicates a programming error.'
		raise Exception(tostring)

	# Fifth, return simple_molecule_graphs
	return simple_molecule_graphs

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------

