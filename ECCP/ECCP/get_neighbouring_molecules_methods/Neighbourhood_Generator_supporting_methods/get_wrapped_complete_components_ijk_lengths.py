"""
get_wrapped_complete_components_ijk_lengths.py, Geoffrey Weal, 20/5/24

This method will wrap your molecule, and provide the wrapped graph. It will also provide the ijk positions that the molecules spill into from the origin cell. 
"""
from copy     import deepcopy
from networkx import connected_components
#from SUMELF import get_cell_corner_points
from SUMELF import get_distance, less_than_or_equal_to_max_bondlength, convert_position_into_ijk_lengths

def get_wrapped_complete_components_ijk_lengths(molecule, molecule_graph, crystal_cell_lattice):
	"""
	This method will wrap your molecule, and provide the wrapped graph.

	It will also provide the ijk positions that the molecules spill into from the origin cell. 

	Parameters
	----------
	molecule : ase.Atoms objects
		This is the molecule
	molecule_graphs : networkx.Graph
		This is the graph of the molecule. 
	crystal_cell_lattice : numpy.array or ase.Cell
		This is the unit cell lattice cell dimensions for the crystal

	Returns
	-------
	unit_cell_displacements : list of (int, int, int)
		These are the ijk coordinates that the molecules occupies if it spills over from the origin cell (when unwrapped).
	"""

	# First, make copies of the molecule and the graph, and wrap the molecule.
	wrapped_molecule = molecule.copy()
	wrapped_molecule.wrap()
	wrapped_molecule_graph = deepcopy(molecule_graph)

	# Second, determine the bonds that have been broken in the molecule
	disconnected_bonds = []
	wrapped_molecule_elements  = wrapped_molecule.get_chemical_symbols()
	wrapped_molecule_positions = wrapped_molecule.get_positions()
	for atom_index in wrapped_molecule_graph.nodes:
		element1  = wrapped_molecule_elements[atom_index]
		position1 = wrapped_molecule_positions[atom_index]
		for neighbouring_index in wrapped_molecule_graph[atom_index]:
			element2  = wrapped_molecule_elements[neighbouring_index]
			position2 = wrapped_molecule_positions[neighbouring_index]
			distance_between_atoms = round(get_distance(position1, position2), 4)
			if not less_than_or_equal_to_max_bondlength(distance_between_atoms, element1, element2):
				disconnected_bonds.append((atom_index, neighbouring_index))

	# Third, remove the bonds from the wrapped molecule's graph
	wrapped_molecule_graph.remove_edges_from(disconnected_bonds)

	# Fourth, determine how many components their are in the unwrapped molecule.
	molecule_components_graphs = list(connected_components(wrapped_molecule_graph))

	# Fifth, obtain all the cell corners for the crystal. # This method doesnt seem to be used?
	#cell_points, unit_cell_displacements = get_cell_corner_points(crystal_cell_lattice, super_cell_reach=1, bottom_of_range=0, get_corrspeonding_ijk_values=True)

	# Fifth, get all the possible positions of the molecule in a cell.
	unit_cell_displacements = []
	for component in molecule_components_graphs:

		# 5.1: Take just the first index in the component list. Any one will do, so take the first.
		first_index = sorted(component)[0]

		# 5.2: get the position of the atom in the original connected molecule and the wrapped molecule.
		position_of_component_in_connected_molecule = molecule[first_index].position
		position_of_component_in_wrapped_molecule   = wrapped_molecule[first_index].position

		# 5.3: Get the translation from the atom in the original connected position to the wrapped position.
		#      * This translation should be alligned to a unit cell. 
		translation = position_of_component_in_wrapped_molecule - position_of_component_in_connected_molecule

		# 5.4: Determine the unit cell that this translated molecule will be assigned to with respect to the original unwrapped molecule. 
		unit_cell_displacement = convert_position_into_ijk_lengths(translation, crystal_cell_lattice, should_be_interger_length=True)

		# 5.5: Add the unit cell displacement to the unit_cell_displacements list
		if not unit_cell_displacement in unit_cell_displacements:
			unit_cell_displacements.append(unit_cell_displacement)

	# Sixth, return the unit_cell_displacements list
	return unit_cell_displacements
