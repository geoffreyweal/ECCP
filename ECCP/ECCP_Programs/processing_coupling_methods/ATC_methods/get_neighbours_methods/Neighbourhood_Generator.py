"""
Neighbourhood_Generator.py, Geoffrey Weal, 3/2/23

This script is designed to generate all the neighbours that are possible between molecules in a crystal. 
"""
import numpy as np
from copy import deepcopy
from tqdm import tqdm

from networkx import connected_components

from SUMELF import get_cell_corner_points, get_distance, less_than_or_equal_to_max_bondlength
from ECCP.ECCP_Programs.processing_coupling_methods.ATC_methods.get_neighbours_methods.supporting_methods.supporting_methods import is_molecule_already_recorded
from SUMELF import Cell_Generator, convert_ijk_to_displacement_vector, convert_position_into_ijk_lengths

class Neighbourhood_Generator:
	"""
	This object is designed to generator all the neighbours that are possible between molecules in a crystal. 

	Parameters
	----------
	molecules : list of ase.Atoms objects.
		These are all the individual molecules identified in the crystal that you want to determine neighbours for.
	molecule_graphs : list of networkx.Graph
		These are the graph that are associated with the molecules in the molecules list. 
	crystal_cell_lattice : numpy.array
		This is the matrix of the unit cell for the crystal.
	"""
	def __init__(self, molecules, molecule_graphs, crystal_cell_lattice, show_progressbar=True):

		# First, save the input variables
		self.molecules = molecules
		self.molecule_graphs = molecule_graphs
		self.crystal_cell_lattice = crystal_cell_lattice
		self.show_progressbar = show_progressbar
		self.origin_cell_point = np.array((0,0,0))

		# Second, obtain the ijk placements of the components of the wrapped molecule. 
		self.all_wrapped_molecules_ijk_placements = []
		for index in range(len(self.molecules)):
			wrapped_molecules_ijk_placements = self.get_wrapped_complete_components_ijk_lengths(self.molecules[index], self.molecule_graphs[index], crystal_cell_lattice)
			self.all_wrapped_molecules_ijk_placements.append(sorted(wrapped_molecules_ijk_placements))

	def get_wrapped_complete_components_ijk_lengths(self, molecule, molecule_graph, crystal_cell_lattice):
		"""
		This method will wrap your molecule, and provide the wrapped graph.

		Parameters
		----------
		molecule : ase.Atoms objects.
			This is the molecule
		molecule_graph : networkx.Graph
			This is the graph of the molecule. 

		Returns
		-------
		wrapped_molecule : ase.Atoms objects.
			This is the wrapped molecule
		wrapped_molecule_graphs : networkx.Graph
			This is the graph of the wrapped molecule. 
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

		# Fifth, obtain all the cell corners for the crystal.
		#cell_points, unit_cell_displacements = get_cell_corner_points(self.crystal_cell_lattice, super_cell_reach=1, bottom_of_range=0, get_corrspeonding_ijk_values=True)

		# Sixth, all the possible positions of the molecule in a cell.
		unit_cell_displacements = []
		for component in molecule_components_graphs:

			# 6.1: Take just the first index in the component list. Any one will do, so take the first.
			first_index = sorted(component)[0]

			# 6.2: get the position of the atom in the original connected molecule and the wrapped molecule.
			position_of_component_in_connected_molecule = molecule[first_index].position
			position_of_component_in_wrapped_molecule   = wrapped_molecule[first_index].position

			# 6.3: Get the translation from the atom in the original connected position to the wrapped position.
			# This translation should be alligned to a unit cell. 
			translation = position_of_component_in_wrapped_molecule - position_of_component_in_connected_molecule

			# 6.4: Determine the unit cell that this translated molecule will be assigned to with respect to the original unwrapped molecule. 
			unit_cell_displacement = convert_position_into_ijk_lengths(translation, self.crystal_cell_lattice, should_be_interger_length=True)

			# 6.5: Add the unit cell displacement to the unit_cell_displacements list
			if not unit_cell_displacement in unit_cell_displacements:
				unit_cell_displacements.append(unit_cell_displacement)

		# Seventh, return the unit_cell_displacements list
		return unit_cell_displacements

	# -----------------------------------------------------------------------------------------------------------------------

	def generator(self):
		"""
		This method will create the pair of molecules that may be neighbouring each other.

		Return
		------
		index1 : int
			This is the index of the first neighbouring molecule you want to examine from the molecules list.
		index2 : int
			This is the index of the second neighbouring molecule you want to examine from the molecules list.
		positions1 : numpy.array
			This is the position of the first neighbouring molecule in the neighbour pair.
		positions2 : numpy.array
			This is the position of the second neighbouring molecule in the neighbour pair.
		displacement : numpy.array
			This is the displacement of the second neighbouring molecule to create your neighbour pair
		unit_cell_displacement : numpy.array
			This is the displacement of the second neighbouring molecule to create your neighbour pair, scaled to the unit cell. 
		"""

		# First, determine the number of overall molecule positions that will need to be examined by this method. 
		counter = 0
		for index1 in range(len(self.molecules)): 
			for index2 in range(index1,len(self.molecules)): 
				counter += 1

		# Second, setup the process bar.
		if self.show_progressbar:
			self.pbar = tqdm(total=counter, unit='task')

		# Third, get the first molecule. 
		for index1 in range(len(self.molecules)): 
			positions1              = self.molecules[index1].get_positions()
			charges1                = self.molecules[index1].get_initial_charges()
			molecule_1_translations = self.all_wrapped_molecules_ijk_placements[index1]

			# Fourth, get the second molecule. This could be the same molecule as index1, but will be displaced to a different position.
			for index2 in range(index1,len(self.molecules)): 
				positions2              = self.molecules[index2].get_positions()
				charges2                = self.molecules[index2].get_initial_charges()
				molecule_2_translations = self.all_wrapped_molecules_ijk_placements[index2]

				# Fifth, setup the already_created_neighbouring_pairs to prevent neighbouring pairs of moelcules from being 
				# recreated if they are separated by the same unit cell ijk lengths.
				already_created_neighbouring_pairs = {}

				# Sixth, get a component of Monomer index1. 
				for indexA in range(len(molecule_1_translations)):
					molecule_1_translation = molecule_1_translations[indexA]

					# Seventh, get a component of Monomer index2 
					for indexB in range(len(molecule_2_translations)):
						molecule_2_translation = molecule_2_translations[indexB]

						# Eighth, create the cell generator object
						cell_generator = Cell_Generator()

						# Ninth, get the displacement of the molecule index2 from its original unit cell (being from indexB-indexA). 
						for unit_cell_displacement in cell_generator.generate_next_ijk_points(): 

							# Tenth, determine the relative unit cell displacement between molecules index1 and index2, given in ijk units.
							relative_unit_cell_displacement_between_mol1_and_mol2 = tuple([(value1+value2-value3) for value1, value2, value3 in zip(unit_cell_displacement, molecule_2_translation, molecule_1_translation)])

							# Eleventh, get the relative displacement of molecules index1 and index2.
							relative_displacement = convert_ijk_to_displacement_vector(relative_unit_cell_displacement_between_mol1_and_mol2, self.crystal_cell_lattice)

							# Twelfth, do not include a neighbour pair of molecules that involves a molecule with itself.
							if (index1 == index2) and (relative_displacement == self.origin_cell_point).all(): 
								continue

							# Thirteenth, only return (yield) info about this neighbouring pair of molecules if it has not been recorded yet.
							if not (relative_unit_cell_displacement_between_mol1_and_mol2 in already_created_neighbouring_pairs.keys()):

								# Fourtheenth, indicate to the terminal what this generator is returning (yielding) to the program using this generator. 
								if self.show_progressbar:
									self.pbar.set_description("Processing M"+str(index1+1)+' M'+str(index2+1)+'; rel. ijk disp.: '+str(''.join(str(relative_unit_cell_displacement_between_mol1_and_mol2).split())))

								# Fifteenth, yield the neighbouring pair information, and obtain feedback if this neighbouring pair of molecules was accepted by the method usig this generator.
								are_molecules_within_max_neighbour_distance = yield (index1, index2, positions1, positions2, charges1, charges2, relative_displacement, relative_unit_cell_displacement_between_mol1_and_mol2)
								yield 'Go to for loop'

								# Sixteenth, record this unit_cell_displacement for molecules index1 and index2, as we have not created this neighbouring pair of molecules yet.
								already_created_neighbouring_pairs[relative_unit_cell_displacement_between_mol1_and_mol2] = are_molecules_within_max_neighbour_distance

							else:

								# Seventeenth, obtain what the result was if this neighbouring pair of molecules has already been recorded.
								are_molecules_within_max_neighbour_distance = already_created_neighbouring_pairs[relative_unit_cell_displacement_between_mol1_and_mol2]

							# Eighteenth, send the result of if a neighbouring pair of molecules was generated from unit_cell_displacement to the cell_generator
							cell_generator.add_result(are_molecules_within_max_neighbour_distance)

				# Nineteenth, update the process bar.
				if self.show_progressbar:
					self.pbar.update(1)

		# Twentieth, close the progress bar.
		if self.show_progressbar:
			self.pbar.close()

# ---------------------------------------------------------------------------------------------------------------------------


