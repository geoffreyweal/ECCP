"""
Neighbourhood_Generator_Multi_CPU.py, Geoffrey Weal, 3/2/23

This script is designed to generate all the neighbours that are possible between molecules in a crystal. 

Used for Multi-CPU nighbourhood methods.
"""
import numpy as np

from ECCP.ECCP.get_neighbouring_molecules_methods.Neighbourhood_Generator_supporting_methods.get_wrapped_complete_components_ijk_lengths import get_wrapped_complete_components_ijk_lengths
from SUMELF import remove_hydrogens
from SUMELF import Cell_Generator, convert_ijk_to_displacement_vector

class Neighbourhood_Generator_Multi_CPU:
	"""
	This object is designed to generator all the neighbours that are possible between molecules in a crystal. 

	Used for Multi-CPU neighbourhood methods.

	Parameters
	----------
	mol_name1 : int
		This is the name of molecule 1.
	mol_name2 : int
		This is the name of molecule 2.
	molecule1 : ase.Atoms object
		This is the ase.Atoms object for molecule 1. 
	molecule2 : ase.Atoms object
		This is the ase.Atoms object for molecule 2. 
	molecule_graph1 : networkx.Graph
		This is the graph that recording the bonding structure in molecule 1. 
	molecule_graph2 : networkx.Graph
		This is the graph that recording the bonding structure in molecule 2. 
	crystal_cell_lattice : numpy.array
		This is the matrix of the unit cell for the crystal.
	"""
	def __init__(self, mol_name1, mol_name2, molecule1, molecule2, molecule_graph1, molecule_graph2, crystal_cell_lattice):

		# First, save the input variables
		self.mol_name1            = mol_name1
		self.mol_name2            = mol_name2
		self.molecule1            = molecule1
		self.molecule2            = molecule2
		self.molecule_graph1      = molecule_graph1
		self.molecule_graph2      = molecule_graph2
		self.crystal_cell_lattice = crystal_cell_lattice
		self.origin_cell_point    = np.array((0,0,0))

		# Second, obtain the ijk placements of the components of the wrapped molecule. 
		self.molecule_1_translations = sorted(get_wrapped_complete_components_ijk_lengths(molecule1, molecule_graph1, crystal_cell_lattice))
		self.molecule_2_translations = sorted(get_wrapped_complete_components_ijk_lengths(molecule2, molecule_graph2, crystal_cell_lattice))

		# Third, remove all hydrogen from the molecule, as we dont want to include hydeogens in our analysis
		self.molecule1_with_no_hydrogens = remove_hydrogens(molecule1)
		self.molecule2_with_no_hydrogens = remove_hydrogens(molecule2)

		# Fourth, get the positions of the molecules with no hydrogens.
		self.positions1 = self.molecule1_with_no_hydrogens.get_positions()
		self.positions2 = self.molecule2_with_no_hydrogens.get_positions()

	# -----------------------------------------------------------------------------------------------------------------------

	def generator(self):
		"""
		This method will create the pair of molecules that may be neighbouring each other.

		Return
		------
		mol_name1 : int
			This is the name of the first neighbouring molecule you want to examine from the molecules list.
		mol_name2 : int
			This is the name of the second neighbouring molecule you want to examine from the molecules list.
		positions1 : numpy.array
			This is the position of the first neighbouring molecule in the neighbour pair.
		positions2 : numpy.array
			This is the position of the second neighbouring molecule in the neighbour pair.
		displacement : numpy.array
			This is the displacement of the second neighbouring molecule to create your neighbour pair
		unit_cell_displacement : numpy.array
			This is the displacement of the second neighbouring molecule to create your neighbour pair, scaled to the unit cell. 
		"""

		# First, setup the already_created_neighbouring_pairs to prevent neighbouring pairs from being 
		#        recreated if they are separated by the same unit cell ijk lengths.
		already_created_neighbouring_pairs = {}

		# Second, get a component of Monomer mol_name1. 
		for indexA in range(len(self.molecule_1_translations)):
			molecule_1_translation = self.molecule_1_translations[indexA]

			# Third, get a component of Monomer mol_name2 
			for indexB in range(len(self.molecule_2_translations)):
				molecule_2_translation = self.molecule_2_translations[indexB]

				# Fourth, create the cell generator object
				cell_generator = Cell_Generator()

				# Fifth, get the displacement of the molecule mol_name2 from its original unit cell (being from indexB-indexA). 
				for unit_cell_displacement in cell_generator.generate_next_ijk_points(): 

					# Sixth, determine the relatice unit cell displacement between molecules mol_name1 and mol_name2, given in ijk units.
					#        * This is obtain by unit_cell_displacement - molecule_1_translation + molecule_2_translation
					relative_unit_cell_displacement_between_mol1_and_mol2 = tuple([(value1+value2-value3) for value1, value2, value3 in zip(unit_cell_displacement, molecule_2_translation, molecule_1_translation)])

					# Seventh, get the relative displacement of molecules mol_name1 and mol_name2.
					relative_displacement = convert_ijk_to_displacement_vector(relative_unit_cell_displacement_between_mol1_and_mol2, self.crystal_cell_lattice)

					# Eighth, do not include a neighbour pair of molecules that involves a molecule with itself.
					if (self.mol_name1 == self.mol_name2) and (relative_displacement == self.origin_cell_point).all(): 
						continue

					# Ninth, only return (yield) info about this neighbouring pair of molecules if it has not been recorded yet.
					if not (relative_unit_cell_displacement_between_mol1_and_mol2 in already_created_neighbouring_pairs.keys()):

						# Tenth, yield the neighbouring pair information, and obtain feedback if this neighbouring pair of molecules was accepted by the method usig this generator.
						are_molecules_within_max_neighbour_distance = yield (self.mol_name1, self.mol_name2, self.positions1, self.positions2, relative_displacement, relative_unit_cell_displacement_between_mol1_and_mol2)
						yield 'Go to get_neighbours method for loop'

						# Eleventh, record this unit_cell_displacement for molecules mol_name1 and mol_name2, as we have not created this neighbouring pair of molecules yet.
						already_created_neighbouring_pairs[relative_unit_cell_displacement_between_mol1_and_mol2] = are_molecules_within_max_neighbour_distance

					else:

						# Twelfth, obtain what the result was if this neighbouring pair of molecules has already been recorded.
						are_molecules_within_max_neighbour_distance = already_created_neighbouring_pairs[relative_unit_cell_displacement_between_mol1_and_mol2]

					# Thirteenth, send the result of if a neighbouring pair of molecules was generated from unit_cell_displacement to the cell_generator
					cell_generator.add_result(are_molecules_within_max_neighbour_distance)

# ---------------------------------------------------------------------------------------------------------------------------

