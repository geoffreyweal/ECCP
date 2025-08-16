"""
Neighbourhood_Generator.py, Geoffrey Weal, 3/2/23

This script is designed to generate all the neighbours that are possible between molecules in a crystal. 
"""
import numpy as np
from tqdm import tqdm

from ECCP.ECCP.get_neighbouring_molecules_methods.Neighbourhood_Generator_supporting_methods.get_wrapped_complete_components_ijk_lengths import get_wrapped_complete_components_ijk_lengths
from SUMELF import remove_hydrogens
from SUMELF import Cell_Generator, convert_ijk_to_displacement_vector

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

		# First, save the input variables.
		self.molecules = molecules
		self.molecule_graphs = molecule_graphs
		self.crystal_cell_lattice = crystal_cell_lattice
		self.show_progressbar = show_progressbar
		self.origin_cell_point = np.array((0,0,0))

		# Second, obtain the ijk placements of the components of the wrapped molecules.
		self.all_wrapped_molecules_ijk_placements = {}
		for mol_name in self.molecules.keys():

			# 2.1: Get the ijk placements of the components of the wrapped molecule.
			wrapped_molecules_ijk_placements = get_wrapped_complete_components_ijk_lengths(self.molecules[mol_name], self.molecule_graphs[mol_name], self.crystal_cell_lattice)

			# 2.2: Sort and add the ijk placements to all_wrapped_molecules_ijk_placements.
			self.all_wrapped_molecules_ijk_placements[mol_name] = sorted(wrapped_molecules_ijk_placements)

		# Third, remove all hydrogen from the molecule, as we dont want to include hydeogens in our analysis.
		self.molecules_with_no_hydrogens = {mol_name: remove_hydrogens(molecule) for mol_name, molecule in self.molecules.items()}

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

		# First, get the sorted keys from molecules_with_no_hydrogens.
		molecule_names = sorted(self.molecules_with_no_hydrogens.keys())

		# Second, determine the number of overall molecule positions that will need to be examined by this method. 
		total = int((len(molecule_names)*(len(molecule_names)+1.0))/2.0)

		# Third, setup the process bar.
		if self.show_progressbar:
			self.pbar = tqdm(total=total, unit='task')

		# Fourth, get the first molecule. 
		for index1 in range(len(molecule_names)): 
			mol_name1               = molecule_names[index1]
			positions1              = self.molecules_with_no_hydrogens[mol_name1].get_positions()
			molecule_1_translations = self.all_wrapped_molecules_ijk_placements[mol_name1]

			# Fifth, get the second molecule. This could be the same molecule as mol_name1, but will be displaced to a different position.
			for index2 in range(index1,len(molecule_names)): 
				mol_name2               = molecule_names[index2]
				positions2              = self.molecules_with_no_hydrogens[mol_name2].get_positions()
				molecule_2_translations = self.all_wrapped_molecules_ijk_placements[mol_name2]

				# Sixth, setup the already_created_neighbouring_pairs to prevent neighbouring pairs of moelcules from being 
				# recreated if they are separated by the same unit cell ijk lengths.
				already_created_neighbouring_pairs = {}

				# Seventh, get a component of Monomer mol_name1 that spills into another unit cell from the origin (molecule_1_translations will include the origin).
				for indexA in range(len(molecule_1_translations)):
					molecule_1_translation = molecule_1_translations[indexA]

					# Eighth, get a component of Monomer mol_name2  that spills into another unit cell from the origin (molecule_1_translations will include the origin).
					for indexB in range(len(molecule_2_translations)):
						molecule_2_translation = molecule_2_translations[indexB]

						# Ninth, create the cell generator object
						cell_generator = Cell_Generator()

						# Tenth, get the displacement of the molecule mol_name2 from its original unit cell (being from indexB-indexA). 
						for unit_cell_displacement in cell_generator.generate_next_ijk_points(): 

							# Eleventh, determine the relative unit cell displacement between molecules mol_name1 and mol_name2, given in ijk units.
							#           * This is obtain by unit_cell_displacement - molecule_1_translation + molecule_2_translation
							relative_unit_cell_displacement_between_mol1_and_mol2 = tuple([(value1+value2-value3) for value1, value2, value3 in zip(unit_cell_displacement, molecule_2_translation, molecule_1_translation)])

							# Twelfth, get the relative displacement of molecules mol_name1 and mol_name2 (in angstroms). 
							relative_displacement = convert_ijk_to_displacement_vector(relative_unit_cell_displacement_between_mol1_and_mol2, self.crystal_cell_lattice)

							# Thirteenth, do not include a neighbour pair of molecules that involves a molecule with itself.
							if (mol_name1 == mol_name2) and (relative_displacement == self.origin_cell_point).all(): 
								continue

							# Fourtheenth, only return (yield) info about this neighbouring pair of molecules if it has not been recorded yet.
							if not (relative_unit_cell_displacement_between_mol1_and_mol2 in already_created_neighbouring_pairs.keys()):

								# Fifteenth, indicate to the terminal what this generator is returning (yielding) to the program using this generator. 
								if self.show_progressbar:
									self.pbar.set_description("Processing M"+str(mol_name1)+' M'+str(mol_name2)+'; rel. ijk disp.: '+str(''.join(str(relative_unit_cell_displacement_between_mol1_and_mol2).split())))

								# Sixteenth, yield the neighbouring pair information, and obtain feedback if this neighbouring pair of molecules was accepted by the method usig this generator.
								are_molecules_within_max_neighbour_distance = yield (mol_name1, mol_name2, positions1, positions2, relative_displacement, relative_unit_cell_displacement_between_mol1_and_mol2)
								yield 'Go to get_neighbours method for loop'

								# Seventeenth, record this unit_cell_displacement for molecules mol_name1 and mol_name2, as we have not created this neighbouring pair of molecules yet.
								already_created_neighbouring_pairs[relative_unit_cell_displacement_between_mol1_and_mol2] = are_molecules_within_max_neighbour_distance

							else:

								# Eighteenth, obtain what the result was if this neighbouring pair of molecules has already been recorded.
								are_molecules_within_max_neighbour_distance = already_created_neighbouring_pairs[relative_unit_cell_displacement_between_mol1_and_mol2]

							# Nineteenth, send the result of if a neighbouring pair of molecules was generated from unit_cell_displacement to the cell_generator
							cell_generator.add_result(are_molecules_within_max_neighbour_distance)

				# Twentieth, update the process bar.
				if self.show_progressbar:
					self.pbar.update(1)

		# Twenty-first, close the progress bar.
		if self.show_progressbar:
			self.pbar.close()

# ---------------------------------------------------------------------------------------------------------------------------


