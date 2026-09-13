"""
remove_solvents_from_molecules_dict.py, Geoffrey Weal, 19/5/24

This method is designed to remove solvents from the molecules and molecule_graphs dictionaries.
"""

def remove_solvents_from_molecules_dict(molecules, molecule_graphs, SolventsList):
	"""
	This method is designed to remove solvents from the molecules and molecule_graphs dictionaries.

	Parameters
	----------
	molecules : dict.
		This dictionary holds the molecules in the crystal.
	molecule_graphs : dict.
		This dictionary holds the graphs of the molecules in the crystal.
	SolventsList : list.
		This is the list of molecules that have been considered as solvents
	"""

	# First, for each solvent in SolventsList:
	for solvent_mol_name in sorted(SolventsList, reverse=True):

		# Second, check if the solvent is in the molecules and molecule_graphs dictionaries.
		solvent_in_molecules       = solvent_mol_name in molecules.keys()
		solvent_in_molecule_graphs = solvent_mol_name in molecule_graphs.keys()

		# Third, if the solvent is in neither dictionary, it has already been removed.
		#        This happens when SolventsList names a molecule that was already taken out of
		#        molecules and molecule_graphs earlier (for example, as a symmetric duplicate).
		#        There is nothing left to do for this solvent, so move on to the next one.
		if (not solvent_in_molecules) and (not solvent_in_molecule_graphs):
			continue

		# Fourth, if the solvent is in the molecule_graphs dictionary but not the molecules
		#         dictionary, the two dictionaries have gone out of step with each other.
		if (not solvent_in_molecules) and solvent_in_molecule_graphs:
			to_string  = 'Error: The solvent is not in the molecule dictionary, but it was found in the graphs dictionary.\n'
			to_string += 'There is a problem, this should not happen.\n'
			to_string += f'SolventsList = {SolventsList}\n'
			to_string += f'solvent_mol_name being currently checked = {solvent_mol_name}\n'
			to_string += f'molecules = {sorted(molecules.keys())}\n'
			to_string += f'molecule_graphs = {sorted(molecule_graphs.keys())}\n'
			raise Exception(to_string)

		# Fifth, if the solvent is in the molecules dictionary but not the molecule_graphs
		#        dictionary, the two dictionaries have gone out of step with each other.
		if solvent_in_molecules and (not solvent_in_molecule_graphs):
			to_string  = 'Error: The solvent is not in the graphs dictionary, but it was found in the molecules dictionary.\n'
			to_string += 'There is a problem, this should not happen.\n'
			to_string += f'SolventsList = {SolventsList}\n'
			to_string += f'solvent_mol_name being currently checked = {solvent_mol_name}\n'
			to_string += f'molecules = {sorted(molecules.keys())}\n'
			to_string += f'molecule_graphs = {sorted(molecule_graphs.keys())}\n'
			raise Exception(to_string)

		# Sixth, the solvent is in both dictionaries, so remove it from both.
		del molecules[solvent_mol_name]
		del molecule_graphs[solvent_mol_name]
