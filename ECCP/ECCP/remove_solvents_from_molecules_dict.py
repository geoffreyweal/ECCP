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

	print('This method has not been tested yet. Make sure it works before removing this checkpoint code.')
	print('This about if SolventsList should contain symmetric molecules that will have been removed from molecules and molecule_graphs already.')
	print('Actually, I think symmetruc molecules should not have been removed at this point. Do a check here')
	import pdb; pdb.set_trace()

	# First, for each solvent in SolventsList:
	for solvent_mol_name in sorted(SolventsList, reverse=True):

		# Second, Check if the solvent if in the molecules and molecule_graphs dictionaries. 
		solvent_in_molecules       = solvent_mol_name not in molecules.keys()
		solvent_in_molecule_graphs = solvent_mol_name not in molecule_graphs.keys()

		# Third, if solvent is not in the molecules dictionary but it is in the molecule_graphs dictionary, there is a problem
		if not solvent_in_molecules and solvent_in_molecule_graphs:
			to_string  = 'Error: The solvent is not in the molecule dictionary, but it was found in the graphs dictionary.\n'
			to_string += 'There is a problem, this should not happen.\n'
			to_string += f'SolventsList = {SolventsList}\n'
			to_string += f'solvent_mol_name being currently checked = {solvent_mol_name}\n'
			to_string += f'molecules = {sorted(molecules.keys())}\n'
			to_string += f'molecule_graphs = {sorted(molecule_graphs.keys())}\n'
			raise Exception(to_string)

		# Fourth, if solvent_mol_name is not in molecules dictionary, there is a problem.
		#         * solvent_mol_name should be in molecules currently. 
		if solvent_mol_name not in molecules.keys():
			to_string  = 'Error: The solvent is not in the molecule dictionary.\n'
			to_string += 'There is a problem, all the mol_names in SolventsList should be in molecules here.\n'
			to_string += f'SolventsList = {SolventsList}\n'
			to_string += f'solvent_mol_name being currently checked = {solvent_mol_name}\n'
			to_string += f'molecules = {sorted(molecules.keys())}\n'
			raise Exception(to_string)

		# Fifth, if solvent is not in the molecule_graphs dictionary at this point, there is a problem
		if not solvent_in_molecule_graphs:
			to_string  = 'Error: The solvent is not in the graphs dictionary, but it was found in the molecules dictionary.\n'
			to_string += 'There is a problem, this should not happen.\n'
			to_string += f'SolventsList = {SolventsList}\n'
			to_string += f'solvent_mol_name being currently checked = {solvent_mol_name}\n'
			to_string += f'molecules = {sorted(molecules.keys())}\n'
			to_string += f'molecule_graphs = {sorted(molecule_graphs.keys())}\n'
			raise Exception(to_string)

		# Sixth, remove solvent_mol_name from molecules dictionary.
		del molecules[solvent_mol_name]

		# Seventh, remove solvent_mol_name from molecule_graphs dictionary.
		del molecule_graphs[solvent_mol_name]
