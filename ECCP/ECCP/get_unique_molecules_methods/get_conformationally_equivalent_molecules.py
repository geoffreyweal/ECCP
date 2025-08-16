"""
get_conformationally_unique_molecules_indices.py, Geoffrey Weal, 21/4/24

This method is designed to determine which molecules are unique and which molecules are equivalent.
"""
from copy import deepcopy

from networkx import is_isomorphic
import networkx.algorithms.isomorphism as iso

nm = iso.categorical_node_match(attr=['E', 'H'], default=[None, 0])
def get_conformationally_equivalent_molecules(structurally_equivalent_molecule_groups, molecules, molecule_graphs):
	"""
	This method is designed to determine which molecules are unique and which molecules are equivalent.

	Note: If there is an issue, it may be based on the is_isomorphic method, as it is not 100% if it determines isomorphicity or not. 

	Parameters
	----------
	structurally_equivalent_molecule_groups : dict. 
		This dictionary contains all the structurally equivalent groups in the format {structurally unique molecule: list of structurally equivalent molecules}
	molecules : list of ase.Atoms
		These are all the molecules in the crystal unit cell.
	molecule_graphs : list of networkx.Graph
		These are all the bonding connectivity graphs of molecules in the crystal.

	Returns
	-------
	conformationally_equivalent_molecule_groups : dict.
		This is a dictionary of the conformationally unique molecule, with a lit of all the other molecules in the crystal that are equivalent to this unique molecule. 
	"""

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# Part 1: Obtain groups of molecules from the structurally unique molecules that are conformationally equivalent to each other. 

	# First, make a dictionary to add conformationally equivalent molecules to
	conformationally_equivalent_molecules = {}

	# Second, we create a list for remembering which molecules have been recorded as conformationally equivalent
	#         so that we do not process them in the third step for loop. 
	conformationally_equivalent_molecules_list = []

	# Third, get the structurally equivalent molecule names, and sort them. 
	structurally_unique_molecules_names = sorted(structurally_equivalent_molecule_groups.keys())

	# Note: We define the conformationally equivalent molecule as the one in a conformationally equivalent group
	#       that has the lowest name. 

	# Fourth, go through the list of structurally_unique_molecules_names
	for i1 in range(len(structurally_unique_molecules_names)):

		# 4.1: Get the molecule name for i1
		mol_name1 = structurally_unique_molecules_names[i1]

		# 4.2: If mol_name1 is in the conformationally_equivalent_molecules_list, 
		#      then it is conformationally equivalent to another molecule, so it has
		#      been processed and we move on.
		if mol_name1 in conformationally_equivalent_molecules_list:
			continue

		# 4.3: If we have not processed mol_name1 yet, it is by definition the conformationally unique
		#      and is the lowest value in its group as we performed the sorting algorithm in step 2, 
		#      so we create a new conformational equivalence group for it. 
		conformationally_equivalent_molecules[mol_name1] = []

		# 4.3: Get the molecule graph associated with mol_name1
		molecule_graph1 = molecule_graphs[mol_name1]

		# Fifth, go through the list of structurally_unique_molecules_names
		for i2 in range(i1+1,len(structurally_unique_molecules_names)):

			# 5.1: Get the mol_name molecule name for i1
			mol_name2 = structurally_unique_molecules_names[i2]

			# 5.2: Get the molecule graph associated with mol_name1
			molecule_graph2 = molecule_graphs[mol_name2]

			# Sixth, check if the two molecules are conformationally the same by 
			#         checking if their molecule connectivity graphs are equivalent.
			if is_isomorphic(molecule_graph1, molecule_graph2, node_match=nm):

				# 6.1: If these two molecules graphs are the same, they are conformationally equivalent.
				#      Add mol_name2 to conformationally_equivalent_molecules
				conformationally_equivalent_molecules.setdefault(mol_name1,[]).append(mol_name2)
				conformationally_equivalent_molecules_list.append(mol_name2)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# Part 2: Now that we have obtain how the structurally unique molecules are conformationally equivalent to each other,
	#         we want to add the structurally equivalent molecules to this dictionary. 
	#         For this exercise, conformationally_equivalent_molecules will turn into conformationally_equivalent_molecule_groups

	# Seventh, We want to make a copy of conformationally_equivalent_molecules, where the values are a list of all conformationally equivalent molecules. 
	conformationally_equivalent_molecule_groups = deepcopy(conformationally_equivalent_molecules)

	# Eighth, for each equivalent conformer, add all the structurally equivalent counterparts from structurally_equivalent_molecule_groups
	for con_unique_mol, con_equivalent_mols in conformationally_equivalent_molecules.items():
		for con_equivalent_mol in con_equivalent_mols:
			conformationally_equivalent_molecule_groups[con_unique_mol] += list(structurally_equivalent_molecule_groups[con_equivalent_mol])

	# Ninth, for each unique conformer, add all the structurally equivalent counterparts from structurally_equivalent_molecule_groups
	for con_unique_mol in conformationally_equivalent_molecules.keys():
		conformationally_equivalent_molecule_groups[con_unique_mol] += list(structurally_equivalent_molecule_groups[con_unique_mol])

	# Tenth, sort all the inputs in conformationally_equivalent_molecule_groups:
	for con_unique_mol in conformationally_equivalent_molecule_groups.keys():
		conformationally_equivalent_molecule_groups[con_unique_mol].sort()

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Eleventh, check that no molecule occurs twice in the conformationally_equivalent_molecule_groups dictionary
	checking_conformationally_equivalent_molecules_list = list(conformationally_equivalent_molecule_groups.keys())
	checking_conformationally_equivalent_molecules_list += [j for sub in conformationally_equivalent_molecule_groups.values() for j in sub]
	if not sorted(checking_conformationally_equivalent_molecules_list) == sorted(set(checking_conformationally_equivalent_molecules_list)):
		toString  = 'Error in def add_str_eqi_mols_to_con_equ_mol_dict:'+'\n'
		toString += 'At least 1 molecule in the molecule conformation dictionary appears twice after being processed by this method'+'\n'
		toString += 'conformationally_equivalent_molecule_groups = '+str(conformationally_equivalent_molecule_groups)
		raise Exception(toString)

	# Twelfth, return conformationally_equivalent_molecule_groups
	return conformationally_equivalent_molecule_groups






