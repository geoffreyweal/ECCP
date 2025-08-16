"""
write_crystal_structures_to_disk.py, Geoffrey Weal, 14/4/22

This script is designed to write all the various crystal files of the crystal to disk.
"""
from copy                                  import deepcopy
from ase.io                                import write
from SUMELF                                import add_graph_to_ASE_Atoms_object
from SUMELF                                import check_molecule_against_file
from ECCP.ECCP.make_human_friendly_crystal import make_human_friendly_crystal

def write_crystal_structures_to_disk(original_crystal, original_molecules, path_to_eccp_folder, crystal_graph=None, molecule_graphs=None, structurally_equivalent_molecules=None, conformationally_equivalent_molecules=None):
	"""
	This method is designed to write all the various crystal files of the crystal to disk.

	Parameters
	----------
	crystal : ase.Atoms 
		This is the crystal.
	molecules : list of ase.Atoms 
		This is a list of all the molecules that make up the crystal.
	path_to_eccp_folder : str. 
		This is the path to the eccp folder.
	crystal_graph : networkx.Graph
		This is the graph for the crystal object. If None given, this graph will not be used to write crystal structures to file. Default: None.
	molecule_graphs : list of networkx.Graph
		This is a list of all the graph of the molecules in the molecules list. If None given, this graph will not be used to write crystal structures to file. Default: None.
	structurally_equivalent_molecules : dict.
		This is a dictionary of the structurally unique molecule, with a lit of all the other molecules in the crystal that are equivalent to this unique molecule. Default: None
	conformationally_equivalent_molecules : dict.
		This is a dictionary of the conformationally unique molecule, with a lit of all the other molecules in the crystal that are equivalent to this unique molecule. Default: None
	"""

	# First, make a copy of the original original_crystal ase.Atoms object
	crystal = original_crystal.copy()

	# Second, if the crystal graph given, attach it to the crystal ase.Atoms object
	if crystal_graph is not None:
		add_graph_to_ASE_Atoms_object(crystal, crystal_graph)

	# Third, set the structurally equivalent molecules in the crystal object
	if structurally_equivalent_molecules is not None:
		structurally_equivalent_molecules = convert_dictionary_to_string(structurally_equivalent_molecules)
		if not structurally_equivalent_molecules == '':
			crystal.info['structurally_equivalent_molecules'] = structurally_equivalent_molecules

	# Fourth, set the conformationally equivalent molecules in the crystal object
	if conformationally_equivalent_molecules is not None:
		conformationally_equivalent_molecules = convert_dictionary_to_string(conformationally_equivalent_molecules)
		if not conformationally_equivalent_molecules == '':
			crystal.info['conformationally_equivalent_molecules'] = conformationally_equivalent_molecules

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Fifth, create the path to the crystal.xyz file in the ECCP_Information folder.
	path_to_crystal_xyz_file = path_to_eccp_folder+'/crystal.xyz'

	# Sixth, check that the crystal object is consistent with that currently on file.
	check_molecule_against_file(crystal, path_to_crystal_xyz_file)

	# Seventh, write the crystal as an xyz file.
	write(path_to_crystal_xyz_file, crystal)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Eighth, make a copy of the original molecules.
	molecules = {mol_name: molecule.copy() for mol_name, molecule in sorted(original_molecules.items(), key=lambda x: x[0])}

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# Ninth, make a small human-friendly to view crystal file.
	#        * This file contains only the atoms that could fit in the unit cell.

	# 9.1: Create a version of the crystal in a human friendly format.
	human_friendly_crystal_small, human_friendly_crystal_small_graph = make_human_friendly_crystal(molecules, super_cell_reach=0, molecule_graphs=molecule_graphs)
	
	# 9.2: Add the graph object of human_friendly_crystal_small to the ase.atoms object.
	add_graph_to_ASE_Atoms_object(human_friendly_crystal_small, human_friendly_crystal_small_graph)
	
	# 9.3: Create the path to the human_friendly_crystal_small.xyz file in the ECCP_Information folder.
	small_crystal_name = path_to_eccp_folder+'/human_friendly_crystal_small.xyz'

	# 9.4: Check that the human_friendly_crystal_small object is consistent with that currently on file.
	check_molecule_against_file(human_friendly_crystal_small, small_crystal_name)

	# 9.5: Write human_friendly_crystal_small to file.
	write(small_crystal_name, human_friendly_crystal_small)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# Tenth, make a larger human-friendly to view crystal file.
	#             This file contains molecules in the origin unit cell, as well as neighbouring unit cells about the origin cell.

	# 10.1: Create a version of the crystal in a human friendly format.
	human_friendly_crystal_large, human_friendly_crystal_large_graph = make_human_friendly_crystal(molecules, super_cell_reach=1, molecule_graphs=molecule_graphs)

	# 10.2: Add the graph object of human_friendly_crystal_large to the ase.atoms object.
	add_graph_to_ASE_Atoms_object(human_friendly_crystal_large, human_friendly_crystal_large_graph)

	# 10.3: Create the path to the human_friendly_crystal_large.xyz file in the ECCP_Information folder.
	large_crystal_name = path_to_eccp_folder+'/human_friendly_crystal_large.xyz'

	# 10.4: Check that the human_friendly_crystal_large object is consistent with that currently on file.
	check_molecule_against_file(human_friendly_crystal_small, small_crystal_name)

	# 10.5: Write human_friendly_crystal_large to file.
	write(large_crystal_name, human_friendly_crystal_large)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

def convert_dictionary_to_string(original_equivalent_molecules):
	"""
	This method is designed to convert the dictionaries for structurally_equivalent_molecules and conformationally_equivalent_molecules from its dictionary form into a string. 

	Parameters
	----------
	original_equivalent_molecules : dict
		This dictionary indicates which unique molecules are equivalent to which other molecules in the crystal. 

	Returns
	-------
	equivalent_molecules_string : str
		This is the equivalent_molecules dictionary in compressed string form
	"""

	# First, make a copy of the equivalent_molecules dictionary
	equivalent_molecules = deepcopy(original_equivalent_molecules)

	# Second, sort all the lists in equivalent_molecules.values()
	for unique_index in equivalent_molecules:
		equivalent_molecules[unique_index].sort()

	# Third, order the dictionary by the key to make it easier to see this in the xyz file. 
	equivalent_molecules = {key: value for (key, value) in sorted(equivalent_molecules.items())}

	# Fourth, convert dictionary into string
	equivalent_molecules_string = []
	for unique_index, list_of_equivalent_mol_indices in sorted(equivalent_molecules.items()):
		to_string = str(unique_index)+':'+','.join([str(value) for value in list_of_equivalent_mol_indices])
		equivalent_molecules_string.append(to_string)
	equivalent_molecules_string = ';'.join(equivalent_molecules_string)

	# Fifth, return equivalent_molecules_string
	if equivalent_molecules_string == '':
		return '-'
	else:
		return equivalent_molecules_string



