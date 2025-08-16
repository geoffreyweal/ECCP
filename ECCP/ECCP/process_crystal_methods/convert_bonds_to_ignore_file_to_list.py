"""
convert_bonds_to_ignore_file_to_list.py, Geoffrey Weal, 28/2/23

This script is designed to turn a file containing the bonds between atoms to ignore in a crystal.

This is particularly useful if you have a polymer crystal or MOF/COF that need to be separated into components.
"""

def convert_bonds_to_ignore_file_to_list(path_to_bonds_to_ignore_file):
	"""
	This method is designed to turn a file containing the bonds between atoms to ignore in a crystal.

	This is particularly useful if you have a polymer crystal or MOF/COF that need to be separated into components.

	Parameters
	----------
	path_to_bonds_to_ignore_file : str.
		This is the path to the text file that contain information about the bonds to ignore. 
	Returns
	-------
	bonds_to_ignore : list
		This is a list of all the bonds to ignore, given as the atom index pair for the bond to ignore. 
	"""

	# First, initialise the list to record which bonds to ignore. 
	bonds_to_ignore = []

	# Second, read the file containing which bonds to ignore
	with open(path_to_bonds_to_ignore_file,'r') as bonds_to_ignoreFILE:

		# 2.1: For each line in the bonds_to_ignore.txt file.
		for line in bonds_to_ignoreFILE:

			# 2.2: Read the line, strip it, and split it into a list separated by spaces in line.
			line  = line.rstrip().split()

			# 2.3: Record atom index 1 in the bond.
			index1 = line[0]

			# 2.4: Record atom index 2 in the bond.
			index2 = line[1]

			# 2.5: Create the bonding pair.
			bonding_pair = tuple(sorted([index1, index2]))

			# Add the bonding pair to the bonds_to_ignore list.
			bonds_to_ignore.append(bonding_pair)

	# Third, return bonds_to_ignore
	return bonds_to_ignore
	