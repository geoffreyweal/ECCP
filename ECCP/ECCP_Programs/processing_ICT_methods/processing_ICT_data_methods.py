'''
Geoffrey Weal, get_matrix_from_file.py, 16/5/22

This program is designed to import the matrix from a file into python in numpy
'''

import os
import numpy as np

def get_folder_list(path_to_folder):
    """
    This method will provide a list of the names of all the subfolders in this folder

    Parameters
    ----------
    path_to_folder : str.
        This is the path to the folder to obtain the names of its subfolders
    """
    return [name for name in os.listdir(path_to_folder) if os.path.isdir(path_to_folder+'/'+name)]

# ---------------------------------------------------------------------------------------------------

def get_matrix_from_file(filename):
	"""
	This method will return the matrix from a text file as a numpy array

	Parameters
	----------
	filename : str.
		This is the path of the file to read.

	Returns
	-------
	matrix : numpy.array
		This is the matrix that is contained in the file.
	"""

	# First, read the rows of the matrix
	matrix = []
	with open(filename,'r') as fileTXT:
		for line in fileTXT:
			row = [float(value) for value in line.rstrip().split()]
			matrix.append(row)

	# Second, get some test variables.	
	is_a_symmetric_matrix = ([len(row) for row in matrix] == list(range(1,len(matrix)+1)))
	length_of_last_row = len(matrix[-1])

	# Third, write the matrix
	if is_a_symmetric_matrix:

		# 3.1: If the matrix is symmetric, then you need to do the following to make the full matrix

		# 3.1.1: Make a lower triangle matrix and add 0s to the ends of it
		for index in range(len(matrix)):
			row = matrix[index]
			new_row = row + [0]*(length_of_last_row-len(row))
			matrix[index] = new_row

		# 3.1.2: Make the lower triangle
		lower_triangle = np.array(matrix)

		# 3.1.3: Make the uper triangle that excludes the diagonal
		upper_triangle = lower_triangle.T
		upper_triangle_plus_1 = np.triu(upper_triangle,+1)

		# 3.1.4: Return the completed symmetric matric
		full_symmetric_matrix = lower_triangle + upper_triangle_plus_1
		return full_symmetric_matrix

	else:
		
		# 3.2: For a 1D matrix or a non-symmetric 2D matrix, just convert your list to a numpy.array
		return np.array(matrix)

def get_MO_orbital_names(filename):
	"""
	This method will obtain the rows in each matrix that go with each atom, just in case of rearrangements in atom indices in the dimer compared to the monomers. 

	Parameters
	----------
	filename : str.
		This is the path of the file to read.

	Returns
	-------
	atom_MOs : dict.
		This contains the indices of the rows that each atom in the monomer/dimer goes with (i.e. the basis sets for each atom)
	"""
	atom_MOs = {}
	with open(filename,'r') as fileTXT:
		counter = 0
		for line in fileTXT:
			if counter == 0:
				atom_index = int(line.rstrip().replace('Atom Index: ',''))
			elif counter == 1:
				MO_coefficient_rows_related_to_atom_index = [int(value) for value in line.rstrip().split()]
				atom_MOs[atom_index] = MO_coefficient_rows_related_to_atom_index
			elif counter == 2:
				counter = -1
			counter += 1
	return atom_MOs

def get_MO_occupancies(filename):
	"""
	This method obtains the occupancies of electrons in each MO from a file

	Parameters
	----------
	filename : str.
		This is the path of the file to read.

	Returns
	-------
	occupancies : list.
		This list indicates the electron occupancy in each MO.
	"""

	with open(filename,'r') as fileTXT:
		occupancies = fileTXT.readline()
		occupancies = [value for value in occupancies.rstrip().split()]
	return occupancies

# ---------------------------------------------------------------------------------------------------

def assign_MO_coefficients_with_atoms(MO_coefficients, MO_orbital_names):
	"""
	This method will assign the rows of the MO_coefficients to each atom in your monomer.

	Parameters
	----------
	MO_coefficients : numpy.array
		This is the matrix that contains the MO coefficients for each MO
	MO_orbital_names : dict.
		These are the rows in the MO_coefficients matrix that go with each atom in the monomer

	Returns
	-------
	assigned_MO_coefficients : list.
		These are the MO coefficients of each row in the MO_coefficients that go with each atom.
	"""
	assigned_MO_coefficients = {}
	for atom_index, rows in MO_orbital_names.items():
		assigned_MO_coefficients[atom_index] = {row_index: MO_coefficients[row_index,:] for row_index in rows}
	return assigned_MO_coefficients

def create_MO_coefficients_matrix(MO_coefficients_data):
	"""
	This method will create the MO_coefficient matrix for the monomer in the atom order required for the dimer.

	Parameters
	----------
	MO_coefficients : list of lists
		This is the matrix that contains the MO coefficients for each MO as a list of lists

	Returns
	-------
	MO_coefficients : numpy.array
		This is the matrix that contains the MO coefficients for each MO as a numpy.array
	"""

	MO_coefficients_matrix = []
	for atom_index in sorted(MO_coefficients_data.keys()):
		for index in sorted(MO_coefficients_data[atom_index].keys()):
			MO_coefficients_matrix.append(MO_coefficients_data[atom_index][index])

	MO_coefficients_matrix = np.vstack(MO_coefficients_matrix)

	return MO_coefficients_matrix

# ---------------------------------------------------------------------------------------------------






