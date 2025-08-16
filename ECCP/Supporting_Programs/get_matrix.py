#!/usr/bin/env python3
'''
Geoffrey Weal, get_matrix.py, 9//22

This program is designed to provide the matrix given into python in numpy
'''

import numpy as np

def get_matrix(filename):
	"""
	This method will return a 
	"""

	matrix = []
	with open(filename,'r') as fileTXT:
		for line in fileTXT:
			row = [float(value) for value in line.rstrip().split()]
			matrix.append(row)

	symmetric_matrix = ([len(row) for row in matrix] == list(range(1,len(matrix)+1)))
	length_of_last_row = len(matrix[-1])

	if symmetric_matrix:
		for index in range(len(matrix)):
			row = matrix[index]
			new_row = row + [0]*(length_of_last_row-len(row))
			matrix[index] = new_row

		lower_triangle = np.array(matrix)

		# If this is the case, you are looking at a symmetric 2D matrix
		upper_triangle = lower_triangle.T
		upper_triangle_plus_1 = np.triu(upper_triangle,+1)
		return (lower_triangle + upper_triangle_plus_1)
	else:
		# For a 1D matrix or a non-symmetric 2D matrix, just convert your list to a numpy.array
		return np.array(matrix)