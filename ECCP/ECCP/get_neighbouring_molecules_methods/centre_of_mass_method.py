"""
centre_of_mass_method.py, Geoffrey Weal, 17/2/22

This script includes methods for obtaining neighbouring pairs of molecules using the centre of mass.
"""
import numpy as np

from SUMELF import get_cell_corner_points, get_distance
from ECCP.ECCP.get_neighbouring_molecules_methods.supporting_methods.supporting_methods import make_neighbouring_pair

def get_neighbours_centre_of_mass_method(molecules, max_neighbour_distance):
	"""
	This method will obtain neighbouring pairs of molecules in the crystal based on the distances between their centre of masses.

	Parameters
	----------
	molecules : list of ase.Atoms objects
		These are all the individual molecules identified in the crystal that you want to determine neighbouring pairs of molecules for.
	max_neighbour_distance : float.
		This is the maximum distance that the centre of mass between two molecules can be to be considered a neighbouring pair of molecules. Given in Å. 

	Returns
	-------
	neighbouring_pairs : list
		This is a list of all the neighbouring pairs of molecules identified.		
	"""

	raise Exception('May be broken, molecules not given by names but by index')

	raise Exception('Check if this method includes multiprocessing')

	# First, get the cell points of the super cell around the origin unit cell with reach 1 
	crystal_cell_lattice = molecules[0].get_cell()
	cell_points, unit_cell_displacements = get_cell_corner_points(crystal_cell_lattice, super_cell_reach=1, get_unit_cell_displacements=True)
	origin_cell_point = np.array((0,0,0))

	# Second, holding the first molecule in place, move the second molecule around the reach 1 super cell and 
	# record which the displacement required by molecule 2 to make a neighbouring pair betweem molecule 1 and molecule 2. 
	neighbouring_pairs = []
	centre_of_masses = [molecules[index].get_center_of_mass() for index in range(len(molecules))]
	for index1 in range(len(centre_of_masses)):
		cm1 = centre_of_masses[index1]
		# index2 starts at index1 rather than index1+1 to capture any neighbouring pair with itself in an adjacent unit cell.
		# We dont include any index2 that are lower than index1 because a (index1,index2) neighbouring pair will be the same 
		# as a (index2,index1) neighbouring pair, even if index1 and index2 are chemically different (i.e. two different chemicals). 
		for index2 in range(index1,len(centre_of_masses)):
			cm2 = centre_of_masses[index2]
			for displacement, unit_cell_displacement in zip(cell_points, unit_cell_displacements):
				# Do not include a neighbouring pair of a molecule with itself in the same position.
				if (index1 == index2) and (displacement == origin_cell_point).all(): 
					continue
				# If the distance between the centre of mass is less than max_neighbour_distance, the two molecules form a neighbouring pair.
				distance = round(get_distance(cm1,cm2 + displacement),4)
				if distance <= max_neighbour_distance:
					neighbouring_pair = (index1, index2, unit_cell_displacement, displacement, distance)
					neighbouring_pairs.append(neighbouring_pair)
	
	# Third, return the neighbouring_pairs list.
	return neighbouring_pairs


