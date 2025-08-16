"""
determine_number_of_permutations_of_dimers.py, Geoffrey Weal, 12/2/24

This method is designed to warn the user if the dimer comparisons proceedure will take a long time. 
"""

def determine_number_of_permutations_of_dimers(equivalent_molecule_atom_indices_comparison):
	"""
	This method is designed to warn the user if the dimer comparisons proceedure will take a long time. 

	Parameters
	----------
	equivalent_molecule_atom_indices_comparison : dict.
		This is all the ways that two molecules in non_hydrogen_graphs can map onto each other. 
	"""

	# First, use equivalent_molecule_atom_indices_comparison to determine the total number of permutations of mapping molecules onto each other that may need to be performed.  
	permutations_of_indices_in_molecules = {key: len(value) for key, value in equivalent_molecule_atom_indices_comparison.items()}

	# Second, determine the total number of permutations of mapping dimers onto each other that may need to be performed. 
	permutations_of_indices_in_dimers = {}
	for (key1a, key1b) in permutations_of_indices_in_molecules.keys():
		no_of_perm_molecule1_d1_to_d2 = permutations_of_indices_in_molecules[(key1a, key1b)]
		for (key2a, key2b) in equivalent_molecule_atom_indices_comparison.keys():
			no_of_perm_molecule2_d1_to_d2 = permutations_of_indices_in_molecules[(key2a, key2b)]
			permutations_of_indices_in_dimers[((key1a+1, key1b+1),(key2a+1, key2b+1))] = no_of_perm_molecule1_d1_to_d2 * no_of_perm_molecule2_d1_to_d2

	# Third, report any issues. Note: 4096 is 64^2
	if any([(pid > 4096) for pid in permutations_of_indices_in_dimers.values()]):
		print('----------------------------------------------')
		print('WARNING FROM INVARIENCE METHOD')
		print('The invarience method works by determining how to rotate and reflect a dimer ontop of another dimer to determine if the two dimers are variants of each other.')
		print('However in order to do this, the method needs to determine which atom index in dimer 1 goes with which atom index in dimer 2.')
		print('It is likely this process will take a long time due to the number of permutations of atom indices of dimer 2 to be rearranged to make a rotational and reflective comparison with dimer 2')
		print()
		print('Permutations of indices in molecules (reasonable max is about 64) (given also as (prime number: abundance, ...):')
		for (mol1_name, mol2_name) in sorted(permutations_of_indices_in_molecules.keys()):
			if mol1_name == mol2_name:
				no_of_perms = permutations_of_indices_in_molecules[ (mol1_name, mol2_name) ]
				print('Molecule '+str(mol1_name)+' with itself: '+str(no_of_perms)+' '+str(get_prime_number_composite(no_of_perms)))
		for (mol1_name, mol2_name) in sorted(permutations_of_indices_in_molecules.keys()):
			if mol1_name >= mol2_name: 
				continue
			no_of_perms = permutations_of_indices_in_molecules[ (mol1_name, mol2_name) ]
			print('Molecule '+str(mol1_name)+' with Molecule '+str(mol2_name)+': '+str(no_of_perms)+' '+str(get_prime_number_composite(no_of_perms)))
		print()
		print('Permutations of indices in dimers (reasonable max is about 4096) given also as (prime number: abundance, ...):')
		for (mol1_name, mol2_name) in sorted(permutations_of_indices_in_dimers.keys()):
			no_of_perms = permutations_of_indices_in_dimers[ (mol1_name, mol2_name) ]
			print('Dimer '+str(mol1_name)+' with Dimer '+str(mol2_name)+': '+str(no_of_perms)+' '+str(get_prime_number_composite(no_of_perms)))
		print('----------------------------------------------')

# ------------------------------------------------------------------------------------------------------------------------------------------

def get_prime_number_composite(value):
	"""
	This method will determine the lowest common denominators for this value from its prime numbers

	Parameters
	----------
	value : int
		This is the value to determine which prime numbers it is made of.

	Returns
	-------
	all_prime_numbers : list
		This is the number of each prime number that makes up value.
	"""

	# First setup the variable
	all_prime_numbers = {}
	current_value = float(value)
	divided_number = 2.0

	# Second, keep this algorithm going until divided_number <= current_value
	while divided_number <= current_value:
		if not (current_value % 1.0 == 0.0):
			exit('problem with get_prime_number_composite')
		while True:
			current_current_number = current_value / divided_number
			if not (current_current_number % 1.0 == 0.0):
				break
			all_prime_numbers[int(divided_number)] = all_prime_numbers.get(divided_number,0) + 1
			current_value = current_current_number
		divided_number += 1.0

	# Third, convert dict into string which is a numerically ordered dict
	all_prime_numbers_str = '('
	for prime_number, abundance in sorted(all_prime_numbers.items()):
		all_prime_numbers_str += str(prime_number)+':'+str(abundance)+', '
	all_prime_numbers_str = all_prime_numbers_str[:-2:] + ')'

	return all_prime_numbers_str

# ------------------------------------------------------------------------------------------------------------------------------------------




