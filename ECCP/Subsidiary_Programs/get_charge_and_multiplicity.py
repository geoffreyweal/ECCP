


def get_charge_and_multiplicity_from_gaussian(optimisation_filepath):
	"""

	"""
	with open(optimisation_filepath,'r') as OPTFILE:
		for line in OPTFILE:
			if ('Charge =' in line) and ('Multiplicity =' in line):
				line = line.rstrip().split()
				charge = int(line[2])
				multiplicity = int(line[5])
				return charge, multiplicity
			if 'Standard basis' in line:
				break
	raise Exception('Error: Could not find the charge or multiplicity in the gaussian output file: '+str(optimisation_filepath))


def get_charge_and_multiplicity_from_orca(optimisation_filepath):
	"""

	"""
	charge = None
	multiplicity = None
	with open(optimisation_filepath,'r') as OPTFILE:
		for line in OPTFILE:
			if ('Total Charge' in line):
				line = line.rstrip().split()
				charge = int(line[4])
				if (charge is not None) and (multiplicity is not None):
					return charge, multiplicity
			if ('Multiplicity' in line):
				line = line.rstrip().split()
				multiplicity = int(line[3])
				if (charge is not None) and (multiplicity is not None):
					return charge, multiplicity
			if 'SCF ITERATIONS' in line:
				break
	raise Exception('Error: Could not find the charge or multiplicity in the gaussian output file: '+str(optimisation_filepath))

