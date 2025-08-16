"""
make_initial_gaussian_RE_optimisation_gjf_file.py, Geoffrey Weal, 8/5/22

This method is designed to create the initial gaussian gjf files for running the RE optimisation calculation
"""
from copy                                                                         import deepcopy
from SUMELF                                                                       import check_molecule_against_file
from ECCP.ECCP.write_molecules_to_disk_methods.write_methods.gaussian_modified_RE import write_gaussian_in as write_gaussian_in_RE

def make_initial_gaussian_RE_optimisation_gjf_file(ground_structure_calc_folder, ground_structure_GS_DFT_main_opt, molecule, perform_TD, molecule_name, gaussian_parameters_GS):
	"""
	This method is designed to create the initial gaussian gjf files for running the RE optimisation calculation
	"""

	# First, make a copy of the molecule for creating files from, and set the 
	#        periodic boundary condition (pbc) for the system as False. 
	molecule_for_input = molecule.copy()
	molecule_for_input.set_pbc(False)

	# Second, create the gaussian .gjf file for optimising the ground structure.
	if ('pre_method' in gaussian_parameters_GS) or ('pre_basis' in gaussian_parameters_GS):

		# Third, Make changes to gaussian_parameters to change gaussian filenames to include "main_opt" in their name.
		gaussian_parameters_GS_copy = deepcopy(gaussian_parameters_GS)
		for gaussian_file in ['chk', 'd2e', 'int', 'rwf', 'skr']:
			if gaussian_file in gaussian_parameters_GS_copy:
				gaussian_filename = gaussian_parameters_GS_copy[gaussian_file]
				gaussian_filename = gaussian_filename.split('.')
				gaussian_filename = gaussian_filename[0]+'_main_preopt.'+gaussian_filename[1]
				gaussian_parameters_GS_copy[gaussian_file] = gaussian_filename

		# Fourth, create the path to the gaussian file.
		full_path_to_gjf_file = ground_structure_calc_folder+'/'+ground_structure_GS_DFT_main_opt.replace('.gjf','_preopt.gjf')

		# Fifth, check if there already exists this molecule on file, check if the molecules are the same.
		check_molecule_against_file(molecule_for_input, full_path_to_gjf_file)

		# Sixth, create the pre-optimisation gjf file
		with open(full_path_to_gjf_file, 'w') as fd:	
			initial_gaussian_parameters_GS = deepcopy(gaussian_parameters_GS_copy)
			if ('pre_method' in gaussian_parameters_GS_copy):
				initial_gaussian_parameters_GS['method'] = initial_gaussian_parameters_GS.pop('pre_method')
				del gaussian_parameters_GS_copy['pre_method']
			if ('pre_basis'  in gaussian_parameters_GS_copy):
				initial_gaussian_parameters_GS['basis']  = initial_gaussian_parameters_GS.pop('pre_basis')
				del gaussian_parameters_GS_copy['pre_basis']
			write_gaussian_in_RE(fd, molecule_for_input, perform_opt=True, perform_CalcAll=False, perform_TD=perform_TD, perform_freq=False, perform_raman=False, perform_density=perform_TD, perform_pop=perform_TD, read_chk_file=False, molecule_name=molecule_name, **initial_gaussian_parameters_GS)
		
		# Seventh, indicate that you are running a preoptimisation first
		run_pre_optimisation = True

	else:

		# Eighth, make changes to gaussian_parameters to change gaussian filenames to include "main_opt" in their name.
		gaussian_parameters_GS_copy = deepcopy(gaussian_parameters_GS)
		for gaussian_file in ['chk', 'd2e', 'int', 'rwf', 'skr']:
			if gaussian_file in gaussian_parameters_GS_copy:
				gaussian_filename = gaussian_parameters_GS_copy[gaussian_file]
				gaussian_filename = gaussian_filename.split('.')
				gaussian_filename = gaussian_filename[0]+'_main_opt.'+gaussian_filename[1]
				gaussian_parameters_GS_copy[gaussian_file] = gaussian_filename

		# Ninth, create the path to the gaussian file.
		full_path_to_gjf_file = ground_structure_calc_folder+'/'+ground_structure_GS_DFT_main_opt

		# Tenth, check if there already exists this molecule on file, check if the molecules are the same.
		check_molecule_against_file(molecule_for_input, full_path_to_gjf_file)

		# Eleventh, create the optimisation gjf file
		with open(full_path_to_gjf_file, 'w') as fd:	
			write_gaussian_in_RE(fd, molecule_for_input, perform_opt=True, perform_CalcAll=False, perform_TD=perform_TD, perform_freq=False, perform_raman=False, perform_density=perform_TD, perform_pop=perform_TD, read_chk_file=False, molecule_name=molecule_name, **gaussian_parameters_GS_copy)
		
		# Twelfth, indicate that you are not running a preoptimisation first
		run_pre_optimisation = False

	# Thirteenth, return run_pre_optimisation
	return run_pre_optimisation

# ------------------------------------------------------------------------------------------------------------------------------
