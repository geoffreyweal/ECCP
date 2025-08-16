"""
prepare_gau_and_submit_information.py, Geoffrey Weal, 19/2/22

This script is designed to hold methods useful for many components of this program.
"""

from copy import deepcopy

# ------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------

def prepare_ATC_calc_params_and_submit_information(atc_file_creation_information):
	"""
	This method is designed to check the gaussian_parameters and submission_information parameters used for performing Atomic Transition Charge (ATC) calculations and turn them into lists for the rest of the program.

	Parameters
	----------
	atc_file_creation_information : list of dict. or None
		This variable contains the Gaussian Parameters and Submission information required for performing Atomic Transition Charge (ATC) calculations in Gaussian and Multiwfn. Set to None if you do not want to perform ATC calculations. Default: None
	
	Attributes
	----------
	gaussian_parameters_for_ATCs : dict./list of dict.
		This contains all the information required by gaussian to perform atomic transition charge (ATC) calculations, including the functional and basis set. 
	submission_information_for_ATCs : dict./list of dict.
		This contains all the information required for the submit.sl file to perform atomic transition charge (ATC) calculations, required for submitting gaussian jobs to slurm. 
	submission_information_for_multiwfn : dict./list of dict.
		This contains all the information required for the multiwfn_submit.sl file, required for submitting multiwfn jobs to slurm to get atc chg files. 

	Returns
	-------
	all_gaussian_parameters_for_ATCs : list
		This list contain all the information required to be included in the Gaussian file(s).
	all_gaussian_parameters_for_multiwfn : list
		This list contain all the information required to be included in the MultiWFN file(s).
	all_submission_information_for_ATCs : list
		This list contain all the information required to be included in the submit.sl script to submit Gaussian jobs.
	all_submission_information_for_multiwfn : list
		This list contain all the information required to be included in the multiwfn_submit.sl script to submit Multiwfn jobs.
	get_molecule_atcs : bool.
		True means the user wants to perform ATC calculations.
	"""

	# First, if atc_file_creation_information is false, set all to None and move on, we dont want to perform gaussian jobs.
	if atc_file_creation_information is None:
		return None, None, None, None, False

	# Second, obtain all the dictionaries from atc_file_creation_information, and set get_molecule_atcs to True
	gaussian_parameters_for_ATCs, submission_information_for_ATCs, submission_information_for_multiwfn = atc_file_creation_information
	get_molecule_atcs = True

	# Third, verify and format the gaussian and submission parameters for ATC calcs.
	all_gaussian_parameters_for_ATCs, all_submission_information_for_ATCs = check_calc_and_submission_variables(gaussian_parameters_for_ATCs, submission_information_for_ATCs, job_type='ATC')
	check_defaults_for_all_calc_parameters(all_gaussian_parameters_for_ATCs)
	
	# Fourth, verify and format the gaussian and submission parameters for multiwfn calcs. 
	all_gaussian_parameters_for_multiwfn, all_submission_information_for_multiwfn = check_calc_and_submission_variables(gaussian_parameters_for_ATCs, submission_information_for_multiwfn, job_type='ATC MultiWFN', no_quantum_chemistry_computing_program_required=True)

	# ---------------------------------------------------------------------------------------------------------------------
	# Fifth, check to make sure their are no issues with list of dicts.
	check_gaussian_list = []
	check_submission_list = []
	if all_gaussian_parameters_for_ATCs is not None:
		check_gaussian_list.append(len(all_gaussian_parameters_for_ATCs))
		check_submission_list.append(len(all_submission_information_for_ATCs))
	if all_gaussian_parameters_for_multiwfn is not None:
		check_gaussian_list.append(len(all_gaussian_parameters_for_multiwfn))
		check_submission_list.append(len(all_submission_information_for_multiwfn))

	if ((all_gaussian_parameters_for_ATCs is not None) or (all_gaussian_parameters_for_multiwfn is not None)) and (not len(set(check_gaussian_list+check_submission_list)) == 1):
		print('Error: The lengths of each gaussian_parameters lists and submission_information lists are not the same')
		print('Number of entries in ')
		print('\t all_gaussian_parameters_for_ATCs: '+str(len(all_gaussian_parameters_for_ATCs)))
		print('\t all_gaussian_parameters_for_multiwfn: '+str(len(all_gaussian_parameters_for_multiwfn)))
		print('\t all_submission_information_for_ATCs: '+str(len(all_submission_information_for_ATCs)))
		print('\t all_submission_information_for_multiwfn: '+str(len(all_submission_information_for_multiwfn)))
		print()
		print('Check your lists:')
		print()
		print('\t all_gaussian_parameters_for_ATCs: '+str((all_gaussian_parameters_for_ATCs)))
		print('\t all_gaussian_parameters_for_multiwfn: '+str((all_gaussian_parameters_for_multiwfn)))
		print('\t all_submission_information_for_ATCs: '+str((all_submission_information_for_ATCs)))
		print('\t all_submission_information_for_multiwfn: '+str((all_submission_information_for_multiwfn)))
		print()
		exit('This program will finish without beginning')
	# ---------------------------------------------------------------------------------------------------------------------

	return deepcopy(all_gaussian_parameters_for_ATCs), deepcopy(all_gaussian_parameters_for_multiwfn), deepcopy(all_submission_information_for_ATCs), deepcopy(all_submission_information_for_multiwfn), get_molecule_atcs

# ------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------

def prepare_RE_calc_params_and_submit_information(re_file_creation_information):
	"""
	This method is designed to check the gaussian_parameters and submission_information parameters used for performing reorganisation energy (RE) calculations and turn them into lists for the rest of the program.

	Parameters
	----------
	re_file_creation_information : list of dict. or None
		This variable contains the Gaussian Parameters and Submission information required for performing reorganisation energy (RE) calculations in Gaussian and Multiwfn. Set to None if you do not want to perform RE calculations. Default: None
	
	Attributes
	----------
	gaussian_parameters_for_REs : dict./list of dict.
		This contains all the information required by gaussian to perform reorganisation energy (RE) calculations, including the functional and basis set. 
	submission_information_for_REs : dict./list of dict.
		This contains all the information required for the submit.sl file to perform reorganisation energy (RE) calculations, required for submitting gaussian jobs to slurm. 

	Returns
	-------
	all_gaussian_parameters_for_REs : list
		This list contain all the information required to be included in the Gaussian file(s) for obtaining reorganisation energies of molecules. 
	all_submission_information_for_REs : list
		This list contain all the information required to be included in the submit.sl script to submit Gaussian jobs for obtaining reorganisation energies of molecules. 
	get_molecule_res : bool.
		True means the user wants to perform RE calculations.
	"""

	# ---------------------------------------------------------------------------------------------------------------------
	# First, re_file_creation_information is None, set all to None and move on, we dont want to perform gaussian jobs
	if re_file_creation_information is None:
		return None, None, False

	# ---------------------------------------------------------------------------------------------------------------------
	# Second, obtain all the dictionaries from re_file_creation_information, and set get_molecule_res to True
	gaussian_parameters_for_REs, submission_information_for_REs = re_file_creation_information
	get_molecule_res = True

	# ---------------------------------------------------------------------------------------------------------------------
	# Third, verify and format the gaussian and submission parameters for reorganisation energy calcs. 
	all_gaussian_parameters_for_REs, all_submission_information_for_REs = check_calc_and_submission_variables(gaussian_parameters_for_REs, submission_information_for_REs, job_type='RE')
	check_defaults_for_all_calc_parameters_for_RE(all_gaussian_parameters_for_REs)
	# ---------------------------------------------------------------------------------------------------------------------

	# Fourth, check to make sure their are no issues with list of dicts.
	check_gaussian_list = []
	check_submission_list = []
	if all_gaussian_parameters_for_REs is not None:
		check_gaussian_list.append(len(all_gaussian_parameters_for_REs))
		check_submission_list.append(len(all_submission_information_for_REs))

	if (all_gaussian_parameters_for_REs is not None) and (not len(set(check_gaussian_list+check_submission_list)) == 1):
		print('Error: The lengths of each gaussian_parameters lists and submission_information lists are not the same')
		print('Number of entries in ')
		print('\t all_gaussian_parameters_for_RE: '+str(len(all_gaussian_parameters_for_REs)))
		print('\t all_submission_information_for_RE: '+str(len(all_submission_information_for_REs)))
		print()
		print('Check your lists:')
		print()
		print('\t all_gaussian_parameters_for_RE: '+str((all_gaussian_parameters_for_REs)))
		print('\t all_submission_information_for_RE: '+str((all_submission_information_for_REs)))
		print()
		exit('This program will finish without beginning')
	# ---------------------------------------------------------------------------------------------------------------------

	# Fifth, return quantities.
	return deepcopy(all_gaussian_parameters_for_REs), deepcopy(all_submission_information_for_REs), get_molecule_res

# ------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------

def prepare_FC_calc_params_and_submit_information(fc_file_creation_information):
	"""
	This method is designed to check the gaussian_parameters and submission_information parameters used for performing Franck-Condon/Huang-Rhys (FC/HR) calculations and turn them into lists for the rest of the program.

	Parameters
	----------
	fc_file_creation_information : list of dict. or None
		This variable contains the Gaussian Parameters and Submission information required for performing Franck-Condon/Huang-Rhys (FC/HR) calculations in Gaussian. Set to None if you do not want to perform FC calculations. Default: None
	
	Attributes
	----------
	gaussian_parameters_for_FCs : dict./list of dict.
		This contains all the information required by gaussian to perform Franck-Condon/Huang-Rhys (FC/HR) calculations, including the functional and basis set. 
	submission_information_for_FCs : dict./list of dict.
		This contains all the information required for the submit.sl file to perform Franck-Condon/Huang-Rhys (FC/HR) calculations, required for submitting gaussian jobs to slurm. 

	Returns
	-------
	all_gaussian_parameters_for_FCs : list
		This list contain all the information required to be included in the Gaussian file(s) for obtaining Franck-Condon/Huang-Rhys (FC/HR) factors of molecules. 
	all_submission_information_for_FCs : list
		This list contain all the information required to be included in the submit.sl script to submit Gaussian jobs for obtaining Franck-Condon/Huang-Rhys (FC/HR) factors of molecules. 
	get_molecule_res : bool.
		True means the user wants to perform FC calculations.
	"""

	# ---------------------------------------------------------------------------------------------------------------------
	# First, re_file_creation_information is None, set all to None and move on, we dont want to perform gaussian jobs
	if fc_file_creation_information is None:
		return None, None, False

	# ---------------------------------------------------------------------------------------------------------------------
	# Second, obtain all the dictionaries from re_file_creation_information, and set get_molecule_res to True
	gaussian_parameters_for_FCs, submission_information_for_FCs = fc_file_creation_information
	get_molecule_fcs = True

	# ---------------------------------------------------------------------------------------------------------------------
	# Third, verify and format the gaussian and submission parameters for reorganisation energy calcs. 
	all_gaussian_parameters_for_FCs, all_submission_information_for_FCs = check_calc_and_submission_variables(gaussian_parameters_for_FCs, submission_information_for_FCs, job_type='FC')
	check_defaults_for_all_calc_parameters_for_RE(all_gaussian_parameters_for_FCs)
	# ---------------------------------------------------------------------------------------------------------------------

	# Fourth, check to make sure their are no issues with list of dicts.
	check_gaussian_list = []
	check_submission_list = []
	if all_gaussian_parameters_for_FCs is not None:
		check_gaussian_list.append(len(all_gaussian_parameters_for_FCs))
		check_submission_list.append(len(all_submission_information_for_FCs))

	if (all_gaussian_parameters_for_FCs is not None) and (not len(set(check_gaussian_list+check_submission_list)) == 1):
		print('Error: The lengths of each gaussian_parameters lists and submission_information lists are not the same')
		print('Number of entries in ')
		print('\t all_gaussian_parameters_for_FC: '+str(len(all_gaussian_parameters_for_FCs)))
		print('\t all_submission_information_for_FC: '+str(len(all_submission_information_for_FCs)))
		print()
		print('Check your lists:')
		print()
		print('\t all_gaussian_parameters_for_FC: '+str((all_gaussian_parameters_for_FCs)))
		print('\t all_submission_information_for_FC: '+str((all_submission_information_for_FCs)))
		print()
		exit('This program will finish without beginning')
	# ---------------------------------------------------------------------------------------------------------------------

	# Fifth, return quantities.
	return deepcopy(all_gaussian_parameters_for_FCs), deepcopy(all_submission_information_for_FCs), get_molecule_fcs

# ------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------

def prepare_dimer_calc_params_and_submit_information(eet_file_creation_information, job_type):
	"""
	This method is designed to check the gaussian_parameters and submission_information parameters used for performing electronic energy transfer (EET) calculations and turn them into lists for the rest of the program.

	Parameters
	----------
	eet_file_creation_information : list of dict. or None
		This variable contains the Gaussian Parameters and Submission information required for performing electronic energy transfer (EET) calculations in Gaussian and Multiwfn. Set to None if you do not want to perform EET calculations. Default: None
	job_type : str.
		This is the job type for the job being checked.

	Attributes
	----------
	gaussian_parameters_for_electronic_energy_transfer : dict./list of dict.
		This contains all the information required by gaussian to perform electronic energy transfer (EET) calculations, including the functional and basis set. 
	submission_information_for_electronic_energy_transfer : dict./list of dict.
		This contains all the information required for the submit.sl file to perform electronic energy transfer (EET) calculations, required for submitting gaussian jobs to slurm. 

	Returns
	-------
	all_gaussian_parameters_for_Dimers : list
		This list contain all the information required to be included in the Gaussian file(s).
	all_submission_information_for_Dimers : list
		This list contain all the information required to be included in the submit.sl script to submit Gaussian jobs.
	get_molecule_eets : bool.
		True means the user wants to perform EET calculations.

	"""

	# First, if eet_file_creation_information is None, set all to None and move on, we dont want to perform gaussian jobs
	if eet_file_creation_information is None:
		return None, None, False

	# Second, obtain all the dictionaries from eet_file_creation_information, and set get_molecule_eets to True
	gaussian_parameters_for_Dimers, submission_information_for_Dimers = eet_file_creation_information
	get_molecule_eets = True

	# Third, verify and format the gaussian and submission parameters for ATC and EET calcs.
	all_gaussian_parameters_for_Dimers, all_submission_information_for_Dimers = check_calc_and_submission_variables(gaussian_parameters_for_Dimers, submission_information_for_Dimers, job_type=job_type)
	check_defaults_for_all_calc_parameters(all_gaussian_parameters_for_Dimers)

	# ---------------------------------------------------------------------------------------------------------------------
	# Fourth, check to make sure their are no issues with list of dicts.
	check_gaussian_list = []
	check_submission_list = []
	if all_gaussian_parameters_for_Dimers is not None:
		check_gaussian_list.append(len(all_gaussian_parameters_for_Dimers))
		check_submission_list.append(len(all_submission_information_for_Dimers))

	if (all_gaussian_parameters_for_Dimers is not None) and (not len(set(check_gaussian_list+check_submission_list)) == 1):
		print('Error: The lengths of each gaussian_parameters lists and submission_information lists are not the same')
		print('Number of entries in ')
		print('\t all_gaussian_parameters_for_Dimers: '+str(len(all_gaussian_parameters_for_Dimers)))
		print('\t all_submission_information_for_Dimers: '+str(len(all_submission_information_for_Dimers)))
		print()
		print('Check your lists:')
		print()
		print('\t all_gaussian_parameters_for_Dimers: '+str((all_gaussian_parameters_for_Dimers)))
		print('\t all_submission_information_for_Dimers: '+str((all_submission_information_for_Dimers)))
		print()
		exit('This program will finish without beginning')
	# ---------------------------------------------------------------------------------------------------------------------

	# Fifth, return quantities.
	return deepcopy(all_gaussian_parameters_for_Dimers), deepcopy(all_submission_information_for_Dimers), get_molecule_eets

# ------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------

def check_calc_and_submission_variables(gaussian_parameters, submission_information, job_type, no_quantum_chemistry_computing_program_required=False):
	"""
	This method will check the Gaussian parameters dictionary against the submission information to make sure they are consistant with each other.

	This method will also turn these dictionaries into list of dictionaries that can be used in the ECCP program.

	This conversion to lists of dictionaries allows you to do things like use the same submission_information dictionary to perform calculations with the same dimersw using a range of different basis set and functionals, given by a number of gaussian_parameters is a list.

	Or, you make be using different submission_information for your different gaussian_parameters.

	Parameters
	----------
	gaussian_parameters : dict. or list of dicts.
		This is a single dictionary or a list of many dictionaries that contain information about making the Gaussian input (.gjf) files.
	submission_information : dict. or list of dicts.
		This is a single dictionary or a list of many dictionaries that contain information about making the submit.sl files for submitting jobs to slurm. 
	job_type : str.
		This is the job type for the job being checked.
	no_quantum_chemistry_computing_program_required : bool.
		This tag indicates if this calculation is being performed with a quantum chemistry computing program or not. Default: False

	Returns
	-------
	all_gaussian_parameters : list of dicts.
		This is a list of dictionaries that contain information about making the Gaussian input (.gjf) files.
	all_submission_information : list of dicts.
		This is a list of dictionaries that contain information about making the submit.sl files for submitting jobs to slurm. 
	"""

	# First, check to see if the gaussian_parameters variable has been given in the correct format
	# While doing these checks, also turn gaussian_parameters into an appropriate list called all_gaussian_parameters.
	an_error = False
	if (gaussian_parameters is None) or (gaussian_parameters == []):
		all_gaussian_parameters = None
		submission_information = None
	elif isinstance(gaussian_parameters,dict):
		all_gaussian_parameters = [gaussian_parameters]
	elif (isinstance(gaussian_parameters,list) or isinstance(gaussian_parameters,tuple)):
		all_gaussian_parameters = gaussian_parameters
		if not all([isinstance(gaussian_parameters,dict) for gaussian_parameters in all_gaussian_parameters]):
			an_error = True
	else:
		an_error = True
	if an_error:
		print('Error: the gaussian_parameters variable must be a dictionary or a list of dictionaries')
		print('gaussian_parameters: '+str(gaussian_parameters))
		print('check this out')
		exit('This program will finish without beginning')

	# Second, check to see if the submission_information variable has been given in the correct format
	# While doing these checks, also turn submission_information into an appropriate list called all_submission_information.
	an_error = False
	if (submission_information is None) or (submission_information == []):
		all_submission_information = None
	elif isinstance(submission_information,dict):
		all_submission_information = [submission_information]*len(all_gaussian_parameters)
	elif (isinstance(submission_information,list) or isinstance(submission_information,tuple)):
		all_submission_information = submission_information
		if not all([isinstance(submission_information,dict) for submission_information in all_submission_information]):
			an_error = True
	else:
		an_error = True
	if an_error:
		print('Error: the submission_information variable must be a dictionary or a list of dictionaries')
		print('submission_information: '+str(submission_information))
		print('check this out')
		exit('This program will finish without beginning')

	# Third, check that Gaussian and ORCA settings are given appropriately
	if not no_quantum_chemistry_computing_program_required:
		for gaussian_parameters, submission_information in zip(all_gaussian_parameters, all_submission_information):
			if 'calc_software' not in gaussian_parameters:
				toString = "Error in "+str(job_type)+" checks: One of more of your entries for all_gaussian_parameters does not contain an entry for 'calc_software', which is needed.\n 'calc_software' must be either 'Gaussian' or 'ORCA'.\n"
				toString += 'Your all_gaussian_parameters list is:\n'
				for gaussian_parameters in all_gaussian_parameters:
					toString += str(gaussian_parameters)+'\n'
				raise Exception(toString)
			if ('gaussian_version' not in submission_information) and ('orca_version' not in submission_information):
				toString = "Error in "+str(job_type)+" checks: One of more of your entries for all_submission_information does not contain an entry for 'gaussian_version' or 'orca_version', which is needed fur running either Gaussian or ORCA calculations respectively.\n"
				toString += 'Your all_submission_information list is:\n'
				for submission_information in all_submission_information:
					toString += str(submission_information)+'\n'
				raise Exception(toString)
			if   gaussian_parameters['calc_software'].lower() == 'gaussian':
				if 'gaussian_version' not in submission_information:
					toString = "Error in "+str(job_type)+" checks: You are requesting to use Gaussian in an entry in gaussian_parameters, however you have not given an input for 'gaussian_version' in submission_information\n"
					toString += 'Your all_gaussian_parameters, all_submission_information list is:\n'
					toString += 'all_gaussian_parameters\t|\tall_submission_information:\n'
					for gaussian_parameters, submission_information in zip(all_gaussian_parameters, all_submission_information):
						toString += str(gaussian_parameters)+'\t|\t'+str(submission_information)+'\n'
					raise Exception(toString)
			elif gaussian_parameters['calc_software'].lower() == 'orca':
				if 'orca_version' not in submission_information:
					toString = "Error in "+str(job_type)+" checks: You are requesting to use ORCA in an entry in gaussian_parameters, however you have not given an input for 'orca_version' in submission_information\n"
					toString += 'Your all_gaussian_parameters, all_submission_information list is:\n'
					toString += 'all_gaussian_parameters\t|\tall_submission_information:\n'
					for gaussian_parameters, submission_information in zip(all_gaussian_parameters, all_submission_information):
						toString += str(gaussian_parameters)+'\t|\t'+str(submission_information)+'\n'
					raise Exception(toString)
			else:
				toString = "Error in "+str(job_type)+" checks: One of more of your entries for all_gaussian_parameters has an incorrect inpit for 'calc_software'.\n 'calc_software' must be either 'Gaussian' or 'ORCA'.\n"
				toString += 'Your all_gaussian_parameters list is:\n'
				for gaussian_parameters in all_gaussian_parameters:
					toString += str(gaussian_parameters)+'\n'
				raise Exception(toString)

	# Fourth, return dictionaries.
	return all_gaussian_parameters, all_submission_information

# ------------------------------------------------------------------------------------------------------------------------------

def check_defaults_for_all_calc_parameters(all_gaussian_parameters):
	"""
	This method will check the gaussian_parameters in all_gaussian_parameters to make sure they contain certain inputs. 

	If they don't, will provide the default values for those inputs.

	Parameters
	----------
	gaussian_parameters : dict. or list of dicts.
		This is a single dictionary or a list of many dictionaries that contain information about making the Gaussian input (.gjf) files.
	"""
	for index in range(len(all_gaussian_parameters)):
		gaussian_parameters = all_gaussian_parameters[index]
		if not 'td_settings' in gaussian_parameters.keys():
			gaussian_parameters['td_settings'] = 'td'
		if not 'obtain_excitation_amplitudes' in gaussian_parameters.keys():
			gaussian_parameters['obtain_excitation_amplitudes'] = False


def check_defaults_for_all_calc_parameters_for_RE(all_gaussian_parameters_for_RE):
	"""
	This method will check the gaussian_parameters in all_gaussian_parameters to make sure they contain certain inputs for performing reorganisation energy (RE) calculations. 

	If they don't, will provide the default values for those inputs.

	Parameters
	----------
	all_gaussian_parameters_for_RE : dict. or list of dicts.
		This is a single dictionary or a list of many dictionaries that contain information about making the Gaussian input (.gjf) files for performing reorganisation energy (RE) calculations. 
	"""
	for index in range(len(all_gaussian_parameters_for_RE)):
		gaussian_parameters = all_gaussian_parameters_for_RE[index]
		if not 'td_settings' in gaussian_parameters.keys():
			gaussian_parameters['td_settings'] = 'td'

# ------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------









