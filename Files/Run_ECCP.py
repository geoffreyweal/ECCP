"""
Run_ECCP.py, Geoffrey Weal, 18/2/22

This script is an example input script for the Electronic Crystal Calculation Prep program.
"""
import os, copy
from ECCP import ECCP

# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
# These are the parameters needed for the dimer method. 
# See https://github.com/geoffreyweal/ECCP for more information about these parameters

# This is the path to the xyz/cif file containing the crystal of interest. 
filepath = 'MUPMOC.xyz'

# This string points to a file that indicates are any bonds you want to ignore in the "filepath" file. If you do not want to give this, set this to None
bonds_to_ignore = None

# This is the method use to reassemble individual molecule from the crystal. 
make_molecule_method = 'component_assembly_approach'
# This dictionary include information about determining which molecules are equivalent. Required if you want to perform ATC calculations on molecules.
molecule_equivalence_method = {'method': 'invariance_method', 'type': 'combination'} 

# This is the method use to obtain dimers between molecules in the system.
make_dimer_method = {'method': 'nearest_atoms_method', 'max_dimer_distance': 8.0}
# This dictionary provides information for determining which dimers are equivalent
dimer_equivalence_method = {'method': 'invariance_method', 'type': 'combination'} 

# This dictionary includes info about how to treat the enivornment surrounding dimers (where applicable).  
environment_settings = {'include_environment_where_possible': False, 'max_environment_distance': 8.0}

# This tag indicates if you want to remove solvents from the crystal. This requires the input file to have a reference to which molecules are solvents called "SolventList"
remove_solvents = False

# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
# The following settings provide details that the ECCP will use to create folders and run the ECCP program

# The following will be added to the end of the ECCP-created folder. If you dont want this folder to have a suffix, set this to ''.
overall_folder_suffix_name = 'Default'

# The number of cpus that will be used to run the ECCP program. 
no_of_cpus = 16

# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
# The following dictionaries are required if you want to perform Gaussian/ORCA and Multiwfn jobs on your molecules and dimers

# --------------------------------------------------------------------------------------------------------------
# The following dictionaries provide the Gaussian/ORCA parameters and Submission information required for performing 
# Atomic Transition Charge (ATC) calculations on your molecules.

# This dictionary will add tags to your Gaussian/ORCA input file.
calculation_parameters_for_atomic_transition_charges = {}
calculation_parameters_for_atomic_transition_charges['calc_software']                = 'Gaussian'
calculation_parameters_for_atomic_transition_charges['mem']                          = '64GB'
calculation_parameters_for_atomic_transition_charges['method']                       = 'wB97XD'
calculation_parameters_for_atomic_transition_charges['basis']                        = '6-31+G(d,p)'
calculation_parameters_for_atomic_transition_charges['td_settings']                  = 'tda(nstates=10)'
calculation_parameters_for_atomic_transition_charges['obtain_excitation_amplitudes'] = False
calculation_parameters_for_atomic_transition_charges['temp_folder_path']             = '/tmp/wealge'
calculation_parameters_for_atomic_transition_charges['extra']                        = '# maxdisk=2TB scf=(xqc,maxcycle=512)'

# This dictionary will add tags to your submit.sl file for performing Gaussian/ORCA Calculations to get the initial ATC wfn files and to perofmr EET calculations
submission_information_for_atomic_transition_charges = {}
submission_information_for_atomic_transition_charges['cpus_per_task']    = 32
submission_information_for_atomic_transition_charges['mem']              = '68GB' # This has been set to calculation_parameters_for_atomic_transition_charges['mem'] + 12 GBs
submission_information_for_atomic_transition_charges['time']             = '10-00:00'
submission_information_for_atomic_transition_charges['partition']        = 'parallel'
submission_information_for_atomic_transition_charges['constraint']       = 'AVX'
submission_information_for_atomic_transition_charges['email']            = 'geoffreywealslurmnotifications@gmail.com'
submission_information_for_atomic_transition_charges['gaussian_version'] = 'gaussian/g16'
submission_information_for_atomic_transition_charges['python_version']   = 'python/3.8.1'

# This tag indicates if you want to obtain Gaussian/ORCA input files of molecules for performing ATC calculations. 
submission_information_for_multiwfn = {}
submission_information_for_multiwfn['cpus_per_task'] = 16
submission_information_for_multiwfn['mem']           = '32GB'
submission_information_for_multiwfn['time']          = '0-03:00'
submission_information_for_multiwfn['partition']     = 'parallel'
submission_information_for_multiwfn['constraint']    = 'AVX'
submission_information_for_multiwfn['email']         = 'geoffreywealslurmnotifications@gmail.com'

# --------------------------------------------------------------------------------------------------------------
# The following dictionaries provide the Gaussian/ORCA parameters and Submission information required for performing 
# reorganisation energy (RE) calculations on your molecules.

# This dictionary will add tags to your Gaussian/ORCA input file.
calculation_parameters_for_reorganisation_energy = {}
calculation_parameters_for_reorganisation_energy['calc_software']    = 'Gaussian'
calculation_parameters_for_reorganisation_energy['mem']              = '256GB'
calculation_parameters_for_reorganisation_energy['method']           = calculation_parameters_for_atomic_transition_charges['method']
calculation_parameters_for_reorganisation_energy['basis']            = calculation_parameters_for_atomic_transition_charges['basis']
calculation_parameters_for_reorganisation_energy['td_settings']      = calculation_parameters_for_atomic_transition_charges['td_settings']
#calculation_parameters_for_reorganisation_energy['temp_folder_path'] = calculation_parameters_for_atomic_transition_charges['temp_folder_path']
calculation_parameters_for_reorganisation_energy['extra']            = '# maxdisk=2TB scf=(xqc,maxcycle=512)'

# This tag indicates if you want to obtain Gaussian/ORCA input files of molecules for performing ATC calculations. 
submission_information_for_reorganisation_energy = {}
submission_information_for_reorganisation_energy['cpus_per_task']    = 128
submission_information_for_reorganisation_energy['mem']              = '264GB' # This has been set to calculation_parameters_for_reorganisation_energy['mem'] + 24 GBs
submission_information_for_reorganisation_energy['time']             = '10-00:00'
submission_information_for_reorganisation_energy['partition']        = 'parallel'
submission_information_for_reorganisation_energy['constraint']       = 'AVX'
submission_information_for_reorganisation_energy['email']            = 'geoffreywealslurmnotifications@gmail.com'
submission_information_for_reorganisation_energy['gaussian_version'] = 'gaussian/g16'
submission_information_for_reorganisation_energy['python_version']   = 'python/3.8.1'

# --------------------------------------------------------------------------------------------------------------
# The following dictionaries provide the Gaussian/ORCA parameters and Submission information required for performing 
# Franck-Condon (FC) (and Huang-Rhys (HR)) calculations on your molecules. 
# NOTE: If you want to use this, you also need to set Gaussian/ORCA parameters and submission information for obtaining reorganisation energies. 

# This dictionary will add tags to your Gaussian/ORCA input file.
calculation_parameters_for_franck_condon_factors = {}
calculation_parameters_for_franck_condon_factors['calc_software'] = 'Gaussian'
calculation_parameters_for_franck_condon_factors['mem']           = '8GB'
calculation_parameters_for_franck_condon_factors['method']        = calculation_parameters_for_atomic_transition_charges['method']
calculation_parameters_for_franck_condon_factors['basis']         = calculation_parameters_for_atomic_transition_charges['basis']
calculation_parameters_for_franck_condon_factors['extra']         = '# maxdisk=2TB scf=(xqc,maxcycle=512)'

# This tag indicates if you want to obtain Gaussian/ORCA input files of molecules for performing FC calculations. 
submission_information_for_franck_condon_factors = {}
submission_information_for_franck_condon_factors['cpus_per_task']    = 2
submission_information_for_franck_condon_factors['mem']              = '16GB' # This has been set to calculation_parameters_for_reorganisation_energy['mem'] + 24 GBs
submission_information_for_franck_condon_factors['time']             = '00-02:00'
submission_information_for_franck_condon_factors['partition']        = 'quicktest'
submission_information_for_franck_condon_factors['constraint']       = 'AVX'
submission_information_for_franck_condon_factors['email']            = 'geoffreywealslurmnotifications@gmail.com'
submission_information_for_franck_condon_factors['gaussian_version'] = 'gaussian/g16'
submission_information_for_franck_condon_factors['python_version']   = 'python/3.8.1'

# --------------------------------------------------------------------------------------------------------------
# The following dictionaries provide the Gaussian/ORCA parameters and Submission information required for performing 
# Electronic Energy Transfer (EET) calculations on your molecules.

# This dictionary will add tags to your Gaussian/ORCA .gjf file
calculation_parameters_for_electronic_energy_transfer = {}
calculation_parameters_for_electronic_energy_transfer['calc_software']                = 'Gaussian'
calculation_parameters_for_electronic_energy_transfer['mem']                          = '64GB'
calculation_parameters_for_electronic_energy_transfer['method']                       = calculation_parameters_for_atomic_transition_charges['method']
calculation_parameters_for_electronic_energy_transfer['basis']                        = calculation_parameters_for_atomic_transition_charges['basis']
calculation_parameters_for_electronic_energy_transfer['td_settings']                  = calculation_parameters_for_atomic_transition_charges['td_settings']
calculation_parameters_for_electronic_energy_transfer['obtain_excitation_amplitudes'] = False
#calculation_parameters_for_electronic_energy_transfer['temp_folder_path']             = calculation_parameters_for_atomic_transition_charges['temp_folder_path']
calculation_parameters_for_electronic_energy_transfer['extra']                        = '# maxdisk=2TB scf=(xqc,maxcycle=512)'

# This dictionary will add tags to your submit.sl file for performing Gaussian/ORCA Calculations to get the initial ATC wfn files and to perform EET calculations
submission_information_for_electronic_energy_transfer = {}
submission_information_for_electronic_energy_transfer['cpus_per_task']    = 32
submission_information_for_electronic_energy_transfer['mem']              = '68GB' # This has been set to calculation_parameters_for_electronic_energy_transfer['mem'] + 12 GBs
submission_information_for_electronic_energy_transfer['time']             = '10-00:00'
submission_information_for_electronic_energy_transfer['partition']        = 'parallel'
submission_information_for_electronic_energy_transfer['constraint']       = 'AVX'
submission_information_for_electronic_energy_transfer['email']            = 'geoffreywealslurmnotifications@gmail.com'
submission_information_for_electronic_energy_transfer['gaussian_version'] = 'gaussian/g16'
submission_information_for_electronic_energy_transfer['python_version']   = 'python/3.8.1'

# --------------------------------------------------------------------------------------------------------------
# The following dictionaries provide the Gaussian/ORCA parameters and Submission information required for obtaining 
# eigendata (such as overlap orbtials and molecular orbital energies and coefficients).

# This dictionary will add tags to your Gaussian/ORCA input file
calculation_parameters_for_eigendata = dict(calculation_parameters_for_electronic_energy_transfer)

# This dictionary will add tags to your submit.sl file for performing Gaussian/ORCA Calculations to obtain eigendata.
submission_information_for_eigendata = dict(submission_information_for_electronic_energy_transfer)

# --------------------------------------------------------------------------------------------------------------
# These tag provide the information to the ECCP about Gaussian/ORCA parameters and submission information for performing ATC, RE, and/or EET calculations and/or obtain eigendata (such as overlap orbtials and molecular orbital energies and coefficients).
# If you dont want to perform a task, set the appropriate take to None. For example: if you don't want to perform ATC calcs, set atc_file_creation_information = None 

# This tag indicates if you want to obtain Gaussian/ORCA input files of molecules for performing ATC calculations. 
atc_file_creation_information = (calculation_parameters_for_atomic_transition_charges, submission_information_for_atomic_transition_charges, submission_information_for_multiwfn)
# This tag indicates if you want to obtain Gaussian/ORCA input files to obtain the disorder energies of the molecules in the crystal.
re_file_creation_information = (calculation_parameters_for_reorganisation_energy, submission_information_for_reorganisation_energy)
# This tag indicates if you want to obtain Gaussian/ORCA input files to obtain the franck-condon factors (and huang-rhys factors) of the molecules in the crystal.
fc_file_creation_information = (calculation_parameters_for_franck_condon_factors, submission_information_for_franck_condon_factors)
# This tag indicates if you want to obtain Gaussian/ORCA input files of dimers for performing EET calculations. 
eet_file_creation_information = (calculation_parameters_for_electronic_energy_transfer, submission_information_for_electronic_energy_transfer)
# This tag indicates if you want to obtain Gaussian/ORCA input files of dimers for obtaining eigendata (such as overlap orbtials and molecular orbital energies and coefficients).
ict_file_creation_information = (calculation_parameters_for_eigendata, submission_information_for_eigendata)

# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
# This will run this method
ECCP(filepath, make_molecule_method=make_molecule_method, molecule_equivalence_method=molecule_equivalence_method, make_dimer_method=make_dimer_method, dimer_equivalence_method=dimer_equivalence_method, environment_settings=environment_settings, remove_solvents=remove_solvents, atc_file_creation_information=atc_file_creation_information, re_file_creation_information=re_file_creation_information, fc_file_creation_information=fc_file_creation_information, eet_file_creation_information=eet_file_creation_information, ict_file_creation_information=ict_file_creation_information, overall_folder_suffix_name=overall_folder_suffix_name, no_of_cpus=no_of_cpus)
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------