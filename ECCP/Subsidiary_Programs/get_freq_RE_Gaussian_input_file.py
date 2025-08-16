#!/usr/bin/env python3
'''
Geoffrey Weal, get_freq_RE_Gaussian_input_file.py, 21/4/22

This program is designed to  
'''
import os, sys, shutil
from ase.io import read
from copy import deepcopy

from ECCP.ECCP_Programs.shared_general_methods.shared_gaussian_methods            import did_gaussian_opt_job_complete
from ECCP.Subsidiary_Programs.can_read_data_from_checkpoint_file                  import can_read_data_from_checkpoint_file
from ECCP.Subsidiary_Programs.get_charge_and_multiplicity                         import get_charge_and_multiplicity_from_gaussian
from ECCP.ECCP.write_molecules_to_disk_methods.write_methods.gaussian_modified_RE import write_gaussian_in as write_gaussian_in_RE

def run_method(dirpath, optimisation_filename, freq_calc_filename, perform_TD, perform_raman, gaussian_parameters={}):
    """
    
    """

    # First, check that the output file has successfully completed
    did_first_opt_complete, has_fully_converged, successful_image_index, total_no_of_images = did_gaussian_opt_job_complete(dirpath+'/'+optimisation_filename, get_most_converged_image=True, get_total_no_of_images=True)
    if not did_first_opt_complete:
        raise Exception('Error in get_freq_RE_Gaussian_input_file method. Your DFT optimisation calculation in the reorganisation energy did not complete successfully. Check this and try again.')

    # Second, if there exists an oldchk file to read from, check that it exists and if it does, create the gjf file to read from it.
    gaussian_parameters_new = deepcopy(gaussian_parameters)
    read_chk_file = False
    if 'oldchk' in gaussian_parameters_new:
        if can_read_data_from_checkpoint_file(gaussian_parameters_new['oldchk']):
            read_chk_file = True
        else:
            del gaussian_parameters_new['oldchk']

    # Third, import the optimsed structure from the output.log file
    print('Extracting image '+str(total_no_of_images + successful_image_index + 1)+' from '+dirpath+'/'+optimisation_filename+' (numbering in file goes from 1 to '+str(total_no_of_images)+')')
    optimised_structure = read(dirpath+'/'+optimisation_filename,index=successful_image_index)
    charge, multiplicity = get_charge_and_multiplicity_from_gaussian(dirpath+'/'+optimisation_filename)
    gaussian_parameters_new['charge'] = charge
    gaussian_parameters_new['mult']   = multiplicity

    # Fourth, make the next input gjf file for performing the single point calculation in Gaussian
    with open(freq_calc_filename,'w') as fd:
        write_gaussian_in_RE (fd, optimised_structure, perform_opt=False, perform_CalcAll=False, perform_TD=perform_TD, perform_freq=True, perform_raman=perform_raman, perform_density=perform_TD, perform_pop=True, read_chk_file=read_chk_file, **gaussian_parameters_new)

# ----------------------------------------------------------

# First, get the dirpath and opt_filename variables
dirpath = os.getcwd()
optimisation_filename = str(sys.argv[1])
freq_calc_filename    = str(sys.argv[2])
perform_TD            = str(sys.argv[3])
perform_raman         = str(sys.argv[4])

if perform_TD == 'True':
    perform_TD = True
elif perform_TD == 'False':
    perform_TD = False
else:
    raise Exception('perform_TD must be either "True" or "False". perform_TD = '+str(perform_TD))

if perform_raman == 'True':
    perform_raman = True
elif perform_raman == 'False':
    perform_raman = False
else:
    raise Exception('perform_raman must be either "True" or "False". perform_raman = '+str(perform_raman))

# Second, get the gaussian_parameters used in the optimisation calculation to input into the single point calculation.
gaussian_parameters = eval(' '.join(sys.argv[5:]))

# Third, run method
run_method(dirpath, optimisation_filename, freq_calc_filename, perform_TD, perform_raman, gaussian_parameters=gaussian_parameters)

print('Successfully converted optimised image of '+str(optimisation_filename)+' into '+str(freq_calc_filename))

# ----------------------------------------------------------
