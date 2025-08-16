#!/usr/bin/env python3
'''
Geoffrey Weal, get_optimisation_RE_Gaussian_input_file.py, 21/4/22

This program is designed to  
'''
import os, sys, shutil
from ase.io import read
from copy import deepcopy

from ECCP.ECCP_Programs.shared_general_methods.shared_gaussian_methods            import did_gaussian_opt_job_complete
from ECCP.Subsidiary_Programs.can_read_data_from_checkpoint_file                  import can_read_data_from_checkpoint_file
from ECCP.Subsidiary_Programs.get_charge_and_multiplicity                         import get_charge_and_multiplicity_from_gaussian
from ECCP.ECCP.write_molecules_to_disk_methods.write_methods.gaussian_modified_RE import write_gaussian_in as write_gaussian_in_RE

def run_method(dirpath, optimisation_filename, single_point_filename, perform_TD, gaussian_parameters={}):
    """
    
    """

    # First, check that the output file has successfully completed
    did_finish_successfully, has_fully_converged, successful_image_index, total_no_of_images = did_gaussian_opt_job_complete(dirpath+'/'+optimisation_filename, get_most_converged_image=True, get_total_no_of_images=True)
    if not did_finish_successfully:
        raise Exception('Error in get_optimisation_RE_Gaussian_input_file method. Your DFT optimisation calculation in the reorganisation energy did not complete successfully. Check this and try again.')

    # Second, we do not want to use the previous chk file for this calculation.
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
    with open(single_point_filename,'w') as fd:
        perform_pop = perform_TD
        write_gaussian_in_RE(fd, optimised_structure, perform_opt=True, perform_CalcAll=False, perform_TD=perform_TD, perform_freq=False, perform_raman=False, perform_density=perform_TD, perform_pop=perform_pop, read_chk_file=False, **gaussian_parameters_new)

# ----------------------------------------------------------

# First, get the dirpath and opt_filename variables
dirpath = os.getcwd()
optimisation_filename = str(sys.argv[1])
single_point_filename = str(sys.argv[2])

# Second, get the gaussian_parameters used in the optimisation calculation to input into the single point calculation.
gaussian_parameters = eval(' '.join(sys.argv[3:]))

# Third, determine if this calculation will be run using TD-DFT
perform_TD = 'td_settings' in gaussian_parameters

# Fourth, run method
run_method(dirpath, optimisation_filename, single_point_filename, perform_TD=perform_TD, gaussian_parameters=gaussian_parameters)

print('Successfully converted optimised image of '+str(optimisation_filename)+' into '+str(single_point_filename))

# ----------------------------------------------------------
