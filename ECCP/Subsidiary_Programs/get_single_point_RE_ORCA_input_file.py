#!/usr/bin/env python3
'''
Geoffrey Weal, get_single_point_RE_ORCA_input_file.py, 21/4/22

This program is designed to  
'''
import os, sys, shutil
from ase.io import read
from copy import deepcopy

from ECCP.ECCP_Programs.shared_general_methods.shared_orca_methods import did_orca_opt_job_complete
from ECCP.Subsidiary_Programs.get_charge_and_multiplicity          import get_charge_and_multiplicity_from_orca
from ECCP.ECCP.write_molecules_to_disk_methods.orca_modified_RE    import write_orca_in_RE

def run_method(dirpath, optimisation_filename, single_point_filename, perform_TD, orca_parameters={}):
    """
    
    """

    # First, check that the output file has successfully completed
    did_first_opt_complete, has_fully_converged, successful_image_index, total_no_of_images = did_orca_opt_job_complete(dirpath+'/'+optimisation_filename, get_most_converged_image=True, get_total_no_of_images=True)
    if not did_first_opt_complete:
        raise Exception('Error in get_single_point_RE_ORCA_input_file method. Your DFT optimisation calculation in the reorganisation energy did not complete successfully. Check this and try again.')

    # Second, we do not want to use the previous chk file for this calculation.
    orca_parameters_new = deepcopy(orca_parameters)

    # Third, import the optimsed structure from the output.log file
    print('Extracting image '+str(total_no_of_images + successful_image_index + 1)+' from '+dirpath+'/'+optimisation_filename+' (numbering in file goes from 1 to '+str(total_no_of_images)+')')
    if False: #os.path.exists(dirpath+'/'+optimisation_filename):
        optimised_structure = read(dirpath+'/'+optimisation_filename,index=successful_image_index)
    else:
        optimised_structure = read(dirpath+'/'+optimisation_filename.replace('.out','_trj.xyz'),index=successful_image_index)
    charge, multiplicity = get_charge_and_multiplicity_from_orca(dirpath+'/'+optimisation_filename)
    orca_parameters_new['charge'] = charge
    orca_parameters_new['mult']   = multiplicity

    # Fourth, make the next input gjf file for performing the single point calculation in ORCA
    with open(single_point_filename,'w') as fd:
        perform_pop = perform_TD
        write_orca_in_RE(fd, optimised_structure, perform_opt=False, perform_CalcAll=False, perform_TD=perform_TD, perform_freq=False, perform_raman=False, perform_pop=perform_pop, **orca_parameters_new)

# ----------------------------------------------------------

# First, get the dirpath and opt_filename variables
dirpath = os.getcwd()
optimisation_filename = str(sys.argv[1])
single_point_filename = str(sys.argv[2])
perform_TD            = str(sys.argv[3]) 

if perform_TD == 'True':
    perform_TD = True
elif perform_TD == 'False':
    perform_TD = False
else:
    raise Exception('perform_TD must be either "True" or "False". perform_TD = '+str(perform_TD))

# Second, get the orca_parameters used in the optimisation calculation to input into the single point calculation.
orca_parameters = eval(' '.join(sys.argv[4:]))

# Third, run method
run_method(dirpath, optimisation_filename, single_point_filename, perform_TD, orca_parameters=orca_parameters)

print('Successfully converted optimised image of '+str(optimisation_filename)+' into '+str(single_point_filename))

# ----------------------------------------------------------
