#!/usr/bin/env python3
'''
Geoffrey Weal, has_optimisation_converged.py, 16/6/23

This program is designed to  
'''
import os, sys
from ECCP.ECCP_Programs.shared_general_methods.shared_gaussian_methods import did_gaussian_opt_job_complete

# First, get the dirpath and opt_filename variables
dirpath = os.getcwd()
optimisation_filename = str(sys.argv[1])

# Second, determine if the calculation converged successfully or not.
did_finish_successfully, has_fully_converged, successful_image_index, total_no_of_images = did_gaussian_opt_job_complete(dirpath+'/'+optimisation_filename, get_most_converged_image=True, get_total_no_of_images=True)

# Third, return the result to the terminal.
print(int(not did_finish_successfully))
