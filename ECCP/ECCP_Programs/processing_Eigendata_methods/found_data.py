'''
Geoffrey Weal, found_data.py, 9/3/22

This script provides methods for determining if the Eigendata process has completed and if an output.log file should contain eigendata.
'''

import os

MO_coefficients_filename  = 'MO_coefficients.txt'
MO_orbital_names_filename = 'MO_orbital_names.txt'
MO_occupancies_filename   = 'MO_occupancies.txt'
fort_7_filename = 'fort.7'

def found_eigendata_files(path_to_eigendata):
    """
    This method is designed to determine if the files that contain eigendata have been obtained.

    The files that this method looks for were originally contained in the output.log file

    Parameters
    ----------
    path_to_eigendata : str.
        This is the path to the eigendata files

    Returns
    -------
    found_MO_coefficients : bool.
        This boolean indicates if the MO_coefficients.txt has been found.
    found_MO_orbital_names : bool.
        This boolean indicates if the MO_orbital_names.txt has been found.
    found_MO_occupancies : bool.
        This boolean indicates if the MO_occupancies.txt has been found.
    """
    found_MO_coefficients  = MO_coefficients_filename  in os.listdir(path_to_eigendata)
    found_MO_orbital_names = MO_orbital_names_filename in os.listdir(path_to_eigendata)
    found_MO_occupancies   = MO_occupancies_filename   in os.listdir(path_to_eigendata)
    return (found_MO_coefficients, found_MO_orbital_names, found_MO_occupancies)

def should_this_calc_contain_eigendata(path_to_eigendata, log_filename):
    """
    This method will look through the output.log file to determine if this calculation is currently or has produced eigendata information. 
    """
    output_has_punch = False
    output_has_pop = False
    with open(path_to_eigendata+'/'+log_filename,'r') as outputLOG:
        for line in outputLOG:
            if (output_has_punch and output_has_pop):
                if 'punch=mo' in line.lower():
                    output_has_punch = True
                if 'pop=full' in line.lower():
                    output_has_pop = True
                if 'Input orientation:' in line:
                    if not (output_has_punch and output_has_pop):
                        return False
            else:
                if '*** Overlap ***' in line:
                    return True
    return False














