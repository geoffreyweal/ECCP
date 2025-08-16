'''
Geoffrey Weal, processing_OPV_Dimer_data.py, 9/3/22

This script contains methods for processing_OPV_Dimer_data.py

'''
import os

# Constants and conversions that are useful for processing data.
planks_constant = 4.135667696 * (10.0 ** -15.0) #eVs
speed_of_light = 2.99792458 * (10.0 ** 8.0) #ms-1
metre_to_centimetre = 100.0 #cm m-1
eV_to_meV = 1000.0 #meV eV-1
eV_to_micro_eV = 1000000.0 #micro eV eV-1
eV_to_inverse_cm = 1/(metre_to_centimetre * speed_of_light * planks_constant) #cm eV-1

# -----------------------------------------------------------------

hartree_to_eV = 27.211386245988
def convert_hartree_to_eV(hartree_energy):
    return hartree_energy * hartree_to_eV

def get_hartree_to_eV_conversion_value():
    return hartree_to_eV

def get_energy_diff(final_energy, initial_energy):
    return final_energy - initial_energy

def get_reorganisation_energy(donor_eGS_gGS_energy, acceptor_eES_gGS_energy, donor_eGS_gES_energy, acceptor_eES_gES_energy):
    # See https://pubs.rsc.org/en/content/articlepdf/2021/tc/d0tc05697a
    reorganisation_energy = (acceptor_eES_gGS_energy - acceptor_eES_gES_energy) + (donor_eGS_gES_energy - donor_eGS_gGS_energy)
    return reorganisation_energy

# -----------------------------------------------------------------

def found_a_re_jobset(directory_path):
    """
    This method will look through the folder for the ground_structure and excited_state folders.

    Parameters
    ----------
    directory_path : str.
        This is the path to look for the ground_structure and excited_structure folders.

    Returns
    -------
    True if this calc is an eet calc, False if not.
    """

    found_ground_structure_folder  = ('ground_structure'  in os.listdir(directory_path)) and os.path.isdir(directory_path+'/'+'ground_structure')
    found_excited_structure_folder = ('excited_structure' in os.listdir(directory_path)) and os.path.isdir(directory_path+'/'+'excited_structure')

    return found_ground_structure_folder, found_excited_structure_folder

# -----------------------------------------------------------------

def is_finished_reading(energy, maximum_force_converged, rms_force_converged, maximum_distance_converged, rms_distance_converged, opt_job):
    if opt_job:
        return all((entry is not None) for entry in [energy, maximum_force_converged, rms_force_converged, maximum_distance_converged, rms_distance_converged])
    else:
        return (energy is not None)

# -----------------------------------------------------------------

def format_worksheet(workbook):
    '''
    This method is designed to 

    Parameters
    ----------
    workbook : xlsxwriter.Workbook
        This is the workbook you want to add formats too. 

    Returns
    -------
    merge_format_mol : xlsxwriter.format
        This the format for the molecule name title
    merge_format_dimer : xlsxwriter.format
        This the format for the dimer name title
    merge_format_barrier : xlsxwriter.format
        This the format for the barrier between molecule data
    table_format : xlsxwriter.format
        This the format for the row and column names in a table of data
    number_format : xlsxwriter.format
        This is the format for the numbers in each cell.
    '''
    
    merge_format_mol = workbook.add_format({
        'bold':     True,
        #'border':   6,
        'align':    'center',
        'valign':   'vcenter',
        'fg_color': '#D7E4BC',
    })

    merge_format_dimer = workbook.add_format({
        'bold':     True,
        #'border':   6,
        'align':    'center',
        'valign':   'vcenter',
        'fg_color': '#D7E4BC',
    })

    merge_format_barrier = workbook.add_format({
        'bold':     True,
        #'border':   6,
        'align':    'center',
        'valign':   'vcenter',
        'fg_color': '#000000',
    })

    table_format = workbook.add_format({
        'bold':     True,
        'align':    'left',
        'valign':   'vcenter',
    })

    number_format = workbook.add_format({
        'num_format': '0'
    })

    return merge_format_mol, merge_format_dimer, merge_format_barrier, table_format, number_format

# -----------------------------------------------------------------

