'''
Geoffrey Weal, processing_OPV_Dimer_data.py, 9/3/22

This script contains methods for processing_OPV_Dimer_data.py

'''
import os

from ECCP.ECCP_Programs.shared_general_methods.shared_general_methods import reverse_readline

# Constants and conversions that are useful for processing data.
planks_constant = 4.135667696 * (10.0 ** -15.0) #eVs
speed_of_light = 2.99792458 * (10.0 ** 8.0) #ms-1
metre_to_centimetre = 100.0 #cm m-1
eV_to_meV = 1000.0 #meV eV-1
eV_to_micro_eV = 1000000.0 #micro eV eV-1
eV_to_inverse_cm = 1/(metre_to_centimetre * speed_of_light * planks_constant) #cm eV-1

# -----------------------------------------------------------------

def is_this_calc_an_eet_calc(path_to_outputLOG):
    """
    This method will look through the output.log file to see if this job was an eet job

    Parameters
    ----------
    root : str.
        This is the path to where the calculation files are located

    Returns
    -------
    True if this calc is an eet calc, False if not.
    """

    looking_at_EET_output = False
    with open(path_to_outputLOG,'r') as outputLOG:
        for line in outputLOG:
            if 'eet' in line:
                looking_at_EET_output = True
                break
            if 'Input orientation:' in line:
                break
    return looking_at_EET_output

# -----------------------------------------------------------------

def get_electronic_coupling_of_lowest_TD_state(log_filepath):
    '''
    This method is designed to obtain the electronic coupling for the lowest TD state as obtained using EET. This is between state 1 in fragment 1 and state 1 in fragement 2. 

    Parameters
    ----------
    log_filepath : str
        This is the path to the log file.

    Returns
    -------
    True if both the input .gjf file and the output .log files are found. 
    '''

    # all values to obtain
    delta_w_value = None
    coulomb_value = None
    exact_exchange_value = None
    exchange_correlation_value = None
    overlap_value = None
    w_avg_value = None
    w_avg_times_Overlap_value = None
    total_coupling_value = None

    for line in reverse_readline(log_filepath):
        if 'Electronic Coupling for Excitation Energy' in line:
            break # if this line is found, all electronic coupling calcs have been found, so dont need to read the file anymore.
        elif line.startswith('   delta-w'):
            line = line.rstrip().replace('=','').split()
            delta_w_value = float(line[1]) # eV
        elif line.startswith('   Coulomb'):
            line = line.rstrip().replace('=','').split()
            coulomb_value = float(line[1]) # eV
        elif line.startswith('   Exact-exchange'):
            line = line.rstrip().replace('=','').split()
            exact_exchange_value = float(line[1]) # eV
        elif line.startswith('   Exchange-correlation'):
            line = line.rstrip().replace('=','').split()
            exchange_correlation_value = float(line[1]) # eV
        elif line.startswith('   w-avg*Overlap'):
            line1 = line.rstrip().replace('   w-avg*Overlap             =','').lstrip().split()
            w_avg_times_Overlap_value = float(line1[0]) # eV
            line2 = line.rstrip().replace('   w-avg*Overlap             =','').lstrip().split('(')[1].replace(')','')
            w_avg_value, overlap_value = line2.split(',')
            w_avg_value   = float(w_avg_value.replace('w-avg=','').replace('eV','')) # eV
            overlap_value = float(overlap_value.replace('Ovlp=','').replace('D','E').replace('eV',''))
        elif line.startswith('   Total coupling'):
            line = line.rstrip().replace('=','').split()
            total_coupling_value = float(line[2]) # eV

    # If this was a Hartree-Fock calculation, there will be no Excharge Correlation component to this analysis
    # In this case, set this to 0.0 eV
    if (exchange_correlation_value is None) and is_a_HF_calculation(log_filepath):
        exchange_correlation_value = 0.0

    # Return all the values obtained for coupling constants. 
    electronic_coupling_datum = [delta_w_value,coulomb_value,exact_exchange_value,exchange_correlation_value,w_avg_times_Overlap_value,w_avg_value,overlap_value,total_coupling_value]

    # Default turn Nones into zero for now
    for index in range(len(electronic_coupling_datum)):
        if electronic_coupling_datum[index] is None:
            electronic_coupling_datum[index] = 0.0

    if any([(value is None) for value in electronic_coupling_datum]):
        return 'error' # something weird has happened, return an error message. 
    else:
        return tuple(electronic_coupling_datum)

def is_a_HF_calculation(log_filepath):
    """
    Determine if this calculation is a Hartree-Fock calculation

    Parameters
    ----------
    log_filepath : str
        This is the path to the log file.

    Returns
    -------
    Returns True if this is a Hartree-Fock calculation, False if not. 
    """

    LOGfile = open(log_filepath,'r')
    for line in LOGfile:
        if line.startswith(' #') and ('HF' in line):
            LOGfile.close()
            return True
        if 'Population analysis using the SCF Density.' in line:
            break
    LOGfile.close()
    return False

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


