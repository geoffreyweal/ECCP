'''
Geoffrey Weal, processing_OPV_Dimer_data.py, 7/3/22

This program is designed to process the data from the Gaussian TD DFT calculations and put them into formats that the user can use for further work. 

'''
import os, time
import numpy as np
from datetime import datetime, timedelta

from SUMELF import remove_folder, make_folder
from ECCP.ECCP_Programs.processing_EET_methods.get_EET_data import get_EET_data
from ECCP.ECCP_Programs.processing_EET_methods.write_data_to_excel import write_data_to_excel

# ---------------------------------------------------------------------

class CLICommand:
    """Will process EET Data into text files and an excel file.
    """

    @staticmethod
    def add_arguments(parser):
        pass

    @staticmethod
    def run(args):
        Run_method()

# ---------------------------------------------------------------------

def Run_method():
    """
    This method is the main method for running this program
    """
    # General variables for processing data.
    overall_path = os.getcwd()
    log_filename = 'output.log'
    eet_data_foldername = 'EET_Data'
    individual_eet_data_foldername = 'Individual_EET_Data'
    for foldername in [eet_data_foldername, individual_eet_data_foldername]:
        remove_folder(foldername)
        make_folder(foldername)

    print('------------------------------------------------')
    dt_string = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    print('Starting the process_EET program at ('+str(dt_string)+')')
    start_time = time.time()
    print('------------------------------------------------')

    # First, obtain the electronic coupling data from the Gaussian output.log files. 
    electronic_coupling_data, issues = get_EET_data(overall_path, log_filename, start_time)

    # Thirteenth, write the EET to an excel file.
    write_data_to_excel(electronic_coupling_data, eet_data_foldername, individual_eet_data_foldername, start_time)

    # Fourteenth, write any issues to the terminal.
    print('------------------------------------------------')
    if len(issues) > 0:
        print('The following Gaussian jobs could not be processed because they have not finished running or did not complete successfully.')
        for issue in issues:
            print(issue)
        print('------------------------------------------------')
    print(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+' - This EET calculations program has finished successfully!')
    print('Total running time (HH:MM:SS): '+str(timedelta(seconds=time.time() - start_time)))
    print('------------------------------------------------')

# ------------------------------------------------------------------------------------------------------------------------------------------------------------



