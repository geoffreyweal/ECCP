'''
Geoffrey Weal, ECCP_processing_Eigendata_data.py, 10/6/22

This program is designed to process the Eigen data from Gaussian output files

'''
import os, time
from datetime import datetime, timedelta
from SUMELF import remove_folder, make_folder

from ECCP.ECCP_Programs.processing_Eigendata_methods.process_Eigendata_to_disk import process_Eigendata_to_disk

# ---------------------------------------------------------------------

class CLICommand:
    """Will process the eigendata and extract from your output.log files and place them in individual txt files.
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
    This method is the main method for running this program. 
    """

    # First, get the general variables for processing data.
    overall_path = os.getcwd()
    log_filename = 'output.log'

    # Second, start the process and start timing
    print('------------------------------------------------')
    dt_string = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    print('Starting the process_ICT program at ('+str(dt_string)+')')
    start_time = time.time()
    print('------------------------------------------------')

    # Third, obtain the eigendata from the Gaussian output.log files. 
    print('Extracting Eigen-data from output.log files to external files')
    print('Note: This program may take some time if you have recorded eigendata, such as orbital overlap matrices, as these matrices can be very large depending on the number of atoms in your dimer.')
    eigendata = {}
    issues = []
    for root, dirs, files in os.walk(overall_path):
        dirs.sort()

        # Fourth, check to see if you have these files in your folder. 
        # If you do, then you are in the right place for obtaining eigendata on your dimer and molecules.
        if not (('Dimer' in dirs) and ('Monomer_1' in dirs) and ('Monomer_2' in dirs)):
            continue

        # Fifth, indicate that a Gaussian job has been found
        print('------------------------------------------------')
        print(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+' - Found a Gaussian job: '+str(root))
        dirs[:]  = []
        files[:] = []

        # Sixth, process the eigendata.
        issue = process_Eigendata_to_disk(root, log_filename, start_time)
        if issue is not None:
            issues.append(issue)
            print(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+' - Their was an issue with this job.')
            print('------------------------------------------------')
            continue

    # Seventh, write any issues to the terminal.
    print('------------------------------------------------')
    if len(issues) > 0:
        print('The following Gaussian jobs could not be processed because they have not finished running or did not complete successfully.')
        for issue in issues:
            print(issue)
        print('------------------------------------------------')
    print(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+' - This ICT calculations program has finished successfully!')
    print('Total running time (HH:MM:SS): '+str(timedelta(seconds=time.time() - start_time)))
    print('------------------------------------------------')

# ----------------------------------------------------------------------------------------------------------------------------------



