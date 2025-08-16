'''
Geoffrey Weal, ECCP_collect_ATC_files.py, 18/2/24

This program is designed to collect the ATC charge files. 

'''
import os, time
from tqdm import tqdm
from shutil import copyfile
from datetime import datetime, timedelta

from SUMELF import remove_folder, make_folder

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

    # First, check that the folders that contains the ATC chg files exist. 
    if (not 'Unique_ATC_Calc_Jobs' in os.listdir('.')) or (not 'All_ATC_Calc_Jobs' in os.listdir('.')):
        toString = 'Error: You need to run this program inside a "ECCP_Data" file where either or both "All_ATC_Calc_Jobs" and "Unique_ATC_Calc_Jobs" exist.'
        exit(toString)

    # Second, general variables for processing data.
    overall_path = os.getcwd()
    chg_filename = 'output.chg'
    atc_data_foldername = 'ATC_charge_files'
    for foldername in [atc_data_foldername]:
        remove_folder(foldername)
        make_folder(foldername)

    # Third, write a message to the user. 
    print('------------------------------------------------')
    dt_string = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    print('Starting the collect_ATC_files program at ('+str(dt_string)+')')
    start_time = time.time()
    print('------------------------------------------------')

    # Fourth, initalise the progress bar so the user can see what is being copied
    pbar = tqdm(desc='chg files')

    # Fifth, copy all chg files from "All_ATC_Calc_Jobs" and "Unique_ATC_Calc_Jobs" into "ATC_charge_files"
    for overall_general_folder, molecule_folder_type in [('Unique_ATC_Calc_Jobs', 'Unique_Molecules'), ('All_ATC_Calc_Jobs', 'All_Molecules')]: 

        # 5.1: Look through each crystal folder in overall_general_folder
        for crystal_folder in sorted(os.listdir(overall_general_folder)):
            if not os.path.isdir(overall_general_folder+'/'+crystal_folder):
                continue

            # 5.2: Look through each molecule folder in overall_general_folder+'/'+crystal_folder
            for molecule_folder in sorted(os.listdir(overall_general_folder+'/'+crystal_folder)):
                if not os.path.isdir(overall_general_folder+'/'+crystal_folder+'/'+molecule_folder):
                    continue

                # 5.3: Record the path to the molecule folder for the current crystal.
                path_to_molcule_folder = overall_general_folder+'/'+crystal_folder+'/'+molecule_folder

                # 5.4: For every folder in path_to_molcule_folder
                for root, dirs, files in os.walk(path_to_molcule_folder):

                    # 5.5: Sort dirs and files
                    dirs.sort(); files.sort()

                    # 5.6: If root contains a .chg file(s), we want to do stuff in this root.
                    for file in files:
                        if os.path.isfile(root+'/'+file) and file.endswith('.chg'):
                            break
                    else:
                        continue

                    # 5.7: Copy the .chg files to 'ATC_charge_files'
                    for file in files:
                        if os.path.isfile(root+'/'+file) and file.endswith('.chg'):

                            # 5.8: Note to the user what chg will be copied
                            pbar.set_description(root+'/'+file); pbar.update(1)

                            # 5.9: Make the folder to place chg files into
                            make_folder(atc_data_foldername+'/'+root)

                            # 5.10: Copy the file into the "ATC_charge_files" folder.  
                            copyfile(root+'/'+file, atc_data_foldername+'/'+root+'/'+file)

                            # 5.11: If the "ECCP_Information" folder exists, try to copy the molecule xyz file into the folder as well.
                            path_to_xyz_file = 'ECCP_Information/'+str(crystal_folder)+'/'+str(molecule_folder_type)+'/'+str(molecule_folder)+'.xyz'
                            if os.path.exists(path_to_xyz_file) and os.path.isfile(path_to_xyz_file):
                                copyfile(path_to_xyz_file, atc_data_foldername+'/'+root+'/'+str(molecule_folder)+'.xyz')

    # Sixth, close the progress bar.
    pbar.close()

    # Seventh, print the running time and final message for the user. 
    print(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+' - This EET calculations program has finished successfully!')
    print('Total running time (HH:MM:SS): '+str(timedelta(seconds=time.time() - start_time)))
    print('------------------------------------------------')

# ------------------------------------------------------------------------------------------------------------------------------------------------------------











