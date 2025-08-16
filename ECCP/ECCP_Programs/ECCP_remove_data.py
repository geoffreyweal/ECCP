'''
Geoffrey Weal, ECCP_tidy_data.py, 6/6/24

This program is designed to remove a data set if you want to remove it from your database. 
'''
import os, shutil
from tqdm import tqdm

class CLICommand:
    """This program is designed to remove a data set if you want to remove it from your database. 
    """

    @staticmethod
    def add_arguments(parser):
        parser.add_argument('--path_to_ECCP_Data_folder', nargs='*', help='This is the path to the ECCP_data folder', default=['ECCP_Data'])
        parser.add_argument('remove_dataset_names', nargs='*', help='These are the names of the dataset you want to remove.')

    @staticmethod
    def run(args):

        # First, obtain the path to the ECCP Data database folder. 
        path_to_ECCP_Data_folder = args.path_to_ECCP_Data_folder
        if len(path_to_ECCP_Data_folder) > 1:
            to_string  = 'Error: You can only enter one path for "path_to_ECCP_Data_folder"\n'
            to_string += f'path_to_ECCP_Data_folder = {path_to_ECCP_Data_folder}\n'
            to_string += 'Check this and rerun this program.'
            raise Exception(to_string)
        path_to_ECCP_Data_folder = path_to_ECCP_Data_folder[0]

        # Second, obtain the list of dataaset names from args
        remove_dataset_names = args.remove_dataset_names

        # Third, run the algorithm.
        Run_method(path_to_ECCP_Data_folder, remove_dataset_names)

def Run_method(path_to_ECCP_Data_folder, remove_dataset_names):
    """
    This program is designed to remove a data set if you want to remove it from your database. 

    Parameters
    ----------
    path_to_ECCP_Data_folder : str.
        This is the path to the ECCP_Data folder.
    remove_dataset_names : list of str.
        These are the names of the crystals you want to remove from the ECCP_Data folder.
    """

    print(f'Removing the following crystals from {path_to_ECCP_Data_folder}: {remove_dataset_names}')

    # First, for each calculation type and ECCP_Information in path_to_ECCP_Data_folder
    pbar1 = tqdm(sorted(os.listdir(path_to_ECCP_Data_folder)), leave=False)
    for foldername in pbar1:

        # Second, write the path to the ECCP folder to check through 
        path_to_ECCP_Data_type = path_to_ECCP_Data_folder+'/'+foldername

        # Third, report when you are up to.
        pbar1.description('Looking through: '+str(path_to_ECCP_Data_type))

        # Fourth, make sure that foldername is a folder.
        if not os.path.isdir(path_to_ECCP_Data_folder+'/'+foldername):
            continue

        # Fifth, Look through the folders in path_to_ECCP_Data_type.
        pbar2 = tqdm(sorted(os.listdir(path_to_ECCP_Data_type)), desc='Reading through: '+str(path_to_ECCP_Data_type), leave=False)
        for crystal_foldername in pbar2:

            # Sixth, this is the path that we want to potentially remove.
            path_to_crystal_data = path_to_ECCP_Data_type+'/'+crystal_foldername

            # Sixth, check that path_to_crystal_data is a path.
            if not os.path.isdir(path_to_crystal_data):
                continue

            # Seventh, Check that crystal_foldername is one of the crystals we want to remove from the ECCP Data database.
            if crystal_foldername in remove_dataset_names:

                # Seventh, if this is the case, remove this crystal data from your ECCP Data database.
                print(f'Removing {crystal_foldername} from {path_to_ECCP_Data_type}')
                shutil.rmtree(path_to_crystal_data)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

