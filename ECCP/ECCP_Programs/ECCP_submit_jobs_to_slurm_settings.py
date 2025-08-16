'''
Geoffrey Weal, Run_Adsorber_submitSL_slurm.py, 16/06/2021

This program is designed to submit all sl files called submit.sl to slurm.
'''
import os
from subprocess import Popen, PIPE

from ECCP.ECCP_Programs.ECCP_submit_jobs_to_slurm_settings_methods.settings_methods import check_submit_settingsTXT, change_settings
from ECCP.ECCP_Programs.ECCP_submit_jobs_to_slurm_settings_methods.print_settings   import print_settings

# Get the path to the settings script.
this_scripts_path = os.path.dirname(os.path.abspath(__file__))
submit_settings_name = 'ECCP_submit_gaussian_jobs_to_slurm_settings_methods/submit_settings.txt'
path_to_settings_txt_file = this_scripts_path+'/'+submit_settings_name

class CLICommand:
    """Change the settings for the ECCP program that submits ECCP jobs to slurm.
    """

    @staticmethod
    def add_arguments(parser):
        parser.add_argument('settings', nargs='*', help='Enter in here a setting for submit that you would like to change.')

    @staticmethod
    def run(args_submit):
        args_submit_settings = args_submit.settings
        check_submit_settingsTXT(path_to_settings_txt_file)
        if len(args_submit_settings) > 0:
            change_settings(args_submit_settings)
        else:
            print_settings()

# =========================================================================================================================================

