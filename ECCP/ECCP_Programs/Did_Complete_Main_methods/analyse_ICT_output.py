'''
analyse_ICT_output.py, Geoffrey Weal, 29/12/22

This method is designed to check if an ICT Gaussian job has completed or not.
'''
import os
from ECCP.ECCP_Programs.shared_general_methods.shared_general_methods import reverse_readline

def analyse_ICT_output(path_to_outputLOG):
    """
    This method is designed to check if an ICT Gaussian job has completed or not.

    Parameters
    ----------
    path_to_outputLOG : str.
        This is the path to the ICT log file. 

    Returns
    -------
    The results of the job: str.
        * 'NBY': Not begun yet.
        * 'NC' : Not complete.
        * 'C'  : Complete.
    """
    pass