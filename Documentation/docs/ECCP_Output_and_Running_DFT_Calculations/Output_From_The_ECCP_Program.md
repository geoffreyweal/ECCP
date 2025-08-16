# Outputs from the Electronic Crystal Calculation Prep program

The Electronic Crystal Calculation Prep program will create three to five folders. These are: 

General files are placed in:

* ``ECCP_Data``: This contains all the information about molecules and dimers in the crystal, as well as the xyz files of molecules and dimers. 

Atomic Transition Charge (ATC) files are placed in:

* ``All_ATC_Calc_Jobs``: This folder contains the ATC Gaussian/ORCA files of all the molecules in the crystal. These files are created in case you need to look through them, but you don't necessarily need to run these in Gaussian/ORCA. 
* ``Unique_ATC_Calc_Jobs``: This folder contains the ATC Gaussian/ORCA files of just the structurally unique molecules in the crystal. These are the files you want to run in Gaussian/ORCA. 

Reorganisation energy (RE) files are placed in:

* ``All_RE_Calc_Jobs``: This folder contains the RE Gaussian/ORCA files of all the molecules in the crystal. These files are created in case you need to look through them, but you don't necessarily need to run these in Gaussian/ORCA. 
* ``Unique_RE_Calc_Jobs``: This folder contains the RE Gaussian/ORCA files of just the conformationally unique molecules in the crystal. These are the files you want to run in Gaussian/ORCA. 

Franck-Condon and Huang-Rhys factors (FC/HR) files are placed in:

* ``All_FC_Calc_Jobs``: This folder contains the FC/HR Gaussian/ORCA files of all the molecules in the crystal. These files are created in case you need to look through them, but you don't necessarily need to run these in Gaussian/ORCA. 
* ``Unique_FC_Calc_Jobs``: This folder contains the FC/HR Gaussian/ORCA files of just the conformationally unique molecules in the crystal. These are the files you want to run in Gaussian/ORCA. 

Electronic Energy Transfer (EET) files are placed in:

* ``All_EET_Calc_Jobs``: This folder contains the EET Gaussian/ORCA files of all the dimers in the crystal. These files are created in case you need to look through them, but you don't necessarily need to run these in Gaussian/ORCA. 
* ``Unique_EET_Calc_Jobs``: This folder contains the EET Gaussian/ORCA files of just the structurally unique dimers in the crystal. These are the files you want to run in Gaussian/ORCA. 

Eigendata files that is used to obtain Intermolecular Charge Transfer (ICT) energies are placed in:

* ``All_Eigendata_Calc_Jobs``: This folder contains the Gaussian/ORCA files of all the dimers in the crystal for obtaining molecular orbtial (MO) energies, coefficients, and overlap matrices for dimers and their associated monomers. These files are created in case you need to look through them, but you don't necessarily need to run these in Gaussian/ORCA. 
* ``Unique_Eigendata_Calc_Jobs``: This folder contains the Gaussian/ORCA files of just the structurally unique dimers in the crystal for obtaining molecular orbtial (MO) energies, coefficients, and overlap matrices for dimers and their associated monomers. These are the files you want to run in Gaussian/ORCA. 


## Files in the ``ECCP_Information`` folder

The folder called ``ECCP_Information`` contains all the information about molecules and dimers in the crystal, as well as the xyz files of molecules and dimers. The files included in this folder include: 

* ``All_Dimer_Information.txt``: This file contains information about all the dimers that were found in the crystal, including the molecules in the dimer and how they are positioned relative to one another in the unit cell. 
* ``Unique_Dimer_Information.txt``: This file includes which equivalent dimers are the same as which unique dimers. This may not contain all the unique dimers in it, as it only assigns equivalent dimers to unique dimers. 
* ``crystal.xyz``: The original crystal in xyz format. Name given will be the same as the crystal file name given to Electronic_Crystal_Calculation_Prep.
* ``human_friendly_crystal_small.xyz``: The crystal in a version that makes it easier to view and understand the crystal packing in the crystal. This is a smaller version of this view. 
* ``human_friendly_crystal_large.xyz``: The crystal in a version that makes it easier to view and understand the crystal packing in the crystal. This is a larger version of this view. 
* ``All_Molecules``: This is a folder containing all the molecules found in the crystal. 
* ``Unique_Molecules``: This is a folder containing all the unique molecules found in the crystal. 
* ``All_Dimers``: This is a folder containing all the dimers found in the crystal. 
* ``Unique_Dimers``: This is a folder containing all the unique dimers found in the crystal. 


## Files in the ATC, EET, and Eigendata folders

These folders contain two files in them. These are:

* An input file  that contain all the information about the molecule or the dimer to run. This is either:

	* Gaussian input file (``.gjf``), or a
	* ORCA input file (``.inp``)

* ``submit.sl``: This is the submission file requires to submit the gaussian job to slurm. 


## Files in the RE folders

These folders contain a few files. 

In the ``ground_structure`` folder will be:

* ``GS_GS_PM6.gjf``: The Gaussian file that contain information about how to optimise your molecule to the ground state structure. This initial calculation will be perform using the PM6 force field. This calculation is performed before the resulting optimised structure is further optimised with your intended functional. Optimising the structure using the PM6 force field initially saves computational time. 
* ``GS_GS_submit.sl``: File to submit the ``GS_GS.gjf`` Gaussian file to slurm.
* ``GS_GS_freq_submit.sl``: File to submit the ``GS_GS_freq.gjf`` Gaussian file to slurm.
* ``GS_ES_submit.sl``: File to submit the ``GS_ES.gjf`` Gaussian file to slurm.

During the course of calculation, the following files will be created:

* ``GS_GS.gjf``: The Gaussian file that contain information about how to optimise your molecule to the ground state structure. This calculation will be performed using the functional you intended to use. 
* ``GS_GS_freq.gjf``: This Gaussian file contains information for performing a frequency calculation upon the ground state. 
* ``GS_ES.gjf``: This Gaussian file contains information for performing a single point calculation of the excited state of your molecule with the ground state structure. 

In the ``excited_structure`` folder will be:

* ``ES_ES_PM6.gjf``: The Gaussian file that contain information about how to optimise your molecule to the excited state structure. This initial calculation will be perform using the PM6 force field. This calculation is performed before the resulting optimised structure is further optimised with your intended functional. Optimising the structure using the PM6 force field initially saves computational time. 
* ``ES_ES_submit.sl``: File to submit the ``ES_ES.gjf`` Gaussian file to slurm.
* ``ES_ES_freq_submit.sl``: File to submit the ``ES_ES_freq.gjf`` Gaussian file to slurm.
* ``ES_GS_submit.sl``: File to submit the ``ES_GS.gjf`` Gaussian file to slurm.

During the course of calculation, another file will be created:

* ``ES_ES.gjf``: The Gaussian file that contain information about how to optimise your molecule to the excited state structure. This calculation will be performed using the functional you intended to use. 
* ``ES_ES_freq.gjf``: This Gaussian file contains information for performing a frequency calculation on the excited state structure. 
* ``ES_GS.gjf``: This Gaussian file contains information for performing a single point calculation of the ground state of your molecule with the excited state structure. 

## Files in the FC folders

These folders contain two files in them. These are:

* A Gaussian input file (``FC.gjf``) that contain all the information and commands needed to obtain Franck-Condon and Huang-Rhys factors. 
* ``FC_submit.sl``: This is the submission file requires to submit the gaussian job to slurm. 

NOTE: During this job, the ``freq_gaussian.chk`` checkpoint files obtained during the reorganisation energies will be moved from the ``ground_structure`` and ``excited_structure`` folders, into the franck-condon folder and renamed as ``GS_gaussian.chk`` and ``ES_gaussian.chk``, respectively. This is performed by the ``copy_checkpoint_files_for_franck_condon_calculation.py`` program that you will see is performed during the ``FC_submit.sl`` job. 










