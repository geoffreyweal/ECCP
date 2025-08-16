# Post-processing programs that are available in the ECCP program

There are a few programs that you can run after you have submitted atomic transition charge (ATC), reorganisation energy (RE), electronic energy transfer (EET), eigendata, and intermolecular charge transfer (ICT) calculations to slurm. These are:

* ``ECCP submit``: For submitting Gaussian jobs to your slurm computer cluster.
* ``ECCP submit_multiwfn``: For submitting Multiwfn jobs once your ATC gaussian jobs have finished. 
* ``ECCP did_complete``: Check if jobs have completed successfully or not.
* ``ECCP reset``: For resuming reorganisatino calculations if they did not completed in the allotted time.
* ``ECCP process_RE``: For obtaining results from the reorganisation energy calculations.
* ``ECCP process_EET``: For obtaining result from the electronic energy transfer calculations.
* ``ECCP process_ICT``: For obtaining result from the intermolecular charge transfer calculations.
* ``ECCP process_Eigendata``: For solely obtaining matrix files from Eigendata jobs. 
* ``ECCP tidy``: For removing unnecessary files that are not needed anymore and take up a lot of hard disk space. 

These programs are described in more detail below.


## ``ECCP submit``: Submitting Gaussian jobs to slurm

It is possible to get the ECCP to submit all the Gaussian jobs to slurm automatically. 

To use ``ECCP submit``, you want to move into the parent folder that contains all the subdirectories of ECCP Gaussian jobs you want to submit to slurm. For example, if you have a bunch of EET calculations you would like to submit, you would want to move into the ``Unique_EET_Gaussian_Jobs`` before running ``ECCP submit``:

```bash
cd path/to/parent/folder
```

Once you are in the parent folder you want to submit ECCP Gaussian jobs for, type ``ECCP submit`` into the terminal. 

```bash
ECCP submit
```

When you run ``ECCP submit`` on this folder in the terminal, ECCP will look through all the subdirectories in search of slurm submit files, specifically those called ``submit.sl`` or one of the reorganisation energy submit files (``eGS_gGS_main_opt_submit.sl`` , ``eGS_gGS_refine_opt_submit.sl``, ``eGS_gGS_freq_submit.sl``, ``eES_gGS_submit.sl``, ``eES_gES_main_opt_submit.sl``, ``eES_gES_refine_opt_submit.sl``, ``eES_gES_freq_submit.sl``, ``eGS_gES_submit.sl``). For reorganisation energy calculations, only the appropriate submit file will be submitted for the ground and excited state calculations. 


## ``ECCP mass_submit_multiwfn``: Submitting MultiWFN ATC jobs manually

When running ATC jobs, the ``ECCP`` program will first run Gaussian on the ``input.gjf`` file (containing the molecules you want to obtain ATC results for). Once this calculation has finished, it will then run the MultiWFN program on the ``.wfn`` file that is created by Gaussian automatically. However, if an issue arises that stops MultiWFN from running, it is still possible to run MultiWFN on slurm by submitting the associated ``submit_multiwfn.sl`` file to slurm: 

```bash
cd path/to/atc/submit_multiwfn/file
sbatch submit_multiwfn.sl
```

This will run the MultiWFN program as long as the ``.wfn`` file exists for that job. 

You can also mass submit all these file by typing ``ECCP mass_submit_multiwfn`` in to the terminal. 

```bash
cd path/to/ATC/jobs
ECCP mass_submit_multiwfn
```

This will search within the subdirectories for ``submit_multiwfn.sl`` files, which it will then submit to slurm. 


## ``ECCP did_complete``: Checking which jobs have completed successfully

You can check which gaussian jobs have finished successfully and which jobs have not finished successfully by typing ``ECCP did_complete`` into the terminal. ``ECCP`` will go through all the folders and will look for ``output.log`` files. It will then look through these files for key words that indicate the job has completed successfully, whether that be ATC, RE, EET, or Eigendata (ICT) calculations. To do this, move into the folder that contains your jobs and  run ``ECCP did_complete`` in the terminal: 

```bash
cd path/to/gaussian/jobs
ECCP did_complete
```


## ``ECCP reset``: Resetting unfinished reorganisation energy calculation for resubmission to slurm

Some of your Gaussian calculations may not complete or, in the case of reorganisation energy calculations, you may need to reset your geometry optimisations so you can restart a Gaussian jobs that may have timed our. It is possible to reset all the files in your ECCP folders (or setup your reorganisation energy calculations to resume geometry optimisation calculations) so that you can submit them all to slurm without manually doing anything using ``ECCP submit``. 

To run this program, move into the panret folder that contains your jobs and run ``ECCP reset`` in the terminal: 

```bash
cd path/to/gaussian/jobs
ECCP reset
```

``ECCP reset`` will move through all the subdirectories for uncompleted jobs to reset. 

!!! tip "IMPORTANT"

	Make sure that there are no jobs currently running in slurm before you use this program. Cancel any running slurm jobs before running ``ECCP reset`` in the parent folder of these currently running jobs. 

!!! note

	This program **WILL NOT** modify any files or folders that have successfully finished geometry optimisations. You can run it from any folder, and it will only modify subdirectories that have not completed. 

For reorganisation energy geometry optimisation jobs, ``ECCP reset`` will look to see if the checkpoint file for the job exists. If it does, it will use this to resume the geometry optimisation, as this contains the initial SCF data which can save time initially for optimisations. If the checkpoint file does not exist or it can not be used for some reason, the checkpoint file will not be used and instead the last geometry the calculation optimisation to will be taken from the Gaussian output ``log`` file. 

Once you have run this program, you can resubmit jobs for resuming reorganisation energy calculations using the ``ECCP submit`` command.


## Processing results of atomic transition charge Gaussian jobs using ``ECCP process_ATC``

To do




## Processing results of reorganisation energy Gaussian jobs using ``ECCP process_RE``

To do




## Processing results of Electronic Energy Transfer jobs using ``ECCP process_EET``

It is possible for this program to process the electronic coupling results from Gaussian jobs and present the results in text and in an excel spreadsheet. If the Gaussian job runs successfully, it will give a ``output.log`` file that contains values for the electronic coupling between the molecules in the dimer using the EET method in Gaussian. To do this, move into your ``Unique_EET_Gaussian_Jobs`` folder and run ``ECCP process_EET`` by typing into the terminal: 

```bash
cd Unique_EET_Gaussian_Jobs
ECCP process_EET
```

This method will go through all the folders in ``Unique_EET_Gaussian_Jobs`` and extract the data from every successfully run ``output.log`` file and 

1. Place this data into an excel file called ``Unique_EET_Gaussian_Jobs.xlsx``, and two text files called ``Unique_EET_Gaussian_Jobs.txt`` and ``Unique_EET_Gaussian_Jobs_wavenumber.txt``. These are located in the ``EET_Data`` folder that has just been created. 
2. Also in the newly created ``EET_Data`` folder are text files of the electronic coupling values for particular functionals and basis sets are also be given in the  ``TXT_of_Func_and_basis_sets_Energy`` and ``TXT_of_Func_and_basis_sets_Wavenumber`` folders. These files contain the electron coupling information in meV and in cm<sup>-1</sup>. 
3. The ``Individual_EET_Data`` contains electronic coupling values for each dimer in individual text files. 


## Processing results of Intermolecular Charge Transfer (ICT) jobs using ``ECCP process_ICT``

**To Do**

It is possible for this program to process the intermolecular charge transfer results from Gaussian jobs and present the results in text and in an excel spreadsheet. If the Gaussian job runs successfully, it will give a ``output.log`` file that matrix data for obtaining intermolecular charge transfer energies, such as hole and electron transfer energies. To do this, move into your ``Unique_Eigendata_Gaussian_Jobs`` folder and run ``ECCP process_ICT`` by typing into the terminal: 

```bash
cd Unique_Eigendata_Gaussian_Jobs
ECCP process_ICT
```

This method will go through all the folders in ``Unique_Eigendata_Gaussian_Jobs`` and extract the data from every successfully run ``output.log`` file and: 

1. Locate and extract the overlap matrix and molecular orbital (MO) energies and coefficients from the ``output.log`` and save these matrices as txt files called:

	* ``orbital_overlap_matrix.txt``: This is the overlap matrix. 
	* ``MO_energies.txt``: These are the energies of the MOs in your monomer/dimer.
	* ``MO_coefficients.txt``: These are the coefficients of the MOs in your monomer/dimer.
	* ``MO_orbital_names.txt``: These are the names and the indices of the MO coefficients that are involved with each atom in your monomer/dimer.
	* ``MO_occupancies.txt``: This file indicates which orbtials are occupied and which orbtials are vacant. 

2. Remove any matrices from the ``output.log``  files. This is necessary as these matrices can make a ``output.log`` incredibly large (GBs is size). This will reduce the size of the ``output.log`` to a few MBs or less. All the necessary matrix data will be located in txt file as mentioned in (1.). 
3. Calculate the hole and electron transfer energies and place these values into an excel file called ``Unique_ICT_Gaussian_Jobs.xlsx``, and two text files called ``Unique_ICT_Gaussian_Jobs.txt`` and ``Unique_ICT_Gaussian_Jobs.txt``. These are located in the ``ICT_Data`` folder that has just been created. 
4. Also in the newly created ``ICT_Data`` folder are text files of the intermolecular charge transfer values for particular functionals and basis sets are also be given in the  ``TXT_of_Func_and_basis_sets_Energy`` and ``TXT_of_Func_and_basis_sets_Wavenumber`` folders. These files contain the intermolecular charge transfer energy values in meV and in cm :sup:`-1`. 
5. The ``Individual_ICT_Data`` contains intermolecular charge transfer energy values for each dimer in individual text files. 


## Processing matrix values from Eigendata jobs using ``ECCP process_Eigendata``

**To Do**

Instead of processing the ICT jobs completely, you may just want to obtain the matrix text files from the Eigendata folder. If you only want to process the matrix data from your successfully run ``output.log`` files, move into your ``Unique_Eigendata_Gaussian_Jobs`` folder and run ``ECCP process_Eigendata`` by typing into the terminal: 

```bash
cd Unique_Eigendata_Gaussian_Jobs
ECCP process_Eigendata
```

This method will go through all the folders in ``Unique_Eigendata_Gaussian_Jobs`` and extract the data from every successfully run ``output.log`` file and: 

1. Locate and extract the overlap matrix and molecular orbital (MO) energies and coefficients from the ``output.log`` and save these matrices as txt files called:
	
	* ``orbital_overlap_matrix.txt``: This is the overlap matrix. 
	* ``MO_energies.txt``: These are the energies of the MOs in your monomer/dimer.
	* ``MO_coefficients.txt``: These are the coefficients of the MOs in your monomer/dimer.
	* ``MO_orbital_names.txt``: These are the names and the indices of the MO coefficients that are involved with each atom in your monomer/dimer.
	* ``MO_occupancies.txt``: This file indicates which orbtials are occupied and which orbtials are vacant. 

	2. Remove any matrices from the ``output.log``  files. This is necessary as these matrices can make a ``output.log`` incredibly large (GBs is size). This will reduce the size of the ``output.log`` to a few MBs or less. All the necessary matrix data will be located in txt file as mentioned in (1.). 


## Remove large and unnecessary files using ``ECCP tidy``

Gaussian produces a number of files that are not needed once the calculation has finished. Many of these calculation are large and do not provide any further necessary information. These include ``.chk``, ``.d2e``, ``.int``, ``.rwf``, and ``.skr`` files. Removing these files saves space and makes data management much easier. 

To use this program, move into the overall path that contains all your Gaussian jobs, and type ``ECCP tidy`` into the terminal:

```bash
cd path/to/gaussian/jobs
ECCP tidy
```

!!! note

	Note that this program will only remove the files of jobs that have successfully finished. ``ECCP`` will check to make sure the job has completed successfully before removing unnecessary files. 









