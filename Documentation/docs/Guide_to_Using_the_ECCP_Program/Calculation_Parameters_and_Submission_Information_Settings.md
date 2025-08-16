# Calculation Parameters and Submission Information Settings

The following provide information about the settings and other advice for the Gaussian/ORCA parameters and submission information dictionaries required for running this program.


## Tags needed for the ``calculation_parameters`` dictionaries

The ``calculation_parameters_for_atomic_transition_charges``, ``calculation_parameters_for_reorganisation_energy``, ``calculation_parameters_for_electronic_energy_transfer``, and ``calculation_parameters_for_intermolecular_charge_transfer`` dictionaries allow you to give the parameters needed for the Gaussian/ORCA input files. These include:

* ``'calc_software'`` (*str.*): This is the program that you would like to use to perform your quantum calculations. This can be either ``'Gaussian'`` or ``'ORCA'``. 
* ``'mem'`` (*str.*): This is the amount of memory that is used by Gaussian. **Gaussian needs this variable. ORCA does not use this variable**.
* ``'method'`` (*str.*): The level of theory you want to use (i.e. the functional). Example: ``'wB97XD'`` for Gaussian, or ``wB97X-D3`` for ORCA. **Note that the input names for the same functionals may be different for Gaussian and ORCA. Check the manuals before running the ECCP program to make sure you are calculating molecules and dimers with the desired functional.**
* ``'basis'`` (*str.*): The basis set you want to use. Example: ``'6-311G(d)'``. **Note that the input names for the same basis sets may be different for Gaussian and ORCA. Check the manuals before running the ECCP program to make sure you are calculating molecules and dimers with the desired basis set.**
* ``'td_settings'`` (*str.*): This allows you to provide the settings for performing TD-DFT calculations. For example, if you want to perform TD-DFT calculations with the Tamm-Dancoff approximation and solve for 3 states, you want to set this to ``'tda(nstates=3)'``. Default: ``'TD'``
* ``'temp_folder_path'`` (*str.*): This is the scratch drive to place files into. This is useful if you want the Gaussian/ORCA temporary files to be saved onto a local disk rather than a NFS if running on a computer cluster. Set to ``None`` if you do not want to give this. Default: ``None``. 

Other parameters that can be given in the ``calculation_parameters_for_atomic_transition_charges``, ``calculation_parameters_for_electronic_energy_transfer``, and ``calculation_parameters_for_intermolecular_charge_transfer`` dictionaries are:

* ``'extra'`` (*str.*): These are any extra tags that you need to add. For example in Gaussian: ``'# maxdisk=2TB scf=(xqc,maxcycle=512) scrf=(iefpcm,read)'``
* ``'addsec'`` (*str.*): These are the lines that you want to write to the ``.gjf`` file that come after the atomic positions. For example, if you added ``'scrf=(iefpcm,read)'`` in your ``'extra'`` entry, you will want to read in the relative electric permittivity after the atomic positions by setting ``calculation_parameters_for_atomic_transition_charges['addsec'] = 'Eps=4.0'``, where in this example we have set the relative electric permittivity to 4.0. **This tag is only needed for Gaussian. This is because all the variable you would write in here for Gaussian you write in** ``calculation_parameters['extra']`` **for ORCA.**

You can also specify where you would like certain files to be placed, particularly if you want to place some of the temporary Gaussian/ORCA files into a SCRATCH directory. If you want all the temporary files to be placed in a general SCRATCH files, provide a entry for ``'scratch_path'`` in ``calculation_parameters`` dictionary: 

* ``'temp_folder_path'`` (*str.*): This is the place that you want to store temporary Gaussian/ORCA files in (the scratch drive). The ECCP will write the script for that appropriate folders are made in your ``'temp_folder_path'`` that temporary Gaussian/ORCA files will be saved to. Checkpoint files will be saved to the same place as the ``.gjf`` file to allow you to easily resubmit incomplete jobs. 
* ``'place_chk_in_temp'`` (*bool.*): If you have given a directory for ``'temp_folder_path'``, you can choose to also place the checkpoint file in the ``'temp_folder_path'`` directory as well if this is set to ``True``. If you want the checkpoint file to be place in the same place as the ``.gjf`` (to allow you to easily continue the calculation if it does not finish in time), set this to ``False``. Default: ``False``

Examples of the ``calculation_parameters_for_atomic_transition_charges``, ``calculation_parameters_for_reorganisation_energy``, ``calculation_parameters_for_electronic_energy_transfer``, and ``calculation_parameters_for_intermolecular_charge_transfer`` dictionaries are given below:


```python title="Examples of calculation_parameters_for_atomic_transition_charges" show_lines="54:63" linenums="54"
--8<-- "docs/Files/Run_ECCP.py"
```

```python title="Examples of calculation_parameters_for_reorganisation_energy" show_lines="89:97" linenums="89"
--8<-- "docs/Files/Run_ECCP.py"
```

```python title="Examples of calculation_parameters_for_electronic_energy_transfer" show_lines="138:147" linenums="138"
--8<-- "docs/Files/Run_ECCP.py"
```

```python title="Examples of calculation_parameters_for_intermolecular_charge_transfer" show_lines="164:165" linenums="164"
--8<-- "docs/Files/Run_ECCP.py"
```


## Tags needed for the ``calculation_parameters_for_franck_condon_factors`` dictionary

THIS WILL PROBABLY CHANGE: UPDATE THIS WHEN THIS PROCESS HAS BEEN SORTED.

The Franck-Condon calculation performed by the Gaussian/ORCA program is a bit different to other calculations, as all the calculations that need to be done are performed by the ``Reorganisation Energy`` component. Here, Gaussian will use the frequency checkpoint files created during the ground and excited state optimisation calculations (performed by the ``Reorganisation Energy`` component). For this reason, we don't need to specify all the inputs as in [Calculation Parameter Settings](Calculation_Parameters_and_Submission_Information_Settings.md#tags-needed-for-the-calculation_parameters_for_franck_condon_factors-dictionary). The inputs needed for ``submission_information_for_franck_condon_factors`` are:

* ``'mem'`` (*str.*): This is the amount of memory that is used by Gaussian. Recommended to set this to ``'8GB'``
* ``'method'`` (*str.*): The level of theory you want to use (i.e. the functional). Example: ``'wB97XD'``
* ``'basis'`` (*str.*): The basis set you want to use. Example: ``'6-311G(d)'``

Note: Recommendations for ``submission_information_for_franck_condon_factors``, ``submission_information_for_franck_condon_factors['cpus_per_task'] = 2`` and ``submission_information_for_franck_condon_factors['mem'] = '16GB'``. 

An example of ``calculation_parameters_for_franck_condon_factors`` is given below

```python title="Example of calculation_parameters_for_franck_condon_factors" show_lines="115:121" linenums="115"
--8<-- "docs/Files/Run_ECCP.py"
```


### Advice for settings for the ``calculation_parameters`` dictionaries

The Gaussian jobs that use the ``calculation_parameters_for_atomic_transition_charges``, ``calculation_parameters_for_reorganisation_energy``, ``calculation_parameters_for_electronic_energy_transfer``, and ``calculation_parameters_for_intermolecular_charge_transfer`` dictionaries perform single point calculations, as well as other calculations for TD-DFT calculations. The Gaussian jobs that use the ``calculation_parameters_for_reorganisation_energy`` dictionary perform geometry optimisation calculations, which are far more computationally intensive to perform. 

The amount of memory that you will need to use for ``calculation_parameters_for_reorganisation_energy`` should be double or more needed for ``calculation_parameters``.

Recommendations for ``calculation_parameters_for_franck_condon_factors`` and ``submission_information_for_franck_condon_factors`` dictionaries: 

* ``calculation_parameters_for_franck_condon_factors['mem'] = '14GB'``
* ``submission_information_for_franck_condon_factors['mem'] = '16GB'``
* ``submission_information_for_franck_condon_factors['cpus_per_task'] = 2``


## Tags needed for the ``submission_information`` dictionaries

The ``submission_information_for_atomic_transition_charges``, ``submission_information_for_multiwfn``, ``submission_information_for_reorganisation_energy``, ``submission_information_for_franck_condon_factors``, ``submission_information_for_electronic_energy_transfer``, and ``dimer_eigendata_file_creation_information`` dictionaries allow you to place the parameters needed for the ``submit.sl`` , ``multiwfn_submit.sl``, and the various reorganisation ``submit.sl`` files. These include:

* ``'cpus_per_task'``: This the the number of CPU's you want to use. This information is also passed on to your Gaussian .gjf file. **THIS VARIABLE IS ONLY USED FOR GAUSSIAN CALCULATIONS**.
* ``'ntasks'``: This the the number of CPU's you want to use. This information is also passed on to your ORCA .inp file. **THIS VARIABLE IS ONLY USED FOR ORCA CALCULATIONS**.
* ``'mem'``: This is the total amount of RAM memory you want to use across all your CPUs. 
* ``'time'``: This is the amount of time you want to run this job for. Written HH:MM:SS
* ``'partition'``: This is the name of the partition you want to run your job on.
* ``'constraint'``: This assigns particular nodes to run Gaussian jobs. This is a variable that is needed at Victoria University of Wellington. See ``https://slurm.schedmd.com/sbatch.html`` for more information about this. This is set to ``'AVX'`` on the Rāpoi computer cluster at Victoria University of Wellington. 
* ``'email'``: This is the email you want to use to notify you about this job
* ``'gaussian_version'``: This is the version of gaussian you want to use. For example: ``'g16'``. This setting is required for ``submission_information`` and ``submission_information_for_reorganisation_energy``, as these processes use Gaussian to perform ATC, RE, and EET calculations. This tag is not needed for the ``submission_information_for_multiwfn`` dictionary. **THIS VARIABLE IS ONLY USED FOR GAUSSIAN CALCULATIONS**.
* ``'ORCA_version'``: This is the version of ORCA you want to use. For example: ``'ORCA/5.0.3'``. This setting is required for ``submission_information`` and ``submission_information_for_reorganisation_energy``, as these processes use ORCA to perform ATC, RE, and EET calculations. This tag is not needed for the ``submission_information_for_multiwfn`` dictionary. **THIS VARIABLE IS ONLY USED FOR GAUSSIAN CALCULATIONS**.
* ``'python_version'``: This is the version of python you want to use. For example: ``'3.8.1'`` for Python 3.8.1. This tag is only needed for the ``submission_information_for_reorganisation_energy`` dictionary, as python is only used here to create single point calculations from the results of the geometrically optimised Gaussian/ORCA calculations. 

For the ``submission_information_for_multiwfn`` dictionary, you only need to provide the ``'cpus_per_task'``/``'ntasks'``, ``'mem'``, ``'time'``, and ``'partition'`` input parameters (and optionally the ``'constraint'`` and ``'email'`` input parameters). 

Examples of the ``submission_information_for_atomic_transition_charges``, ``submission_information_for_multiwfn``, ``submission_information_for_reorganisation_energy``, ``submission_information_for_franck_condon_factors``, ``submission_information_for_electronic_energy_transfer``, and ``dimer_eigendata_file_creation_information`` dictionaries are given below:

```python title="Examples of submission_information_for_atomic_transition_charges and submission_information_for_multiwfn" show_lines="65:83" linenums="65"
--8<-- "docs/Files/Run_ECCP.py"
```

```python title="Examples of submission_information_for_reorganisation_energy" show_lines="99:108" linenums="99"
--8<-- "docs/Files/Run_ECCP.py"
```

```python title="Examples of submission_information_for_franck_condon_factors" show_lines="123:132" linenums="123"
--8<-- "docs/Files/Run_ECCP.py"
```

```python title="Examples of submission_information_for_electronic_energy_transfer" show_lines="149:158" linenums="149"
--8<-- "docs/Files/Run_ECCP.py"
```

```python title="Examples of dimer_eigendata_file_creation_information" show_lines="167:168" linenums="167"
--8<-- "docs/Files/Run_ECCP.py"
```


### Advice for settings for the ``submission_information`` dictionaries

The Gaussian/ORCA calculations that use the ``submission_information_for_atomic_transition_charges``, ``submission_information_for_electronic_energy_transfer``, and ``dimer_eigendata_file_creation_information`` dictionaries perform single point calculations, as well as other calculations for TD-DFT calculations. The Gaussian/ORCA jobs that use the ``submission_information_for_reorganisation_energy`` dictionary perform geometry optimisation calculations, which are far more computationally intensive to perform. 

The amount of cpus, memory and computational time that you should need to use for ``submission_information_for_reorganisation_energy`` should be double or more needed for ``submission_information_for_atomic_transition_charges`` and ``submission_information_for_electronic_energy_transfer`` dictionaries. 

The Multiwfn program to perform the ATC calculations is not as computationally expensive as Gaussian, and likely will finish within a fraction of time needed for the Gaussian calculation. For this reason, the amount of cpus, memory and computational time ``submission_information_for_multiwfn`` can be the same or less than for ```submission_information_for_atomic_transition_charges``, and ``submission_information_for_electronic_energy_transfer`` dictionaries. 

Recommendations for ``calculation_parameters_for_franck_condon_factors`` and ``submission_information_for_franck_condon_factors`` dictionaries: 

* ``calculation_parameters_for_franck_condon_factors['mem'] = '14GB'``
* ``submission_information_for_franck_condon_factors['mem'] = '16GB'``
* ``submission_information_for_franck_condon_factors['cpus_per_task'] = 2``



