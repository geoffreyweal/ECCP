## How To Use The ECCP Program

The Electronic Crystal Calculation Prep (ECCP) program is run using a script called ``Run_ECCP.py``. This contains information about all the crystals you want to obtain molecules and dimers for, as well as the parameters required to create atomic transition charge (ATC), reorganisation energy (RE), and electronic energy transfer (EET) ``Gaussian``/``ORCA`` input files for those molecules and dimers. 

An example of this script is shown below. General recommendation for settings are given in [Recommended Settings](Recommended_Settings.md#Recommended-Settings).

```python title="Run_ECCP.py" linenums="1"
--8<-- "docs/Files/Run_ECCP.py"
```


The first set of parameters allow you to give the path to the crystal file you want to run with the ECCP, set which methods you want to use to extract molecules and dimers from the crystal structure, and determine if you would like to include an environment and solvents. 

* ``filepath`` (*str.*): The path to the crystal file of organic photovoltaic (OPV) crystal or chemical system crystal you would like to analyse with this program.
* ``bonds_to_ignore`` (*str*): This is the path to a text file (generally called ``bonds_to_ignore.txt``) that contains the list of bonded atom index pairs that you want to ignore as bonds.
	* You can determine the indices of the atoms in the bonds you want to ignore using a GUI program, such as ``ase gui``. 

!!! example

	An example of this ``bonds_to_ignore.txt`` is shown below:

	```txt title="bonds_to_ignore.txt"
	0 1
	4 5
	7 11
	23 7
	```

* ``make_molecule_method`` (*str.*): This is the method you would like to use to reconnect a molecule from a crystal that has been cut due to the molecule spanning across a unit cell boundary. See [Methods available for ``make_molecule_method``](Molecule_and_Dimer_Methods.md#methods-available-for-make_molecule_method) for more information about the ``make_molecule_method`` methods available. 
* ``molecule_equivalence_method`` (*dict.*): This dictionary indicates the method that you want to use to determine which molecules in the crystal are unique (and which ones are equivalent to each other). See [Methods available for ``molecule_equivalence_method``](Molecule_and_Dimer_Methods.md#methods-available-for-molecule_equivalence_method) for more information about the ``molecule_equivalence_method`` methods available. 

* ``make_dimer_method`` (*dict.*): This dictionary contains the information required for determining how dimers are determined/obtained by this program. See [Methods available for ``make_dimer_method``](Molecule_and_Dimer_Methods.md#methods-available-for-make_dimer_method) for more information about the ``make_dimer_method`` methods available. 
* ``dimer_equivalence_method`` (*dict.*): This dictionary contains information required for determining which dimers are equivalent and which are unique. See [Methods available for ``dimer_equivalence_method``](Molecule_and_Dimer_Methods.md#methods-available-for-dimer_equivalence_method) for more information about the ``dimer_equivalence_method`` methods available. 

* ``environment_settings`` (*dict.*): This dictionary contains information about if you would like to include the environment about molecules and dimers in your ``Gaussian``/``ORCA`` calculations.  See [Environment Settings](Environment_Settings.md#environment-settings) for more information about the ``environment_settings`` settings. 
* ``remove_solvents`` (*bool.*): If ``True``, will include solvents as monomers in your dimers. If ``False``, will not include solvents as monomers in your dimers.


The second set of parameters allows you to determine the names of folders to save ECCP files to and other general housekeeping parameters. 

* ``overall_folder_suffix_name`` (*str.*): This is the suffix that you can add to the ``'ECCP_Data'`` folder name if you need to distinguish it in any way. Folder created will be called ``'ECCP_Data_XXX'``, where ``XXX`` is the suffix name. If you don't need to add a suffix, set this to ``overall_folder_suffix_name = ''``

* ``no_of_cpus`` (*int.*): This is the number of CPUs that you would like ECCP to run. 

	* **NOTE**: This is different to the number of CPUs you would like to be use in your ``Gaussian``/``ORCA`` calculations (see below). This variable is purely the number of CPUs that ECCP uses to run.

!!! warning

	You can only use multiple cpus for running the ECCP program on Linux. Mac seems to have a problem when you try to run ECCP on multiple cpus. 

	* If you are running the ECCP program pon a Mac, set ``no_of_cpus`` to 1.  


The third set of parameters involves indicating which types of ``Gaussian``/``ORCA`` files you would like to create for the molecules and dimers obtained with the ECCP program.

* Molecule-based calculations:

	* ``atc_file_creation_information``: This variable indicates if you want to write atomic transition charge (ATC) ``Gaussian``/``ORCA`` input files for the molecules found in this crystal. If you would like to perform ATC calculations, see [Dictionaries needed for obtaining atomic transition charge (ATC) input files for Monomers](Using_The_ECCP_Program.md#dictionaries-needed-for-obtaining-atomic-transition-charge-atc-input-files-for-monomers) below. If you do not want to perform ATC calculations, set this variable to ``None``. 
	* ``re_file_creation_information``: This variable indicates if you want to write reorganisation energy (RE) ``Gaussian``/``ORCA`` input files for the molecules found in this crystal. If you would like to perform RE calculations, see [Dictionaries needed for obtaining atomic transition charge (ATC) input files for Monomers](Using_The_ECCP_Program.md#dictionaries-needed-for-obtaining-structural-and-reorganisation-energies-re-input-files-for-monomers) below. If you do not want to perform RE calculations, set this variable to ``None``. 
	* ``fc_file_creation_information``: This variable indicates if you want to write franck-condon (FC) and huang-rhys (HR) ``Gaussian``/``ORCA`` input files for the molecules found in this crystal. If you would like to perform FC/HR calculations, see [Dictionaries needed for obtaining atomic transition charge (ATC) input files for Monomers](Using_The_ECCP_Program.md#dictionaries-needed-for-obtaining-franck-condonhuang-rhys-factor-fchr-input-files-for-monomers) below. If you do not want to perform FC/HR calculation, set this variable to ``None``. Note: If you do want to obtain FC/HR factors, you also need to perform the calculations needed for obtaining reoriganisation energies (RE) as well. 

* Dimer-based calculations:

	* ``eet_file_creation_information``: This variable indicates if you want to write electronic energy transfer (EET) ``Gaussian``/``ORCA`` input files for the dimers found in this crystal. If you would like to perform EET calculations, see [Dictionaries needed for obtaining atomic transition charge (ATC) input files for Monomers](Using_The_ECCP_Program.md#dictionaries-needed-for-obtaining-electronic-energy-transfer-eet-input-files-for-dimers) below. If you do not want to perform EET calculations, set this variable to ``None``. 
	* ``dimer_eigendata_file_creation_information``: NOT IMPLEMENTED YET -> This variable indicates if you want to write intermolecular charge transfer (ICT) ``Gaussian``/``ORCA`` input files for the dimers found in this crystal. ``Gaussian``/``ORCA`` files will be written to extract the eigendata matrices (such as overlap intergral matrices, molecular orbital (MO) energies, and MO coefficients) for the monomers and dimers found in this crystal. If you would like to perform ICT calculations, see [Dictionaries needed for obtaining atomic transition charge (ATC) input files for Monomers](Using_The_ECCP_Program.md#dictionaries-needed-for-obtaining-intermolecular-charge-transfer-ict-input-files-for-dimers) below. If you do not want to perform ICT calculations, set this variable to ``None``. 


## Dictionaries needed for obtaining atomic transition charge (ATC) input files for Monomers

If you would like to obtain input files for performing atomic transition charge (ATC) calculations, you will want to provide three dictionaries for the ``get_molecule_atcs`` variable. These are:

* ``calculation_parameters_for_atomic_transition_charges`` (*dict.*): This dictionary contains information required for the ``Gaussian``/``ORCA`` files for performing ATC calculations.
* ``submission_information_for_atomic_transition_charges`` (*dict.*): This dictionary contains information required for making the ``submit.sl`` file for submitting ATC jobs to slurm. 
* ``submission_information_for_multiwfn`` (*dict.*): This dictionary contains information required for the ``multiwfn_submit.sl`` file. This submit script will run the Multiwfn component of the calculation, where the ``.wfn`` created by ``Gaussian``/``ORCA`` is used to perform the ATC calculation to obtain the ``.chg`` file of the molecule. This ``.chg`` file contains the atomic transition charges for the molecule. 

If you would like to perform ATC calculations for the molcules, you want to set ``atc_file_creation_information`` as:

```python
atc_file_creation_information = (calculation_parameters_for_atomic_transition_charges, submission_information_for_atomic_transition_charges, submission_information_for_multiwfn)
```

See [Calculation Parameter Settings](Calculation_Parameters_and_Submission_Information_Settings.md#tags-needed-for-the-calculation_parameters_for_franck_condon_factors-dictionary) and [Submission Information Settings](Calculation_Parameters_and_Submission_Information_Settings.md#tags-needed-for-the-submission_information-dictionaries) for more information about what settings you want to include in these dictionaries. 

If you do not want to perform ATC calculations for the molcules, you want to set ``atc_file_creation_information = None``:

```python
atc_file_creation_information = None
```

Examples of the ``calculation_parameters_for_atomic_transition_charges``, ``submission_information_for_atomic_transition_charges``, and ``submission_information_for_multiwfn`` dictionaries as given below:

```python title="Examples of dictionaries needed for obtaining atomic transition charges (ATC) input files" show_lines="54:86" linenums="54"
--8<-- "docs/Files/Run_ECCP.py"
```


## Dictionaries needed for obtaining structural and reorganisation energies (RE) input files for Monomers

If you would like to obtain the ``Gaussian``/``ORCA`` input files required for running structural optimisations and obtain ground and excited state energies, as well as to obtain reorganisation energies (RE), you will want to provide two dictionaries for the ``re_file_creation_information`` variable. These are:

* ``calculation_parameters_for_reorganisation_energy`` (*dict.*): This dictionary contains information required for the ``Gaussian``/``ORCA`` files for performing RE calculations.
* ``submission_information_for_reorganisation_energy`` (*dict.*): This dictionary contains information required for making the ``submit.sl`` file for submitting RE jobs to slurm. 

If you would like to perform RE calculations for the dimers, you want to set ``re_file_creation_information`` as:

```python
re_file_creation_information = (calculation_parameters_for_reorganisation_energy, submission_information_for_reorganisation_energy)
```

See [Calculation Parameter Settings](Calculation_Parameters_and_Submission_Information_Settings.md#tags-needed-for-the-calculation_parameters_for_franck_condon_factors-dictionary) and [Submission Information Settings](Calculation_Parameters_and_Submission_Information_Settings.md#tags-needed-for-the-submission_information-dictionaries) for more information about what settings you want to include in these dictionaries. 

If you do not want to perform RE calculations for the molcules, you want to set ``re_file_creation_information = None``:

```python
re_file_creation_information = None
```

Examples of the ``calculation_parameters_for_reorganisation_energy`` and ``submission_information_for_reorganisation_energy`` dictionaries as given below:

```python title="Examples of dictionaries needed for obtaining structural and reorganisation energies (RE) input files" show_lines="89:111" linenums="89"
--8<-- "docs/Files/Run_ECCP.py"
```


## Dictionaries needed for obtaining franck-condon/huang-rhys factor (FC/HR) input files for Monomers

If you would like to obtain input files for performing franck-condon/huang-rhys factor (FC/HR) calculations, **you will want to provide the two dictionaries for the** ``re_file_creation_information`` **variable (see above)** as well as two dictionaries for the ``fc_file_creation_information`` variable. These are:

* ``calculation_parameters_for_franck_condon_factors`` (*dict.*): This dictionary contains information required for the ``Gaussian``/``ORCA`` files for performing FC/HR calculations.
* ``submission_information_for_franck_condon_factors`` (*dict.*): This dictionary contains information required for making the ``submit.sl`` file for submitting FC/HR jobs to slurm. 

If you would like to perform FC/HR calculations for the dimers, you want to set ``fc_file_creation_information`` as:

```python
fc_file_creation_information = (calculation_parameters_for_franck_condon_factors, submission_information_for_franck_condon_factors)
```

See [Calculation Parameter Settings](Calculation_Parameters_and_Submission_Information_Settings.md#tags-needed-for-the-calculation_parameters_for_franck_condon_factors-dictionary) and [Submission Information Settings](Calculation_Parameters_and_Submission_Information_Settings.md#tags-needed-for-the-submission_information-dictionaries) for more information about what settings you want to include in these dictionaries. 

If you do not want to perform FC/HR calculations for the molcules, you want to set ``fc_file_creation_information = None``:

```python
fc_file_creation_information = None
```

Examples of the ``calculation_parameters_for_franck_condon_factors`` and ``submission_information_for_franck_condon_factors`` dictionaries as given below:

```python title="Examples of dictionaries needed for obtaining franck-condon/huang-rhys factor (FC/HR) input files" show_lines="114:135" linenums="114"
--8<-- "docs/Files/Run_ECCP.py"
```


## Dictionaries needed for obtaining electronic energy transfer (EET) input files for Dimers

If you would like to obtain input files for performing electronic energy transfer (EET) calculations, you will want to provide two dictionaries for the ``eet_file_creation_information`` variable. These are:

* ``calculation_parameters_for_electronic_energy_transfer`` (*dict.*): This dictionary contains information required for the ``Gaussian``/``ORCA`` files for performing EET calculations.
* ``submission_information_for_electronic_energy_transfer`` (*dict.*): This dictionary contains information required for making the ``submit.sl`` file for submitting EET jobs to slurm. 

If you would like to perform EET calculations for the dimers, you want to set ``eet_file_creation_information`` as:

```python
eet_file_creation_information = (calculation_parameters_for_electronic_energy_transfer, submission_information_for_electronic_energy_transfer)
```

See [Calculation Parameter Settings](Calculation_Parameters_and_Submission_Information_Settings.md#tags-needed-for-the-calculation_parameters_for_franck_condon_factors-dictionary) and [Submission Information Settings](Calculation_Parameters_and_Submission_Information_Settings.md#tags-needed-for-the-submission_information-dictionaries) for more information about what settings you want to include in these dictionaries. 

If you do not want to perform EET calculations for the dimers, you want to set ``eet_file_creation_information = None``:

```python
eet_file_creation_information = None
```

Examples of the ``calculation_parameters_for_electronic_energy_transfer`` and ``submission_information_for_electronic_energy_transfer`` dictionaries as given below:

```python title="Examples of dictionaries needed for obtaining electronic energy transfer (EET) input files" show_lines="138:161" linenums="138"
--8<-- "docs/Files/Run_ECCP.py"
```


## Dictionaries needed for obtaining intermolecular charge transfer (ICT) input files for Dimers

**Not Functioning Yet**

Intermolecular charge transfer coupling values can be obtained for charge transfer between monomers (in a dimer) by obtaining the eigendata from ``Gaussian``/``ORCA`` claculations upon the dimer and the two monomers that are apart of the dimer. This includes obtaining overlap intergral matrices, molecular orbital (MO) energies, and MO coefficients for the monomers and dimers found in this crystal.

If you would like to obtain input files for performing intermolecular charge transfer (ICT) calculations, you will want to provide two dictionaries for the ``dimer_eigendata_file_creation_information`` variable. These are:

* ``calculation_parameters_for_eigendata`` (*dict.*): This dictionary contains information required for the ``Gaussian``/``ORCA`` files for obtaining eigendata from your monomers and dimers, required for obtaining ICT coupling values.
* ``submission_information_for_eigendata`` (*dict.*): This dictionary contains information required for making the ``submit.sl`` file for submitting eigendata gathering jobs to slurm. 

If you would like to perform ICT calculations for the dimers, you want to set ``ict_file_creation_information`` as:

```python
ict_file_creation_information = (calculation_parameters_for_intermolecular_charge_transfer, submission_information_for_intermolecular_charge_transfer)
```

See [Calculation Parameter Settings](Calculation_Parameters_and_Submission_Information_Settings.md#tags-needed-for-the-calculation_parameters_for_franck_condon_factors-dictionary) and [Submission Information Settings](Calculation_Parameters_and_Submission_Information_Settings.md#tags-needed-for-the-submission_information-dictionaries) for more information about what settings you want to include in these dictionaries. 

If you do not want to perform ICT calculations for the dimers, you want to set ``ict_file_creation_information = None``:

```python
ict_file_creation_information = None
```

NOTE: Both ``calculation_parameters_for_intermolecular_charge_transfer`` and ``submission_information_for_intermolecular_charge_transfer`` are usually the same as  ``calculation_parameters_for_electronic_energy_transfer`` and ``submission_information_for_electronic_energy_transfer``, respectively. See below for examples of ``calculation_parameters_for_intermolecular_charge_transfer`` and ``submission_information_for_intermolecular_charge_transfer``: 

```python title="Examples of dictionaries needed for obtaining intermolecular charge transfer (ICT) input files" show_lines="164:171" linenums="164"
--8<-- "docs/Files/Run_ECCP.py"
```


## Examples of Input Files

The folder called ``Examples`` contains all the files that are needed to run the ECCP program. This includes examples of the ``Run_ECCP.py`` file. [Click here to access these files on Github](https://github.com/geoffreyweal/ECCP/tree/main/Examples).


