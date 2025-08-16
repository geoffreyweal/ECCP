# Environment Settings

**NOTE: NOT IMPLEMENTED YET.**

In some cases, you may like to include the neighbouring molecules around each molecule and dimer in Gaussian/ORCA calculations. This may be desired if you want to obtain coupling values that are influenced by the local electric field surrounding each molecule and dimer. For example, ADD EXAMPLE HERE.

There are a few parameters that are included in this dictinoary:

* ``'include_environment_in_molecule_calcs'`` (*bool.*): Set this to ``True`` if you would like include neighbours surrounding each molecule in the crystal in your Gaussian/ORCA calculations. 
* ``'include_environment_in_dimer_calcs'`` (*bool.*): Set this to ``True`` if you would like include neighbours surrounding each dimer in the crystal system in your Gaussian/ORCA calculations. 

If you set either ``'include_environment_in_molecule_calcs'`` or ``'include_environment_in_dimer_calcs'`` as true, you will need to include in your dictionary the following variable:

* ``'max_environment_distance'`` (*float*): This is the maximum distance between any atoms in each molecules from each other for those two molecules to be considered neighbours.

NOTE: The Nearest Atoms method is used to determine neighbouring molecules in the crystal. 

