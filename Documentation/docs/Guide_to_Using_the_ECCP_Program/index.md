# Guide To Using The ECCP Program

The ECCP program is designed allow the user to create ``Gaussian``/``ORCA`` files for calculating electronic properties of molecules in crystals, including exciton and charge diffusion. 

!!! info

	The ECCP program works by separating the crystal into molecules, combining them into dimers of neighbouring molecules, and creating DFT files for calculating electronic energy coupling between molecules, as well as reorganisation energies of molecules caused by a molecule being excited from the ground to excited state. 

The ECCP program works as following; given an input crystal file for a molecular crystal, the ECCP program will: 

1. Determine all the individual molecules in the crystal.
2. Remove solvents if desired.
3. Determine which molecules neighbour each other, including those that neighbour each other across different cells. 
4. Save all the molecules to disk, including ``Gaussian``/``ORCA`` files for obtaining structural energies, reorganisation energies (RE) and atomic transition charges (ATC). 
5. If desired, determine which of those molecules are structurally and conformationally equivalent to each other. 

	* Save which molecules are structurally equivalent and conformationally equivalent into a text file.
	* Create ``Gaussian``/``ORCA`` files for obtaining structural and reorganisation energies (with conformationally equivalent molecules), as well as atomic transition charges (with structurally equivalent molecules).

6. Obtain all the dimers within the crystal.
7. Save all the dimers to disk, including ``Gaussian``/``ORCA`` files for obtaining electronic energy transfer (EET) and intermolecular charge transfers (ICT) coupling values. 
8. If desired, determine which of those dimers are structurally equivalent to each other. 

	* Save those structurally equivalent dimers into a text file.
	* Create ``Gaussian``/``ORCA`` files for obtaining electronic energy transfer (EET) and intermolecular charge transfers (ICT) coupling values with those structurally equivalent dimers.
	
9. Write the summaries of all the molecules and dimers extracted from the crystal. 

#### Notes about structurally and conformationally equivalent molecules

* Structurally equivalent means two molecules are exactly the same, include the position of all non-hydrogen atoms
* Conformationally equivalent means two molecules are the same, but the exact position of atoms does not necessarily need to be the same. 