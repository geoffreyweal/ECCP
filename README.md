# The Electronic Crystal Calculation Prep (ECCP) Program

Authors: Dr. Geoffrey Weal<sup>\*,†</sup>, Dr. Chayanit Wechwithayakhlung<sup>†</sup>, Dr. Josh Sutton<sup>\*</sup>, Dr. Daniel Packwood<sup>†</sup>, Dr. Paul Hume<sup>\*</sup>, Prof. Justin Hodgkiss<sup>\*</sup>

<sup>\*</sup> Victoria University of Wellington, Wellington, New Zealand; The MacDiarmid Institute for Advanced Materials and Nanotechnology, Wellington, New Zealand. 

<sup>†</sup> Institute for Integrated Cell-Material Sciences (iCeMS), Kyoto University, Kyoto, Japan.

Group pages: https://people.wgtn.ac.nz/paul.hume/grants, https://www.packwood.icems.kyoto-u.ac.jp/, https://people.wgtn.ac.nz/justin.hodgkiss/grants


## What is the Electronic Crystal Calculation Prep (ECCP) Program

The Electronic Crystal Calculation Prep (ECCP) program is designed create the files needed to perform DFT calculations upon the molecules in molecular crystals. These DFT calculations allow the user to understand excited state properties, reorganisation energy, electronic coupling, and exciton and charge diffusion.

The Electronic Crystal Calculation Prep (ECCP) program is designed create the ``Gaussian``/``ORCA`` files needed to perform DFT calculations upon the molecules in molecular crystals. The DFT calculations allow the user to calculate the electronic properties of molecules, such as excited state properties and reorganisation energies, as well as electronic intermolecular properties such as electronic coupling, and exciton and charge diffusion. 

The ECCP program works as follows: 

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

## Installation

It is recommended to read the installation page before using the ECCP program. See [Installation: Setting Up ECCP and Pre-Requisites Packages](https://geoffreyweal.github.io/ECCP/Installation) for more information. Note that you can install ECCP through ``pip3`` and ``conda``. 

## Guide To Using ECCP

The ECCP program is one in a series of programs that are designed to be used in the workflow shown below. After you have installed ECCP, see [How To Use The ECCP Program](https://geoffreyweal.github.io/ECCP/Guide_to_Using_the_ECCP_Program/Using_The_ECCP_Program) to learn about how to use this program. 

## The Grand Scheme

The ECCP program is used as part of a grand scheme for calculating the excited-state electronic properties of molecules in a crystal. This includes simulations of exciton and charge diffusion through crystal structures, in particular for organic molecules (but not limited to them). This scheme is shown below, along with where the ECCP program is used in this scheme. 

<img alt="Schematic of Grand Scheme" src="Documentation/docs/Shared_Images/Grand_Scheme/Grand_Scheme.png" />

## Websites and Github Repositories for All Associated Programs

### Instructional Websites

* ACSD: https://geoffreyweal.github.io/ACSD
* ReCrystals: https://geoffreyweal.github.io/ReCrystals
* RSGC: https://geoffreyweal.github.io/RSGC
* ReJig: https://geoffreyweal.github.io/ReJig
* ECCP: https://geoffreyweal.github.io/ECCP
* EKMC: https://geoffreyweal.github.io/EKMC
* SORE: https://geoffreyweal.github.io/SORE
* SUMELF: https://geoffreyweal.github.io/SUMELF

### Github Repositories

* ACSD: https://github.com/geoffreyweal/ACSD
* ReCrystals: https://github.com/geoffreyweal/ReCrystals
* RSGC: https://github.com/geoffreyweal/RSGC
* ReJig: https://github.com/geoffreyweal/ReJig
* ECCP: https://github.com/geoffreyweal/ECCP
* EKMC: https://github.com/geoffreyweal/EKMC
* SORE: https://github.com/geoffreyweal/SORE
* SUMELF: https://github.com/geoffreyweal/SUMELF

