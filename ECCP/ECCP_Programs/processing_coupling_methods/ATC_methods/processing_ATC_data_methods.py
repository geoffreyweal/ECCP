'''
Geoffrey Weal, processing_OPV_Dimer_data.py, 9/3/22

This script contains methods for processing_OPV_Dimer_data.py

'''
import os
from tqdm import tqdm
from ase.io import read

from SUMELF import obtain_graph
from ECCP.ECCP_Programs.processing_coupling_methods.ATC_methods.invariance_method import assigned_ATCs_to_molecules_invariance_method

def is_this_calc_an_atc_calc(path_to_root):
    """
    This method will look through the output.log file to see if this job was an eet job

    Parameters
    ----------
    root : str.
        This is the path to where the calculation files are located

    Returns
    -------
    True if this calc is an eet calc, False if not.
    """

    return ('output.log' in os.listdir(path_to_root)) and ('output.chg' in os.listdir(path_to_root))
    
# -----------------------------------------------------------------

def get_structurally_equivalent_molecules(path_to_txtfile, molecule_nos):
    """
    This method is designed to provide a dictionary to indicate which molecules are structurally equivalent.

    Parameters
    ----------
    path_to_txtfile : str.
        This is the path to the Structurally_Unique_Molecule_Information.txt file.

    Returns
    -------
    structurally_equivalent_molecules : dict.
        This dictionary indicates which molecules are sturcturally equivalent, as given in the Structurally_Unique_Molecule_Information.txt file.
    molecule_nos : list
        These are the numbers/names of the molecules that are structurally unique
    """

    # First, initialise the structurally_equivalent_molecules dictionary with structurally unique molecules from molecule_nos.
    structurally_equivalent_molecules = {no: no for no in molecule_nos}

    # Second, look through Structurally_Unique_Molecule_Information.txt and record structurally equivalent molecules.
    with open(path_to_txtfile+'/Structurally_Unique_Molecule_Information.txt') as SUMITXT:
        SUMITXT.readline()
        for line in SUMITXT:
            equivalent, _, unique = line.rstrip().split()
            equivalent = int(equivalent); unique = int(unique)
            if equivalent in structurally_equivalent_molecules.keys():
                import pdb; pdb.set_trace()
                raise Exception('huh?')
            structurally_equivalent_molecules[equivalent] = unique

    # Third, return structurally_equivalent_molecules
    return structurally_equivalent_molecules

def obtain_all_molecules_with_ATC_charges_in_crystal(path_to_molecules, ATC_coupling_data_for_crystal, structurally_equivalent_molecules):
    """
    This method is designed to attach the ATC charges to all the molecules in the crystal from the unique molecules.
    """

    # First, initiate the dictionary that will hold all the molecules in the crystal with associated ATC charges. 
    molecules_in_crystal        = {a_functional_basis_set: {} for a_functional_basis_set in ATC_coupling_data_for_crystal.keys()}
    molecules_in_crystal_graphs = {}

    # Second, obtain one of the functional + basis set to contain structures for later on in this algorithm.
    a_functional_basis_set = list(ATC_coupling_data_for_crystal.keys())[0]

    # Third, find each molecule in the "All_Molecules" folder, add ATC charges, and add it to the molecules_in_crystal dictionary.
    print('Adding ATC charges to all molecules in the crystal.')
    pbar = tqdm(sorted([file for file in os.listdir(path_to_molecules+'/'+'All_Molecules') if file.endswith('.xyz')]))
    for file in pbar:

        # 3.1: Write description
        pbar.set_description(file)

        # 3.2: Read the molecule of interest from file.    
        molecule = read(path_to_molecules+'/'+'All_Molecules'+'/'+file)
        molecule_graph = obtain_graph(molecule)

        # 3.3: Obtain the molecule number for this molecule.
        molecule_no = int(file.replace('molecule_','').replace('.xyz',''))

        if molecule_no in molecules_in_crystal_graphs:
            raise Exception('huh?')
        molecules_in_crystal_graphs[molecule_no] = molecule_graph

        # 3.4: Obtain the unique molecule that is structurally equivalent to this molecule.
        unique_molecule_no = structurally_equivalent_molecules[molecule_no]
        
        # 3.5: If not (unique_molecule_no == molecule_no), obtain which atoms in the ATC go with each atom in the molecule.  
        if not (unique_molecule_no == molecule_no):
            unique_molecule = ATC_coupling_data_for_crystal[a_functional_basis_set][unique_molecule_no]
            unique_molecule_graph = obtain_graph(unique_molecule)
            molecule_to_atc = assigned_ATCs_to_molecules_invariance_method([unique_molecule], [unique_molecule_graph], [molecule], [molecule_graph], max_disparity=None)
            molecule_to_atc = molecule_to_atc[0][2]

        # 3.6: For each functional and basis set:
        for a_functional_basis_set in ATC_coupling_data_for_crystal.keys():

            # 3.5.1: Obtain the unique molecule with ATC charges.
            unique_molecule = ATC_coupling_data_for_crystal[a_functional_basis_set][unique_molecule_no]

            # 3.5.2: Obtain the ATC charges for each atom.
            if unique_molecule_no == molecule_no:
                ATC_charges_to_add = unique_molecule.get_initial_charges()
            else:
                ATC_charges_to_add = [unique_molecule[ATC_index].charge for ATC_index in molecule_to_atc]

            # 3.5.3: Add the ATC charges to the current molecule. 
            for index, ATC_charge in zip(range(len(molecule)), ATC_charges_to_add):
                molecule[index].charge = ATC_charge

            # 3.5.4: Record the ATC charge atoms object to the molecules_in_crystal dictionary. 
            molecules_in_crystal[a_functional_basis_set][molecule_no] = molecule
            
    # Fourth, return molecules_in_crystal
    return molecules_in_crystal, molecules_in_crystal_graphs

