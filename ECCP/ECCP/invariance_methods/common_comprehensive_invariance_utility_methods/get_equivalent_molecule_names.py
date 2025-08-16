"""
get_equivalent_molecule_names.py, Geoffrey Weal, 7/2/24

This script is designed to obtain all the names between equivalent moleules. 
"""
import sys
from tqdm import tqdm

import multiprocessing as mp
from tqdm.contrib.concurrent import process_map

from SUMELF import GraphMatcher

def get_equivalent_molecule_names(non_hydrogen_graphs, include_comparisons_with_itself=False, no_of_cpus=1):
    """
    This method is designed to obtain all the names between equivalent moleules. 

    Parameters
    ----------
    non_hydrogen_graphs : dict of networkx.Graphs
        This dict contains all the graphs of the molecules in the crystal, where hydrogens have been removed and added to attached atoms as node features.
    include_comparisons_with_itself : bool.
        This boolean indicates if the user want to give the ways that a molecule could map onto itself. 
            * For determining equivalent molecules, we dont want to do this, so set to False
            * For determining equivalent dimers, we do want to do this, so set this to True.
    no_of_cpus : int.
        This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.

    Returns
    -------
    equivalent_molecule_names : dict.
        This is all the ways that two molecules in non_hydrogen_graphs can map onto each other. 
    """

    # First, obtain a list of equivalent_molecule_names
    # * This list is formated (mol_name1, mol_name2): All the ways that mol_name1 maps onto mol_name2
    #
    if no_of_cpus == 1: # If the user want to perform this task with 1 cpu, use standard simple approach.
        list_of_equivalent_molecule_names = compare__non_hydrogen_graphs__with_one_cpu(non_hydrogen_graphs, include_comparisons_with_itself)
    else: 
        list_of_equivalent_molecule_names = compare__non_hydrogen_graphs__multiple_cpu(non_hydrogen_graphs, include_comparisons_with_itself, no_of_cpus)

    # Second, initialise the equivalent_molecule_names dictionary that will hold all the ways that two molecules in non_hydrogen_graphs can map on to each other. 
    equivalent_molecule_names = {}

    # Third, move data from list_of_equivalent_molecule_names onto equivalent_molecule_names in dictionary format rather than list format.
    for (mol_name1, mol_name2), all_unique_matches in list_of_equivalent_molecule_names:

        # 3.1: Check that none of the entries in list_of_equivalent_molecule_names has mol_name1
        if include_comparisons_with_itself: 
            if mol_name1 > mol_name2: # We have included a comparision of the molecule with itself. This is used for determining symmetric dimers. 
                import pdb; pdb.set_trace()
                raise Exception('Error: mol_name 1 is bigger than mol_name 2. This may indicate a programming error. Check this.')
        else:
            if mol_name1 >= mol_name2: # We have not included a comparision of the molecule with itself. This is used for determining symmetric molecules. 
                raise Exception('Error: mol_name 1 is bigger than or equal to mol_name 2. This may indicate a programming error. Check this.')

        # 3.2: Add all_unique_matches to (mol_name1,mol_name2) in equivalent_molecule_names
        equivalent_molecule_names[(mol_name1,mol_name2)] = all_unique_matches

        # 3.3: If mol_name1 is not equal to mol_name2 (because mol_name1 is larger than mol_name2), then add the reverse connections to equivalent_molecule_names for mapping molecule 2 onto molecule 1
        if mol_name1 != mol_name2:
            equivalent_molecule_names[(mol_name2, mol_name1)] = [{v: k for k, v in a_unique_match.items()} for a_unique_match in all_unique_matches]

    # Fourth, return equivalent_molecule_names
    return equivalent_molecule_names

# ----------------------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------

def compare__non_hydrogen_graphs__with_one_cpu(non_hydrogen_graphs, include_comparisons_with_itself):
    """
    This method is desiged to compare the molecules in the non_hydrogen_graphs and obtain all the ways that molecules can be mapped onto each other.

    THis method is designed to use only 1 cpu

    Parameters
    ----------
    non_hydrogen_graphs : dict of networkx.Graphs
        This dict contains all the graphs of the molecule, where hydrogens have been removed and added to attached atoms as node features
    include_comparisons_with_itself : bool.
        This boolean determines if the user want to give the ways that a molecule could map onto itself. 
            * For determining equivalent molecules, we dont want to do this, so set to False
            * For determining equivalent dimers, we do want to do this, so set this to True.

    Returns
    -------
    list_of_equivalent_molecule_names : list
        This is the list to store all the names that map molecule 1 onto molecule 2. 
    """

    # First, initialise the list to hold all the ways that molecules in non_hydrogen_graphs can map on to eachother.
    list_of_equivalent_molecule_names = []

    # Second, create a tqdm instance to notify the user what it is doing
    no_of_neighbourhood_sets = int((len(non_hydrogen_graphs) * (len(non_hydrogen_graphs) - 1)) / 2)
    pbar = tqdm(total=no_of_neighbourhood_sets,unit='molecule pair')

    # Third, compare all molecules together in non_hydrogen_graphs
    for input_data in get_inputs(non_hydrogen_graphs, include_comparisons_with_itself, list_of_equivalent_molecule_names):

        # 3.1: Write output to user.
        mol_name1 = input_data[0]
        mol_name2 = input_data[1]
        if mol_name1 == mol_name2:
            pbar.set_description('Comparing molecule '+str(mol_name1)+' with itself')
        else:
            pbar.set_description('Comparing molecules '+str(mol_name1)+' and '+str(mol_name2))

        # 3.2: Compare mol_name1 and mol_name2 in non_hydrogen_graphs together.
        compare_two_molecules_single_process(input_data)

        # 3.3: Update pbar, as one job is completed.
        pbar.update(1)

    # Fourth, close pbar.
    pbar.close()

    # Fifth, return list_of_equivalent_molecule_names
    return list_of_equivalent_molecule_names

def compare__non_hydrogen_graphs__multiple_cpu(non_hydrogen_graphs, include_comparisons_with_itself, no_of_cpus):
    """
    This method is desiged to compare the molecules in the non_hydrogen_graphs and obtain all the ways that molecules can be mapped onto each other.

    This method is designed to use multiple cpus. 

    Parameters
    ----------
    non_hydrogen_graphs : dict of networkx.Graphs
        This dict contains all the graphs of the molecule, where hydrogens have been removed and added to attached atoms as node features
    include_comparisons_with_itself : bool.
        This boolean indicates if the user want to give the ways that a molecule could map onto itself. 
            * For determining equivalent molecules, we dont want to do this, so set to False
            * For determining equivalent dimers, we do want to do this, so set this to True.
    no_of_cpus : int.
        This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.

    Returns
    -------
    list_of_equivalent_molecule_names : list
        This is the list to store all the names that map molecule 1 onto molecule 2. 
    """

    # First, create the manager to save lists to
    with mp.Manager() as manager:

        # Second, create the list to hold all the ways that molecules in non_hydrogen_graphs can map on to eachother.
        list_of_equivalent_molecule_names = manager.list()

        # Third, run the multiprocessing jobs.
        print('Obtaining neighbourhoods between molecules (Please wait until after 100%, as the process will still be running.)', file=sys.stderr)
        no_of_neighbourhood_sets = int((len(non_hydrogen_graphs) * (len(non_hydrogen_graphs) - 1)) / 2)
        #process_map(compare_two_molecules_single_process, get_inputs(non_hydrogen_graphs, include_comparisons_with_itself, list_of_equivalent_molecule_names), total=no_of_neighbourhood_sets, desc='Obtaining neighbouring pairs of molecules', unit='calc', max_workers=no_of_cpus)
        pool = mp.Pool(no_of_cpus)
        pool.map_async(compare_two_molecules_single_process, tqdm(get_inputs(non_hydrogen_graphs, include_comparisons_with_itself, list_of_equivalent_molecule_names), total=no_of_neighbourhood_sets, desc='Obtaining neighbouring pairs of molecules', unit='calc'))
        pool.close()
        pool.join()
        
        # Fourth, save the list from a multiprocessing Manager list to a regular list
        list_of_equivalent_molecule_names = list(list_of_equivalent_molecule_names)

    # Fifth, return list_of_equivalent_molecule_names
    return list_of_equivalent_molecule_names

# ----------------------------------------------------------------------------------------------------------------------------------------------------------

def compare_two_molecules_single_process(input_data):
    """
    This method is designed to perform a single 
    """

    # First, unpack the data from input_data
    mol_name1, mol_name2, molecule1_graph, molecule2_graph, list_of_equivalent_molecule_names = input_data

    # Second, create an instance of GraphMatcher
    GM = GraphMatcher(molecule1_graph, molecule2_graph) # Turn mol1 into mol2

    # Third, determine all the ways that molecule 1 could be mapped onto molecule 2 to (potentially) give the same molecule. 
    all_unique_matches = GM.get_all_unique_isomorphic_graphs()

    # Fourth, record all these matches in 
    list_of_equivalent_molecule_names.append(((mol_name1, mol_name2), all_unique_matches))

def get_inputs(non_hydrogen_graphs, include_comparisons_with_itself, list_of_equivalent_molecule_names):
    """
    This generator is designed to return all the input methods required for compare_two_molecules_single_cpu

    Parameters
    ----------
    non_hydrogen_graphs : dict of networkx.Graphs
        This dict contains all the graphs of the molecule, where hydrogens have been removed and added to attached atoms as node features.
    include_comparisons_with_itself : bool.
        This boolean indicates it the user want to give the ways that a molecule could map onto itself. 
            * For determining equivalent molecules, we dont want to do this, so set to False
            * For determining equivalent dimers, we do want to do this, so set this to True.
    list_of_equivalent_molecule_names : list
        This is the list to store all the names that map molecule 1 onto molecule 2. 

    Returns
    -------
    mol_name1 : int
        This is the mol_name of molecule 1 in non_hydrogen_graphs
    mol_name2 : int
        This is the mol_name of molecule 2 in non_hydrogen_graphs
    non_hydrogen_graphs : dict of networkx.Graph
        This is the dict of all the graphs of all molecules in the crystal, where hydrogens have been removed and added to attached atoms as node features.
    list_of_equivalent_molecule_names : list
        This is the list to store all the names that map molecule 1 onto molecule 2. 
    """

    # First, get the names of the molecules in 
    mol_names = sorted(non_hydrogen_graphs.keys())

    # Second, provide the additional index increment if include_comparisons_with_itself == True.
    index_increment = 0 if include_comparisons_with_itself else 1

    # Third, for each molecule in non_hydrogen_graphs.
    for index1 in range(len(mol_names)):

        # Fourth, obtain the name of the first molecule to eexamine.
        mol_name1 = mol_names[index1]

        # Fifth, make a pointer to the non-hydrogen graph of molecule 1.
        molecule1_graph = non_hydrogen_graphs[mol_name1]

        # Sixth, for every other molecule in non_hydrogen_graphs.
        for index2 in range(index1+index_increment, len(mol_names)):

            # Seventh, obtain the name of the first molecule to eexamine.
            mol_name2 = mol_names[index2]

            # Eighth, make a pointer to the non-hydrogen graph of molecule 2.
            molecule2_graph = non_hydrogen_graphs[mol_name2]

            # Ninth, yield input variables for compare_two_molecules_single_process
            yield (mol_name1, mol_name2, molecule1_graph, molecule2_graph, list_of_equivalent_molecule_names)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------





