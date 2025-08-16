"""
are_environments_equivalent.py, Geoffrey Weal, 3/3/22

This script is designed to determine if the environment about two molecules or dimers are structurally equivalent or not.
"""
from ase import Atoms
from SUMELF import get_distance

def are_environments_equivalent(rotation_reflection_matrix, molecule1_com, molecule2_com, names_being_compared, neighbouring_molecules_about_molecules, non_hydrogen_molecules):
    """
    This method is designed to determine if the environment about two molecules or dimers are structurally equivalent or not.

    Parameters
    ----------
    rotation_reflection_matrix : 3x3 np.array :
        This matrix indicates how to reflect and rotate molecule 2 so that it fit ontop of molecule 1 (after molecules have been centred to their centre of masses).
    molecule1_com : 3x1 np.array 
        This is the centre of mass of molecule 1.
    molecule2_com : 3x1 np.array 
        This is the centre of mass of molecule 2.

    names_being_compared : list
        These are the names of the molecules being compared.
    neighbouring_molecules_about_molecules : dict.
        This is the information about the molecules that surround (in the vicinity of) each moleule in the crystal.
    non_hydrogen_molecules : list of ase.Atoms
        This is the list of molecules in the crystal, not including hydrogens.

    Returns
    -------
    True if the two molecules/dimers environments are structurally equivalent, False if they are not
    """

    # Preamble: If neighbouring_molecules_about_molecules == {}, we are not worrying about the environment, so return True
    if len(neighbouring_molecules_about_molecules) == 0:
        return True

    # First, extract the names of structurally equivalent molecules that we want to assess.
    mol_name1, mol_name2 = names_being_compared

    # Second, get the environment around molecule 1, having also been moved to molecule 1's centre of mass
    environment_about_molecule1 = get_environment(mol_name1, neighbouring_molecules_about_molecules, non_hydrogen_molecules)
    environment_about_molecule1.set_positions(environment_about_molecule1.get_positions() - molecule1_com)

    # Third, get the environment around molecule 2, having also been moved to molecule 2's centre of mass and rotated and reflected so that molecule 2 lies on top of molecule 1
    environment_about_molecule2 = get_environment(mol_name2, neighbouring_molecules_about_molecules, non_hydrogen_molecules)
    environment_about_molecule2.set_positions(environment_about_molecule2.get_positions() - molecule2_com)
    environment_about_molecule2.set_positions(environment_about_molecule2.get_positions() @ rotation_reflection_matrix.T)

    # Fourth, obtain the chemical elements of atoms in each environment. 
    symbols_of_environment_about_molecule1 = environment_about_molecule1.get_chemical_symbols()
    symbols_of_environment_about_molecule2 = environment_about_molecule2.get_chemical_symbols()

    # Fifth, if the elements do not contain the same elements, the two environments are different
    if not sorted(symbols_of_environment_about_molecule1) == sorted(symbols_of_environment_about_molecule2):
        return False

    # Sixth, get the positions of all the atoms in each molecule. 
    positions_of_environment_about_molecule1 = environment_about_molecule1.get_positions()
    positions_of_environment_about_molecule2 = environment_about_molecule2.get_positions()

    #if not isinstance(mol_name1,int):
    #    import pdb; pdb.set_trace()

    # Seventh, obtain the names of atoms in environment 2 that are closest to atoms in environment 1
    likely_assignment_of_environs = {}
    #all_shortest_distances_between_atoms_in_environ = {}
    for atom_index_em1, element_em1, position_em1 in zip(range(len(symbols_of_environment_about_molecule1)), symbols_of_environment_about_molecule1, positions_of_environment_about_molecule1):
        shortest_index_em2 = -1
        shortest_distance  = float('inf')
        for atom_index_em2, element_em2, position_em2 in zip(range(len(symbols_of_environment_about_molecule2)), symbols_of_environment_about_molecule2, positions_of_environment_about_molecule2):
            if not element_em1 == element_em2:
                continue
            distance = get_distance(position_em1, position_em2)
            if distance < shortest_distance:
                shortest_index_em2 = atom_index_em2
                shortest_distance  = distance
        if shortest_index_em2 == -1:
            raise Exception('huh? could not assign an atom in environment 2 to atom index '+str(atom_index_em1)+' in environment 1.')
        if not (shortest_distance < 0.1):
            return False
        likely_assignment_of_environs[atom_index_em1] = atom_index_em2
        #all_shortest_distances_between_atoms_in_environ[(atom_index_em1, atom_index_em2)] = shortest_distance

     # Eighth, if any two names in likely_assignment_of_environs values are the same, an atom in environment 2 has been assigned to multiple atoms in environment 1, so return False
    if not sorted(likely_assignment_of_environs.values()) == sorted(set(likely_assignment_of_environs.values())):
        return False

    # Ninth, if all distances between identically assigned atoms in each environment have been assigned are small enough, then they are the same.
    #if not all([(shortest_distance < 0.1) for shortest_distance in all_shortest_distances_between_atoms_in_environ.values()]):
    #    return False

    # Tenth, if here, the environments are the same, return True
    return True

# =====================================================================================================================================================

def get_environment(mol_name, neighbouring_molecules_about_molecules, molecules):
    """
    This method is designed to obtain the environment surrounding molecules[mol_name] in the origin unit cell as an ASE Atoms object.

    Parameters
    ----------
    mol_name : int
        This is the molecule that we want to obtain the environment surrounding it (for molecules[mol_name] in the origin unit cell).
    neighbouring_molecules_about_molecules : dict.
        This is the information about the molecules that surround (in the vicinity of) each moleule in the crystal.
    molecules : list of ase.Atoms
        This is the list of molecules in the crystal, which include hydrogens attached to them.

    Returns
    -------
    environment_about_molecule : ase.Atoms
        This is the environment surrounding molecules[mol_name] in the origin unit cell.
    """

    # First, initialise the environment_about_molecule ase.Atoms object
    environment_about_molecule = Atoms()

    # Second, obtain the instance that describes all the neighbouring molecules around the molecule of interest.
    neighbouring_molecules_about_molecule_temp = neighbouring_molecules_about_molecules[mol_name]

    # Third, if neighbouring_molecules_about_molecule_temp is a list, modify it into a dict. in the correct format. 
    if isinstance(neighbouring_molecules_about_molecule_temp,list):
        neighbouring_molecules_about_molecule = {}
        for input_from in neighbouring_molecules_about_molecule_temp:
            neighbour_mol_name                  = input_from[0]
            neighbour_unit_cell_displacement = input_from[1]
            neighbour_displacement           = input_from[2]
            key = (neighbour_mol_name, neighbour_unit_cell_displacement)
            if key in neighbouring_molecules_about_molecule:
                raise Exception('Huh?')
            neighbouring_molecules_about_molecule[key] = neighbour_displacement
    else:
        neighbouring_molecules_about_molecule = neighbouring_molecules_about_molecule_temp

    # Fourth, add molecules that surround the molecule of interest to environment_about_molecule
    for (neighbour_mol_name, neighbour_unit_cell_displacement), neighbour_displacement in sorted(neighbouring_molecules_about_molecule.items()):

        # 4.1: Obtain the surrounding molecule, and translate it based on neighbour_displacement
        neighbour_molecule = molecules[neighbour_mol_name].copy()
        neighbour_molecule.set_positions(neighbour_molecule.get_positions() + neighbour_displacement)

        # 4.2: Add the neighbouring molecule to environment_about_molecule
        environment_about_molecule += neighbour_molecule

    # Fifth, return environment_about_molecule
    return environment_about_molecule

# =====================================================================================================================================================

