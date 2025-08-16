"""
are_environments_equivalent.py, Geoffrey Weal, 3/3/22

This script is designed to determine if the environment about two dimers are structurally equivalent or not.
"""
from ase import Atoms
from SUMELF import get_distance

from ECCP.ECCP.get_unique_molecules_methods.set_of_invariance_methods.are_environments_equivalent import are_environments_equivalent as are_environments_equivalent_for_molecules

def are_environments_equivalent(rotation_reflection_matrix, dimer1_com, dimer2_com, dimer_indices_being_compared, neighbouring_molecules_about_dimers, non_hydrogen_molecules, no_of_cpus=1):
    """
    This method is designed to determine if the environment about two molecules or dimers are structurally equivalent or not.

    Parameters
    ----------
    rotation_reflection_matrix : 3x3 np.array :
        This matrix indicates how to reflect and rotate dimer 2 so that it fit ontop of dimer 1 (after molecules have been centred to their centre of masses).
    dimer1_com : 3x1 np.array 
        This is the centre of mass of dimer 1.
    dimer2_com : 3x1 np.array 
        This is the centre of mass of dimer 2.

    dimer_indices_being_compared : list
        These are the indices/names of the dimers being compared.
    neighbouring_molecules_about_dimers : dict.
        This is the information about the molecules that surround (in the vicinity of) each dimer in the crystal.
    non_hydrogen_molecules : list of ase.Atoms
        This is the list of molecules in the crystal, not including hydrogens.

    Returns
    -------
    True if the two molecules/dimers environments are structurally equivalent, False if they are not
    """

    return are_environments_equivalent_for_molecules(rotation_reflection_matrix, dimer1_com, dimer2_com, dimer_indices_being_compared, neighbouring_molecules_about_dimers, non_hydrogen_molecules, no_of_cpus=no_of_cpus)
