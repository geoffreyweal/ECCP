"""
get_coulomb_energy.py, Geoffrey Weal, 20/3/22

get_coulomb_energy contains a method designed to obtain the coulomb energy between two molecules in eV
"""
import numpy as np

from SUMELF import get_distance

coulomb_constant = 14.3996454784878811182 # units of eV·Å·e−2
def get_coulomb_energy(positions_1, positions_2, charges_1, charges_2, displacement_vector, relative_permittivity=1.0):
    """
    get_coulomb_energy contains a method designed to obtain the coulomb energy between two molecules in eV.

    Parameters
    ----------
    molecule1 : ase.Atoms
        This is the atoms object for molecule 1.
    molecule2 : ase.Atoms
        This is the atoms object for molecule 2.
    displacement_vector : numpy.array
        This is the vector that describes the position of the unit cell that molecule 2 is in compared to the ijk=000 centre unit cell.

    Returns
    -------
    energy : float
        This is the coulomb energy between molecule 1 and molecule 2 (in a unit cell described by displacement_vector). This value is given in eV.
    """

    # First, move molecule 2 into the ijk unit cell described by the displacement_vector vector.
    positions_2_with_disp_vec = positions_2 + displacement_vector

    # Second, get the Coulomb energy between molecule 1 and molecule 2 (in a unit cell described by displacement_vector)
    # The total_coulomb_value is given in e^2/Å
    total_coulomb_value = 0.0 
    for charge_1, position_1 in zip(charges_1,positions_1):
        for charge_2, position_2 in zip(charges_2,positions_2_with_disp_vec):
            distance = get_distance(position_1,position_2)
            coulomb_value = (charge_1*charge_2)/distance
            total_coulomb_value += coulomb_value

    # Third, multiply total_coulomb_value by the coulomb constant, which is given in eV·Å·e−2. This will give your energy in eV
    energy = (coulomb_constant/relative_permittivity) * total_coulomb_value

    # Fourth, return the ATC-based energy value between molecule1 and molecule2 with a difference in cell of displacement_vector
    return energy
