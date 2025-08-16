from io import StringIO
from ase.io import read
from ase.utils import reader

# Made from NWChem interface

no_of_decimal_places = 8
default_maxcore = 1000
def write_orca_in_FC(fd, atoms, params):
    """Function to write ORCA input file
    """

    charge = sum(atoms.get_initial_charges())
    mult   = sum(atoms.get_initial_magnetic_moments())+1

    orcablocks.append('%SCF MaxIter 1000 END')

    fd.write("! %s \n" % params['orcasimpleinput'])
    for orcablock in params['orcablocks']:
        fd.write("%s \n" % orcablock)

    fd.write('\n*xyz')
    fd.write(" %d" % charge)
    fd.write(" %d \n" % mult)
    for atom in atoms:
        symbol = atom.symbol + '   '
        fd.write(symbol + str(atom.position[0]) + ' ' + str(atom.position[1]) + ' ' + str(atom.position[2]) + '\n')
    fd.write('*\n')
