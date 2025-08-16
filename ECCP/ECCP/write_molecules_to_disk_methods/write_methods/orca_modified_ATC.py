from io import StringIO
from ase.io import read #, write
from ase.utils import reader

# Made from NWChem interface

default_maxcore = 1000
def write_orca_in_ATC(fd, atoms, environment_about_molecule, molecule_name, **orca_parameters):
    """Function to write ORCA input file
    """

    # First, create the orca_params that holds the settings for ORCA.
    if environment_about_molecule is None: 
        orcasimpleinput = orca_parameters['method']+' '+orca_parameters['basis']
        orcablocks      = [('%TDDFT '+orca_parameters['td_settings']+' END'), ('%PAL NPROCS '+str(orca_parameters['NPROCS'])+' END')]
    else:
        orcasimpleinput = 'QM/QM2 '+orca_parameters['method']+' '+orca_parameters['basis']
        orcablocks      = [('%QMMM'), ('\tQM2CUSTOMMETHOD "'+orca_parameters['method']+' '+orca_parameters['basis']+'"'), ('\tQMATOMS {0:'+str(len(atoms)-1)+'} END'), ('END'), ('%TDDFT '+orca_parameters['td_settings']+' END'), ('%PAL NPROCS '+str(orca_parameters['NPROCS'])+' END')]

    maxcore = orca_parameters['maxcore'] if ('maxcore' in orca_parameters) else default_maxcore
    orcablocks.append('%maxcore '+str(maxcore))

    orcablocks.append('%SCF MaxIter 1000 END')

    # Second, create the params dictionary.
    params = {'orcasimpleinput': orcasimpleinput, 'orcablocks': orcablocks}

    # Third, create the system that combines the molecule and the environment if given
    if environment_about_molecule is not None:
        total_system = atoms.copy() + environment_about_molecule.copy()
    else:
        total_system = atoms.copy()

    # Fouth, obtain the total charge and multiplicity for this system
    charge = sum(total_system.get_initial_charges())
    mult   = sum(total_system.get_initial_magnetic_moments()) + 1

    # Fifth, write the header of the input file
    fd.write("! %s \n" % params['orcasimpleinput'])
    for orcablock in params['orcablocks']:
        fd.write("%s \n" % orcablock)

    # Sixth, write the atoms and their positions into the input file.
    fd.write('\n*xyzfile')
    fd.write(" %d" % charge)
    fd.write(" %d" % mult)
    #for atom in total_system:
    #    symbol = atom.symbol + '   '
    #    #fd.write(symbol + str(atom.position[0]) + ' ' + str(atom.position[1]) + ' ' + str(atom.position[2]) + '\n')
    #    fd.write('{:<10s}{:20.10f}{:20.10f}{:20.10f}{}'.format(symbol, *atom.position)+'\n')
    fd.write(f" {molecule_name}.xyz")
    fd.write('*\n')
