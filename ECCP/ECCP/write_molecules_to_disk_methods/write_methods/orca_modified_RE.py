from io import StringIO
from ase.io import read
from ase.utils import reader

# Made from NWChem interface

default_maxcore = 1000
def write_orca_in_RE(fd, atoms, molecule_xyz_name, perform_opt=True, perform_CalcAll=False, perform_TD=False, perform_freq=False, perform_raman=False, perform_pop=False, **orca_parameters):
    """Function to write ORCA input file
    """

    orcasimpleinput = orca_parameters['method']+' '+orca_parameters['basis']
    orcablocks      = [('%PAL NPROCS '+str(orca_parameters['NPROCS'])+' END')]

    if perform_opt:
        orcasimpleinput += ' OPT'

    maxcore = orca_parameters['maxcore'] if ('maxcore' in orca_parameters) else default_maxcore
    orcablocks.append('%maxcore '+str(maxcore))

    if perform_CalcAll:
        orcablocks.append('%GEOM Calc_Hess true Recalc_Hess 1 END')

    if perform_TD:
        orcablocks.append('%TDDFT '+orca_parameters['td_settings']+' END')

    if perform_raman:
        orcasimpleinput += ' NUMFREQ'
        orcablocks.append('%ELPROP POLAR 1 END')
    elif perform_freq:
        orcasimpleinput += ' FREQ'

    orcablocks.append('%SCF MaxIter 1000 END')
    
    params = {'orcasimpleinput': orcasimpleinput, 'orcablocks': orcablocks}

    charge = sum(atoms.get_initial_charges()) # orca_parameters['charge'] #
    mult   = sum(atoms.get_initial_magnetic_moments())+1 # orca_parameters['mult']

    fd.write("! %s \n" % params['orcasimpleinput'])
    for orcablock in params['orcablocks']:
        fd.write("%s \n" % orcablock)

    fd.write('\n*xyzfile')
    fd.write(" %d" % charge)
    fd.write(" %d" % mult)
    #for atom in atoms:
    #    symbol = atom.symbol + '   '
    #    fd.write(symbol + str(atom.position[0]) + ' ' + str(atom.position[1]) + ' ' + str(atom.position[2]) + '\n')
    fd.write(f' {molecule_xyz_name}')
    fd.write('*\n')

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


