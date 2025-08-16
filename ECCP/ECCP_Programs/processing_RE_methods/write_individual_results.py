'''
Geoffrey Weal, write_individual_results.py, 10/6/22

This program is designed to write txt files that contains the coupling energies for the dimers for a crystal.
'''

from ECCP.ECCP_Programs.processing_RE_methods.processing_RE_data_methods import get_reorganisation_energy, get_energy_diff, convert_hartree_to_eV, eV_to_meV

def write_individual_results(organised_data, individual_re_data_foldername):
    """
    This method will write txt files that contains the coupling energies for the dimers for a crystal.

    Parameters
    ----------
    organised_data : dict.
        This contains all the EET data to save to the excel worksheet and to disk.
    individual_re_data_foldername : str.
        This is the folder to save this data in to.
    """

    # energies
    filenames_created = []
    for crystal_name in organised_data.keys():
        molecule_names = sorted(organised_data[crystal_name].keys(), key=lambda x: int(x.split('_')[1].replace('S','')))
        # For each dimer.
        for molecule_name in molecule_names:
            # For each functional and basis set.
            for dft_details in organised_data[crystal_name][molecule_name].keys():
                # Get the filepath to save data to.
                filename = str(crystal_name)+'_'+str(dft_details)+"_energies.txt"
                path_to_file = individual_re_data_foldername+'/'+filename
                # Initialise the file if it has not been made yet
                if not (filename in filenames_created):
                    with open(path_to_file, "w") as RE_file:
                        RE_file.write('Molecule Name\tE_GS(GS) (Ha)\tE_GS(ES) (Ha)\tE_ES(GS) (Ha)\tE_ES(ES) (Ha)\tBand Gap (E_ES(ES) - E_GS(GS)) (eV)\n')
                    filenames_created.append(filename)
                root, eGS_gGS_energy, eES_gGS_energy, eGS_gES_energy, eES_gES_energy, negative_eGS_gGS_freqs, negative_eES_gES_freqs = organised_data[crystal_name][molecule_name][dft_details]
                # Get the band gap
                band_gap = convert_hartree_to_eV(get_energy_diff(eES_gES_energy, eGS_gGS_energy))
                with open(path_to_file, "a") as RE_file:
                    RE_file.write(str(molecule_name)+'\t'+str(eGS_gGS_energy)+'\t'+str(eGS_gES_energy)+'\t'+str(eES_gGS_energy)+'\t'+str(eES_gES_energy)+'\t'+str(band_gap)+'\n')

    # Reorganisation Energies
    filenames_created = []
    for crystal_name in organised_data.keys():
        molecule_names = sorted(organised_data[crystal_name].keys(), key=lambda x: int(x.split('_')[1].replace('S','')))
        # For each dimer.
        for molecule_name in molecule_names:
            # For each functional and basis set.
            for dft_details1 in organised_data[crystal_name][molecule_name].keys():

                # First, get the excited structure geometry state of mol1 (where exciton is coming from) (donor)
                root1, donor_eGS_gGS_energy, donor_eES_gGS_energy, donor_eGS_gES_energy, donor_eES_gES_energy, negative_eGS_gGS_freqs, negative_eES_gES_freqs = organised_data[crystal_name][molecule_name][dft_details1]

                for molecule_name2 in molecule_names:
                    # For each functional and basis set.
                    for dft_details2 in organised_data[crystal_name][molecule_name].keys():

                        if not (dft_details1 == dft_details2):
                            continue

                        # Second, get the ground structure geometry state of mol2 (where exciton is going to) (acceptor)
                        root2, acceptor_eGS_gGS_energy, acceptor_eES_gGS_energy, acceptor_eGS_gES_energy, acceptor_eES_gES_energy, negative_eGS_gGS_freqs, negative_eES_gES_freqs = organised_data[crystal_name][molecule_name2][dft_details2]

                        # Get the filepath to save data to.
                        filename = str(crystal_name)+'_'+str(dft_details1)+"_reorg_energies.txt"
                        path_to_file = individual_re_data_foldername+'/'+filename
                        # Initialise the file if it has not been made yet
                        if not (filename in filenames_created):
                            with open(path_to_file, "w") as RE_file:
                                RE_file.write('Donor Molecule\tAcceptor Molecule\tReorganisation Energy (meV)\n')
                            filenames_created.append(filename)
                        # Save coupling energy to file, in meV, 1 -> 2
                        # Get the reorganisation energy in meV
                        reorganisation_energy = convert_hartree_to_eV(get_reorganisation_energy(donor_eGS_gGS_energy, acceptor_eES_gGS_energy, donor_eGS_gES_energy, acceptor_eES_gES_energy)) * eV_to_meV
                        with open(path_to_file, "a") as RE_file:
                            RE_file.write(str(molecule_name)+'\t'+str(molecule_name2)+'\t'+str(reorganisation_energy)+'\n')
