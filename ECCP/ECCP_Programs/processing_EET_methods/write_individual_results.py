'''
Geoffrey Weal, write_individual_results.py, 10/6/22

This program is designed to write txt files that contains the coupling energies for the dimers for a crystal.
'''

def write_individual_results(organised_data, individual_eet_data_foldername):
    """
    This method will write txt files that contains the coupling energies for the dimers for a crystal.

    Parameters
    ----------
    organised_data : dict.
        This contains all the EET data to save to the excel worksheet and to disk.
    individual_eet_data_foldername : str.
        This is the folder to save this data in to.
    """
    filenames_created = []
    # For each crystal.
    for crystal_name in organised_data.keys():
        dimer_names = sorted(organised_data[crystal_name].keys(), key=lambda x: int(x.split('_')[0].replace('Dimer','')))
        # For each dimer.
        for dimer_name in dimer_names:
            # For each functional and basis set.
            for dft_details in organised_data[crystal_name][dimer_name].keys():
                # Get the filepath to save data to.
                filename = str(crystal_name)+'_'+str(dft_details)+".txt"
                path_to_file = individual_eet_data_foldername+'/'+filename
                # Initialise the file if it has not been made yet
                if not (filename in filenames_created):
                    with open(path_to_file, "w") as EET_file:
                        EET_file.write('Dimer Name\tTotal coupling (meV)\n')
                    filenames_created.append(filename)
                # Save coupling energy to file, in meV. 
                coupling_energy = organised_data[crystal_name][dimer_name][dft_details][1][-1] * 1000.0 # Data is given in eV originally, we want to give it in meV
                with open(path_to_file, "a") as EET_file:
                    EET_file.write(str(dimer_name)+'\t'+str(coupling_energy)+'\n')
