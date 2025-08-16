'''
Geoffrey Weal, functional_and_basis_set_storage.py, 13/4/22

This script is designed to hold data for functionals and basis sets.

'''

class functional_and_basis_set_storage:
    """
    This class is designed to hold the data with reference to functionals and basis sets.
    """
    def __init__(self):
        self.data = {}

    def add(self, dft_details, crystal_name, dimer_name, electronic_coupling_data):
        """
        This method is designed to input data into the storage.

        Parameters
        ----------
        dft_details : str.
            This is the dft details for this calculation, including the functional and basis set name
        crystal_name : str.
            This is the name of the crystal file
        dimer_name : str.
            This is the name of the dimer.
        electronic_coupling_data : tuple
            This contains all the coupling data from the output.log file.
        """
        if dft_details not in self.data:
            self.data[dft_details] = {}
        self.data[dft_details].setdefault(crystal_name,[]).append((dimer_name, electronic_coupling_data))

    def keys(self):
        """
        This method will give the list of functionals and basis sets in this storage unit
        
        Returns
        -------
        The keys of the self.data dictionary.
        """
        return sorted(list(self.data.keys()))

    def items(self):
        """
        This method will give the list of functionals and basis sets in this storage unit, as well as the results each.

        Returns
        -------
        The items of the self.data dictionary. This includes keys and values.
        """
        return sorted(list(self.data.items()), key=lambda x:x[0])

    def get(self, dft_details):
        """
        This method will return the data for a specific functional and basis set

        Parameters
        ----------
        dft_details : str.
            This is the dft details for this calculation, including the functional and basis set name

        Returns
        -------
        self.data[dft_details] : dict.
            This is the data for all crystals for this functional and basis set
        """
        return self.data[dft_details]

    def __repr__(self):
        """
        This is the way that this class will be represented in python

        Returns
        -------
        to_string : str.
            This cleanly shows what information is stored in this object.
        """
        to_string  = ''
        to_string += '--------------------------------\n'
        to_string += 'Data in Storage\n'
        to_string += '--------------------------------\n'
        for dft_details, all_data in self.data.items():
            to_string += 'Functional+Basis Set: '+str(dft_details)+'\n'
            to_string += 'Molecule name: Dimer names\n'
            for crystal_name, dimer_data in all_data.items():
                to_string += str(crystal_name)+': '+str([dimer_name for dimer_name, data in dimer_data])+'\n'
            to_string += '--------------------------------\n'
        return to_string

    def no_of_crystals(self, dft_details):
        """
        This is the amount of crystals in the database for a functional.

        Parameters
        ----------
        dft_details : str.
            This is the dft details for this calculation, including the functional and basis set name

        Returns
        -------
        len(self.data[dft_details]) : dict.
            This is the number of crystals in this database.
        """
        return len(self.data[dft_details])











