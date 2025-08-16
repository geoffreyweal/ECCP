'''
Geoffrey Weal, processing_OPV_Dimer_data.py, 9/3/22

This script provides the methods required for processing output.log files for overlap matrix, the coefficients of the MOs, and the energies of the MOs.

'''
import os

# ----------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------

def get_matrix_data(line, max_row, max_col, current_max_col, current_cols, data_heap, not_square=False):
    """
    This method will extract the matrix data from the Gaussian output.log file.

    Parameters
    ----------
    line : str
        This is the line from the output.log file to read.
    max_row : int
        This is the maximum row size of the matrix.
    max_col : int
        This is the maximum column size of the matrix.
    current_max_col : int
        This is the current column that is being read from the Gaussian output.log file. This labels the maximum column index that has been obtained so far.
    current_cols : list of ints
        This is the current columns being read from the Gaussian output.log file.
    data_heap : dict.
        This dictionary contains all the information recorded about the matrix from the Gaussian output.log file. 
    not_square : bool
        This boolean indices if this matrix is not square. 

    Returns
    -------
    max_row : int
        This is the maximum row for this matrix.
    max_col : int
        This is the maximum column for this matrix if it is a square matrix.
    current_max_col : int
        This is the current maximum row that is currently being read from the Gaussian output.log file.
    current_cols : list of ints
        This is the current columns being read from the Gaussian output.log file.
    end_of_matrix_in_file : bool.
        True if you think you still need to read the Gaussian output.log file, False if you have read the full matrix from the Gaussian output.log file.
    """

    # First, if you see lots of spaces, the column numbers are being given
    if line.startswith('             '):
        current_cols = [int(col) for col in line.rstrip().split()]
        if current_cols[-1] > current_max_col:
            current_max_col = current_cols[-1]
        return max_row, max_col, current_max_col, current_cols, True

    # Second, split the line into a list of str.
    line = line.rstrip().split()

    # Third, if this matrix is not square, this will tell us if the matrix has been completed
    if not_square:
        if not line[0].isdigit():
            return max_row, max_col, current_max_col, current_cols, False 

    # Fourth, if the above is not true, we are looking at the values across a row.
    # record the information about value in this line for the matrix
    current_row = int(line.pop(0)) # first entry is the row number
    # The rest of the values are the values in the matrix
    values = []
    for value in line:
        value = ('' if value.startswith('-') else ' ') + value.replace('D','E')
        values.append(value)

    # Fifth, add column values for this row into the data heap 
    for index in range(len(values)):
        current_col = current_cols[index]
        if data_heap is not None:
            value = values[index]
            try:
                float(value)
            except:
                raise Exception('huh?')
            data_heap[(current_row,current_col)] = value

    # Sixth, determine if we are at the end of the matrix.
    if not not_square:
        if current_cols[0] == 1: #If the first column is labelled 1, we have only just started reading the current matrix from the Gaussian output.log file.
            max_row = current_row
            max_col = max_row
        else:
            if (current_row == max_row) and (current_col == max_col):
                if not (current_max_col == max_col):
                    raise Exception('Issue here. current_max_col: '+str(current_max_col)+'. max_col: '+str(max_col))
                return max_row, max_col, current_max_col, current_cols, False

    return max_row, max_col, current_max_col, current_cols, True

# ----------------------------------------------------------------------------------------------------------------------------------

def get_eigenvalue_and_MO_coefficients_data(line, found_eigenvalue, recorded_all_MO_names, recorded_all_electron_positions, current_symbol, current_index, orbital_energies_data_heap, MO_coefficients_data_heap, MO_coefficients_orbital_names_heap, MO_occupancies):
    """
    Get the names of the MO orbtials for MO coefficients, as well as the occupancies of the MOs.

    Parameters
    ----------
    line : str
        This is the line from the output.log file to read.
    found_eigenvalue : bool.
        This bool will indicates if you are currently reading in eigenvalues.
    recorded_all_MO_names : bool.
        This bool indicates if you are reading in the names of the MO coefficients.
    recorded_all_electron_positions : lis
        This bool. indicates if the occupancy of MOs has been recorded.
    current_symbol : str.
        This is the current symbol of the atom being investigated
    current_index : int
        This is the index of the current atom in the molecule.
    orbital_energies_data_list : list
        This list is where to write all the energies of MOs (eigenvalues)
    MO_coefficients_orbital_names_heap : dict.
        This is the heap containing all the infromation about the coefficients of the MOs (eigenvectors)
    MO_occupancies : list
        This is a list of all the MO that are occupied (T) and vacant (F). 

    Returns
    -------
    end_of_matrix_in_file : bool.
        True if you think you still need to read the Gaussian output.log file, False if you have read the full matrix from the Gaussian output.log file.
    found_eigenvalue : bool.
        True if you are looking at eigenvalues, False if not
    recorded_all_MO_names : bool.
        This tag indicates if all the MO orbtial names for the molecule/system have been recorded from the output.log file
    recorded_all_electron_positions : bool.
        This tag indicates if all the MO occupancies have been recorded or not yet. 
    current_symbol : str.
        This is the chemical symbol (element) of the current atom being recorded.
    current_index : int.
        This is the index of the current atom being recorded.
    """

    # First, if you find this, time to end this reading of the output.log file.
    if line.startswith(' Mulliken charges:'):
        return False, found_eigenvalue, recorded_all_MO_names, recorded_all_electron_positions, current_symbol, current_index

    # Second, these tags indicate if you are recording occupancies or if you are done.
    if   recorded_all_electron_positions:
        return True, found_eigenvalue, recorded_all_MO_names, recorded_all_electron_positions, current_symbol, current_index
    elif line.startswith('     Density Matrix:'):
        recorded_all_electron_positions = True
        return True, found_eigenvalue, recorded_all_MO_names, recorded_all_electron_positions, current_symbol, current_index

    # Third, this tag will help to indicate if all the names of the MO orbtials have been recorded.
    if   'Eigenvalues' in line:
        if found_eigenvalue:
            recorded_all_MO_names = True
        found_eigenvalue = True
        for energy in line.strip().replace('Eigenvalues --', '').split():
            orbital_energies_data_heap[len(orbital_energies_data_heap)] = float(energy)
        return True, found_eigenvalue, recorded_all_MO_names, recorded_all_electron_positions, current_symbol, current_index

    # Fourth, if you are looking at the MO coefficients, see which orbtials are vacant and which are occupied.
    if line.startswith('                          '):
        occupancy_information = line.rstrip().split()
        if any([value.isdigit() for value in occupancy_information]):
            return True, found_eigenvalue, recorded_all_MO_names, recorded_all_electron_positions, current_symbol, current_index
        MO_occupancy = [str(value) for value in occupancy_information]
        MO_occupancies += MO_occupancy
        return True, found_eigenvalue, recorded_all_MO_names, recorded_all_electron_positions, current_symbol, current_index

    # Fifth, if the MO names have already been fully recorded, return at this point.
    if recorded_all_MO_names:
        return True, found_eigenvalue, recorded_all_MO_names, recorded_all_electron_positions, current_symbol, current_index

    # Sixth, if you have not obtained the names of the MO's yet, obtain this info from the line.
    MO_name_information = line.rstrip().split()[:-5]
    if len(MO_name_information) == 4:
        current_index  = int(MO_name_information[1])-1
        current_symbol = MO_name_information[2]
        del MO_name_information[2]
        del MO_name_information[1]
    MO_index, MO_name_on_atom = MO_name_information
    MO_index = int(MO_index)-1
    MO_coefficients_orbital_names_heap.setdefault(current_index,[]).append((MO_index, MO_name_on_atom))

    # Seventh, obtain the MO coefficients from the output.log file for this line.
    MO_coefficients_information = line.rstrip().split()[-5:]
    first_orbital_level = len(MO_occupancies)-5+1
    for MO_coefficient_index in range(len(MO_coefficients_information)):
        MO_coefficient = MO_coefficients_information[MO_coefficient_index]
        orbital_level = first_orbital_level + MO_coefficient_index
        basis_function = MO_index + 1
        MO_coefficients_data_heap[(basis_function,orbital_level)] = float(MO_coefficient)

    # Eighth, return quantities.
    return True, found_eigenvalue, recorded_all_MO_names, recorded_all_electron_positions, current_symbol, current_index

# ----------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------

def process_MO_data_from_fort7_file(filename):
    """
    This method is designed to extract the MO energies and coefficients (eigenvalues and eigenvectors) from the fort.7 file. 

    This values are to a greater significant figure than found in the output.log file.

    Parameters
    ----------
    filename : str
        This is the path to the fort.7 file.

    Returns
    -------
    orbital_energies_data_heap : dict.
        This dictionary contains all the information about the MO energies.
    MO_coefficients_data_heap : dict.
        This dictionary contains all the information about the MO coefficients.
    """

    # First, initalise the data heaps
    orbital_energies_data_heap = {}
    MO_coefficients_data_heap = {}

    # Second, read data from the fort.7 file.
    with open(filename,'r') as outputLOG:
        outputLOG.readline()
        orbital_level = None
        for line in outputLOG:
            line.rstrip()
            
            # 2.1: Get orbital level and orbital energy
            if 'MO' in line:
                line = line.split()
                orbital_level = int(line[0])
                orbital_energy = line[-1].replace('OE=','').replace('D','E')
                orbital_energies_data_heap[orbital_level] = orbital_energy
                basis_function = 1
                continue
            
            # 2.2: Get MO coefficients
            line = line.replace('D+','E+')
            line = line.replace('D-','EN')
            line = line.replace('-',' -')
            line = line.replace('N','-')
            
            # 2.3: Get the MO coefficients for the line in the state it needs to be in.
            MO_coefficients = []
            for value in line.split():
                value = ('' if value.startswith('-') else ' ') + value.lstrip()
                MO_coefficients.append(value)
            
            # 2.4: Add the MO coefficients to the heap.
            for MO_coefficient in MO_coefficients:
                # This heap is ordered by (row_position, column_position)
                MO_coefficients_data_heap[(basis_function,orbital_level)] = MO_coefficient
                basis_function += 1

    # Third, return data heaps
    return orbital_energies_data_heap, MO_coefficients_data_heap

def remove_fort7_file(path_to_log_file):
    """
    This method is designed to remove the fort.7 file from the folder once it has been used. 

    This is because he have extracted the data from this file, and this file can be quite large.
    
    Parameters
    ----------
    path_to_log_file : str
        This is the path to the fort.7 file.
    """
    if os.path.exists(path_to_log_file+'/fort.7'):
        os.remove(path_to_log_file+'/fort.7')

# ----------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------

def remove_eigenfile_data_from_outputLOG_file(filename, lines_to_remove):
    """
    This method is designed to remove the lines from the output.log file. 

    These lines will include the matrix data that is large and takes up a lot of space. 

    This matrix data will have been extracted into other files from previous methods. 

    Parameters
    ----------
    filename : str
        This is the path to the output.log file.
    lines_to_remove : list
        These are the start and ends of sections of the output.log file to remove.
    """

    # First, write a new file called filename+'.new', which does not contain the components of the file given by lines_to_remove
    with open(filename+'.new','w') as outputLOGnew:
        with open(filename,'r') as outputLOG:
            line_counter = 0
            remove_index = 0
            to_read = True
            for line in outputLOG:
                line_counter += 1
                if   remove_index == len(lines_to_remove):
                    pass
                elif line_counter == lines_to_remove[remove_index]:
                    remove_index += 1
                    to_read = not to_read
                if to_read:
                    outputLOGnew.write(line)

    # Second, remove the original file and replace it with filename+'.new'
    os.remove(filename)
    os.rename(filename+'.new',filename)

# ----------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------

def write_orbital_names(MO_coefficients_orbital_names_heap, filename):
    """
    This method is designed to write the names of the orbitals of atoms in the molcule/system to disk.

    Parameters
    ----------
    MO_coefficients_orbital_names_heap : dict
        This is the data heap that contains all the information about the MO orbital names. 
    filename : str
        This is the path to save this data to.
    """
    with open(filename, 'w') as filenameTXT:
        for atom_index, orbital_information in sorted(MO_coefficients_orbital_names_heap.items(), key=lambda x:x[0]):
            filenameTXT.write('Atom Index: '+str(atom_index)+'\n')
            MO_indices = []
            MO_names = []
            for MO_index, MO_name in sorted(orbital_information, key=lambda x:x[0]):
                MO_indices.append(str(MO_index))
                MO_names.append(str(MO_name))
            filenameTXT.write(' '.join(MO_indices)+'\n')
            filenameTXT.write(' '.join(MO_names)+'\n')

def write_MO_occupancies(MO_occupancies, filename):
    """
    This method is designed to write the MO occupancies.

    Parameters
    ----------
    MO_occupancies : dict
        This is the data heap that contains all the information about the MO occupancies. 
    filename : str
        This is the path to save this data to.
    """
    with open(filename, 'w') as filenameTXT:
        filenameTXT.write(' '.join([str(value) for value in MO_occupancies]))

# ----------------------------------------------------------------------------------------------------------------------------------

def write_1D_matrix(data_heap, filename):
    """
    This method is designed to save a 1D matrix as simple as possible to a txt file.

    Parameters
    ----------
    data_heap : dict
        This is the data heap that contains all the information about the 1D matrix
    filename : str
        This is the path to save this data to.
    """
    # First, get the max row for this matrix.
    max_row = 0
    for row in data_heap.keys():
        if row > max_row:
            max_row = row
    
    # Second, save the data to disk.
    with open(filename, 'w') as filenameTXT:
        for row in range(1,max_row+1):
            value = data_heap[row]
            filenameTXT.write(str(value))
            filenameTXT.write('\n')

def write_2D_matrix(data_heap, filename, symmetric_matrix=False):
    """
    This method is designed to save a 2D matrix as simple as possible to a txt file.

    Parameters
    ----------
    data_heap : dict
        This is the data heap that contains all the information about the 1D matrix.
    filename : str
        This is the path to save this data to.
    symmetric_matrix : bool.
        True if this is a symmetric matrix, False if not.
    """
    # First, get the max row and max column for this matrix.
    max_row = 0
    max_col = 0
    for (row, col) in data_heap.keys():
        if row > max_row:
            max_row = row
        if col > max_col:
            max_col = col
    
    # Second, save the data to disk.
    filenameTXT = open(filename, 'w')
    if symmetric_matrix:
        write_symmetric_2D_matrix(max_row, max_col, data_heap, filenameTXT)
    else:
        write_non_symmetric_2D_matrix(max_row, max_col, data_heap, filenameTXT)
    filenameTXT.close()

def write_symmetric_2D_matrix(max_row, max_col, data_heap, filenameTXT):
    """
    This method will write a symmetric 2D matrix.

    Parameters
    ----------
    max_row : int
        This is the maximum row in this matrix
    max_col : int
        This is the maximum column in this matrix
    data_heap : dict
        This is the data heap that contains all the information about the 1D matrix.
    filenameTXT : TXTFILE
        This is the text file to write data to.
    """
    for row in range(1,max_row+1):
        for col in range(1,row+1):
            index1 = row if (row >= col) else col
            index2 = col if (row >= col) else row
            value = data_heap[(index1,index2)]
            filenameTXT.write(str(value))
            if col < max_col:
                filenameTXT.write('\t')
        filenameTXT.write('\n')

def write_non_symmetric_2D_matrix(max_row, max_col, data_heap, filenameTXT):
    """
    This method will write a non-symmetric 2D matrix.

    Parameters
    ----------
    max_row : int
        This is the maximum row in this matrix
    max_col : int
        This is the maximum column in this matrix
    data_heap : dict
        This is the data heap that contains all the information about the 1D matrix.
    filenameTXT : TXTFILE
        This is the text file to write data to.
    """
    for row in range(1,max_row+1):
        for col in range(1,max_col+1):
            index1 = row
            index2 = col
            value = data_heap[(index1,index2)]
            filenameTXT.write(str(value))
            if col < max_col:
                filenameTXT.write('\t')
        filenameTXT.write('\n')

# ----------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------











