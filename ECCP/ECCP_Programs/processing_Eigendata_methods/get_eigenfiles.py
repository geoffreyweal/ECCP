'''
Geoffrey Weal, processing_OPV_Dimer_data.py, 9/3/22

This script provides the methods required for processing output.log files for overlap matrix, the coefficients of the MOs, and the energies of the MOs.

'''
import os
from datetime import datetime, timedelta

from ECCP.ECCP_Programs.shared_general_methods.get_eigenfiles_methods import get_matrix_data, get_eigenvalue_and_MO_coefficients_data, remove_eigenfile_data_from_outputLOG_file
from ECCP.ECCP_Programs.shared_general_methods.get_eigenfiles_methods import process_MO_data_from_fort7_file, remove_fort7_file
from ECCP.ECCP_Programs.shared_general_methods.get_eigenfiles_methods import write_1D_matrix, write_2D_matrix, write_orbital_names, write_MO_occupancies

def get_eigenfiles(path_to_log_file, log_filename, remove_eigendata_from_outputLOG_file=False, get_MO_data_from_fort7_file=True):
    """
    This method is designed to extract eigen-information from the output.log and fort.7 file. This includes orbital overlaps, MO energies and MO coefficients.

    Parameters
    ----------
    path_to_log_file : str
        This is the path to the log file.
    log_filename : str
        This is the name of the log file (likely called output.log).
    remove_eigendata_from_outputLOG_file : bool.
        If true, remove eigendata, such as orbital overlaps, MO energies and MO coefficients, from the output file once you have created text files of this. This is useful to turn a file that is GBs in size into KBs.
    """

    # First, if the following text files have already been created, we probably dont need to do this again, especially since the output.log file may not contain this information anymore.
    main_2D_matrices = ['orbital_overlap_matrix.txt', 'MO_orbital_names.txt', 'MO_occupancies.txt', 'MO_energies.txt', 'MO_coefficients.txt'] 
    if check_files_have_been_made(path_to_log_file, main_2D_matrices):
        print('Will not obtain eigenfiles, have already obtained all the necessary eigendata in text files.')
        return

    # Second, obtain the matrix data from the output.log file and save it to txt files
    lines_to_remove = get_eigenfiles_from_outputLOG_file(path_to_log_file, log_filename, get_MO_data_from_fort7_file=get_MO_data_from_fort7_file)

    # Third, save matrix files from fort.7 file (if desired). 
    if get_MO_data_from_fort7_file:
        orbital_energies_data_heap, MO_coefficients_data_heap = process_MO_data_from_fort7_file(path_to_log_file+'/fort.7')
        if not orbital_energies_data_heap == {}:
            write_1D_matrix(orbital_energies_data_heap, filename=path_to_log_file+'/'+'MO_energies.txt')
            print('Made MO_energies.txt file')
        if not MO_coefficients_data_heap == {}:
            write_2D_matrix(MO_coefficients_data_heap, filename=path_to_log_file+'/'+'MO_coefficients.txt', symmetric_matrix=False)
            print('Made MO_coefficients.txt file')
    
    # Fourth, remove the matrices from the output.log file, and remove the fort.7 file.
    if remove_eigendata_from_outputLOG_file and check_files_have_been_made(path_to_log_file, main_2D_matrices):
        message = 'matrix data from output.log '+('and fort.7' if get_MO_data_from_fort7_file else '')+' file'+('s' if get_MO_data_from_fort7_file else '')+'.'
        print('Removing '+str(message))
        remove_eigenfile_data_from_outputLOG_file(path_to_log_file+'/'+log_filename, lines_to_remove)
        if get_MO_data_from_fort7_file:
            remove_fort7_file(path_to_log_file)
        print('Removed '+str(message))

def check_files_have_been_made(path_to_log_file, main_2D_matrices):
    """
    This method will check that all the txt files made by this program exist before removals occur. 
    """
    for check_file in main_2D_matrices:
        if not (check_file in os.listdir(path_to_log_file)):
            return False
    return True


# ----------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------

def get_eigenfiles_from_outputLOG_file(path_to_log_file, filename, get_MO_data_from_fort7_file=True):
    """
    This method is designed to obtain the the overlap matrix, the coefficients of the MOs, and the energies of the MOs for a single molecule or dimer in a non-EET calculation.

    Parameters
    ----------
    path_to_log_file : str.
        This is the path to the output.log file.
    filename : str.
        This is the filename of the output.log file.
    get_MO_data_from_fort7_file : bool.
        This tag indicates if MO energies and coefficients will be obtained from the fort.7 file or not. 

    Returns
    -------
    True if any of these dicts or lists have data in them, else return False
    """

    # First, create all the dictionaries for throwing data from the output file into
    printed_data = ''

    lines_to_remove = []
    indicating_lines = ['*** Overlap ***', '*** Kinetic Energy ***', '***** Potential Energy *****', '****** Core Hamiltonian ******', 'Orthogonalized basis functions:', '< mu | del r + r del | nu >', 'Molecular Orbital Coefficients:']

    # Second, look through the output.log file for info about the 
    found_overlap = False
    found_kinetic_energy = False
    found_potential_energy = False
    found_core_hamiltonian = False
    found_orthogonalized_basis_functions = False
    found_last_matrix_overlap = False
    found_eigenvalue_and_MO_coefficients = False

    # These variables will allow only the first molecule in the EET (the Dimer) to be read in. Any other monomers in an EET, such as the monomers will not be read.
    made_orbital_overlap_matrix_txt_file = False
    made_MO_files = False

    # ----
    # These booleans are for the get_eigenvalue_and_MO_coefficients_data method.
    found_eigenvalue = False
    recorded_all_MO_names = False
    recorded_all_electron_positions = False
    current_symbol = None
    current_index = None
    # ----

    just_switched = False
    #end_counter = None
    with open(path_to_log_file+'/'+filename,'r') as outputLOG:
        
        line_counter = 0

        for line in outputLOG:

            line_counter += 1
            # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
            # 2.1: Determine if we want to record the upcoming data
            for index in range(len(indicating_lines)):

                indicating_line = indicating_lines[index]

                if indicating_line in line:
                    max_row = 0
                    max_col = 0
                    current_max_col = 0
                    current_cols = []

                    if   index == 0:
                        print(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+': * Found data for a molecule in this output.log file.')
                        if not made_orbital_overlap_matrix_txt_file:
                            print(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+': * Found the Orbtial Overlap matrix. Will '+('not' if made_orbital_overlap_matrix_txt_file else '')+'read it in.')
                        found_overlap = True
                        lines_to_remove.append(line_counter)
                        # Reset heaps
                        overlap_data_heap = {}
                        kinetic_energy_data_heap = {}
                        potential_energy_data_heap = {}
                        core_hamiltonian_data_heap = {}
                        orthogonalized_basis_functions_data_heap = {}
                        orbital_energies_data_heap = {}
                        MO_coefficients_data_heap = {}
                        MO_coefficients_orbital_names_heap = {}
                        MO_occupancies = []

                    elif index == 1:
                        print(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+': * Found the Kinetic energy matrix. Will not read it in.')
                        #found_kinetic_energy = True

                    elif index == 2:
                        print(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+': * Found the Potential energy matrix. Will not read it in.')
                        #found_potential_energy = True

                    elif index == 3:
                        print(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+': * Found the Core Hamiltonian energy matrix. Will not read it in.')
                        #found_core_hamiltonian = True

                    elif index == 4:
                        print(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+': * Found the orthogonalized basis functions. Will not read it in.')
                        #found_orthogonalized_basis_functions = True

                    elif index == 5:
                        found_last_matrix_overlap = True

                    elif index == 6:
                        print(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+': * Found Molecular Orbital data for a molecule in this output.log file. Will '+('not' if made_orbital_overlap_matrix_txt_file else '')+'read it in.')
                        found_eigenvalue_and_MO_coefficients = True
                        lines_to_remove.append(end_counter)
                        lines_to_remove.append(line_counter)

                    just_switched = True
                    break

            if just_switched:
                just_switched = False
                continue

            # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
            # 2.1: If we are reading information about the overlap matrix, this is what to do.
            if   found_overlap:
                if not made_orbital_overlap_matrix_txt_file:
                    max_row, max_col, current_max_col, current_cols, found_overlap = get_matrix_data(line, max_row, max_col, current_max_col, current_cols, overlap_data_heap)
                    if not found_overlap:
                        # Save 2D matrix files.
                        print(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+': * Making orbital_overlap_matrix.txt file')
                        write_2D_matrix(overlap_data_heap, filename=path_to_log_file+'/'+'orbital_overlap_matrix.txt', symmetric_matrix=True)
                        print(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+': * Made orbital_overlap_matrix.txt file.')
                        del overlap_data_heap
                        made_orbital_overlap_matrix_txt_file = True
                else:
                    found_overlap = False

            elif found_kinetic_energy:
                pass #max_row, max_col, current_max_col, current_cols, found_kinetic_energy = get_matrix_data(line, max_row, max_col, current_max_col, current_cols, kinetic_energy_data_heap)
                #write_2D_matrix(kinetic_energy_data_heap, filename=path_to_log_file+'/'+str(prefix)+'_kinetic_energy_matrix.txt', symmetric_matrix=True)

            elif found_potential_energy:
                pass #max_row, max_col, current_max_col, current_cols, found_potential_energy = get_matrix_data(line, max_row, max_col, current_max_col, current_cols, potential_energy_data_heap)
                #write_2D_matrix(potential_energy_data_heap, filename=path_to_log_file+'/'+str(prefix)+'_potential_energy_matrix.txt', symmetric_matrix=True)

            elif found_core_hamiltonian:
                pass #max_row, max_col, current_max_col, current_cols, found_core_hamiltonian = get_matrix_data(line, max_row, max_col, current_max_col, current_cols, core_hamiltonian_data_heap)
                #write_2D_matrix(core_hamiltonian_data_heap, filename=path_to_log_file+'/'+str(prefix)+'_core_hamiltonian_matrix.txt', symmetric_matrix=True)

            elif found_orthogonalized_basis_functions:
                pass #max_row, max_col, current_max_col, current_cols, found_orthogonalized_basis_functions = get_matrix_data(line, max_row, max_col, current_max_col, current_cols, orthogonalized_basis_functions_data_heap, not_square=True)
                #write_2D_matrix(orthogonalized_basis_functions_data_heap, filename=path_to_log_file+'/'+str(prefix)+'_orthogonalized_basis_functions_matrix.txt', symmetric_matrix=False)

            elif found_last_matrix_overlap:
                max_row, max_col, current_max_col, current_cols, found_last_matrix_overlap = get_matrix_data(line, max_row, max_col, current_max_col, current_cols, None, not_square=True)
                if not found_last_matrix_overlap:
                    end_counter = line_counter

            elif found_eigenvalue_and_MO_coefficients:
                found_eigenvalue_and_MO_coefficients, found_eigenvalue, recorded_all_MO_names, recorded_all_electron_positions, current_symbol, current_index = get_eigenvalue_and_MO_coefficients_data(line, found_eigenvalue, recorded_all_MO_names, recorded_all_electron_positions, current_symbol, current_index, orbital_energies_data_heap, MO_coefficients_data_heap, MO_coefficients_orbital_names_heap, MO_occupancies)
                
                if not found_eigenvalue_and_MO_coefficients:
                    lines_to_remove.append(line_counter)

                    # Save MO files
                    if not made_MO_files:
                        print(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+': * Making Molecular Orbtial (MO) files')
                        write_orbital_names(MO_coefficients_orbital_names_heap, filename=path_to_log_file+'/'+'MO_orbital_names.txt')
                        print(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+': * Made MO_orbital_names.txt file')
                        write_MO_occupancies(MO_occupancies, filename=path_to_log_file+'/'+'MO_occupancies.txt')
                        print(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+': * Made MO_occupancies.txt file')
                        del MO_coefficients_orbital_names_heap
                        del MO_occupancies

                        if not get_MO_data_from_fort7_file:
                            write_1D_matrix(orbital_energies_data_heap, filename=path_to_log_file+'/'+'MO_energies.txt')
                            print(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+': * Made MO_energies.txt file')
                            write_2D_matrix(MO_coefficients_data_heap, filename=path_to_log_file+'/'+'MO_coefficients.txt', symmetric_matrix=False)
                            print(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+': * Made MO_coefficients.txt file')
                            del orbital_energies_data_heap
                            del MO_coefficients_data_heap

                        made_MO_files = True
            # ----------------------------------------------------------------------------------------------------------------------------------------------------------------

    if not lines_to_remove == sorted(lines_to_remove):
        raise Exception('Huh?')

    return lines_to_remove

# ----------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------

