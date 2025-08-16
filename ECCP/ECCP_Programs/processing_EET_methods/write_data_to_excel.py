'''
Geoffrey Weal, write_data_to_excel.py, 10/6/22

This program is designed to create an excel spreadsheet that contains all the electronic information about the dimer from the eigendata, such as hole and electron tranfer energy between monomers.
'''

import os, time
from datetime import datetime, timedelta
from xlsxwriter import Workbook

from ECCP.ECCP_Programs.processing_EET_methods.functional_and_basis_set_storage import functional_and_basis_set_storage
from ECCP.ECCP_Programs.processing_EET_methods.processing_EET_data_methods import format_worksheet, eV_to_inverse_cm, eV_to_meV

from ECCP.ECCP_Programs.processing_EET_methods.write_individual_results import write_individual_results

electronic_coupling_names = ('delta-w','Coulomb','Exact-exchange','Exchange-correlation','w-avg*Overlap','w-avg','Overlap','Total coupling')
def write_data_to_excel(electronic_coupling_data, eet_data_foldername, individual_eet_data_foldername, start_time):
    """
    This method will create an excel spreadsheet that contains all the electronic information about the dimer from the eigendata, such as hole and electron tranfer energy between monomers.

    Parameters
    ----------
    electronic_coupling_data : list
        This is a list of all electronic energy transfer values. 
    eet_data_foldername : str.
        This is the folder to place the excel files in, as well as other files created. 
    individual_eet_data_foldername : str.
        This is the folder to place txt file about each individual dimer in. 
    start_time : float
        This is the start time of this program
    """

    # As long as their is completed Gaussian data, record it in Excel and other formats. 
    if len(electronic_coupling_data) > 0:

        print('------------------------------------------------')
        print(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+' - Making Excel Spreadsheet and text files with EET Data.')

        # Next, reorganise the data into another format. look at this later. 
        organised_data = {}
        for data_information, electronic_coupling_details in electronic_coupling_data.items():
            crystal_name, dimer_name, dft_details = data_information
            if not (crystal_name in organised_data):
                organised_data[crystal_name] = {}
            if not (dimer_name in organised_data[crystal_name]):
                organised_data[crystal_name][dimer_name] = {}
            organised_data[crystal_name][dimer_name][dft_details] = electronic_coupling_details

        # -------------------------------------------------------------
        # Print data for eet_data_foldername

        print('Writing data to excel and text files.')

        # General variables for writing data to excel and text files. 
        path_to_excel_file = 'EET_Data'

        # Set up the workbook for recording data to excel spreadsheet.
        workbook = Workbook(eet_data_foldername+'/'+path_to_excel_file+'.xlsx')

        # Set up the textfile for reporting data to.
        dataTXT_energy     = open(eet_data_foldername+'/'+path_to_excel_file+'_energy.txt', 'w')
        dataTXT_wavenumber = open(eet_data_foldername+'/'+path_to_excel_file+'_wavenumber.txt', 'w')
        dataTXT_energy    .write('Molecule Name\tDimer Name\tFunctional And Basis Set\t')
        dataTXT_wavenumber.write('Molecule Name\tDimer Name\tFunctional And Basis Set\t')
        for index in range(len(electronic_coupling_names)):
            electronic_coupling_name = electronic_coupling_names[index]
            if index == 6:
                dataTXT_energy    .write(str(electronic_coupling_name)+'\t')
                dataTXT_wavenumber.write(str(electronic_coupling_name)+'\t')
            else:
                dataTXT_energy    .write(str(electronic_coupling_name)+' (meV)\t')
                dataTXT_wavenumber.write(str(electronic_coupling_name)+' (cm-1)\t')
        dataTXT_energy    .write('\n')
        dataTXT_wavenumber.write('\n')

        # Dictionaries to hp;d information that will be useful for making sheets based on data of different functionals and basis sets
        functional_and_basis_set_data = functional_and_basis_set_storage()

        # Write data to excel and text files. 
        write_main_sheet_data(organised_data, workbook, dataTXT_energy, dataTXT_wavenumber, functional_and_basis_set_data)

        # Write data about individual functionals and basis sets to excel and text files. 
        write_functional_and_basis_set_sheet_data(functional_and_basis_set_data, workbook, path_to_excel_file, eet_data_foldername)

        # Close files
        workbook.close()
        dataTXT_energy.close()
        dataTXT_wavenumber.close()

        # -------------------------------------------------------------
        # Print data for individual_eet_data_foldername

        write_individual_results(organised_data, individual_eet_data_foldername)

        # -------------------------------------------------------------

    else:
        print('There were no complete Gaussian jobs found.')

# ---------------------------------------------------------------------

def write_main_sheet_data(organised_data, workbook, dataTXT_energy, dataTXT_wavenumber, functional_and_basis_set_data):
    """
    This method will write the overall information to the excel file.

    Parameters
    ----------
    organised_data : dict.
        This contains all the EET data to save to the excel worksheet and to disk.
    workbook : xlsxwriter.Workbook
        This is the excel spreadsheet to create.
    dataTXT_energy : TXTfile
        This is the text file to save data to in meV.
    dataTXT_wavenumber : TXTfile
        This is the text file to save data to in wavenumbers.
    functional_and_basis_set_data : ECCP.ECCP_Programs.processing_EET_methods.functional_and_basis_set_storage.functional_and_basis_set_storage
        This is an object to save EET data to for each functional and basis set. 
    """
    # First, set up the sheets that will contain the main information.
    worksheet_energy = workbook.add_worksheet('Data_Energy')
    worksheet_wavenumber = workbook.add_worksheet('Data_Wavenumber')
    merge_format_mol_energy, merge_format_dimer_energy, merge_format_barrier_energy, table_format_energy, number_format_energy = format_worksheet(workbook)
    merge_format_mol_wavenumber, merge_format_dimer_wavenumber, merge_format_barrier_wavenumber, table_format_wavenumber, number_format_wavenumber = format_worksheet(workbook)

    # Second, write data to excel and text files. 
    counter_row = counter_row_initial = 0; counter_col = 0
    column_divide_difference = 11
    largest_row_length = 0
    for crystal_name, crystal_data in sorted(organised_data.items(), key=lambda x: [value.lower() for value in x[0].split('_')]):
        # Write the name of the molecule name
        worksheet_energy    .merge_range(counter_row,counter_col+1,counter_row+1,counter_col+1+7,str(crystal_name),merge_format_mol_energy)
        worksheet_wavenumber.merge_range(counter_row,counter_col+1,counter_row+1,counter_col+1+7,str(crystal_name),merge_format_mol_wavenumber)
        counter_row += 3
        for dimer_name, dimer_data in sorted(crystal_data.items(), key=lambda x: int(x[0].split('_')[0].replace('Dimer',''))):
            # Write the data about each dimer found for each molecule.
            dimer_data = sorted(dimer_data.items())
            worksheet_energy    .merge_range(counter_row,counter_col+1,counter_row,counter_col+1+7,str(dimer_name),merge_format_dimer_energy)
            worksheet_wavenumber.merge_range(counter_row,counter_col+1,counter_row,counter_col+1+7,str(dimer_name),merge_format_dimer_wavenumber)
            counter_row += 1
            # Write the titles of each energy component of hamiltonian.
            for index_col in range(len(electronic_coupling_names)):
                electronic_coupling_name = electronic_coupling_names[index_col]
                if index_col == 6:
                    worksheet_energy    .write(counter_row,counter_col+index_col+1,str(electronic_coupling_name))
                    worksheet_wavenumber.write(counter_row,counter_col+index_col+1,str(electronic_coupling_name))
                else:
                    worksheet_energy    .write(counter_row,counter_col+index_col+1,str(electronic_coupling_name)+' (meV)')
                    worksheet_wavenumber.write(counter_row,counter_col+index_col+1,str(electronic_coupling_name)+' (cm-1)')
            counter_row += 1
            # Write each component of the data for each dimer.
            # Each data in dimer_data is the result for a functional and basis set. 
            for index_row in range(len(dimer_data)):
                dft_details, (root, electronic_coupling_data) = dimer_data[index_row]
                # Write data to functional_and_basis_set_data dunctionary
                functional_and_basis_set_data.add(dft_details,crystal_name,dimer_name,electronic_coupling_data)
                # Go back to writing data into main sheets
                worksheet_energy    .write(counter_row+index_row,counter_col,str(dft_details))
                worksheet_wavenumber.write(counter_row+index_row,counter_col,str(dft_details))
                dataTXT_energy      .write(str(crystal_name)+'\t'+str(dimer_name)+'\t'+str(dft_details)+'\t')
                dataTXT_wavenumber  .write(str(crystal_name)+'\t'+str(dimer_name)+'\t'+str(dft_details)+'\t')
                for index_col in range(len(electronic_coupling_data)):
                    electronic_coupling_energy     = electronic_coupling_data[index_col]
                    electronic_coupling_wavenumber = electronic_coupling_data[index_col]
                    if not (index_col == 6):
                        electronic_coupling_energy     *= eV_to_meV # meV
                        electronic_coupling_wavenumber *= eV_to_inverse_cm # cm-1
                    worksheet_energy    .write_number(counter_row+index_row,counter_col+index_col+1,electronic_coupling_energy)
                    worksheet_wavenumber.write_number(counter_row+index_row,counter_col+index_col+1,electronic_coupling_wavenumber)
                    dataTXT_energy      .write(str(electronic_coupling_energy)+'\t')
                    dataTXT_wavenumber  .write(str(electronic_coupling_wavenumber)+'\t')
                dataTXT_energy      .write('\n')
                dataTXT_wavenumber  .write('\n')
            counter_row += (index_row+1+2)
        if counter_row > largest_row_length:
            largest_row_length = counter_row
        counter_row = counter_row_initial
        counter_col += column_divide_difference

    # Third, colour in dividing cells
    for counter_col in range(column_divide_difference-1,len(organised_data)*column_divide_difference+1,column_divide_difference):
        worksheet_energy    .merge_range(counter_row_initial,counter_col,counter_row_initial+largest_row_length,counter_col,'',merge_format_barrier_energy)
        worksheet_wavenumber.merge_range(counter_row_initial,counter_col,counter_row_initial+largest_row_length,counter_col,'',merge_format_barrier_wavenumber)

# ---------------------------------------------------------------------

def write_functional_and_basis_set_sheet_data(functional_and_basis_set_data, workbook, path_to_excel_file, eet_data_foldername):
    """
    It is also useful to separate the date based on basis set and functional

    This code will create sheets that contain the information for each basis set and functional used.

    This method will write the excel sheets for each individual functional and basis set.

    Parameters
    ----------
    functional_and_basis_set_data : ECCP.ECCP_Programs.processing_EET_methods.functional_and_basis_set_storage.functional_and_basis_set_storage
        This is an object to save EET data to for each functional and basis set. 
    workbook : xlsxwriter.Workbook
        This is the excel spreadsheet to create.
    path_to_excel_file : str.
        This is the path to the excel file
    eet_data_foldername : str.
        This is the path to the eet data folder.
    """

    # First, remove any of the previous folders below.
    energy_textfile_folder_name = 'TXT_of_Func_and_basis_sets_Energy'
    wavenumber_textfile_folder_name = 'TXT_of_Func_and_basis_sets_Wavenumber'
    for folder in [energy_textfile_folder_name, wavenumber_textfile_folder_name]:
        if os.path.exists(eet_data_foldername+'/'+folder):
            shutil.rmtree(eet_data_foldername+'/'+folder)
        os.makedirs(eet_data_foldername+'/'+folder)

    # Second, obtain the format for the spreadsheet.
    merge_format_mol_energy, merge_format_dimer_energy, merge_format_barrier_energy, table_format_energy, number_format_energy = format_worksheet(workbook)
    merge_format_mol_wavenumber, merge_format_dimer_wavenumber, merge_format_barrier_wavenumber, table_format_wavenumber, number_format_wavenumber = format_worksheet(workbook)

    # Third, write data to excel and text files. 
    counter_row = counter_row_initial = 0; 
    column_divide_difference = 11
    largest_row_length = 0
    index_row = 0
    for dft_details, all_data in sorted(functional_and_basis_set_data.items()):

        counter_col = 0

        # Set up the textfile for reporting data to.
        dataTXT_energy     = open(eet_data_foldername+'/'+energy_textfile_folder_name+'/'+path_to_excel_file+'-'+dft_details+'.txt', 'w')
        dataTXT_wavenumber = open(eet_data_foldername+'/'+wavenumber_textfile_folder_name+'/'+path_to_excel_file+'-'+dft_details+'_wavenumber.txt', 'w')
        dataTXT_energy    .write('Molecule Name\tDimer Name\t')
        dataTXT_wavenumber.write('Molecule Name\tDimer Name\t')
        for index in range(len(electronic_coupling_names)):
            electronic_coupling_name = electronic_coupling_names[index]
            if index == 6:
                dataTXT_energy    .write(str(electronic_coupling_name)+'\t')
                dataTXT_wavenumber.write(str(electronic_coupling_name)+'\t')
            else:
                dataTXT_energy    .write(str(electronic_coupling_name)+' (meV)\t')
                dataTXT_wavenumber.write(str(electronic_coupling_name)+' (cm-1)\t')
        dataTXT_energy    .write('\n')
        dataTXT_wavenumber.write('\n')

        # Set up the sheet for the functional and basis set.
        worksheet_energy     = workbook.add_worksheet(dft_details+'_E')
        worksheet_wavenumber = workbook.add_worksheet(dft_details+'_W')
        for crystal_name, dimer_data in sorted(all_data.items(), key=lambda x: [value.lower() for value in x[0].split('_')]):
            # Write the name of the molecule name
            worksheet_energy    .merge_range(counter_row,counter_col+1,counter_row+1,counter_col+1+7,str(crystal_name),merge_format_mol_energy)
            worksheet_wavenumber.merge_range(counter_row,counter_col+1,counter_row+1,counter_col+1+7,str(crystal_name),merge_format_mol_wavenumber)
            counter_row += 3

            # This will give the titles for each component of the electronic coupling analysis from Gaussian
            for index_col in range(len(electronic_coupling_names)):
                electronic_coupling_name = electronic_coupling_names[index_col]
                if index_col == 6:
                    worksheet_energy    .write(counter_row,counter_col+index_col+1,str(electronic_coupling_name))
                    worksheet_wavenumber.write(counter_row,counter_col+index_col+1,str(electronic_coupling_name))
                else:
                    worksheet_energy    .write(counter_row,counter_col+index_col+1,str(electronic_coupling_name)+' (meV)')
                    worksheet_wavenumber.write(counter_row,counter_col+index_col+1,str(electronic_coupling_name)+' (cm-1)')
            counter_row += 1

            # Get the electronic coupling data for the dimer in question
            for dimer_name, electronic_coupling_data in dimer_data:
                worksheet_energy    .write(counter_row+index_row,counter_col,str(dimer_name))
                worksheet_wavenumber.write(counter_row+index_row,counter_col,str(dimer_name))
                dataTXT_energy      .write(str(crystal_name)+'\t'+str(dimer_name)+'\t')
                dataTXT_wavenumber  .write(str(crystal_name)+'\t'+str(dimer_name)+'\t')

                # Write each component from the electronic coupling analysis
                for index_col in range(len(electronic_coupling_data)):    
                    electronic_coupling_energy     = electronic_coupling_data[index_col]
                    electronic_coupling_wavenumber = electronic_coupling_data[index_col]
                    if not (index_col == 6):
                        electronic_coupling_energy     *= eV_to_meV # meV
                        electronic_coupling_wavenumber *= eV_to_inverse_cm # cm-1
                    worksheet_energy    .write_number(counter_row+index_row,counter_col+index_col+1,electronic_coupling_energy)
                    worksheet_wavenumber.write_number(counter_row+index_row,counter_col+index_col+1,electronic_coupling_wavenumber)
                    dataTXT_energy      .write(str(electronic_coupling_energy)+'\t')
                    dataTXT_wavenumber  .write(str(electronic_coupling_wavenumber)+'\t')

                dataTXT_energy      .write('\n')
                dataTXT_wavenumber  .write('\n')

                counter_row += 1

            if counter_row > largest_row_length:
                largest_row_length = counter_row
            counter_row = counter_row_initial
            counter_col += column_divide_difference

    # Fourth, colour in dividing cells
    max_no_of_crystals = max([functional_and_basis_set_data.no_of_crystals(dft_details) for dft_details in functional_and_basis_set_data.keys()])
    for counter_col in range(column_divide_difference-1,max_no_of_crystals*column_divide_difference+1,column_divide_difference):
        worksheet_energy    .merge_range(counter_row_initial,counter_col,counter_row_initial+largest_row_length,counter_col,'',merge_format_barrier_energy)
        worksheet_wavenumber.merge_range(counter_row_initial,counter_col,counter_row_initial+largest_row_length,counter_col,'',merge_format_barrier_wavenumber)

# ---------------------------------------------------------------------



