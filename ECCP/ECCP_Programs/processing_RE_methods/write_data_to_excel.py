'''
Geoffrey Weal, write_data_to_excel.py, 10/6/22

This program is designed to create an excel spreadsheet that contains all the electronic information about the dimer from the eigendata, such as hole and electron tranfer energy between monomers.
'''

import os, time
from datetime import datetime, timedelta
from xlsxwriter import Workbook

from ECCP.ECCP_Programs.processing_RE_methods.functional_and_basis_set_storage import functional_and_basis_set_storage
from ECCP.ECCP_Programs.processing_RE_methods.processing_RE_data_methods import format_worksheet, get_reorganisation_energy, get_energy_diff, convert_hartree_to_eV, get_hartree_to_eV_conversion_value, eV_to_meV

from ECCP.ECCP_Programs.processing_RE_methods.write_individual_results import write_individual_results

def write_data_to_excel(reorganisation_energy_data, re_data_foldername, individual_re_data_foldername, start_time):
    """
    This method will create an excel spreadsheet that contains all the reorganisation energies of molecules in the crystal. 

    Parameters
    ----------
    reorganisation_energy_data : list
        This is a list of all reorganisation energy data of molecules in the crystal. 
    re_data_foldername : str.
        This is the folder to place the excel files in, as well as other files created. 
    individual_re_data_foldername : str.
        This is the folder to place txt file about each individual dimer in. 
    start_time : float
        This is the start time of this program
    """

    # As long as their is completed Gaussian data, record it in Excel and other formats. 
    if len(reorganisation_energy_data) > 0:

        print('------------------------------------------------')
        print(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+' - Making Excel Spreadsheet and text files with RE Data.')

        # Next, reorganise the data into another format. look at this later. 
        organised_data = {}
        for data_information, reorganisation_energy_details in reorganisation_energy_data.items():
            crystal_name, molecule_name, dft_details = data_information
            if not (crystal_name in organised_data):
                organised_data[crystal_name] = {}
            if not (molecule_name in organised_data[crystal_name]):
                organised_data[crystal_name][molecule_name] = {}
            organised_data[crystal_name][molecule_name][dft_details] = reorganisation_energy_details

        # -------------------------------------------------------------
        # Print data for eet_data_foldername

        print('Writing data to excel and text files.')

        # General variables for writing data to excel and text files. 
        path_to_excel_file = 'RE_Data'

        # Set up the workbook for recording data to excel spreadsheet.
        workbook = Workbook(re_data_foldername+'/'+path_to_excel_file+'.xlsx')

        # Set up the textfile for reporting data to.
        dataTXT_energy = open(re_data_foldername+'/'+path_to_excel_file+'_energy.txt', 'w')
        dataTXT_energy.write('Crystal Name\tMolecule Name\tFunctional And Basis Set\t')
        dataTXT_energy.write('E_GS(GS) (Ha)\tE_GS(ES) (Ha)\tE_ES(GS) (Ha)\tE_ES(ES) (Ha)\t Band Gap (E_ES(ES) - E_GS(GS)) (eV)')
        dataTXT_energy.write('\n')

        # Dictionaries to hp;d information that will be useful for making sheets based on data of different functionals and basis sets
        functional_and_basis_set_data = functional_and_basis_set_storage()

        # Write data to excel and text files. 
        write_main_sheet_data(organised_data, workbook, dataTXT_energy, functional_and_basis_set_data)

        # Write data about individual functionals and basis sets to excel and text files. 
        write_functional_and_basis_set_sheet_data(functional_and_basis_set_data, workbook, path_to_excel_file, re_data_foldername)

        # Close files
        workbook.close()
        dataTXT_energy.close()

        # -------------------------------------------------------------
        # Print data for individual_eet_data_foldername
        write_individual_results(organised_data, individual_re_data_foldername)
        # -------------------------------------------------------------

    else:
        print('There were no complete Gaussian jobs found.')

# ---------------------------------------------------------------------

def write_main_sheet_data(organised_data, workbook, dataTXT_energy, functional_and_basis_set_data):
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
    functional_and_basis_set_data : ECCP.ECCP_Programs.processing_EET_methods.functional_and_basis_set_storage.functional_and_basis_set_storage
        This is an object to save EET data to for each functional and basis set. 
    """
    # First, set up the sheets that will contain the main information.
    worksheet_energy = workbook.add_worksheet('Data_Energy')
    merge_format_mol_energy, merge_format_dimer_energy, merge_format_barrier_energy, table_format_energy, number_format_energy = format_worksheet(workbook)
    normal = workbook.add_format({})
    subscript = workbook.add_format({'font_script': 2})
    worksheet_energy.write(2,0,'NOTE: E(eV) = E(Hartree) * '+str(get_hartree_to_eV_conversion_value()))

    # Second, write data to excel and text files. 
    counter_row = counter_row_initial = 0; counter_col = 0
    column_divide_difference = 11
    largest_row_length = 0
    for crystal_name, crystal_data in sorted(organised_data.items(), key=lambda x: [value.lower() for value in x[0].split('_')]):
        # Write the name of the molecule name
        worksheet_energy.merge_range(counter_row,counter_col+1,counter_row+1,counter_col+1+7,str(crystal_name),merge_format_mol_energy)
        counter_row += 3
        for mol_name, mol_data in sorted(crystal_data.items(), key=lambda x: int(x[0].split('_')[1].replace('S',''))):
            # Write the data about each molecule found for each crystal.
            mol_data = sorted(mol_data.items())
            worksheet_energy.merge_range(counter_row,counter_col+1,counter_row,counter_col+1+7,str(mol_name),merge_format_dimer_energy)
            counter_row += 1
            # Write the names of the energies of various electronic and structural states here
            worksheet_energy.write_rich_string(counter_row,counter_col+0+1,'E',subscript,'GS',normal,'(GS) (Ha)')
            worksheet_energy.write_rich_string(counter_row,counter_col+1+1,'E',subscript,'GS',normal,'(ES) (Ha)')
            worksheet_energy.write_rich_string(counter_row,counter_col+2+1,'E',subscript,'ES',normal,'(GS) (Ha)')
            worksheet_energy.write_rich_string(counter_row,counter_col+3+1,'E',subscript,'ES',normal,'(ES) (Ha)')
            worksheet_energy.write_rich_string(counter_row,counter_col+4+1,'Band Gap (E',subscript,'ES',normal,'(ES) - E',subscript,'GS',normal,'(GS)) (eV)')
            worksheet_energy.write_string(counter_row,counter_col+5+1,'Negative Frequencies GS_GS')
            worksheet_energy.write_string(counter_row,counter_col+6+1,'Negative Frequencies ES_ES')
            counter_row += 1
            # Write each component of the data for each dimer.
            # Each data in mol_data is the result for a functional and basis set. 
            for index_row in range(len(mol_data)):
                dft_details, electronic_data = mol_data[index_row]
                root, eGS_gGS_energy, eES_gGS_energy, eGS_gES_energy, eES_gES_energy, negative_eGS_gGS_freqs, negative_eES_gES_freqs = electronic_data
                # Get the band gap
                band_gap = convert_hartree_to_eV(get_energy_diff(eES_gES_energy, eGS_gGS_energy))
                # Write data to functional_and_basis_set_data dunctionary
                functional_and_basis_set_data.add(dft_details,crystal_name,mol_name,electronic_data)
                # Go back to writing data into main sheets
                worksheet_energy.write(counter_row+index_row,counter_col,str(dft_details))
                dataTXT_energy  .write(str(crystal_name)+'\t'+str(mol_name)+'\t'+str(dft_details)+'\t')
                # record reorganisation energy data
                info_to_record = [eGS_gGS_energy, eGS_gES_energy, eES_gGS_energy, eES_gES_energy, band_gap]
                for index_col in range(len(info_to_record)):
                    reorganisation_energy_value = info_to_record[index_col]
                    #if (index_col in (len(info_to_record)-2, len(info_to_record)-1)):
                    #    reorganisation_energy_value *= eV_to_meV # meV
                    worksheet_energy.write_number(counter_row+index_row,counter_col+index_col+1,reorganisation_energy_value)
                    dataTXT_energy  .write(str(reorganisation_energy_value)+'\t')
                dataTXT_energy.write('\n')
                # Give a note if there are negative frequencies
                if len(negative_eGS_gGS_freqs) > 0:
                    worksheet_energy.write_string(counter_row+index_row,counter_col+len(info_to_record)+1,str(negative_eGS_gGS_freqs).replace('[','').replace(']',''))
                if len(negative_eES_gES_freqs) > 0:
                    worksheet_energy.write_string(counter_row+index_row,counter_col+len(info_to_record)+2,str(negative_eES_gES_freqs).replace('[','').replace(']',''))
            counter_row += (index_row+1+2)
        if counter_row > largest_row_length:
            largest_row_length = counter_row
        counter_row = counter_row_initial
        counter_col += column_divide_difference

    # Third, colour in dividing cells
    for counter_col in range(column_divide_difference-1,len(organised_data)*column_divide_difference+1,column_divide_difference):
        worksheet_energy.merge_range(counter_row_initial,counter_col,counter_row_initial+largest_row_length,counter_col,'',merge_format_barrier_energy)

# ---------------------------------------------------------------------

def write_functional_and_basis_set_sheet_data(functional_and_basis_set_data, workbook, path_to_excel_file, re_data_foldername):
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
    re_data_foldername : str.
        This is the path to the eet data folder.
    """

    # First, remove any of the previous folders below.
    energy_textfile_folder_name = 'TXT_of_Func_and_basis_sets_Energy'
    if os.path.exists(re_data_foldername+'/'+energy_textfile_folder_name):
        shutil.rmtree(re_data_foldername+'/'+energy_textfile_folder_name)
    os.makedirs(re_data_foldername+'/'+energy_textfile_folder_name)

    # Second, obtain the format for the spreadsheet.
    merge_format_mol_energy, merge_format_dimer_energy, merge_format_barrier_energy, table_format_energy, number_format_energy = format_worksheet(workbook)
    normal = workbook.add_format({})
    subscript = workbook.add_format({'font_script': 2})

    # Third, write data to excel and text files. 
    counter_row = counter_row_initial = 0; 
    column_divide_difference = 11
    largest_row_length = 0
    index_row = 0

    for dft_details, all_data in sorted(functional_and_basis_set_data.items()):

        counter_col = 0

        # Set up the textfile for reporting data to.
        dataTXT_Egap = open(re_data_foldername+'/'+energy_textfile_folder_name+'/'+path_to_excel_file+'-'+dft_details+'_energies.txt', 'w')
        dataTXT_Egap.write('Crystal Name\tMolecule Name\tE_GS(GS) (Ha)\tE_GS(ES) (Ha)\tE_ES(GS) (Ha)\tE_ES(ES) (Ha)\tBand Gap (E_ES(ES) - E_GS(GS)) (eV)')
        dataTXT_Egap.write('\n')

        dataTXT_reorg_energy = open(re_data_foldername+'/'+energy_textfile_folder_name+'/'+path_to_excel_file+'-'+dft_details+'_reorg_energies.txt', 'w')
        dataTXT_reorg_energy.write('Crystal Name\tDonor Molecule\tAcceptor Molecule\tReorganisation Energy (meV)')
        dataTXT_reorg_energy.write('\n')

        # Set up the sheet for the functional and basis set.
        worksheet_energy     = workbook.add_worksheet(dft_details+'_E')
        #worksheet_energy.write_rich_string(2,0,'NOTE 1: RE = (E',subscript,'ES',normal,'(GS) - E',subscript,'ES',normal,'(ES)) + (E',subscript,'GS',normal,'(ES) - E',subscript,'GS',normal,'(GS)) [https://doi.org/10.1039/D0TC05697A]')
        #worksheet_energy.write(2,8,'NOTE 2: E(eV) = E(Hartree) * '+str(get_hartree_to_eV_conversion_value()))
        worksheet_energy.write(2,0,'NOTE: E(eV) = E(Hartree) * '+str(get_hartree_to_eV_conversion_value()))
        for crystal_name, crystal_data in sorted(all_data.items(), key=lambda x: [value.lower() for value in x[0].split('_')]):
            # Write the name of the molecule name
            worksheet_energy    .merge_range(counter_row,counter_col+1,counter_row+1,counter_col+1+7,str(crystal_name),merge_format_mol_energy)
            counter_row += 3

            # This will give the titles for each component of the reorganisation energy data from Gaussian
            worksheet_energy.write_rich_string(counter_row,counter_col+0+1,'E',subscript,'GS',normal,'(GS) (Ha)')
            worksheet_energy.write_rich_string(counter_row,counter_col+1+1,'E',subscript,'GS',normal,'(ES) (Ha)')
            worksheet_energy.write_rich_string(counter_row,counter_col+2+1,'E',subscript,'ES',normal,'(GS) (Ha)')
            worksheet_energy.write_rich_string(counter_row,counter_col+3+1,'E',subscript,'ES',normal,'(ES) (Ha)')
            worksheet_energy.write_rich_string(counter_row,counter_col+4+1,'Band Gap (E',subscript,'ES',normal,'(ES) - E',subscript,'GS',normal,'(GS)) (eV)')
            worksheet_energy.write_rich_string(counter_row,counter_col+5+1,'Negative Frequencies GS_GS')
            worksheet_energy.write_rich_string(counter_row,counter_col+6+1,'Negative Frequencies ES_ES')
            counter_row += 1

            # Get the energy gap data for the molecule in question
            for molecule_name, reorganisation_energy_data in crystal_data:
                worksheet_energy.write(counter_row+index_row,counter_col,str(molecule_name))
                dataTXT_Egap  .write(str(crystal_name)+'\t'+str(molecule_name)+'\t')

                # Write each component from the electronic coupling analysis
                root, eGS_gGS_energy, eES_gGS_energy, eGS_gES_energy, eES_gES_energy, negative_eGS_gGS_freqs, negative_eES_gES_freqs = reorganisation_energy_data

                # Get the band gap
                band_gap = convert_hartree_to_eV(get_energy_diff(eES_gES_energy, eGS_gGS_energy))

                reorganisation_energy_data_new = (eGS_gGS_energy, eGS_gES_energy, eES_gGS_energy, eES_gES_energy, band_gap)
                for index_col in range(len(reorganisation_energy_data_new)):    
                    electronic_coupling_energy = reorganisation_energy_data_new[index_col]
                    #if (index_col == (len(reorganisation_energy_data_new)-1)):
                    #    electronic_coupling_energy *= eV_to_meV # meV
                    worksheet_energy.write_number(counter_row+index_row,counter_col+index_col+1,electronic_coupling_energy)
                    dataTXT_Egap  .write(str(electronic_coupling_energy)+'\t')
                # Give a note if there are negative frequencies
                if len(negative_eGS_gGS_freqs) > 0:
                    worksheet_energy.write_rich_string(counter_row+index_row,counter_col+len(reorganisation_energy_data_new)+1,str(negative_eGS_gGS_freqs).replace('[','').replace(']',''))
                if len(negative_eES_gES_freqs) > 0:
                    worksheet_energy.write_rich_string(counter_row+index_row,counter_col+len(reorganisation_energy_data_new)+2,str(negative_eES_gES_freqs).replace('[','').replace(']',''))
                dataTXT_Egap.write('\n')
                counter_row += 1

            # Get the reorganisation energies between molecules in the crystal.
            counter_row += 1
            counter_row += 1
            worksheet_energy.write(counter_row+index_row,counter_col,'Reorganisation Energies between Molecules in Crystal')
            counter_row += 1
            worksheet_energy.write(counter_row+index_row,counter_col,'Exciton')
            counter_row += 1
            worksheet_energy.write(counter_row+index_row,counter_col+0,'Donor')
            worksheet_energy.write(counter_row+index_row,counter_col+1,'Acceptor')
            worksheet_energy.write(counter_row+index_row,counter_col+2,'Reorganisation (meV)')
            counter_row += 1
            for mol1_index in range(len(crystal_data)):
                molecule_name1, reorganisation_energy_data1 = crystal_data[mol1_index]
                # First, get the excited structure geometry state of mol1 (where exciton is coming from) (donor)
                root1, donor_eGS_gGS_energy, donor_eES_gGS_energy, donor_eGS_gES_energy, donor_eES_gES_energy, negative_eGS_gGS_freqs, negative_eES_gES_freqs = reorganisation_energy_data1
                for mol2_index in range(len(crystal_data)):
                    molecule_name2, reorganisation_energy_data2 = crystal_data[mol2_index]
                    # Second, get the ground structure geometry state of mol2 (where exciton is going to) (acceptor)
                    root2, acceptor_eGS_gGS_energy, acceptor_eES_gGS_energy, acceptor_eGS_gES_energy, acceptor_eES_gES_energy, negative_eGS_gGS_freqs, negative_eES_gES_freqs = reorganisation_energy_data2
                    worksheet_energy.write(counter_row+index_row,counter_col+0,str(molecule_name1))
                    worksheet_energy.write(counter_row+index_row,counter_col+1,str(molecule_name2))
                    # Get the reorganisation energy in meV
                    reorganisation_energy = convert_hartree_to_eV(get_reorganisation_energy(donor_eGS_gGS_energy, acceptor_eES_gGS_energy, donor_eGS_gES_energy, acceptor_eES_gES_energy)) * eV_to_meV
                    worksheet_energy.write(counter_row+index_row,counter_col+2,str(reorganisation_energy))
                    dataTXT_reorg_energy.write(str(molecule_name1)+'\t'+str(molecule_name2)+'\t'+str(reorganisation_energy)+'\n')
                    counter_row += 1

            if counter_row > largest_row_length:
                largest_row_length = counter_row
            counter_row = counter_row_initial
            counter_col += column_divide_difference

    # Fourth, colour in dividing cells
    max_no_of_crystals = max([functional_and_basis_set_data.no_of_crystals(dft_details) for dft_details in functional_and_basis_set_data.keys()])
    for counter_col in range(column_divide_difference-1,max_no_of_crystals*column_divide_difference+1,column_divide_difference):
        worksheet_energy.merge_range(counter_row_initial,counter_col,counter_row_initial+largest_row_length,counter_col,'',merge_format_barrier_energy)

# ---------------------------------------------------------------------



