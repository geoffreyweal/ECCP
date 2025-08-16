#!/usr/bin/env python3
'''
Geoffrey Weal, copy_gaussian_temp_files.py, 15/3/22

This program is designed to  
'''

'''
# Determine if you want your Gaussian jobs to be submitted in parallel or run one after the other.
submit_in_parallel = True
'''

import os, sys, time
from shutil import copyfile

# ------------------------------------------------------------------------------------

def get_lastline(log_filepath):
    '''
    This method is designed to obtain the last line in a file efficiently. 

    Parameters
    ----------
    log_filepath : str
        This is the path to the log file.

    Returns
    -------
    The last line in a file.
    '''
    with open(log_filepath, 'rb') as f:
        try:  # catch OSError in case of a one line file 
            f.seek(-2, os.SEEK_END)
            while f.read(1) != b'\n':
                f.seek(-2, os.SEEK_CUR)
        except OSError:
            f.seek(0)
        last_line = f.readline().decode()
    return last_line


def get_input_gaussian_file(funct_and_basis_name):

    GJF_files = [file for file in os.listdir(funct_and_basis_name) if (file.endswith('.gjf') and os.path.isfile(funct_and_basis_name+'/'+file))]
    if len(GJF_files) > 1:
        exit('Too many gjf giles in: '+str(funct_and_basis_name))
    elif len(GJF_files) == 0:
        exit('No gjf files in: '+str(funct_and_basis_name))
    GJF_file = GJF_files[0]

    return GJF_file

def get_temp_file_names(GJF_file):

    if not os.path.isfile(GJF_file):
        exit('error in def get_temp_file_names, in copy_gaussian_temp_files.py')

    suffix_paths = {tempfile: None for tempfile in ['chk','rwf','int','d2e','skr']}  
    path_to_gjf_folder = os.path.dirname(os.path.abspath(GJF_file))
    with open(GJF_file,'r') as outputLOG:
        for line in outputLOG:
            if not (line.startswith('%') or line.startswith('#')):
                break
            for suffix in suffix_paths.keys():
                start_of_line = '%'+suffix+'='
                if (start_of_line) in line:
                    file = line.rstrip().replace(start_of_line,'')
                    given_filepath = os.path.dirname(file)
                    general_filepath = path_to_gjf_folder.replace(given_filepath,'').replace(path_to_gjf_folder,'')
                    if general_filepath == '':
                        suffix_paths[suffix] = file
                    else:
                        suffix_paths[suffix] = general_filepath+'/'+file
    return suffix_paths


month_str_to_number = {'Jan': '01', 'Feb': '02', 'Mar': '03', 'Apr': '04', 'May': '05', 'Jun': '06', 'Jul': '07', 'Aug': '08', 'Sep': '09', 'Oct': '10', 'Nov': '11', 'Dec': '12'}
def did_gaussian_job_complete(log_filepath):
    '''
    This method is designed to quickly determine if a Gaussian job has completed

    Parameters
    ----------
    log_filepath : str
        This is the path to the log file.

    Returns
    -------
    True if the log file indicates the Gaussian log file indicates that their has been normal termination, otherwise return False
    '''
    last_line = get_lastline(log_filepath)
    did_complete = ('Normal termination of Gaussian' in last_line)
    if not did_complete:
        return False, None, None
    for day in ['Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun']:
        if day in last_line:
            break
    else:
        print('Error')
        exit('finishing without completing')

    last_line = last_line.rstrip().split(day)[1].replace('.','')
    month, day_no, time, year = last_line.split()
    year = int(year)
    month = int(month_str_to_number[month])
    day_no = int(day_no)
    hour, minute, second = time.split(':')
    hour = int(hour)
    minute = int(minute)
    second = int(second)

    path_to_log_filepath = os.path.dirname(log_filepath)
    GJF_file = get_input_gaussian_file(path_to_log_filepath)
    temp_filepaths = get_temp_file_names(path_to_log_filepath+'/'+GJF_file)

    return True, (year, month, day_no, hour, minute, second), temp_filepaths

def determine_completed_jobs(path):
    completed_jobs = []
    functionals_and_basis_sets = os.listdir(path)
    for functional_and_basis_set in functionals_and_basis_sets:
        path_to_outputLOG_file = functional_and_basis_set+'/output.log'
        if not os.path.exists(path_to_outputLOG_file):
            continue
        did_complete, time_finished, temp_filepaths = did_gaussian_job_complete(path_to_outputLOG_file)
        if did_complete:
            completed_jobs.append((functional_and_basis_set, time_finished, temp_filepaths))
    return completed_jobs

# ------------------------------------------------------------------------------------

def first_time_more_recent_than_second_times(full_time1, full_time2):
    def comparison(time1, time2):
        if   time1 > time2:
            return True
        elif time1 < time2:
            return False
        else:
            return None

    year1, month1, day_no1, hour1, minute1, second1 = full_time1
    year2, month2, day_no2, hour2, minute2, second2 = full_time2

    year_comparison = comparison(year1, year2)
    if year_comparison is not None:
        return year_comparison
    month_comparison = comparison(month1, month2)
    if month_comparison is not None:
        return month_comparison
    day_comparison = comparison(day_no1, day_no2)
    if day_comparison is not None:
        return day_comparison
    hour_comparison = comparison(hour1, hour2)
    if hour_comparison is not None:
        return hour_comparison
    minute_comparison = comparison(minute1, minute2)
    if minute_comparison is not None:
        return minute_comparison
    second_comparison = comparison(second1, second2)
    if second_comparison is not None:
        return second_comparison
    return False

def get_most_recent_temp_filepath(completed_jobs):
    most_recently_completed_filepath = None
    most_recently_completed_temp_filepaths = None
    most_recently_completed_filepath_time = (0,0,0,0,0,0)
    for functional_and_basis_set, time_finished, temp_filepaths in completed_jobs:
        if first_time_more_recent_than_second_times(time_finished, most_recently_completed_filepath_time):
            most_recently_completed_filepath = functional_and_basis_set
            most_recently_completed_temp_filepaths = temp_filepaths
            most_recently_completed_filepath_time = time_finished
    return most_recently_completed_filepath, most_recently_completed_temp_filepaths

# ------------------------------------------------------------------------------------

def write_restart_lines_Link_0_section(where_to_place_previous_temp_files_into):
    toString = ''
    for suffix in ['chk','rwf']: # ,'int','d2e','skr']
        tempfilepath = where_to_place_previous_temp_files_into[suffix]
        toString += r'%old'+str(suffix)+'='+str(tempfilepath)+'\n'
    return toString

def write_restart_lines_Route_section(restart_line):
    return str(restart_line)+'\n'

def copy_gaussian_temp_files(most_recent_temp_filepaths, where_to_place_previous_temp_files_into):
    for suffix in ['chk','rwf']: #,'int','d2e','skr']:
        most_recent_temp_filepath = most_recent_temp_filepaths[suffix]
        funct_and_basis_name_next = where_to_place_previous_temp_files_into[suffix]
        copyfile(most_recent_temp_filepath,funct_and_basis_name_next)

def add_restart_line_to_inputGJF_file(funct_and_basis_name_next, most_recent_temp_filepaths):

    GJF_file = get_input_gaussian_file(funct_and_basis_name_next)

    newGJFfilename = funct_and_basis_name_next+'/'+GJF_file+'.temp'
    oldGJFfilename = funct_and_basis_name_next+'/'+GJF_file

    where_to_place_previous_temp_files_into = get_temp_file_names(oldGJFfilename)
    for key in where_to_place_previous_temp_files_into.keys():
        where_to_place_previous_temp_files_into[key] += '.prev'
    
    copy_gaussian_temp_files(most_recent_temp_filepaths, where_to_place_previous_temp_files_into)

    found_restart_line_already = False
    restart_line = '# Geom=Checkpoint Guess=Read'
    with open(newGJFfilename,'w') as newGJFfile:
        with open(oldGJFfilename,'r') as oldGJFfile:
            start_of_input_parameters = False
            end_of_input_parameters   = False
            newGJFfile.write(write_restart_lines_Link_0_section(where_to_place_previous_temp_files_into))
            for line in oldGJFfile:
                if restart_line in line:
                    found_restart_line_already = True
                    break
                if end_of_input_parameters:
                    pass
                elif line.startswith('%') or line.startswith('#'):
                    start_of_input_parameters = True
                elif start_of_input_parameters and (not (line.startswith('%') or line.startswith('#'))):
                    end_of_input_parameters = True
                    newGJFfile.write(write_restart_lines_Route_section(restart_line))
                newGJFfile.write(line)

    if found_restart_line_already:
        os.remove(newGJFfilename)
    else:
        os.remove(oldGJFfilename)
        os.rename(newGJFfilename, oldGJFfilename)

# ------------------------------------------------------------------------------------

current_path = os.getcwd()
funct_and_basis_name_next = str(sys.argv[1])
completed_jobs = determine_completed_jobs(current_path)
if not len(completed_jobs) == 0:
    most_recent_temp_filepath, most_recent_temp_filepaths = get_most_recent_temp_filepath(completed_jobs)
    add_restart_line_to_inputGJF_file(funct_and_basis_name_next, most_recent_temp_filepaths)

