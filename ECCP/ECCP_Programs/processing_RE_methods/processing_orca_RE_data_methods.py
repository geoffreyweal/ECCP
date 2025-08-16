'''
Geoffrey Weal, processing_OPV_Dimer_data.py, 9/3/22

This script contains methods for processing_OPV_Dimer_data.py

'''

from ECCP.ECCP_Programs.shared_general_methods.shared_general_methods    import reverse_readline
from ECCP.ECCP_Programs.processing_RE_methods.processing_RE_data_methods import is_finished_reading
from ECCP.ECCP_Programs.shared_general_methods.shared_gaussian_methods   import did_gaussian_job_complete

# -----------------------------------------------------------------

def get_energy_from_opt_job(log_filepath, energy_state, opt_job=False):
    """
    This method will obtain the optimised energy for the job if the job converged.
    """

    energy = None
    maximum_force_converged    = None
    rms_force_converged        = None
    maximum_distance_converged = None
    rms_distance_converged     = None
    
    for line in reverse_readline(log_filepath):
        line = line.rstrip()

        if 'FINAL SINGLE POINT ENERGY:' in line:
            # Gaussian energy is found here:
            line = line.split()
            energy = float(line[4]) # Hartrees
            if is_finished_reading(energy, maximum_force_converged, rms_force_converged, maximum_distance_converged, rms_distance_converged, opt_job):
                break

        if opt_job:
            if ('THE OPTIMIZATION HAS CONVERGED' in line):
                maximum_force_converged = True
                rms_force_converged = True
                if is_finished_reading(energy, maximum_force_converged, rms_force_converged, maximum_distance_converged, rms_distance_converged, opt_job):
                    break
            if 'MAX gradient' in line:
                line = line.split()
                maximum_force_converged = line[-1]
                maximum_force_converged = True if (maximum_force_converged == 'YES') else False
                if is_finished_reading(energy, maximum_force_converged, rms_force_converged, maximum_distance_converged, rms_distance_converged, opt_job):
                    break
            if 'RMS gradient' in line:
                line = line.split()
                rms_force_converged = line[-1]
                rms_force_converged = True if (rms_force_converged == 'YES') else False
                if is_finished_reading(energy, maximum_force_converged, rms_force_converged, maximum_distance_converged, rms_distance_converged, opt_job):
                    break
            if 'MAX step' in line:
                line = line.split()
                maximum_distance_converged = line[-1]
                maximum_distance_converged = True if (maximum_distance_converged == 'YES') else False
                if is_finished_reading(energy, maximum_force_converged, rms_force_converged, maximum_distance_converged, rms_distance_converged, opt_job):
                    break
            if 'RMS step' in line:
                line = line.split()
                rms_distance_converged = line[-1]
                rms_distance_converged = True if (rms_distance_converged == 'YES') else False
                if is_finished_reading(energy, maximum_force_converged, rms_force_converged, maximum_distance_converged, rms_distance_converged, opt_job):
                    break

    if opt_job and any(entry is None for entry in [energy, maximum_force_converged, rms_force_converged, maximum_distance_converged, rms_distance_converged]):
        toString  = 'Error: '+str(log_filepath)+'\n'
        toString += 'energy: '+str(energy)+'\n'
        toString += 'maximum_force_converged: '+str(maximum_force_converged)+'\n'
        toString += 'rms_force_converged: '+str(rms_force_converged)+'\n'
        toString += 'maximum_distance_converged: '+str(maximum_distance_converged)+'\n'
        toString += 'rms_distance_converged: '+str(rms_distance_converged)+'\n'
        raise Exception(toString)

    if opt_job:
        return energy, maximum_force_converged, rms_force_converged, maximum_distance_converged, rms_distance_converged
    else:
        return energy

def get_frequencies_from_freq_job(log_filepath):
    """
    This method will check to see if the gaussian job has completed successfully. 

    Parameters
    ----------
    log_filepath : str.
        This is the path to the Gaussian job log file.
    """
    # First, determine if the job finished successfully
    #        This is to check that the force convergence using the full Hessian
    #        rather than the approximate from the optimisation is all good.
    did_orca_job_terminate_normally = did_orca_job_complete(freq_log_filepath) 
    if not did_orca_job_terminate_normally:
        return did_orca_job_terminate_normally, None, None

    # Second, obtain the frequencies from the output file.
    frequency_start_line = 'VIBRATIONAL FREQUENCIES'
    frequency_end_line   = 'NORMAL MODES'
    start_recording = False
    all_frequencies = []
    with open(freq_log_filepath) as logFILE:
        for line in logFILE:
            # 2.1: Check if to start frequency finding or finishing looking for frequencies
            if frequency_start_line in line:
                start_recording = True
                continue
            if frequency_end_line in line:
                break
            if not start_recording:
                continue
            # 2.2: Look and record freqencies from file.
            line = line.rstrip().split()
            if len(line) == 0:
                continue
            elif line[0].replace(':','').isdigit():
                frequency = float(line[1])
                all_frequencies.append(frequency)
    all_frequencies.sort()

    # Third, return frequencies found.
    return all_frequencies

# -----------------------------------------------------------------





