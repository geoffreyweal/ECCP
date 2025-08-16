'''
Geoffrey Weal, obtain_gaussian_RE_data.py, 15/7/23

This program is designed to obtain the reorganisation energy data from Gaussian calculations. 

'''
import time
from datetime import datetime, timedelta

from ECCP.ECCP_Programs.processing_RE_methods.processing_orca_RE_data_methods import get_energy_from_opt_job, get_frequencies_from_freq_job
from ECCP.ECCP_Programs.shared_general_methods.shared_orca_methods            import did_orca_opt_job_complete, did_orca_job_complete, orca_temp_files_to_remove

def obtain_orca_RE_data(root, reorganisation_energy_data, start_time, ground_structure_foldername, excited_structure_foldername, lower_limit_negative_frequency, issues):
    """
    This method is designed to obtain the reorganisation energy data from Gaussian calculations. 
    """

    print('------------------------------------------------')
    print(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+': Found a Reorganisation Energy Gaussian jobset.')

    # Paths to reorganisation energy out files.
    eGS_gGS_energy_outpath = root+'/'+ground_structure_foldername +'/'+'eGS_gGS_main_opt.out'
    eGS_gGS_freq_outpath   = root+'/'+ground_structure_foldername +'/'+'eGS_gGS_freq.out'
    eES_gGS_energy_outpath = root+'/'+ground_structure_foldername +'/'+'eES_gGS.out'

    eES_gES_energy_outpath = root+'/'+excited_structure_foldername+'/'+'eES_gES_main_opt.out'
    eES_gES_freq_outpath   = root+'/'+excited_structure_foldername+'/'+'eES_gES_freq.out'
    eGS_gES_energy_outpath = root+'/'+excited_structure_foldername+'/'+'eGS_gES.out'

    # Determine if all the reorganisation energy job completed or not. 
    got_eGS_gGS_energy = did_orca_opt_job_complete(eGS_gGS_energy_outpath)[0]
    got_eGS_gGS_freq   = did_orca_job_complete(eGS_gGS_freq_outpath)
    got_eES_gGS_energy = did_orca_job_complete(eES_gGS_energy_outpath)

    got_eES_gES_energy = did_orca_opt_job_complete(eES_gES_energy_outpath)[0]
    got_eES_gES_freq   = did_orca_job_complete(eES_gES_freq_outpath)
    got_eGS_gES_energy = did_orca_job_complete(eGS_gES_energy_outpath)

    if (got_eGS_gGS_energy and got_eGS_gGS_freq and got_eES_gGS_energy) and (got_eES_gES_energy and got_eES_gES_freq and got_eGS_gES_energy):

        print(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+': Processing: '+str(root))

        # Get data from out files
        eGS_gGS_energy, eGS_gGS_maximum_force_converged, eGS_gGS_rms_force_converged, eGS_gGS_maximum_distance_converged, eGS_gGS_rms_distance_converged = get_energy_from_refine_opt_job(eGS_gGS_energy_outpath,'GS',opt_job=True)
        eGS_gGS_freqs  = get_frequencies_from_freq_job(eGS_gGS_freq_outpath)
        eES_gGS_energy = get_energy_from_refine_opt_job(eES_gGS_energy_outpath,'ES',opt_job=False)

        eES_gES_energy, eES_gES_maximum_force_converged, eES_gES_rms_force_converged, eES_gES_maximum_distance_converged, eES_gES_rms_distance_converged = get_energy_from_refine_opt_job(eES_gES_energy_outpath,'ES',opt_job=True)
        eES_gES_freqs  = get_frequencies_from_freq_job(eES_gES_freq_outpath)
        eGS_gES_energy = get_energy_from_refine_opt_job(eGS_gES_energy_outpath,'GS',opt_job=False)

        # Check results
        eGS_gGS_energy_convergence_error = not all([result for result in [eGS_gGS_maximum_force_converged, eGS_gGS_rms_force_converged]]) # , eGS_gGS_maximum_distance_converged, eGS_gGS_rms_distance_converged]])
        eES_gES_energy_convergence_error = not all([result for result in [eES_gES_maximum_force_converged, eES_gES_rms_force_converged]]) # , eES_gES_maximum_distance_converged, eES_gES_rms_distance_converged]])

        acceptable_negative_eGS_gGS_freq = (eGS_gGS_freqs is None) or (isinstance(eGS_gGS_freqs,list) and any([(freq < lower_limit_negative_frequency) for freq in eGS_gGS_freqs]))
        acceptable_negative_eES_gES_freq = (eES_gES_freqs is None) or (isinstance(eES_gES_freqs,list) and any([(freq < lower_limit_negative_frequency) for freq in eES_gES_freqs]))

        # Determine if there were any errors
        eGS_gGS_error = (eGS_gGS_energy is None) or eGS_gGS_energy_convergence_error or acceptable_negative_eGS_gGS_freq
        eES_gGS_error = (eES_gGS_energy is None)
        eES_gES_error = (eES_gES_energy is None) or eES_gES_energy_convergence_error or acceptable_negative_eES_gES_freq
        eGS_gES_error = (eGS_gES_energy is None)

        # If error is returned, report this Gaussian job as an issue. 
        if eGS_gGS_error or eES_gGS_error or eES_gES_error or eGS_gES_error:
            print(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+': Their was an issue with this job. Path: '+str(root))
            eGS_gGS_energy_data = (eGS_gGS_energy, eGS_gGS_maximum_force_converged, eGS_gGS_rms_force_converged, eGS_gGS_maximum_distance_converged, eGS_gGS_rms_distance_converged)
            eES_gES_energy_data = (eES_gES_energy, eES_gES_maximum_force_converged, eES_gES_rms_force_converged, eES_gES_maximum_distance_converged, eES_gES_rms_distance_converged)
            only_neg_eGS_gGS_freqs = [eGS_gGS_freq for eGS_gGS_freq in eGS_gGS_freqs if eGS_gGS_freq < lower_limit_negative_frequency]
            only_neg_eES_gES_freqs = [eES_gES_freq for eES_gES_freq in eES_gES_freqs if eES_gES_freq < lower_limit_negative_frequency]
            issues.append((root, (eGS_gGS_energy_data, only_neg_eGS_gGS_freqs, eES_gGS_energy, eES_gES_energy_data, only_neg_eES_gES_freqs, eGS_gES_energy)))
            
        else:

            negative_eGS_gGS_freqs = [freq for freq in eGS_gGS_freqs if (freq < 0.0)]
            negative_eES_gES_freqs = [freq for freq in eES_gES_freqs if (freq < 0.0)]

            # Record data, and remove temp Gaussian files.
            calculation_details = tuple(root.split('/')[-3:]) # This tuple contains (crystal_name, Dimer_name, Functional_and_Basis_Set_name)
            
            # Record data
            reorganisation_energy_data[calculation_details] = (root, eGS_gGS_energy, eES_gGS_energy, eGS_gES_energy, eES_gES_energy, negative_eGS_gGS_freqs, negative_eES_gES_freqs)
            print('Processed RE data in (HH:MM:SS): '+str(timedelta(seconds=time.time() - start_time)))
            
            # Remove unnecessary files.
            #orca_temp_files_to_remove(root, files, remove_fort7_file=True)
        
    else:

        # Report this Gaussian job.
        print(str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))+': Their was an issue with this job. Path: '+str(root))
        issues.append((root, (got_eGS_gGS_energy, got_eGS_gGS_freq, got_eES_gGS_energy, got_eES_gES_energy, got_eES_gES_freq, got_eGS_gES_energy)))
