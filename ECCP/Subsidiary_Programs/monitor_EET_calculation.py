#!/usr/bin/env python3
"""
monitor_EET_calculation.py, Geoffrey Weal, 14/8/2024

This program is designed to read the output file from Gaussian every-so-often and check if the Frag 2 State 1 <=> Frag 1 State 1 EET calculation has completed. 

Once Gaussian has done this EET calculation, this program will write a message to the end of the Gaussian output fiole saying we have finished doing what we want to do in Gaussian and cancel the job.
"""

import os, sys
from time import time, sleep

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

def gaussian_has_completed_FRAG2STATE1_to_FRAG1STATE1_EET_calc(path_to_gaussian_log_filename):

	# 
	EET_Begun                       = False
	performed_F2S1_to_F1S1_EET_calc = False
	got_delta_w                     = False
	got_Coulomb                     = False
	got_Exact_exchange              = False
	got_Exchange_correlation        = False
	got_w_avg_times_Overlap         = False
	got_Total_coupling              = False

	if not os.path.exists(path_to_gaussian_log_filename):
		return 'Not Begun EET Calc'

	#
	with open(path_to_gaussian_log_filename, 'r') as Gaussian_Log_File:

		for line in Gaussian_Log_File:

			line = line.rstrip()

			if 'Electronic Coupling for Excitation Energy Tranfer' in line:
				EET_Begun = True
				continue

			if EET_Begun:

				if ('Frag=  2 State=  1' in line) and ('Frag=  1 State=  1' in line):
					performed_F2S1_to_F1S1_EET_calc = True
					continue

				if ('delta-w                   =' in line):
					got_delta_w = True
					continue

				if ('Coulomb                   =' in line):
					got_Coulomb = True
					continue

				if ('Exact-exchange            =' in line):
					got_Exact_exchange = True
					continue

				if ('Exchange-correlation      =' in line):
					got_Exchange_correlation = True
					continue

				if ('w-avg*Overlap             =' in line) and ('w-avg=' in line) and ('Ovlp=' in line):
					got_w_avg_times_Overlap = True
					continue

				if ('Total coupling            =' in line):
					got_Total_coupling = True
					continue

	likely_calculated_job = [EET_Begun, performed_F2S1_to_F1S1_EET_calc]
	definitely_calculated_job = likely_calculated_job + [got_delta_w, got_Coulomb, got_Exact_exchange, got_Exchange_correlation, got_w_avg_times_Overlap, got_Total_coupling]

	if all(definitely_calculated_job): 
		return 'Done'
	elif all(likely_calculated_job):
		return 'Check later'
	else:
		return 'Not Begun EET Calc'



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# First, get the dirpath, gaussian_log_filename, and job_id variables
dirpath = os.getcwd()
gaussian_log_filename = str(sys.argv[1])
job_id = str(sys.argv[2])
no_of_cpus = int(sys.argv[3])

# Second, set other variables for this script.
check_log_file_interval_seconds = 1 * 60 # 20 minutes

start_time = time()

# Third, 
while True:

	command = gaussian_has_completed_FRAG2STATE1_to_FRAG1STATE1_EET_calc(dirpath+'/'+gaussian_log_filename)
		
	if command == 'Done':
		break

	elif command == 'Check later':
		sleep(60)

	elif command == 'Not Begun EET Calc':

		# 3.1: 
		sleep(check_log_file_interval_seconds)

	else:

		raise Exception('huh?')

start_end = time()

elapsed_time = start_end - start_time

minutes, seconds = divmod(elapsed_time, 60)
hours,   minutes = divmod(minutes, 60)
days,    hours   = divmod(hours, 60)

computer_time = elapsed_time * no_of_cpus

ct_minutes, ct_seconds = divmod(computer_time, 60)
ct_hours,   ct_minutes = divmod(ct_minutes, 60)
ct_days,    ct_hours   = divmod(ct_hours, 60)

with open('output.log', 'a') as outputLOG:
	outputLOG.write(' \n')
	outputLOG.write(' Elapsed time:       '+str(ct_days)+' days '+str(ct_hours)+' hours '+str(ct_minutes)+' minutes '+str(ct_seconds)+' seconds.\n')
	outputLOG.write(' Elapsed time:       '+str(days)+' days '+str(hours)+' hours '+str(minutes)+' minutes '+str(seconds)+' seconds.\n')
	outputLOG.write(' \n')
	outputLOG.write(' Normal termination of Gaussian 16 at \n')



# - - - - - - - - - - - - -

os.remove('gaussian.chk')
os.remove('gaussian.skr')
os.remove('gaussian.rwf')
os.remove('gaussian.int')
os.remove('gaussian.d2e')


os.system('scancel '+job_id)




