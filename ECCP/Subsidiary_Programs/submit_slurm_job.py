#!/usr/bin/env python3
'''
Geoffrey Weal, submit_slurm_job.py, 21/4/22

This program is designed to submit a slurm submission script to slurm
'''
import sys, subprocess

submitted_job_filename = sys.argv[1]
print('Submitting '+str(submitted_job_filename)+' to slurm')
proc = subprocess.Popen(['sbatch', submitted_job_filename], shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
outs, errs = proc.communicate()
print(outs.decode("utf-8"))
print(errs.decode("utf-8"), file=sys.stderr)
print('Submitted '+str(submitted_job_filename)+' to slurm')