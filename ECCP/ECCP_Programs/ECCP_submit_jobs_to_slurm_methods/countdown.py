"""
countdown.py, Geoffrey Weal, 4/1/2023

This method will provide a waiting countdown for the user
"""
import sys
from time import sleep

def countdown(t):
    """
    This method will wait with a countdown until the time has run out

    Parameters
    ----------
    t : int
        This is the number of seconds to wait for.
    """
    print('Will wait for ' + str(float(t)/60.0) + ' minutes, then will resume Slurm submissions.\n')
    while t:
        mins, secs = divmod(t, 60)
        timeformat = str(mins) + ':' + str(secs)
        sys.stdout.write("\r                                                                                   ")
        sys.stdout.flush()
        sys.stdout.write("\rCountdown: " + str(timeformat))
        sys.stdout.flush()
        sleep(1)
        t -= 1
    print('Resuming Pan Submissions.\n')