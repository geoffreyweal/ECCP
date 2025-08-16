"""
introduction_message.py, Geoffrey Weal, 17/2/22

introduction_message is a script that give the introduction message for this program
"""

from ECCP import __version__

from math import floor, ceil

def get_version():
	"""
	This method will return the version of this program
	"""
	return str(__version__)

def introduction_message():
	"""
	This is the introdiction message for this program.
	"""
	version = get_version()
	program_string = 'The Electronic Crystal Calculation Prep Program'
	version_string = 'Version: '+str(version)
	no_of_spaces = len(program_string) - len(version_string)
	spaces_LHS = floor(float(no_of_spaces)/2.0)
	spaces_RHS = ceil(float(no_of_spaces)/2.0)
	full_version_string = '|'+'Version: '+str(version)+'|'
	print('.---------------------------------------------------------.')
	print('|                                                         |')
	print('|     The Electronic Crystal Calculation Prep Program     |')
	print('|                                                         |')
	print('|     '+' '*spaces_LHS+version_string+' '*spaces_RHS+'     |')
	print('|                                                         |')
	print('.---------------------------------------------------------.')