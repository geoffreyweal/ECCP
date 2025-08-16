'''
Geoffrey Weal, processing_OPV_Dimer_data.py, 9/3/22

This script contains methods for processing_OPV_Dimer_data.py

'''
import os

# -----------------------------------------------------------------

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

# -----------------------------------------------------------------

def reverse_readline(filename, buf_size=8192):
    """
    A generator that returns the lines of a file in reverse order

    Parameters
    ----------
    filename : str
        This is the path to the file you want to read.
    buf_size : int
        This is the buffer size to read in.

    Returns
    -------
    Returns each line in the file in reverse order.
    """
    
    with open(filename) as fh:
        segment = None
        offset = 0
        fh.seek(0, os.SEEK_END)
        file_size = remaining_size = fh.tell()
        while remaining_size > 0:
            offset = min(file_size, offset + buf_size)
            fh.seek(file_size - offset)
            buffer = fh.read(min(remaining_size, buf_size))
            remaining_size -= buf_size
            lines = buffer.split('\n')
            # The first line of the buffer is probably not a complete line so
            # we'll save it and append it to the last line of the next buffer
            # we read
            if segment is not None:
                # If the previous chunk starts right from the beginning of line
                # do not concat the segment to the last line of new chunk.
                # Instead, yield the segment first 
                if buffer[-1] != '\n':
                    lines[-1] += segment
                else:
                    yield segment
            segment = lines[0]
            for index in range(len(lines) - 1, 0, -1):
                if lines[index]:
                    yield lines[index]
        # Don't yield None if the file was empty
        if segment is not None:
            yield segment

# -----------------------------------------------------------------




