# Delete all files of type 'file_type', eg '.txt'

import os
import glob

def clean(file_type):
    files = '*%s'% file_type
    filenames = glob.glob(files)
    for fle in filenames:
        try:
            os.remove(fle)
        except OSError:
            pass