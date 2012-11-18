# Delete all files of type 'file_type', eg '.txt'

# import os
# import glob
import subprocess

def clean(file_type):
    files_rm = 'rm *%s'% file_type
    subprocess.call(files_rm, shell = True)
    # files = '*%s'% file_type
    # filenames = glob.glob(files)
    # print filenames
    # for fle in filenames:
    #     try:
    #         os.remove(fle)
    #     except OSError:
    #         pass
