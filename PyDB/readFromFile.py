import sys, warnings

def readFromFile(Filename=None):
    #This code opens the file that is given

    if Filename == '':
        warnings.warn('Empty string was given as filename!')
        sys.exit(0) #Check if string is empty
        
    with open(Filename, 'r') as f:
        FILE = f.readlines()

    return(FILE)