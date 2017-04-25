import sys, warnings, urllib.request

def readPDB(ID=None):
    #This code opens the file that has a given PDB ID using internet connection.
    if ID == '':
        warnings.warn('Empty string was given as PDB ID!')
        sys.exit(0)   #Check if string is empty
    else:
        req = urllib.request.Request('https://files.rcsb.org/download/' + 
                                     str(ID).upper() + '.pdb')
        try:
            with urllib.request.urlopen(req) as response:
                html = response.read()
        except urllib.error.HTTPError as e:
            warnings.warn('\n{}\n{}'.format(e.code, e.read()))
            sys.exit(0)
    
    UTF = html.decode('utf-8')
    FILE = str(UTF)[:-1].split('\n')
    return(FILE)