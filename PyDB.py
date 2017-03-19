#!/usr/bin/env python
#This is a package that is designed to read PDB file format.
__author__ = "Csongor Mátyás"

import sys, warnings

##############################################################################

class PDBfile: #PDF file class, each PDB file will be a PDB object
    HEADER = None   #Mandatory
    TITLE  = None   #Mandatory
    COMPND = None   #Mandatory
    SOURCE = None   #Mandatory
    KEYWDS = None   #Mandatory
    EXPDTA = None   #Mandatory
    AUTHOR = None   #Mandatory
    REVDAT = None   #Mandatory
    REMARK_2 = None #Mandatory
    REMARK_3 = None #Mandatory
    SEQRES = None   #Mandatory if ATOM records exist.
    CRYST1 = None   #Mandatory
    ORIGX1 = None   #Mandatory
    ORIGX2 = None   #Mandatory
    ORIGX3 = None   #Mandatory
    SCALE1 = None   #Mandatory
    SCALE2 = None   #Mandatory
    SCALE3 = None   #Mandatory
    MASTER = None   #Mandatory
    END    = None   #Mandatory

##############################################################################

def readPDBfile(ID=None, Extension=None, Filename=None, Verbose=False,
                IDUpper=False, IDLower=False, ExtUpper=False, ExtLower=False):
    #This code opens the file that has a given PDB ID . Extension, or Filename
    #ID or Extension can be changed to Upper or Lower case by changing IDUpper,
    #IDLower, ExtUpper, ExtLower flags to True
    #Verbose=True will print some more data

    if Filename == "":
        print("Empty string was given as filename.")
        sys.exit(0) #Check if string is empty
        
    if ID == "":
        print("Empty string was given as PDB ID")
        sys.exit(0) #Check if string is empty

    if IDUpper == True and IDLower == True:
        print("Conflict between IDUpper and IDLower")
        sys.exit(0) #Check if both flags are raised
    
    if ExtUpper == True and ExtLower == True:
        print("Conflict between ExtUpper and ExtLower")
        sys.exit(0) #Check if both flags are raised

    if IDUpper == True:
        ID = str(ID).upper() #Turn ID to uppercase
        
    if IDLower == True:
        ID = str(ID).lower() #Turn ID to lowercase
        
    if ExtUpper == True:
        Extension = str(Extension).upper()  #Turn Extension to uppercase
        
    if ExtLower == True:
        Extension = str(Extension).lower()  #Turn Extension to lowercase

    if Filename == None:    #Check if filename was given
        if ID == None:      #Check if ID was given
            print("There was no PDB ID or Filename given.")
            sys.exit(0)
        else:
            if Extension == None:   #Check if Extension was given, could be empty
                warnings.warn("Extension was not given; .pdb as default will be used.")
                #Warn if Extension was not given but still use default as .pdb
                ToOpen = str(ID) + ".pdb"               #Create filename
            else:
                ToOpen = str(ID) + "." + str(Extension) #Create filename
    else:
        ToOpen = str(Filename)
    
    if Verbose == True:
        print("Opening {} file.".format(ToOpen))

    File = open(ToOpen, "r")
    
    return(File)
    
##############################################################################

def parsePDBfile(File):
    PDB = PDBfile()
    count = 0
    for line in File:
        if line[0:6].strip() == 'HEADER':
            if PDB.HEADER == None:
                PDB.HEADER = line[6:]
            else:
                print("There is more then one HEADER line.")
        
    return(PDB)

##############################################################################

if __name__ == "__main__":
    print("PyDB is a package that is designed to read PDB file format.")
    File = readPDBfile(ID = '1a0k', Extension = 'pdb')
    PDB = parsePDBfile(File)
    print(PDB.HEADER)
