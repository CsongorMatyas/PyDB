#!/usr/bin/env python
#This is a package that is designed to read PDB file format.
__author__ = 'Csongor Mátyás'

import sys, warnings

##############################################################################

class PDBfile: #PDF file class, each PDB file will be a PDB object
    HEADER = None   #Mandatory 1x1 - First line of the entry, contains PDB ID code,
                                       #classification, and date of deposition.
    TITLE  = None   #Mandatory 1xN - Description of the experiment represented in the entry.
    COMPND = None   #Mandatory 1xN - Description of macromolecular contents of the entry.
    SOURCE = None   #Mandatory 1xN - Biological source of macromolecules in the entry.
    KEYWDS = None   #Mandatory 1xN - List of keywords describing the macromolecule.
    EXPDTA = None   #Mandatory 1xN - Experimental technique used for the structure determination.
    NUMMDL = None   #Optional  1x1 - Number of models.
    AUTHOR = None   #Mandatory 1xN - List of contributors.
    REVDAT = None   #Mandatory Mx1 - Revision date and related information.
    JRNL   = None   #Mandatory 1xN - Literature citation that defines the coordinate set.
    REMARK_2 = None #Mandatory - 
    REMARK_3 = None #Mandatory - 
    SEQRES = None   #Mandatory MxN - Primary sequence of backbone residues.
    CRYST1 = None   #Mandatory 1x1 - Unit cell parameters, space group, and Z.
    ORIGX1 = None   #Mandatory 1x1 - Transformation from orthogonal coordinates
    ORIGX2 = None   #Mandatory 1x1 - to the submitted coordinates.
    ORIGX3 = None   #Mandatory 1x1
    SCALE1 = None   #Mandatory 1x1 - Transformation from orthogonal coordinates
    SCALE2 = None   #Mandatory 1x1 - to fractional crystallographic coordinates.
    SCALE3 = None   #Mandatory 1x1
    MASTER = None   #Mandatory 1x1 - Control record for bookkeeping.
    END    = None   #Mandatory 1x1 - Last record in the file.

##############################################################################

def readPDBfile(ID=None, Extension=None, Filename=None, Verbose=False,
                IDUpper=False, IDLower=False, ExtUpper=False, ExtLower=False):
    #This code opens the file that has a given PDB ID . Extension, or Filename
    #ID or Extension can be changed to Upper or Lower case by changing IDUpper,
    #IDLower, ExtUpper, ExtLower flags to True
    #Verbose=True will print some more data

    if Filename == '':
        print('Empty string was given as filename.')
        sys.exit(0) #Check if string is empty
        
    if ID == '':
        print('Empty string was given as PDB ID')
        sys.exit(0) #Check if string is empty

    if IDUpper == True and IDLower == True:
        print('Conflict between IDUpper and IDLower')
        sys.exit(0) #Check if both flags are raised
    
    if ExtUpper == True and ExtLower == True:
        print('Conflict between ExtUpper and ExtLower')
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
            print('There was no PDB ID or Filename given.')
            sys.exit(0)
        else:
            if Extension == None:   #Check if Extension was given, could be empty
                warnings.warn('Extension was not given; .pdb as default will be used.')
                #Warn if Extension was not given but still use default as .pdb
                ToOpen = str(ID) + '.pdb'               #Create filename
            else:
                ToOpen = str(ID) + '.' + str(Extension) #Create filename
    else:
        ToOpen = str(Filename)
    
    if Verbose == True:
        print('Opening {} file.'.format(ToOpen))

    File = open(ToOpen, 'r')
    
    return(File)
    
##############################################################################

def parsePDBfile(File):
    PDB = PDBfile()
    
    line_number = 0
    
    for line in File:
        if   len(line) < 7:
            RECORD = line.strip()
        else:
            RECORD = line[0:7].strip()
            
        if   RECORD == '':
            warnings.warn('Empty RECORD entry in line {}!!!'.format(line_number))
        
        elif RECORD == 'HEADER':
            if PDB.HEADER == None:
                PDB.HEADER = (line[7:], line_number)
            else:
                print('There is more then one HEADER line.')
                
        elif RECORD == 'TITLE':
            if PDB.TITLE == None:
                PDB.TITLE  = (line[7:], (line_number, line_number))
            elif PDB.TITLE[1][1] != line_number - 1:
                warnings.warn('Another TITLE in line {}.'.format(line_number))
                PDB.TITLE  = (PDB.TITLE[0] + '\nAnother Title\n' + line[7:],
                              (PDB.TITLE[1][0], line_number))
            else:
                PDB.TITLE  = (PDB.TITLE[0] + line[7:],
                              (PDB.TITLE[1][0], line_number))
        
        elif RECORD == 'COMPND':
            if PDB.COMPND == None:
                PDB.COMPND  = (line[7:], (line_number, line_number))
            elif PDB.COMPND[1][1] != line_number - 1:
                warnings.warn('Another COMPND in line {}.'.format(line_number))
                PDB.COMPND = (PDB.COMPND[0] + '\nAnother COMPND\n' + line[7:],
                              (PDB.COMPND[1][0], line_number))
            else:
                PDB.COMPND  = (PDB.COMPND[0] + line[7:],
                              (PDB.COMPND[1][0], line_number))
                
        elif RECORD == 'SOURCE':
            if PDB.SOURCE == None:
                PDB.SOURCE  = (line[7:], (line_number, line_number))
            elif PDB.SOURCE[1][1] != line_number - 1:
                warnings.warn('Another SOURCE in line {}.'.format(line_number))
                PDB.SOURCE = (PDB.SOURCE[0] + '\nAnother SOURCE\n' + line[7:],
                              (PDB.SOURCE[1][0], line_number))
            else:
                PDB.SOURCE  = (PDB.SOURCE[0] + line[7:],
                              (PDB.SOURCE[1][0], line_number))
        
        elif RECORD == 'KEYWDS':
            if PDB.KEYWDS == None:
                PDB.KEYWDS  = (line[7:], (line_number, line_number))
            elif PDB.KEYWDS[1][1] != line_number - 1:
                warnings.warn('Another KEYWDS in line {}.'.format(line_number))
                PDB.KEYWDS = (PDB.KEYWDS[0] + '\nAnother KEYWDS\n' + line[7:],
                              (PDB.KEYWDS[1][0], line_number))
            else:
                PDB.KEYWDS  = (PDB.KEYWDS[0] + line[7:],
                              (PDB.KEYWDS[1][0], line_number))
        
        elif RECORD == 'EXPDTA':
            if PDB.EXPDTA == None:
                PDB.EXPDTA  = (line[7:], (line_number, line_number))
            elif PDB.EXPDTA[1][1] != line_number - 1:
                warnings.warn('Another EXPDTA in line {}.'.format(line_number))
                PDB.EXPDTA = (PDB.EXPDTA[0] + '\nAnother EXPDTA\n' + line[7:],
                              (PDB.EXPDTA[1][0], line_number))
            else:
                PDB.EXPDTA  = (PDB.EXPDTA[0] + line[7:],
                              (PDB.EXPDTA[1][0], line_number))
        
        elif RECORD == 'AUTHOR':
            if PDB.AUTHOR == None:
                PDB.AUTHOR  = (line[7:], (line_number, line_number))
            elif PDB.AUTHOR[1][1] != line_number - 1:
                warnings.warn('Another AUTHOR in line {}.'.format(line_number))
                PDB.AUTHOR = (PDB.AUTHOR[0] + '\nAnother AUTHOR\n' + line[7:],
                              (PDB.AUTHOR[1][0], line_number))
            else:
                PDB.AUTHOR  = (PDB.AUTHOR[0] + line[7:],
                              (PDB.AUTHOR[1][0], line_number))
        
        elif RECORD == 'REVDAT':
            if PDB.REVDAT == None:
                PDB.REVDAT  = (line[7:], (line_number, line_number))
            elif PDB.REVDAT[1][1] != line_number - 1:
                warnings.warn('Another REVDAT in line {}.'.format(line_number))
                PDB.REVDAT = (PDB.REVDAT[0] + '\nAnother REVDAT\n' + line[7:],
                              (PDB.REVDAT[1][0], line_number))
            else:
                PDB.REVDAT = (PDB.REVDAT[0] + line[7:],
                              (PDB.REVDAT[1][0], line_number))
        
        elif RECORD == 'JRNL':
            if PDB.JRNL == None:
                PDB.JRNL  = (line[7:], (line_number, line_number))
            elif PDB.JRNL[1][1] != line_number - 1:
                warnings.warn('Another JRNL in line {}.'.format(line_number))
                PDB.JRNL = (PDB.JRNL[0] + '\nAnother JRNL\n' + line[7:],
                              (PDB.JRNL[1][0], line_number))
            else:
                PDB.JRNL = (PDB.JRNL[0] + line[7:],
                              (PDB.JRNL[1][0], line_number))
       
        
        line_number += 1
    return(PDB)

##############################################################################

if __name__ == '__main__':
    print('PyDB is a package that is designed to read PDB file format.')
    File = readPDBfile(ID = '1a0k', Extension = 'pdb')
    PDB = parsePDBfile(File)
    print(PDB.HEADER)
    print(PDB.TITLE )
    print(PDB.COMPND)
    print(PDB.SOURCE)
    print(PDB.KEYWDS)
    print(PDB.EXPDTA)
    print(PDB.AUTHOR)
    print(PDB.REVDAT)
    print(PDB.JRNL  )
