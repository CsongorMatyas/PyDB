#!/usr/bin/env python
#This is a package that is designed to read PDB file format.
__author__ = 'Csongor Mátyás'

import sys, warnings

##############################################################################

class PDBfile: #PDF file class, each PDB file will be a PDB object
    HEADER = None   #Mandatory 1x1 - First line of the entry, contains PDB ID code,
                                       #classification, and date of deposition.
    OBSLTE = None   #Optional  1xN - Statement that the entry has been removed from
                    #distribution and list of the ID code(s) which replaced it.
    TITLE  = None   #Mandatory 1xN - Description of the experiment represented in the entry.
    SPLIT  = None   #Optional  1xN - List of PDB entries that compose a larger
                    #macromolecular complexes.
    CAVEAT = None   #Optional  1xN - Severe error indicator.    
    COMPND = None   #Mandatory 1xN - Description of macromolecular contents of the entry.
    SOURCE = None   #Mandatory 1xN - Biological source of macromolecules in the entry.
    KEYWDS = None   #Mandatory 1xN - List of keywords describing the macromolecule.
    EXPDTA = None   #Mandatory 1xN - Experimental technique used for the structure determination.
    NUMMDL = None   #Optional  1x1 - Number of models.
    MDLTYP = None   #Optional  1xN - Contains additional annotation pertinent
                    #to the coordinates presented in the entry.
    AUTHOR = None   #Mandatory 1xN - List of contributors.
    REVDAT = None   #Mandatory Mx1 - Revision date and related information.
    SPRSDE = None   #Optional  1xN - List of entries obsoleted from public
                    #release and replaced by current entry.
    JRNL   = None   #Mandatory 1xN - Literature citation that defines the coordinate set.
    REMARK_0 = None #Optional  1xN - Re-refinement notice
    REMARK_1 = None #Optional  1xN - Related publications
    REMARK_2 = None #Mandatory 1xN - Resolution
    REMARK_3 = None #Mandatory 1xN - Final refinement information
    REMARK_4 = None #Optional  1xN - Format
    REMARK_5 = None #Optional  1xN - Obsolete Statement
    REMARK_6 = None #Optional  1xN - 6-99 free text annotation
    REMARK_100=None #Optional  1xN - Deposition or Processing Site
    REMARK_200=None #Optional  1xN - X-ray Diffraction Experimental Details
    REMARK_205=None #Optional  1xN - Fiber Diffraction, Fiber Sample Experiment Details
    REMARK_210=None #Optional  1xN - +215/217  NMR Experiment Details
    REMARK_230=None #Optional  1xN - Neutron Diffraction Experiment Details
    REMARK_240=None #Optional  1xN - Electron Crystallography Experiment Details
    REMARK_245=None #Optional  1xN - Electron Microscopy Experiment Details
    REMARK_247=None #Optional  1xN - Electron Microscopy details
    REMARK_250=None #Optional  1xN - Other Type of Experiment Details
    REMARK_265=None #Optional  1xN - Solution Scattering Experiment Details
    REMARK_280=None #Optional  1xN - Crystal
    REMARK_285=None #Optional  1xN - CRYST1
    REMARK_290=None #Optional  1xN - Crystallographic Symmetry
    REMARK_300=None #Optional  1xN - Biomolecule
    REMARK_350=None #Optional  1xN - Generating the Biomolecule
    REMARK_375=None #Optional  1xN - Special Position
    REMARK_400=None #Optional  1xN - Compound
    REMARK_450=None #Optional  1xN - Source
    REMARK_465=None #Optional  1xN - Missing residues
    REMARK_470=None #Optional  1xN - Missing Atom(s)
    REMARK_475=None #Optional  1xN - Residues modeled with zero occupancy
    REMARK_480=None #Optional  1xN - Polymer atoms modeled with zero occupancy
    REMARK_500=None #Optional  1xN - Geometry and Stereochemistry
    REMARK_525=None #Optional  1xN - Distant Solvent Atoms
    REMARK_600=None #Optional  1xN - Heterogen
    REMARK_610=None #Optional  1xN - Non-polymer residues with missing atoms
    REMARK_615=None #Optional  1xN - Non-polymer residues containing atoms with zero occupancy
    REMARK_620=None #Optional  1xN - Metal coordination
    REMARK_630=None #Optional  1xN - Inhibitor Description
    REMARK_650=None #Optional  1xN - Helix
    REMARK_700=None #Optional  1xN - Sheet
    REMARK_800=None #Optional  1xN - Important Sites
    REMARK_999=None #Optional  1xN - Sequence
    DBREF  = None   #Optional  Mx1 - Reference to the entry in the sequence database(s).
    SEQADV = None   #Optional  Mx1 - Identification of conflicts between PDB
                    #and the named sequence database.
    SEQRES = None   #Mandatory MxN - Primary sequence of backbone residues.
    MODRES = None   #Optional  Mx1 - Identification of modifications to standard residues.
    HET    = None   #Optional  Mx1 - Identification of non-standard groups (heterogens).
    HETNAM = None   #Optional  MxN - Compound name of the heterogens. (JUPAC?)
    HETSYN = None   #Optional  MxN - Synonymous compound names for heterogens. (nonJUPAC?)
    FORMUL = None   #Mandatory MxN - Chemical formula of non-standard groups.
    HELIX  = None   #Optional  Mx1 - Identification of helical substructures.
    SHEET  = None   #Optional  Mx1 - Identification of sheet substructures.
    SSBOND = None   #Optional  Mx1 - Identification of disulfide bonds.
    LINK   = None   #Optional  Mx1 - Identification of inter-residue bonds.
    CISPEP = None   #Optional  Mx1 - Identification of peptide residues in cis conformation.
    SITE   = None   #Optional  MxN - Identification of groups comprising important entity sites.
    CRYST1 = None   #Mandatory 1x1 - Unit cell parameters, space group, and Z.
    ORIGX1 = None   #Mandatory 1x1 - Transformation from orthogonal coordinates
    ORIGX2 = None   #Mandatory 1x1 - to the submitted coordinates.
    ORIGX3 = None   #Mandatory 1x1
    SCALE1 = None   #Mandatory 1x1 - Transformation from orthogonal coordinates
    SCALE2 = None   #Mandatory 1x1 - to fractional crystallographic coordinates.
    SCALE3 = None   #Mandatory 1x1
    MTRIX1 = None   #Optional  MxN - Transformations expressing non-crystallographic
    MTRIX2 = None   #Optional  MxN   symmetry.
    MTRIX3 = None   #Optional  MxN
    MODEL  = None   #Optional ?MxN - Specification of model number for multiple
                    #structures in a single coordinate entry.
    ATOM   = None   #Mandatory Mx1 - Atomic coordinate records for standard groups.
    ANISOU = None   #Optional  Mx1 - Anisotropic temperature factors.
    TER    = None   #Mandatory Mx1 - Chain terminator.
    HETATM = None   #Optional  Mx1 - Atomic coordinate records for heterogens.
    ENDMDL = None   #Optional ?MxN - End-of-model record for multiple structures
                    #in a single coordinate entry.
    CONECT = None   #Optional  Mx1 - Connectivity records.
    MASTER = None   #Mandatory 1x1 - Control record for bookkeeping.
    END    = None   #Mandatory 1x1 - Last record in the file.
    JUNK   = None   #PyDB property, if it's not None, there were lines starting
                    #with a keyword other than the ones allowed by PDB file format.

##############################################################################

class Helix:    #Each HELIX within a PDB object will be a Helix object
    def __init__(self, TEXT):
        
        # (#) Serial number of the helix. This starts at 1 and increases incrementally.
        if TEXT[7:10].strip() == '':
            self.serNum = None
        else:
            self.serNum = int(TEXT[7:10].strip())
        
        # ('H'#) Helix identifier. In addition to a serial number, each helix is given an alphanumeric character helix identifier.
        if TEXT[11:14].strip() == '':
            self.helixID = None
        else:
            self.helixID = TEXT[11:14].strip()
        
        # ('') Name of the initial residue. 
        if TEXT[15:18].strip() == '':
            self.initResName = None
        else:
            self.initResName = TEXT[15:18].strip()
            
        # (A-Z) Chain identifier for the chain containing this helix.
        if TEXT[19:20].strip() == '':
            self.initChainID = None
        else:
            self.initChainID = TEXT[19:20].strip()
        
        # (#) Sequence number of the initial residue. (SEQRES)
        if TEXT[21:25].strip() == '':
            self.initSeqNum = None
        else:
            self.initSeqNum = int(TEXT[21:25].strip())
        
        # Insertion code of the initial residue.
        if TEXT[25:26].strip() == '':
            self.initICode = None
        else:
            self.initICode = TEXT[25:26].strip()
            
        # ('') Name of the terminal residue of the helix.
        if TEXT[27:30].strip() == '':
            self.endResName = None
        else:
            self.endResName = TEXT[27:30].strip()
            
        # (A-Z) Chain identifier for the chain containing this helix. (WHY? Must be the same chain as initial, since the helix cannot be broken)
        if TEXT[31:32].strip() == '':
            self.endChainID = None
        else:
            self.endChainID = TEXT[31:32].strip()
            
        # (#) Sequence number of the terminal residue.
        if TEXT[33:37].strip() == '':
            self.endSeqNum = None
        else:
            self.endSeqNum = int(TEXT[33:37].strip())
        
        # Insertion code of the terminal residue.
        if TEXT[37:38].strip() == '':
            self.endICode = None
        else:
            self.endICode = TEXT[37:38].strip()
            
        # (#) Helix Class (1-10)
        if TEXT[38:40].strip() == '':
            self.helixClass = None
        else:
            self.helixClass = int(TEXT[38:40].strip())
            
        # ('') Comment about this helix.
        if TEXT[40:70].strip() == '':
            self.comment = None
        else:
            self.comment = TEXT[40:70].strip()
        
        # (#) Length of this helix.
        if TEXT[71:76].strip() == '':
            self.length = None
        else:
            self.length = int(TEXT[71:76].strip())

##############################################################################

class Sheet:    #Each SHEET within a PDB object will be a Sheet object
    def __init__(self, TEXT):
        
        # (#) Strand number which starts at 1 for each strand within a sheet and increases by one.
        if TEXT[7:10].strip() == '':
            self.strand = None
        else:
            self.strand = int(TEXT[7:10].strip())
        
        # (A-Z) Sheet identifier.
        if TEXT[11:14].strip() == '':
            self.sheetID = None
        else:
            self.sheetID = TEXT[11:14].strip()
        
        # (#) Number of strands in sheet. 
        if TEXT[14:16].strip() == '':
            self.numStrands = None
        else:
            self.numStrands = int(TEXT[14:16].strip())
            
        # ('') Residue name of initial residue.
        if TEXT[17:20].strip() == '':
            self.initResName = None
        else:
            self.initResName = TEXT[17:20].strip()
        
        # (A-Z) Chain identifier of initial residue in strand.
        if TEXT[21:22].strip() == '':
            self.initChainID = None
        else:
            self.initChainID = TEXT[21:22].strip()
        
        # (#) Sequence number of initial residue in strand.
        if TEXT[22:26].strip() == '':
            self.initSeqNum = None
        else:
            self.initSeqNum = int(TEXT[22:26].strip())
        
        # Insertion code of initial residue in strand.
        if TEXT[26:27].strip() == '':
            self.initICode = None
        else:
            self.initICode = TEXT[26:27].strip()
            
        # ('') Residue name of terminal residue.
        if TEXT[28:32].strip() == '':
            self.endResName = None
        else:
            self.endResName = TEXT[28:31].strip()
            
        # (A-Z) Chain identifier of terminal residue.
        if TEXT[32:33].strip() == '':
            self.endChainID = None
        else:
            self.endChainID = TEXT[32:33].strip()
        
        # (#) Sequence number of terminal residue.
        if TEXT[33:37].strip() == '':
            self.endSeqNum = None
        else:
            self.endSeqNum = int(TEXT[33:37].strip())
        
        # Insertion code of the terminal residue.
        if TEXT[37:38].strip() == '':
            self.endICode = None
        else:
            self.endICode = TEXT[37:38].strip()
            
        # (#) Sense of strand with respect to previous strand in the sheet. 
        # 0 if first strand, 1 if parallel, and -1 if anti-parallel.
        if TEXT[38:40].strip() == '':
            self.sense = None
        else:
            self.sense = int(TEXT[38:40].strip())
        
        # ('') Registration. Atom name in current strand.
        if TEXT[41:45].strip() == '':
            self.curAtom = None
        else:
            self.curAtom = TEXT[41:45].strip()
        
        # ('') Registration. Residue name in current strand
        if TEXT[45:48].strip() == '':
            self.curResName = None
        else:
            self.curResName = TEXT[45:48].strip()
            
        # (A-Z) Registration. Chain identifier in current strand.
        if TEXT[49:50].strip() == '':
            self.curChainId = None
        else:
            self.curChainId = TEXT[49:50].strip()
        
        # (#) Registration. Residue sequence number in current strand.
        if TEXT[50:54].strip() == '':
            self.curSeqNum = None       #curSeqNum better than curSeqRes
        else:
            self.curSeqNum = int(TEXT[50:54].strip())
        
        # Registration. Insertion code in current strand.
        if TEXT[54:55].strip() == '':
            self.curICode = None
        else:
            self.curICode = TEXT[54:55].strip()
        
        # ('') Registration. Atom name in previous strand.
        if TEXT[56:60].strip() == '':
            self.prevAtom = None
        else:
            self.prevAtom = TEXT[56:60].strip()
        
        # ('') Registration. Residue name in previous strand.
        if TEXT[60:63].strip() == '':
            self.prevResName = None
        else:
            self.prevResName = TEXT[60:63].strip()
            
        # (A-Z) Registration. Chain identifier in previous strand.
        if TEXT[64:65].strip() == '':
            self.prevChainId = None
        else:
            self.prevChainId = TEXT[64:65].strip()
        
        # (#) Registration. Residue sequence number in previous strand.
        if TEXT[65:69].strip() == '':
            self.prevSeqNum = None       #curSeqNum better than curSeqRes
        else:
            self.prevSeqNum = int(TEXT[65:69].strip())
        
        # Registration. Insertion code in previous strand.
        if TEXT[69:70].strip() == '':
            self.prevICode = None
        else:
            self.prevICode = TEXT[69:70].strip()

##############################################################################

class Atom:    #Each ATOM within a PDB object will be a Atom object
    def __init__(self, TEXT):
        
        # (#) Atom serial number.
        if TEXT[6:11].strip() == '':
            self.serNum = None
        else:
            self.serNum = int(TEXT[6:11].strip())
        
        # ('') Atom name.
        if TEXT[12:16].strip() == '':
            self.name = None
        else:
            self.name = TEXT[12:16].strip()
        
        # ('') Alternate location indicator. 
        if TEXT[16:17].strip() == '':
            self.altLoc = None
        else:
            self.altLoc = TEXT[16:17].strip()
            
        # ('') Residue name.
        if TEXT[17:20].strip() == '':
            self.resName = None
        else:
            self.resName = TEXT[17:20].strip()
        
        # Chain identifier.
        if TEXT[21:22].strip() == '':
            self.chainID = None
        else:
            self.chainID = TEXT[21:22].strip()
        
        # Residue sequence number.
        if TEXT[22:26].strip() == '':
            self.resSeqNum = None
        else:
            self.resSeqNum = int(TEXT[22:26].strip())
            
        # Code for insertion of residues.
        if TEXT[26:27].strip() == '':
            self.iCode = None
        else:
            self.iCode = TEXT[26:27].strip()
            
        # (#) Orthogonal coordinates for X in Angstroms.
        if TEXT[30:38].strip() == '':
            self.x = None
        else:
            self.x = float(TEXT[30:38].strip())
            
        # (#) Orthogonal coordinates for Y in Angstroms.
        if TEXT[38:46].strip() == '':
            self.y = None
        else:
            self.y = float(TEXT[38:46].strip())
        
        # (#) Orthogonal coordinates for Z in Angstroms.
        if TEXT[46:54].strip() == '':
            self.z = None
        else:
            self.z = float(TEXT[46:54].strip())
            
        # (#) Occupancy.
        if TEXT[54:60].strip() == '':
            self.occupancy = None
        else:
            self.occupancy = float(TEXT[54:60].strip())
            
        # (#) Temperature factor.
        if TEXT[60:66].strip() == '':
            self.tempFactor = None
        else:
            self.tempFactor = float(TEXT[60:66].strip())
        
        # (#) Element symbol, right-justified.
        if TEXT[76:78].strip() == '':
            self.element = None
        else:
            self.element = TEXT[76:78].strip()
        
        # (#) Charge on the atom.
        if TEXT[78:80].strip() == '':
            self.charge = None
        else:
            self.charge = float(TEXT[78:80].strip())

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
            warnings.warn('Empty RECORD entry at line {}!!!'.format(line_number))
        
        elif RECORD == 'HEADER':#1x1
            TEXT = line[7:-1].strip()
            if PDB.HEADER == None:
                PDB.HEADER = [TEXT, line_number]
            else:
                warnings.warn('Another HEADER section found at line {}!!!'.format(line_number))
                sys.exit(0)
        
        elif RECORD == 'OBSLTE':#1xN
            TEXT = line[7:-1].strip()
            if PDB.OBSLTE == None:
                PDB.OBSLTE = [[[],[]]]
                PDB.OBSLTE[0][0].append(TEXT)
                PDB.OBSLTE[0][1].append(line_number)
                PDB.OBSLTE[0][1].append(line_number)
            else:
                if PDB.OBSLTE[-1][1][1] == line_number - 1:
                    PDB.OBSLTE[-1][0].append(TEXT)
                    PDB.OBSLTE[-1][1][1] = line_number
                else:
                    warnings.warn('Another OBSLTE section found at line {}.'.format(line_number))
                    PDB.OBSLTE.append([[],[]])
                    PDB.OBSLTE[-1][0].append(TEXT)
                    PDB.OBSLTE[-1][1].append(line_number)
                    PDB.OBSLTE[-1][1].append(line_number)
        
        elif RECORD == 'TITLE':#1xN
            TEXT = line[7:-1].strip()
            if PDB.TITLE == None:
                PDB.TITLE = [[[],[]]]
                PDB.TITLE[0][0].append(TEXT)
                PDB.TITLE[0][1].append(line_number)
                PDB.TITLE[0][1].append(line_number)
            else:
                if PDB.TITLE[-1][1][1] == line_number - 1:
                    PDB.TITLE[-1][0].append(TEXT)
                    PDB.TITLE[-1][1][1] = line_number
                else:
                    warnings.warn('Another TITLE section found at line {}.'.format(line_number))
                    PDB.TITLE.append([[],[]])
                    PDB.TITLE[-1][0].append(TEXT)
                    PDB.TITLE[-1][1].append(line_number)
                    PDB.TITLE[-1][1].append(line_number)
        
        elif RECORD == 'SPLIT':#1xN
            TEXT = line[7:-1].strip()
            if PDB.SPLIT == None:
                PDB.SPLIT = [[[],[]]]
                PDB.SPLIT[0][0].append(TEXT)
                PDB.SPLIT[0][1].append(line_number)
                PDB.SPLIT[0][1].append(line_number)
            else:
                if PDB.SPLIT[-1][1][1] == line_number - 1:
                    PDB.SPLIT[-1][0].append(TEXT)
                    PDB.SPLIT[-1][1][1] = line_number
                else:
                    warnings.warn('Another SPLIT section found at line {}.'.format(line_number))
                    PDB.SPLIT.append([[],[]])
                    PDB.SPLIT[-1][0].append(TEXT)
                    PDB.SPLIT[-1][1].append(line_number)
                    PDB.SPLIT[-1][1].append(line_number)
        
        elif RECORD == 'CAVEAT':#1xN
            TEXT = line[7:-1].strip()
            if PDB.CAVEAT == None:
                PDB.CAVEAT = [[[],[]]]
                PDB.CAVEAT[0][0].append(TEXT)
                PDB.CAVEAT[0][1].append(line_number)
                PDB.CAVEAT[0][1].append(line_number)
            else:
                if PDB.CAVEAT[-1][1][1] == line_number - 1:
                    PDB.CAVEAT[-1][0].append(TEXT)
                    PDB.CAVEAT[-1][1][1] = line_number
                else:
                    warnings.warn('Another CAVEAT section found at line {}.'.format(line_number))
                    PDB.CAVEAT.append([[],[]])
                    PDB.CAVEAT[-1][0].append(TEXT)
                    PDB.CAVEAT[-1][1].append(line_number)
                    PDB.CAVEAT[-1][1].append(line_number)
        
        elif RECORD == 'COMPND':#1xN
            TEXT = line[7:-1].strip()
            if PDB.COMPND == None:
                PDB.COMPND = [[[],[]]]
                PDB.COMPND[0][0].append(TEXT)
                PDB.COMPND[0][1].append(line_number)
                PDB.COMPND[0][1].append(line_number)
            else:
                if PDB.COMPND[-1][1][1] == line_number - 1:
                    PDB.COMPND[-1][0].append(TEXT)
                    PDB.COMPND[-1][1][1] = line_number
                else:
                    warnings.warn('Another COMPND section found at line {}.'.format(line_number))
                    PDB.COMPND.append([[],[]])
                    PDB.COMPND[-1][0].append(TEXT)
                    PDB.COMPND[-1][1].append(line_number)
                    PDB.COMPND[-1][1].append(line_number)
        
        elif RECORD == 'SOURCE':#1xN
            TEXT = line[7:-1].strip()
            if PDB.SOURCE == None:
                PDB.SOURCE = [[[],[]]]
                PDB.SOURCE[0][0].append(TEXT)
                PDB.SOURCE[0][1].append(line_number)
                PDB.SOURCE[0][1].append(line_number)
            else:
                if PDB.SOURCE[-1][1][1] == line_number - 1:
                    PDB.SOURCE[-1][0].append(TEXT)
                    PDB.SOURCE[-1][1][1] = line_number
                else:
                    warnings.warn('Another SOURCE section found at line {}.'.format(line_number))
                    PDB.SOURCE.append([[],[]])
                    PDB.SOURCE[-1][0].append(TEXT)
                    PDB.SOURCE[-1][1].append(line_number)
                    PDB.SOURCE[-1][1].append(line_number)
        
        elif RECORD == 'KEYWDS':#1xN
            TEXT = line[7:-1].strip()
            if PDB.KEYWDS == None:
                PDB.KEYWDS = [[[],[]]]
                PDB.KEYWDS[0][0].append(TEXT)
                PDB.KEYWDS[0][1].append(line_number)
                PDB.KEYWDS[0][1].append(line_number)
            else:
                if PDB.KEYWDS[-1][1][1] == line_number - 1:
                    PDB.KEYWDS[-1][0].append(TEXT)
                    PDB.KEYWDS[-1][1][1] = line_number
                else:
                    warnings.warn('Another KEYWDS section found at line {}.'.format(line_number))
                    PDB.KEYWDS.append([[],[]])
                    PDB.KEYWDS[-1][0].append(TEXT)
                    PDB.KEYWDS[-1][1].append(line_number)
                    PDB.KEYWDS[-1][1].append(line_number)
        
        elif RECORD == 'EXPDTA':#1xN
            TEXT = line[7:-1].strip()
            if PDB.EXPDTA == None:
                PDB.EXPDTA = [[[],[]]]
                PDB.EXPDTA[0][0].append(TEXT)
                PDB.EXPDTA[0][1].append(line_number)
                PDB.EXPDTA[0][1].append(line_number)
            else:
                if PDB.EXPDTA[-1][1][1] == line_number - 1:
                    PDB.EXPDTA[-1][0].append(TEXT)
                    PDB.EXPDTA[-1][1][1] = line_number
                else:
                    warnings.warn('Another EXPDTA section found at line {}.'.format(line_number))
                    PDB.EXPDTA.append([[],[]])
                    PDB.EXPDTA[-1][0].append(TEXT)
                    PDB.EXPDTA[-1][1].append(line_number)
                    PDB.EXPDTA[-1][1].append(line_number)
        
        elif RECORD == 'NUMMDL':#1xN
            TEXT = line[7:-1].strip()
            if PDB.NUMMDL == None:
                PDB.NUMMDL = [[[],[]]]
                PDB.NUMMDL[0][0].append(TEXT)
                PDB.NUMMDL[0][1].append(line_number)
                PDB.NUMMDL[0][1].append(line_number)
            else:
                if PDB.NUMMDL[-1][1][1] == line_number - 1:
                    PDB.NUMMDL[-1][0].append(TEXT)
                    PDB.NUMMDL[-1][1][1] = line_number
                else:
                    warnings.warn('Another NUMMDL section found at line {}.'.format(line_number))
                    PDB.NUMMDL.append([[],[]])
                    PDB.NUMMDL[-1][0].append(TEXT)
                    PDB.NUMMDL[-1][1].append(line_number)
                    PDB.NUMMDL[-1][1].append(line_number)
        
        elif RECORD == 'MDLTYP':#1xN
            TEXT = line[7:-1].strip()
            if PDB.MDLTYP == None:
                PDB.MDLTYP = [[[],[]]]
                PDB.MDLTYP[0][0].append(TEXT)
                PDB.MDLTYP[0][1].append(line_number)
                PDB.MDLTYP[0][1].append(line_number)
            else:
                if PDB.MDLTYP[-1][1][1] == line_number - 1:
                    PDB.MDLTYP[-1][0].append(TEXT)
                    PDB.MDLTYP[-1][1][1] = line_number
                else:
                    warnings.warn('Another MDLTYP section found at line {}.'.format(line_number))
                    PDB.MDLTYP.append([[],[]])
                    PDB.MDLTYP[-1][0].append(TEXT)
                    PDB.MDLTYP[-1][1].append(line_number)
                    PDB.MDLTYP[-1][1].append(line_number)
        
        elif RECORD == 'AUTHOR':#1xN
            TEXT = line[7:-1].strip()
            if PDB.AUTHOR == None:
                PDB.AUTHOR = [[[],[]]]
                PDB.AUTHOR[0][0].append(TEXT)
                PDB.AUTHOR[0][1].append(line_number)
                PDB.AUTHOR[0][1].append(line_number)
            else:
                if PDB.AUTHOR[-1][1][1] == line_number - 1:
                    PDB.AUTHOR[-1][0].append(TEXT)
                    PDB.AUTHOR[-1][1][1] = line_number
                else:
                    warnings.warn('Another AUTHOR section found at line {}.'.format(line_number))
                    PDB.AUTHOR.append([[],[]])
                    PDB.AUTHOR[-1][0].append(TEXT)
                    PDB.AUTHOR[-1][1].append(line_number)
                    PDB.AUTHOR[-1][1].append(line_number)
        
        elif RECORD == 'REVDAT':#Mx1
            TEXT = line[7:-1].strip()
            if PDB.REVDAT == None:
                PDB.REVDAT = [[],[]]
                PDB.REVDAT[0].append(TEXT)
                PDB.REVDAT[1].append(line_number)
            else:
                PDB.REVDAT[0].append(TEXT)
                PDB.REVDAT[1].append(line_number)
                
        elif RECORD == 'SPRSDE':#1xN
            TEXT = line[7:-1].strip()
            if PDB.SPRSDE == None:
                PDB.SPRSDE = [[[],[]]]
                PDB.SPRSDE[0][0].append(TEXT)
                PDB.SPRSDE[0][1].append(line_number)
                PDB.SPRSDE[0][1].append(line_number)
            else:
                if PDB.SPRSDE[-1][1][1] == line_number - 1:
                    PDB.SPRSDE[-1][0].append(TEXT)
                    PDB.SPRSDE[-1][1][1] = line_number
                else:
                    warnings.warn('Another SPRSDE section found at line {}.'.format(line_number))
                    PDB.SPRSDE.append([[],[]])
                    PDB.SPRSDE[-1][0].append(TEXT)
                    PDB.SPRSDE[-1][1].append(line_number)
                    PDB.SPRSDE[-1][1].append(line_number)
        
        elif RECORD == 'JRNL':#1xN
            TEXT = line[7:-1].strip()
            if PDB.JRNL == None:
                PDB.JRNL = [[[],[]]]
                PDB.JRNL[0][0].append(TEXT)
                PDB.JRNL[0][1].append(line_number)
                PDB.JRNL[0][1].append(line_number)
            else:
                if PDB.JRNL[-1][1][1] == line_number - 1:
                    PDB.JRNL[-1][0].append(TEXT)
                    PDB.JRNL[-1][1][1] = line_number
                else:
                    warnings.warn('Another JRNL section found at line {}.'.format(line_number))
                    PDB.JRNL.append([[],[]])
                    PDB.JRNL[-1][0].append(TEXT)
                    PDB.JRNL[-1][1].append(line_number)
                    PDB.JRNL[-1][1].append(line_number)
        
        elif RECORD == 'REMARK':
            RNUM = int(line[7:10].strip())
            if RNUM == 0:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_0 == None:
                    PDB.REMARK_0 = [[[],[]]]
                    PDB.REMARK_0[0][0].append(TEXT)
                    PDB.REMARK_0[0][1].append(line_number)
                    PDB.REMARK_0[0][1].append(line_number)
                else:
                    if PDB.REMARK_0[-1][1][1] == line_number - 1:
                        PDB.REMARK_0[-1][0].append(TEXT)
                        PDB.REMARK_0[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_0 section found at line {}.'.format(line_number))
                        PDB.REMARK_0.append([[],[]])
                        PDB.REMARK_0[-1][0].append(TEXT)
                        PDB.REMARK_0[-1][1].append(line_number)
                        PDB.REMARK_0[-1][1].append(line_number)
                        
            elif RNUM == 1:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_1 == None:
                    PDB.REMARK_1 = [[[],[]]]
                    PDB.REMARK_1[0][0].append(TEXT)
                    PDB.REMARK_1[0][1].append(line_number)
                    PDB.REMARK_1[0][1].append(line_number)
                else:
                    if PDB.REMARK_1[-1][1][1] == line_number - 1:
                        PDB.REMARK_1[-1][0].append(TEXT)
                        PDB.REMARK_1[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_1 section found at line {}.'.format(line_number))
                        PDB.REMARK_1.append([[],[]])
                        PDB.REMARK_1[-1][0].append(TEXT)
                        PDB.REMARK_1[-1][1].append(line_number)
                        PDB.REMARK_1[-1][1].append(line_number)
                        
            elif RNUM == 2:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_2 == None:
                    PDB.REMARK_2 = [[[],[]]]
                    PDB.REMARK_2[0][0].append(TEXT)
                    PDB.REMARK_2[0][1].append(line_number)
                    PDB.REMARK_2[0][1].append(line_number)
                else:
                    if PDB.REMARK_2[-1][1][1] == line_number - 1:
                        PDB.REMARK_2[-1][0].append(TEXT)
                        PDB.REMARK_2[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_2 section found at line {}.'.format(line_number))
                        PDB.REMARK_2.append([[],[]])
                        PDB.REMARK_2[-1][0].append(TEXT)
                        PDB.REMARK_2[-1][1].append(line_number)
                        PDB.REMARK_2[-1][1].append(line_number)
                        
            elif RNUM == 3:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_3 == None:
                    PDB.REMARK_3 = [[[],[]]]
                    PDB.REMARK_3[0][0].append(TEXT)
                    PDB.REMARK_3[0][1].append(line_number)
                    PDB.REMARK_3[0][1].append(line_number)
                else:
                    if PDB.REMARK_3[-1][1][1] == line_number - 1:
                        PDB.REMARK_3[-1][0].append(TEXT)
                        PDB.REMARK_3[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_3 section found at line {}.'.format(line_number))
                        PDB.REMARK_3.append([[],[]])
                        PDB.REMARK_3[-1][0].append(TEXT)
                        PDB.REMARK_3[-1][1].append(line_number)
                        PDB.REMARK_3[-1][1].append(line_number)
                        
            elif RNUM == 4:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_4 == None:
                    PDB.REMARK_4 = [[[],[]]]
                    PDB.REMARK_4[0][0].append(TEXT)
                    PDB.REMARK_4[0][1].append(line_number)
                    PDB.REMARK_4[0][1].append(line_number)
                else:
                    if PDB.REMARK_4[-1][1][1] == line_number - 1:
                        PDB.REMARK_4[-1][0].append(TEXT)
                        PDB.REMARK_4[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_4 section found at line {}.'.format(line_number))
                        PDB.REMARK_4.append([[],[]])
                        PDB.REMARK_4[-1][0].append(TEXT)
                        PDB.REMARK_4[-1][1].append(line_number)
                        PDB.REMARK_4[-1][1].append(line_number)
                        
            elif RNUM == 5:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_5 == None:
                    PDB.REMARK_5 = [[[],[]]]
                    PDB.REMARK_5[0][0].append(TEXT)
                    PDB.REMARK_5[0][1].append(line_number)
                    PDB.REMARK_5[0][1].append(line_number)
                else:
                    if PDB.REMARK_5[-1][1][1] == line_number - 1:
                        PDB.REMARK_5[-1][0].append(TEXT)
                        PDB.REMARK_5[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_5 section found at line {}.'.format(line_number))
                        PDB.REMARK_5.append([[],[]])
                        PDB.REMARK_5[-1][0].append(TEXT)
                        PDB.REMARK_5[-1][1].append(line_number)
                        PDB.REMARK_5[-1][1].append(line_number)
                        
            elif RNUM < 100:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_6 == None:
                    PDB.REMARK_6 = [[[],[]]]
                    PDB.REMARK_6[0][0].append(TEXT)
                    PDB.REMARK_6[0][1].append(line_number)
                    PDB.REMARK_6[0][1].append(line_number)
                else:
                    if PDB.REMARK_6[-1][1][1] == line_number - 1:
                        PDB.REMARK_6[-1][0].append(TEXT)
                        PDB.REMARK_6[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_6 section found at line {}.'.format(line_number))
                        PDB.REMARK_6.append([[],[]])
                        PDB.REMARK_6[-1][0].append(TEXT)
                        PDB.REMARK_6[-1][1].append(line_number)
                        PDB.REMARK_6[-1][1].append(line_number)
                        
            elif RNUM == 100:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_100 == None:
                    PDB.REMARK_100 = [[[],[]]]
                    PDB.REMARK_100[0][0].append(TEXT)
                    PDB.REMARK_100[0][1].append(line_number)
                    PDB.REMARK_100[0][1].append(line_number)
                else:
                    if PDB.REMARK_100[-1][1][1] == line_number - 1:
                        PDB.REMARK_100[-1][0].append(TEXT)
                        PDB.REMARK_100[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_100 section found at line {}.'.format(line_number))
                        PDB.REMARK_100.append([[],[]])
                        PDB.REMARK_100[-1][0].append(TEXT)
                        PDB.REMARK_100[-1][1].append(line_number)
                        PDB.REMARK_100[-1][1].append(line_number)
                        
            elif RNUM == 200:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_200 == None:
                    PDB.REMARK_200 = [[[],[]]]
                    PDB.REMARK_200[0][0].append(TEXT)
                    PDB.REMARK_200[0][1].append(line_number)
                    PDB.REMARK_200[0][1].append(line_number)
                else:
                    if PDB.REMARK_200[-1][1][1] == line_number - 1:
                        PDB.REMARK_200[-1][0].append(TEXT)
                        PDB.REMARK_200[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_200 section found at line {}.'.format(line_number))
                        PDB.REMARK_200.append([[],[]])
                        PDB.REMARK_200[-1][0].append(TEXT)
                        PDB.REMARK_200[-1][1].append(line_number)
                        PDB.REMARK_200[-1][1].append(line_number)
                        
            elif RNUM == 205:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_205 == None:
                    PDB.REMARK_205 = [[[],[]]]
                    PDB.REMARK_205[0][0].append(TEXT)
                    PDB.REMARK_205[0][1].append(line_number)
                    PDB.REMARK_205[0][1].append(line_number)
                else:
                    if PDB.REMARK_205[-1][1][1] == line_number - 1:
                        PDB.REMARK_205[-1][0].append(TEXT)
                        PDB.REMARK_205[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_205 section found at line {}.'.format(line_number))
                        PDB.REMARK_205.append([[],[]])
                        PDB.REMARK_205[-1][0].append(TEXT)
                        PDB.REMARK_205[-1][1].append(line_number)
                        PDB.REMARK_205[-1][1].append(line_number)
                        
            elif RNUM < 218:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_210 == None:
                    PDB.REMARK_210 = [[[],[]]]
                    PDB.REMARK_210[0][0].append(TEXT)
                    PDB.REMARK_210[0][1].append(line_number)
                    PDB.REMARK_210[0][1].append(line_number)
                else:
                    if PDB.REMARK_210[-1][1][1] == line_number - 1:
                        PDB.REMARK_210[-1][0].append(TEXT)
                        PDB.REMARK_210[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_210 section found at line {}.'.format(line_number))
                        PDB.REMARK_210.append([[],[]])
                        PDB.REMARK_210[-1][0].append(TEXT)
                        PDB.REMARK_210[-1][1].append(line_number)
                        PDB.REMARK_210[-1][1].append(line_number)
                        
            elif RNUM == 230:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_230 == None:
                    PDB.REMARK_230 = [[[],[]]]
                    PDB.REMARK_230[0][0].append(TEXT)
                    PDB.REMARK_230[0][1].append(line_number)
                    PDB.REMARK_230[0][1].append(line_number)
                else:
                    if PDB.REMARK_230[-1][1][1] == line_number - 1:
                        PDB.REMARK_230[-1][0].append(TEXT)
                        PDB.REMARK_230[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_230 section found at line {}.'.format(line_number))
                        PDB.REMARK_230.append([[],[]])
                        PDB.REMARK_230[-1][0].append(TEXT)
                        PDB.REMARK_230[-1][1].append(line_number)
                        PDB.REMARK_230[-1][1].append(line_number)
                        
            elif RNUM == 240:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_240 == None:
                    PDB.REMARK_240 = [[[],[]]]
                    PDB.REMARK_240[0][0].append(TEXT)
                    PDB.REMARK_240[0][1].append(line_number)
                    PDB.REMARK_240[0][1].append(line_number)
                else:
                    if PDB.REMARK_240[-1][1][1] == line_number - 1:
                        PDB.REMARK_240[-1][0].append(TEXT)
                        PDB.REMARK_240[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_240 section found at line {}.'.format(line_number))
                        PDB.REMARK_240.append([[],[]])
                        PDB.REMARK_240[-1][0].append(TEXT)
                        PDB.REMARK_240[-1][1].append(line_number)
                        PDB.REMARK_240[-1][1].append(line_number)
                        
            elif RNUM == 245:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_245 == None:
                    PDB.REMARK_245 = [[[],[]]]
                    PDB.REMARK_245[0][0].append(TEXT)
                    PDB.REMARK_245[0][1].append(line_number)
                    PDB.REMARK_245[0][1].append(line_number)
                else:
                    if PDB.REMARK_245[-1][1][1] == line_number - 1:
                        PDB.REMARK_245[-1][0].append(TEXT)
                        PDB.REMARK_245[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_245 section found at line {}.'.format(line_number))
                        PDB.REMARK_245.append([[],[]])
                        PDB.REMARK_245[-1][0].append(TEXT)
                        PDB.REMARK_245[-1][1].append(line_number)
                        PDB.REMARK_245[-1][1].append(line_number)
                        
            elif RNUM == 247:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_247 == None:
                    PDB.REMARK_247 = [[[],[]]]
                    PDB.REMARK_247[0][0].append(TEXT)
                    PDB.REMARK_247[0][1].append(line_number)
                    PDB.REMARK_247[0][1].append(line_number)
                else:
                    if PDB.REMARK_247[-1][1][1] == line_number - 1:
                        PDB.REMARK_247[-1][0].append(TEXT)
                        PDB.REMARK_247[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_247 section found at line {}.'.format(line_number))
                        PDB.REMARK_247.append([[],[]])
                        PDB.REMARK_247[-1][0].append(TEXT)
                        PDB.REMARK_247[-1][1].append(line_number)
                        PDB.REMARK_247[-1][1].append(line_number)
                        
            elif RNUM == 250:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_250 == None:
                    PDB.REMARK_250 = [[[],[]]]
                    PDB.REMARK_250[0][0].append(TEXT)
                    PDB.REMARK_250[0][1].append(line_number)
                    PDB.REMARK_250[0][1].append(line_number)
                else:
                    if PDB.REMARK_250[-1][1][1] == line_number - 1:
                        PDB.REMARK_250[-1][0].append(TEXT)
                        PDB.REMARK_250[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_250 section found at line {}.'.format(line_number))
                        PDB.REMARK_250.append([[],[]])
                        PDB.REMARK_250[-1][0].append(TEXT)
                        PDB.REMARK_250[-1][1].append(line_number)
                        PDB.REMARK_250[-1][1].append(line_number)
                        
            elif RNUM == 265:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_265 == None:
                    PDB.REMARK_265 = [[[],[]]]
                    PDB.REMARK_265[0][0].append(TEXT)
                    PDB.REMARK_265[0][1].append(line_number)
                    PDB.REMARK_265[0][1].append(line_number)
                else:
                    if PDB.REMARK_265[-1][1][1] == line_number - 1:
                        PDB.REMARK_265[-1][0].append(TEXT)
                        PDB.REMARK_265[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_265 section found at line {}.'.format(line_number))
                        PDB.REMARK_265.append([[],[]])
                        PDB.REMARK_265[-1][0].append(TEXT)
                        PDB.REMARK_265[-1][1].append(line_number)
                        PDB.REMARK_265[-1][1].append(line_number)
                        
            elif RNUM == 280:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_280 == None:
                    PDB.REMARK_280 = [[[],[]]]
                    PDB.REMARK_280[0][0].append(TEXT)
                    PDB.REMARK_280[0][1].append(line_number)
                    PDB.REMARK_280[0][1].append(line_number)
                else:
                    if PDB.REMARK_280[-1][1][1] == line_number - 1:
                        PDB.REMARK_280[-1][0].append(TEXT)
                        PDB.REMARK_280[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_280 section found at line {}.'.format(line_number))
                        PDB.REMARK_280.append([[],[]])
                        PDB.REMARK_280[-1][0].append(TEXT)
                        PDB.REMARK_280[-1][1].append(line_number)
                        PDB.REMARK_280[-1][1].append(line_number)
                        
            elif RNUM == 285:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_285 == None:
                    PDB.REMARK_285 = [[[],[]]]
                    PDB.REMARK_285[0][0].append(TEXT)
                    PDB.REMARK_285[0][1].append(line_number)
                    PDB.REMARK_285[0][1].append(line_number)
                else:
                    if PDB.REMARK_285[-1][1][1] == line_number - 1:
                        PDB.REMARK_285[-1][0].append(TEXT)
                        PDB.REMARK_285[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_285 section found at line {}.'.format(line_number))
                        PDB.REMARK_285.append([[],[]])
                        PDB.REMARK_285[-1][0].append(TEXT)
                        PDB.REMARK_285[-1][1].append(line_number)
                        PDB.REMARK_285[-1][1].append(line_number)
                        
            elif RNUM == 290:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_290 == None:
                    PDB.REMARK_290 = [[[],[]]]
                    PDB.REMARK_290[0][0].append(TEXT)
                    PDB.REMARK_290[0][1].append(line_number)
                    PDB.REMARK_290[0][1].append(line_number)
                else:
                    if PDB.REMARK_290[-1][1][1] == line_number - 1:
                        PDB.REMARK_290[-1][0].append(TEXT)
                        PDB.REMARK_290[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_290 section found at line {}.'.format(line_number))
                        PDB.REMARK_290.append([[],[]])
                        PDB.REMARK_290[-1][0].append(TEXT)
                        PDB.REMARK_290[-1][1].append(line_number)
                        PDB.REMARK_290[-1][1].append(line_number)
                        
            elif RNUM == 300:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_300 == None:
                    PDB.REMARK_300 = [[[],[]]]
                    PDB.REMARK_300[0][0].append(TEXT)
                    PDB.REMARK_300[0][1].append(line_number)
                    PDB.REMARK_300[0][1].append(line_number)
                else:
                    if PDB.REMARK_300[-1][1][1] == line_number - 1:
                        PDB.REMARK_300[-1][0].append(TEXT)
                        PDB.REMARK_300[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_300 section found at line {}.'.format(line_number))
                        PDB.REMARK_300.append([[],[]])
                        PDB.REMARK_300[-1][0].append(TEXT)
                        PDB.REMARK_300[-1][1].append(line_number)
                        PDB.REMARK_300[-1][1].append(line_number)
                        
            elif RNUM == 350:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_350 == None:
                    PDB.REMARK_350 = [[[],[]]]
                    PDB.REMARK_350[0][0].append(TEXT)
                    PDB.REMARK_350[0][1].append(line_number)
                    PDB.REMARK_350[0][1].append(line_number)
                else:
                    if PDB.REMARK_350[-1][1][1] == line_number - 1:
                        PDB.REMARK_350[-1][0].append(TEXT)
                        PDB.REMARK_350[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_350 section found at line {}.'.format(line_number))
                        PDB.REMARK_350.append([[],[]])
                        PDB.REMARK_350[-1][0].append(TEXT)
                        PDB.REMARK_350[-1][1].append(line_number)
                        PDB.REMARK_350[-1][1].append(line_number)
                        
            elif RNUM == 375:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_375 == None:
                    PDB.REMARK_375 = [[[],[]]]
                    PDB.REMARK_375[0][0].append(TEXT)
                    PDB.REMARK_375[0][1].append(line_number)
                    PDB.REMARK_375[0][1].append(line_number)
                else:
                    if PDB.REMARK_375[-1][1][1] == line_number - 1:
                        PDB.REMARK_375[-1][0].append(TEXT)
                        PDB.REMARK_375[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_375 section found at line {}.'.format(line_number))
                        PDB.REMARK_375.append([[],[]])
                        PDB.REMARK_375[-1][0].append(TEXT)
                        PDB.REMARK_375[-1][1].append(line_number)
                        PDB.REMARK_375[-1][1].append(line_number)
                        
            elif RNUM == 400:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_400 == None:
                    PDB.REMARK_400 = [[[],[]]]
                    PDB.REMARK_400[0][0].append(TEXT)
                    PDB.REMARK_400[0][1].append(line_number)
                    PDB.REMARK_400[0][1].append(line_number)
                else:
                    if PDB.REMARK_400[-1][1][1] == line_number - 1:
                        PDB.REMARK_400[-1][0].append(TEXT)
                        PDB.REMARK_400[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_400 section found at line {}.'.format(line_number))
                        PDB.REMARK_400.append([[],[]])
                        PDB.REMARK_400[-1][0].append(TEXT)
                        PDB.REMARK_400[-1][1].append(line_number)
                        PDB.REMARK_400[-1][1].append(line_number)
                        
            elif RNUM == 450:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_450 == None:
                    PDB.REMARK_450 = [[[],[]]]
                    PDB.REMARK_450[0][0].append(TEXT)
                    PDB.REMARK_450[0][1].append(line_number)
                    PDB.REMARK_450[0][1].append(line_number)
                else:
                    if PDB.REMARK_450[-1][1][1] == line_number - 1:
                        PDB.REMARK_450[-1][0].append(TEXT)
                        PDB.REMARK_450[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_450 section found at line {}.'.format(line_number))
                        PDB.REMARK_450.append([[],[]])
                        PDB.REMARK_450[-1][0].append(TEXT)
                        PDB.REMARK_450[-1][1].append(line_number)
                        PDB.REMARK_450[-1][1].append(line_number)
                        
            elif RNUM == 465:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_465 == None:
                    PDB.REMARK_465 = [[[],[]]]
                    PDB.REMARK_465[0][0].append(TEXT)
                    PDB.REMARK_465[0][1].append(line_number)
                    PDB.REMARK_465[0][1].append(line_number)
                else:
                    if PDB.REMARK_465[-1][1][1] == line_number - 1:
                        PDB.REMARK_465[-1][0].append(TEXT)
                        PDB.REMARK_465[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_465 section found at line {}.'.format(line_number))
                        PDB.REMARK_465.append([[],[]])
                        PDB.REMARK_465[-1][0].append(TEXT)
                        PDB.REMARK_465[-1][1].append(line_number)
                        PDB.REMARK_465[-1][1].append(line_number)
                        
            elif RNUM == 470:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_470 == None:
                    PDB.REMARK_470 = [[[],[]]]
                    PDB.REMARK_470[0][0].append(TEXT)
                    PDB.REMARK_470[0][1].append(line_number)
                    PDB.REMARK_470[0][1].append(line_number)
                else:
                    if PDB.REMARK_470[-1][1][1] == line_number - 1:
                        PDB.REMARK_470[-1][0].append(TEXT)
                        PDB.REMARK_470[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_470 section found at line {}.'.format(line_number))
                        PDB.REMARK_470.append([[],[]])
                        PDB.REMARK_470[-1][0].append(TEXT)
                        PDB.REMARK_470[-1][1].append(line_number)
                        PDB.REMARK_470[-1][1].append(line_number)
                        
            elif RNUM == 475:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_475 == None:
                    PDB.REMARK_475 = [[[],[]]]
                    PDB.REMARK_475[0][0].append(TEXT)
                    PDB.REMARK_475[0][1].append(line_number)
                    PDB.REMARK_475[0][1].append(line_number)
                else:
                    if PDB.REMARK_475[-1][1][1] == line_number - 1:
                        PDB.REMARK_475[-1][0].append(TEXT)
                        PDB.REMARK_475[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_475 section found at line {}.'.format(line_number))
                        PDB.REMARK_475.append([[],[]])
                        PDB.REMARK_475[-1][0].append(TEXT)
                        PDB.REMARK_475[-1][1].append(line_number)
                        PDB.REMARK_475[-1][1].append(line_number)
                        
            elif RNUM == 480:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_480 == None:
                    PDB.REMARK_480 = [[[],[]]]
                    PDB.REMARK_480[0][0].append(TEXT)
                    PDB.REMARK_480[0][1].append(line_number)
                    PDB.REMARK_480[0][1].append(line_number)
                else:
                    if PDB.REMARK_480[-1][1][1] == line_number - 1:
                        PDB.REMARK_480[-1][0].append(TEXT)
                        PDB.REMARK_480[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_400 section found at line {}.'.format(line_number))
                        PDB.REMARK_480.append([[],[]])
                        PDB.REMARK_480[-1][0].append(TEXT)
                        PDB.REMARK_480[-1][1].append(line_number)
                        PDB.REMARK_480[-1][1].append(line_number)
                        
            elif RNUM == 500:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_500 == None:
                    PDB.REMARK_500 = [[[],[]]]
                    PDB.REMARK_500[0][0].append(TEXT)
                    PDB.REMARK_500[0][1].append(line_number)
                    PDB.REMARK_500[0][1].append(line_number)
                else:
                    if PDB.REMARK_500[-1][1][1] == line_number - 1:
                        PDB.REMARK_500[-1][0].append(TEXT)
                        PDB.REMARK_500[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_500 section found at line {}.'.format(line_number))
                        PDB.REMARK_500.append([[],[]])
                        PDB.REMARK_500[-1][0].append(TEXT)
                        PDB.REMARK_500[-1][1].append(line_number)
                        PDB.REMARK_500[-1][1].append(line_number)
                        
            elif RNUM == 525:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_525 == None:
                    PDB.REMARK_525 = [[[],[]]]
                    PDB.REMARK_525[0][0].append(TEXT)
                    PDB.REMARK_525[0][1].append(line_number)
                    PDB.REMARK_525[0][1].append(line_number)
                else:
                    if PDB.REMARK_525[-1][1][1] == line_number - 1:
                        PDB.REMARK_525[-1][0].append(TEXT)
                        PDB.REMARK_525[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_525 section found at line {}.'.format(line_number))
                        PDB.REMARK_525.append([[],[]])
                        PDB.REMARK_525[-1][0].append(TEXT)
                        PDB.REMARK_525[-1][1].append(line_number)
                        PDB.REMARK_525[-1][1].append(line_number)
                        
            elif RNUM == 600:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_600 == None:
                    PDB.REMARK_600 = [[[],[]]]
                    PDB.REMARK_600[0][0].append(TEXT)
                    PDB.REMARK_600[0][1].append(line_number)
                    PDB.REMARK_600[0][1].append(line_number)
                else:
                    if PDB.REMARK_600[-1][1][1] == line_number - 1:
                        PDB.REMARK_600[-1][0].append(TEXT)
                        PDB.REMARK_600[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_600 section found at line {}.'.format(line_number))
                        PDB.REMARK_600.append([[],[]])
                        PDB.REMARK_600[-1][0].append(TEXT)
                        PDB.REMARK_600[-1][1].append(line_number)
                        PDB.REMARK_600[-1][1].append(line_number)
                        
            elif RNUM == 610:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_610 == None:
                    PDB.REMARK_610 = [[[],[]]]
                    PDB.REMARK_610[0][0].append(TEXT)
                    PDB.REMARK_610[0][1].append(line_number)
                    PDB.REMARK_610[0][1].append(line_number)
                else:
                    if PDB.REMARK_610[-1][1][1] == line_number - 1:
                        PDB.REMARK_610[-1][0].append(TEXT)
                        PDB.REMARK_610[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_610 section found at line {}.'.format(line_number))
                        PDB.REMARK_610.append([[],[]])
                        PDB.REMARK_610[-1][0].append(TEXT)
                        PDB.REMARK_610[-1][1].append(line_number)
                        PDB.REMARK_610[-1][1].append(line_number)
                        
            elif RNUM == 615:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_615 == None:
                    PDB.REMARK_615 = [[[],[]]]
                    PDB.REMARK_615[0][0].append(TEXT)
                    PDB.REMARK_615[0][1].append(line_number)
                    PDB.REMARK_615[0][1].append(line_number)
                else:
                    if PDB.REMARK_615[-1][1][1] == line_number - 1:
                        PDB.REMARK_615[-1][0].append(TEXT)
                        PDB.REMARK_615[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_615 section found at line {}.'.format(line_number))
                        PDB.REMARK_615.append([[],[]])
                        PDB.REMARK_615[-1][0].append(TEXT)
                        PDB.REMARK_615[-1][1].append(line_number)
                        PDB.REMARK_615[-1][1].append(line_number)
                        
            elif RNUM == 620:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_620 == None:
                    PDB.REMARK_620 = [[[],[]]]
                    PDB.REMARK_620[0][0].append(TEXT)
                    PDB.REMARK_620[0][1].append(line_number)
                    PDB.REMARK_620[0][1].append(line_number)
                else:
                    if PDB.REMARK_620[-1][1][1] == line_number - 1:
                        PDB.REMARK_620[-1][0].append(TEXT)
                        PDB.REMARK_620[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_620 section found at line {}.'.format(line_number))
                        PDB.REMARK_620.append([[],[]])
                        PDB.REMARK_620[-1][0].append(TEXT)
                        PDB.REMARK_620[-1][1].append(line_number)
                        PDB.REMARK_620[-1][1].append(line_number)
                        
            elif RNUM == 630:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_630 == None:
                    PDB.REMARK_630 = [[[],[]]]
                    PDB.REMARK_630[0][0].append(TEXT)
                    PDB.REMARK_630[0][1].append(line_number)
                    PDB.REMARK_630[0][1].append(line_number)
                else:
                    if PDB.REMARK_630[-1][1][1] == line_number - 1:
                        PDB.REMARK_630[-1][0].append(TEXT)
                        PDB.REMARK_630[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_630 section found at line {}.'.format(line_number))
                        PDB.REMARK_630.append([[],[]])
                        PDB.REMARK_630[-1][0].append(TEXT)
                        PDB.REMARK_630[-1][1].append(line_number)
                        PDB.REMARK_630[-1][1].append(line_number)
                        
            elif RNUM == 650:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_650 == None:
                    PDB.REMARK_650 = [[[],[]]]
                    PDB.REMARK_650[0][0].append(TEXT)
                    PDB.REMARK_650[0][1].append(line_number)
                    PDB.REMARK_650[0][1].append(line_number)
                else:
                    if PDB.REMARK_650[-1][1][1] == line_number - 1:
                        PDB.REMARK_650[-1][0].append(TEXT)
                        PDB.REMARK_650[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_650 section found at line {}.'.format(line_number))
                        PDB.REMARK_650.append([[],[]])
                        PDB.REMARK_650[-1][0].append(TEXT)
                        PDB.REMARK_650[-1][1].append(line_number)
                        PDB.REMARK_650[-1][1].append(line_number)
                        
            elif RNUM == 700:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_700 == None:
                    PDB.REMARK_700 = [[[],[]]]
                    PDB.REMARK_700[0][0].append(TEXT)
                    PDB.REMARK_700[0][1].append(line_number)
                    PDB.REMARK_700[0][1].append(line_number)
                else:
                    if PDB.REMARK_700[-1][1][1] == line_number - 1:
                        PDB.REMARK_700[-1][0].append(TEXT)
                        PDB.REMARK_700[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_700 section found at line {}.'.format(line_number))
                        PDB.REMARK_700.append([[],[]])
                        PDB.REMARK_700[-1][0].append(TEXT)
                        PDB.REMARK_700[-1][1].append(line_number)
                        PDB.REMARK_700[-1][1].append(line_number)
                        
            elif RNUM == 800:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_800 == None:
                    PDB.REMARK_800 = [[[],[]]]
                    PDB.REMARK_800[0][0].append(TEXT)
                    PDB.REMARK_800[0][1].append(line_number)
                    PDB.REMARK_800[0][1].append(line_number)
                else:
                    if PDB.REMARK_800[-1][1][1] == line_number - 1:
                        PDB.REMARK_800[-1][0].append(TEXT)
                        PDB.REMARK_800[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_800 section found at line {}.'.format(line_number))
                        PDB.REMARK_800.append([[],[]])
                        PDB.REMARK_800[-1][0].append(TEXT)
                        PDB.REMARK_800[-1][1].append(line_number)
                        PDB.REMARK_800[-1][1].append(line_number)
                        
            elif RNUM == 999:
                TEXT = line[10:-1].strip()
                if PDB.REMARK_999 == None:
                    PDB.REMARK_999 = [[[],[]]]
                    PDB.REMARK_999[0][0].append(TEXT)
                    PDB.REMARK_999[0][1].append(line_number)
                    PDB.REMARK_999[0][1].append(line_number)
                else:
                    if PDB.REMARK_999[-1][1][1] == line_number - 1:
                        PDB.REMARK_999[-1][0].append(TEXT)
                        PDB.REMARK_999[-1][1][1] = line_number
                    else:
                        warnings.warn('Another REMARK_999 section found at line {}.'.format(line_number))
                        PDB.REMARK_999.append([[],[]])
                        PDB.REMARK_999[-1][0].append(TEXT)
                        PDB.REMARK_999[-1][1].append(line_number)
                        PDB.REMARK_999[-1][1].append(line_number)
                        
        elif RECORD == 'DBREF':#Mx1
            TEXT = line[7:-1].strip()
            if PDB.DBREF == None:
                PDB.DBREF  = [[],[]]
                PDB.DBREF[0].append(TEXT)
                PDB.DBREF[1].append(line_number)
            else:
                PDB.DBREF[0].append(TEXT)
                PDB.DBREF[1].append(line_number)
        
        elif RECORD == 'SEQADV':#Mx1
            TEXT = line[7:-1].strip()
            if PDB.SEQADV == None:
                PDB.SEQADV  = [[],[]]
                PDB.SEQADV[0].append(TEXT)
                PDB.SEQADV[1].append(line_number)
            else:
                PDB.SEQADV[0].append(TEXT)
                PDB.SEQADV[1].append(line_number)
        
        elif RECORD == 'SEQRES':#MxN
            TEXT = line[7:-1].strip()
            if PDB.SEQRES == None:
                PDB.SEQRES = [[[],[]]]
                PDB.SEQRES[0][0].append(TEXT)
                PDB.SEQRES[0][1].append(line_number)
                PDB.SEQRES[0][1].append(line_number)
            else:
                if PDB.SEQRES[-1][1][1] == line_number - 1:
                    PDB.SEQRES[-1][0].append(TEXT)
                    PDB.SEQRES[-1][1][1] = line_number
                else:
                    PDB.SEQRES.append([[],[]])
                    PDB.SEQRES[-1][0].append(TEXT)
                    PDB.SEQRES[-1][1].append(line_number)
                    PDB.SEQRES[-1][1].append(line_number)
        
        elif RECORD == 'MODRES':#Mx1
            TEXT = line[7:-1].strip()
            if PDB.MODRES == None:
                PDB.MODRES  = [[],[]]
                PDB.MODRES[0].append(TEXT)
                PDB.MODRES[1].append(line_number)
            else:
                PDB.MODRES[0].append(TEXT)
                PDB.MODRES[1].append(line_number)
        
        elif RECORD == 'HET':#Mx1
            TEXT = line[7:-1].strip()
            if PDB.HET == None:
                PDB.HET  = [[],[]]
                PDB.HET[0].append(TEXT)
                PDB.HET[1].append(line_number)
            else:
                PDB.HET[0].append(TEXT)
                PDB.HET[1].append(line_number)
        
        elif RECORD == 'HETNAM':#MxN
            TEXT = line[7:-1].strip()
            if PDB.HETNAM == None:
                PDB.HETNAM = [[[],[]]]
                PDB.HETNAM[0][0].append(TEXT)
                PDB.HETNAM[0][1].append(line_number)
                PDB.HETNAM[0][1].append(line_number)
            else:
                if PDB.HETNAM[-1][1][1] == line_number - 1:
                    PDB.HETNAM[-1][0].append(TEXT)
                    PDB.HETNAM[-1][1][1] = line_number
                else:
                    PDB.HETNAM.append([[],[]])
                    PDB.HETNAM[-1][0].append(TEXT)
                    PDB.HETNAM[-1][1].append(line_number)
                    PDB.HETNAM[-1][1].append(line_number)
        
        elif RECORD == 'HETSYN':#MxN
            TEXT = line[7:-1].strip()
            if PDB.HETSYN == None:
                PDB.HETSYN = [[[],[]]]
                PDB.HETSYN[0][0].append(TEXT)
                PDB.HETSYN[0][1].append(line_number)
                PDB.HETSYN[0][1].append(line_number)
            else:
                if PDB.HETSYN[-1][1][1] == line_number - 1:
                    PDB.HETSYN[-1][0].append(TEXT)
                    PDB.HETSYN[-1][1][1] = line_number
                else:
                    PDB.HETSYN.append([[],[]])
                    PDB.HETSYN[-1][0].append(TEXT)
                    PDB.HETSYN[-1][1].append(line_number)
                    PDB.HETSYN[-1][1].append(line_number)
        
        elif RECORD == 'FORMUL':#MxN
            TEXT = line[7:-1].strip()
            if PDB.FORMUL == None:
                PDB.FORMUL = [[[],[]]]
                PDB.FORMUL[0][0].append(TEXT)
                PDB.FORMUL[0][1].append(line_number)
                PDB.FORMUL[0][1].append(line_number)
            else:
                if PDB.FORMUL[-1][1][1] == line_number - 1:
                    PDB.FORMUL[-1][0].append(TEXT)
                    PDB.FORMUL[-1][1][1] = line_number
                else:
                    PDB.FORMUL.append([[],[]])
                    PDB.FORMUL[-1][0].append(TEXT)
                    PDB.FORMUL[-1][1].append(line_number)
                    PDB.FORMUL[-1][1].append(line_number)
        
        elif RECORD == 'HELIX':#Mx1
            TEXT = line#[7:-1].strip()
            if PDB.HELIX == None:
                PDB.HELIX  = [[],[]]
                PDB.HELIX[0].append(TEXT)
                PDB.HELIX[1].append(line_number)
            else:
                PDB.HELIX[0].append(TEXT)
                PDB.HELIX[1].append(line_number)
        
        elif RECORD == 'SHEET':#Mx1
            TEXT = line#[7:-1].strip()
            if PDB.SHEET == None:
                PDB.SHEET  = [[],[]]
                PDB.SHEET[0].append(TEXT)
                PDB.SHEET[1].append(line_number)
            else:
                PDB.SHEET[0].append(TEXT)
                PDB.SHEET[1].append(line_number)
        
        elif RECORD == 'SSBOND':#Mx1
            TEXT = line[7:-1].strip()
            if PDB.SSBOND == None:
                PDB.SSBOND  = [[],[]]
                PDB.SSBOND[0].append(TEXT)
                PDB.SSBOND[1].append(line_number)
            else:
                PDB.SSBOND[0].append(TEXT)
                PDB.SSBOND[1].append(line_number)
        
        elif RECORD == 'LINK':#Mx1
            TEXT = line[7:-1].strip()
            if PDB.LINK == None:
                PDB.LINK  = [[],[]]
                PDB.LINK[0].append(TEXT)
                PDB.LINK[1].append(line_number)
            else:
                PDB.LINK[0].append(TEXT)
                PDB.LINK[1].append(line_number)
        
        elif RECORD == 'CISPEP':#Mx1
            TEXT = line[7:-1].strip()
            if PDB.CISPEP == None:
                PDB.CISPEP  = [[],[]]
                PDB.CISPEP[0].append(TEXT)
                PDB.CISPEP[1].append(line_number)
            else:
                PDB.CISPEP[0].append(TEXT)
                PDB.CISPEP[1].append(line_number)
        
        elif RECORD == 'SITE':#MxN
            TEXT = line[7:-1].strip()
            if PDB.SITE == None:
                PDB.SITE = [[[],[]]]
                PDB.SITE[0][0].append(TEXT)
                PDB.SITE[0][1].append(line_number)
                PDB.SITE[0][1].append(line_number)
            else:
                if PDB.SITE[-1][1][1] == line_number - 1:
                    PDB.SITE[-1][0].append(TEXT)
                    PDB.SITE[-1][1][1] = line_number
                else:
                    PDB.SITE.append([[],[]])
                    PDB.SITE[-1][0].append(TEXT)
                    PDB.SITE[-1][1].append(line_number)
                    PDB.SITE[-1][1].append(line_number)
        
        elif RECORD == 'CRYST1':#1x1
            TEXT = line[7:-1].strip()
            if PDB.CRYST1 == None:
                PDB.CRYST1 = [TEXT, line_number]
            else:
                warnings.warn('Another CRYST1 section found at line {}!!!'.format(line_number))
                sys.exit(0)
                
        elif RECORD == 'ORIGX1':#1x1
            TEXT = line[7:-1].strip()
            if PDB.ORIGX1 == None:
                PDB.ORIGX1 = [TEXT, line_number]
            else:
                warnings.warn('Another ORIGX1 section found at line {}!!!'.format(line_number))
                sys.exit(0)
                
        elif RECORD == 'ORIGX2':#1x1
            TEXT = line[7:-1].strip()
            if PDB.ORIGX2 == None:
                PDB.ORIGX2 = [TEXT, line_number]
            else:
                warnings.warn('Another ORIGX2 section found at line {}!!!'.format(line_number))
                sys.exit(0)
                
        elif RECORD == 'ORIGX3':#1x1
            TEXT = line[7:-1].strip()
            if PDB.ORIGX3 == None:
                PDB.ORIGX3 = [TEXT, line_number]
            else:
                warnings.warn('Another ORIGX3 section found at line {}!!!'.format(line_number))
                sys.exit(0)
                
        elif RECORD == 'SCALE1':#1x1
            TEXT = line[7:-1].strip()
            if PDB.SCALE1 == None:
                PDB.SCALE1 = [TEXT, line_number]
            else:
                warnings.warn('Another SCALE1 section found at line {}!!!'.format(line_number))
                sys.exit(0)
                
        elif RECORD == 'SCALE2':#1x1
            TEXT = line[7:-1].strip()
            if PDB.SCALE2 == None:
                PDB.SCALE2 = [TEXT, line_number]
            else:
                warnings.warn('Another SCALE2 section found at line {}!!!'.format(line_number))
                sys.exit(0)
                
        elif RECORD == 'SCALE3':#1x1
            TEXT = line[7:-1].strip()
            if PDB.SCALE3 == None:
                PDB.SCALE3 = [TEXT, line_number]
            else:
                warnings.warn('Another SCALE3 section found at line {}!!!'.format(line_number))
                sys.exit(0)
                
        elif RECORD == 'MTRIX1':#MxN
            TEXT = line[7:-1].strip()
            if PDB.MTRIX1 == None:
                PDB.MTRIX1 = [[[],[]]]
                PDB.MTRIX1[0][0].append(TEXT)
                PDB.MTRIX1[0][1].append(line_number)
                PDB.MTRIX1[0][1].append(line_number)
            else:
                if PDB.MTRIX1[-1][1][1] == line_number - 1:
                    PDB.MTRIX1[-1][0].append(TEXT)
                    PDB.MTRIX1[-1][1][1] = line_number
                else:
                    PDB.MTRIX1.append([[],[]])
                    PDB.MTRIX1[-1][0].append(TEXT)
                    PDB.MTRIX1[-1][1].append(line_number)
                    PDB.MTRIX1[-1][1].append(line_number)
        
        elif RECORD == 'MTRIX2':#MxN
            TEXT = line[7:-1].strip()
            if PDB.MTRIX2 == None:
                PDB.MTRIX2 = [[[],[]]]
                PDB.MTRIX2[0][0].append(TEXT)
                PDB.MTRIX2[0][1].append(line_number)
                PDB.MTRIX2[0][1].append(line_number)
            else:
                if PDB.MTRIX2[-1][1][1] == line_number - 1:
                    PDB.MTRIX2[-1][0].append(TEXT)
                    PDB.MTRIX2[-1][1][1] = line_number
                else:
                    PDB.MTRIX2.append([[],[]])
                    PDB.MTRIX2[-1][0].append(TEXT)
                    PDB.MTRIX2[-1][1].append(line_number)
                    PDB.MTRIX2[-1][1].append(line_number)
        
        elif RECORD == 'MTRIX3':#MxN
            TEXT = line[7:-1].strip()
            if PDB.MTRIX3 == None:
                PDB.MTRIX3 = [[[],[]]]
                PDB.MTRIX3[0][0].append(TEXT)
                PDB.MTRIX3[0][1].append(line_number)
                PDB.MTRIX3[0][1].append(line_number)
            else:
                if PDB.MTRIX3[-1][1][1] == line_number - 1:
                    PDB.MTRIX3[-1][0].append(TEXT)
                    PDB.MTRIX3[-1][1][1] = line_number
                else:
                    PDB.MTRIX3.append([[],[]])
                    PDB.MTRIX3[-1][0].append(TEXT)
                    PDB.MTRIX3[-1][1].append(line_number)
                    PDB.MTRIX3[-1][1].append(line_number)
        
        elif RECORD == 'MODEL':#MxN
            TEXT = line[7:-1].strip()
            if PDB.MODEL == None:
                PDB.MODEL = [[[],[]]]
                PDB.MODEL[0][0].append(TEXT)
                PDB.MODEL[0][1].append(line_number)
                PDB.MODEL[0][1].append(line_number)
            else:
                if PDB.MODEL[-1][1][1] == line_number - 1:
                    PDB.MODEL[-1][0].append(TEXT)
                    PDB.MODEL[-1][1][1] = line_number
                else:
                    PDB.MODEL.append([[],[]])
                    PDB.MODEL[-1][0].append(TEXT)
                    PDB.MODEL[-1][1].append(line_number)
                    PDB.MODEL[-1][1].append(line_number)
        
        elif RECORD == 'ATOM':#Mx1
            TEXT = line
            if PDB.ATOM == None:
                PDB.ATOM  = [[],[]]
                PDB.ATOM[0].append(TEXT)
                PDB.ATOM[1].append(line_number)
            else:
                PDB.ATOM[0].append(TEXT)
                PDB.ATOM[1].append(line_number)
        
        elif RECORD == 'ANISOU':#Mx1
            TEXT = line[7:-1].strip()
            if PDB.ANISOU == None:
                PDB.ANISOU  = [[],[]]
                PDB.ANISOU[0].append(TEXT)
                PDB.ANISOU[1].append(line_number)
            else:
                PDB.ANISOU[0].append(TEXT)
                PDB.ANISOU[1].append(line_number)
        
        elif RECORD == 'TER':#Mx1
            TEXT = line[7:-1].strip()
            if PDB.TER == None:
                PDB.TER  = [[],[]]
                PDB.TER[0].append(TEXT)
                PDB.TER[1].append(line_number)
            else:
                PDB.TER[0].append(TEXT)
                PDB.TER[1].append(line_number)    
        
        elif RECORD == 'HETATM':#Mx1
            TEXT = line[7:-1].strip()
            if PDB.HETATM == None:
                PDB.HETATM  = [[],[]]
                PDB.HETATM[0].append(TEXT)
                PDB.HETATM[1].append(line_number)
            else:
                PDB.HETATM[0].append(TEXT)
                PDB.HETATM[1].append(line_number)    
        
        elif RECORD == 'ENDMDL':#MxN
            TEXT = line[7:-1].strip()
            if PDB.ENDMDL == None:
                PDB.ENDMDL = [[[],[]]]
                PDB.ENDMDL[0][0].append(TEXT)
                PDB.ENDMDL[0][1].append(line_number)
                PDB.ENDMDL[0][1].append(line_number)
            else:
                if PDB.ENDMDL[-1][1][1] == line_number - 1:
                    PDB.ENDMDL[-1][0].append(TEXT)
                    PDB.ENDMDL[-1][1][1] = line_number
                else:
                    PDB.ENDMDL.append([[],[]])
                    PDB.ENDMDL[-1][0].append(TEXT)
                    PDB.ENDMDL[-1][1].append(line_number)
                    PDB.ENDMDL[-1][1].append(line_number)
        
        elif RECORD == 'CONECT':#Mx1
            TEXT = line[7:-1].strip()
            if PDB.CONECT == None:
                PDB.CONECT  = [[],[]]
                PDB.CONECT[0].append(TEXT)
                PDB.CONECT[1].append(line_number)
            else:
                PDB.CONECT[0].append(TEXT)
                PDB.CONECT[1].append(line_number)    
        
        elif RECORD == 'MASTER':#1x1
            TEXT = line[7:-1].strip()
            if PDB.MASTER == None:
                PDB.MASTER = [TEXT, line_number]
            else:
                warnings.warn('Another MASTER section found at line {}!!!'.format(line_number))
                sys.exit(0)
        
        elif RECORD == 'END':#1x1
            if PDB.END == None:
                PDB.END = ['Done', line_number]
            else:
                warnings.warn('Another END section found at line {}!!!'.format(line_number))
                sys.exit(0)
        
        else:
            TEXT = line[:-1]
            if PDB.JUNK == None:
                PDB.JUNK  = [[],[]]
                PDB.JUNK[0].append(TEXT)
                PDB.JUNK[1].append(line_number)
                warnings.warn('Unrecognized section found at line {}.'.format(line_number))
            else:
                PDB.JUNK[0].append(TEXT)
                PDB.JUNK[1].append(line_number)
                warnings.warn('Unrecognized section found at line {}.'.format(line_number))
        
        line_number += 1
        
    #Any tests would go here that involve testing after reading and before analysing
    
    #Analysis starts here - This part is not complete yet
    
    if PDB.HELIX == None:
        pass
    else:
        HELIX = []
        for line in PDB.HELIX[0]:
            HELIX.append(Helix(line))
        PDB.HELIX[0] = HELIX

    if PDB.SHEET == None:
        pass
    else:
        SHEET = []
        for line in PDB.SHEET[0]:
            SHEET.append(Sheet(line))
        PDB.SHEET[0] = SHEET
    
    if PDB.ATOM == None:
        pass
    else:
        ATOM = []
        for line in PDB.ATOM[0]:
            ATOM.append(Atom(line))
        PDB.ATOM[0] = ATOM
    
    
    
    return(PDB)



##############################################################################

if __name__ == '__main__':
    print('PyDB is a package that is designed to read PDB file format.')
    File = readPDBfile(ID = '1a0k', Extension = 'pdb')
    PDB = parsePDBfile(File)
    #print(PDB.HEADER)
    #print(PDB.TITLE )
    #print(PDB.COMPND)
    #print(PDB.SOURCE)
    #print(PDB.KEYWDS)
    #print(PDB.EXPDTA)
    #print(PDB.NUMMDL)
    #print(PDB.AUTHOR)
    #print(PDB.REVDAT)
    #print(PDB.JRNL  )
    #print(PDB.REMARK_2)
    #print(PDB.REMARK_3)
    #print(PDB.DBREF )
    #print(PDB.SEQRES)
    #print(PDB.FORMUL)
    #print(PDB.HELIX )
    #print(PDB.SHEET )
    #print(PDB.CISPEP)
    #print(PDB.SITE  )
    #print(PDB.CRYST1)
    #print(PDB.ORIGX1)
    #print(PDB.ORIGX2)
    #print(PDB.ORIGX3)
    #print(PDB.SCALE1)
    #print(PDB.SCALE2)
    #print(PDB.SCALE3)
    #print(PDB.ATOM  )
    #print(PDB.TER   )
    #print(PDB.HETATM)
    #print(PDB.MASTER)
    #print(PDB.END   )
    print(PDB.JUNK  )
    
    