#!/usr/bin/env python
#This is a package that is designed to read PDB file format 
#and convert it to a python object.
__author__ = 'Csongor Mátyás'

import sys, warnings
from datetime import datetime

##############################################################################
#CLASSES######################################################################
##############################################################################

class PDBfile: #PDF file class, each PDB file will be a PDB object
    HEADER = None   #Mandatory 1x1 - First line of the entry, contains PDB ID code,
                                       #classification, and date of deposition.
    depDate= None   #Date of deposition - from HEADER line 
    idCode = None   #Unique PDB ID code - from HEADER line
    OBSLTE = None   #Optional  1xN - Statement that the entry has been removed from
                    #distribution and list of the ID code(s) which replaced it.
    TITLE  = None   #Mandatory 1xN - Description of the experiment represented in the entry.
    SPLIT  = None   #Optional  1xN - List of PDB entries that compose a larger
                    #macromolecular complexes.
    CAVEAT = None   #Optional  1xN - Severe error indicator.    
    COMPND = None   #Mandatory 1xN - Description of macromolecular contents of the entry.
    CHAINS = None   #                From COMPND, contains list of chains
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
    DBREF1 = None   #Optional  Mx1 - Same as DBREF but new version, goes with DBREF2
    DBREF2 = None   #Optional  Mx1 - goes with DBREF1
    SEQADV = None   #Optional  Mx1 - Identification of conflicts between PDB
                    #and the named sequence database.
    SEQRES = None   #Mandatory MxN - Primary sequence of backbone residues. As dictionary
    SEQRESlen=None  #                Number of amino acids in each chain as dictionary
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
class Dbref:    #Each DBREF or DBREF1/2 within a PDB object will be a Dbref object
    def __init__(self, LIST):
        if LIST[0] == 0:
            TEXT = LIST[1]
            
            # ('') ID code of this entry.
            if TEXT[7:11].strip() == '':
                self.idCode = None
            else:
                self.idCode = TEXT[7:11].strip()
        
            # ('') Chain identifier.
            if TEXT[12:13].strip() == '':
                self.chainID = None
            else:
                self.chainID = TEXT[12:13].strip()
            
            # (#) Initial sequence number of the PDB sequence segment.
            if TEXT[14:18].strip() == '':
                self.seqBegin = None
            else:
                self.seqBegin = int(TEXT[14:18].strip())
            
            # ('') Initial insertion code of the PDB sequence segment.
            if TEXT[18:19].strip() == '':
                self.insertBegin = None
            else:
                self.insertBegin = TEXT[18:19].strip()
            
            # (#) Ending sequence number of the PDB sequence segment.
            if TEXT[20:24].strip() == '':
                self.seqEnd = None
            else:
                self.seqEnd = int(TEXT[20:24].strip())
            
            # ('') Ending insertion code of the PDB sequence segment.
            if TEXT[24:25].strip() == '':
                self.insertEnd = None
            else:
                self.insertEnd = TEXT[24:25].strip()
            
            # ('') Sequence database name.
            if TEXT[26:32].strip() == '':
                self.database = None
            else:
                self.database = TEXT[26:32].strip()
            
            # ('') Sequence database accession code.
            if TEXT[33:41].strip() == '':
                self.dbAccession = None
            else:
                self.dbAccession = TEXT[33:41].strip()
            
            # ('') Sequence database identification code.
            if TEXT[42:54].strip() == '':
                self.dbIdCode = None
            else:
                self.dbIdCode = TEXT[42:54].strip()
            
            # (#) Initial sequence number of the database seqment.
            if TEXT[55:60].strip() == '':
                self.dbseqBegin = None
            else:
                self.dbseqBegin = int(TEXT[55:60].strip())
            
            # ('') Insertion code of initial residue of the segment, if PDB is the reference.
            if TEXT[60:61].strip() == '':
                self.dbinsBeg = None
            else:
                self.dbinsBeg = TEXT[60:61].strip()
            
            # (#) Ending sequence number of the database segment.
            if TEXT[62:67].strip() == '':
                self.dbseqEnd = None
            else:
                self.dbseqEnd = int(TEXT[62:67].strip())
            
            # ('') Insertion code of the ending residue of the segment, if PDB is the reference.
            if TEXT[67:68].strip() == '':
                self.dbinsEnd = None
            else:
                self.dbinsEnd = TEXT[67:68].strip()
                
        elif LIST[0] == 1:
            TEXT = LIST[1]
            
            # ('') ID code of this entry.
            if TEXT[7:11].strip() == '':
                self.idCode = None
            else:
                self.idCode = TEXT[7:11].strip()
        
            # ('') Chain identifier.
            if TEXT[12:13].strip() == '':
                self.chainID = None
            else:
                self.chainID = TEXT[12:13].strip()
            
            # (#) Initial sequence number of the PDB sequence segment.
            if TEXT[14:18].strip() == '':
                self.seqBegin = None
            else:
                self.seqBegin = int(TEXT[14:18].strip())
            
            # ('') Initial insertion code of the PDB sequence segment.
            if TEXT[18:19].strip() == '':
                self.insertBegin = None
            else:
                self.insertBegin = TEXT[18:19].strip()
            
            # (#) Ending sequence number of the PDB sequence segment.
            if TEXT[20:24].strip() == '':
                self.seqEnd = None
            else:
                self.seqEnd = int(TEXT[20:24].strip())
            
            # ('') Ending insertion code of the PDB sequence segment.
            if TEXT[24:25].strip() == '':
                self.insertEnd = None
            else:
                self.insertEnd = TEXT[24:25].strip()
            
            # ('') Sequence database name.
            if TEXT[26:32].strip() == '':
                self.database = None
            else:
                self.database = TEXT[26:32].strip()
                        
            # ('') Sequence database identification code.
            if TEXT[47:67].strip() == '':
                self.dbIdCode = None
            else:
                self.dbIdCode = TEXT[47:67].strip()
            
            TEXT = LIST[2]
            
            # ('') ID code of this entry.
            if TEXT[7:11].strip() == self.idCode:
                pass
            else:
                warnings.warn('idCode of DBREF1 is not equal to DBREF2 {} != {}.'.format(self.idCode, TEXT[7:11].strip()))
        
            # ('') Chain identifier.
            if TEXT[12:13].strip() == self.chainID:
                pass
            else:
                warnings.warn('chainID of DBREF1 is not equal to DBREF2 {} != {}.'.format(self.chainID, TEXT[12:13].strip()))
            
            # ('') Sequence database accession code.
            if TEXT[18:40].strip() == '':
                self.dbAccession = None
            else:
                self.dbAccession = TEXT[18:40].strip()
            
            # (#) Initial sequence number of the database seqment.
            if TEXT[45:55].strip() == '':
                self.dbseqBegin = None
            else:
                self.dbseqBegin = int(TEXT[45:55].strip())
            
            # ('') Insertion code of initial residue of the segment, if PDB is the reference.
            self.dbinsBeg = None

            # (#) Ending sequence number of the database segment.
            if TEXT[57:67].strip() == '':
                self.dbseqEnd = None
            else:
                self.dbseqEnd = int(TEXT[57:67].strip())
            
            # ('') Insertion code of the ending residue of the segment, if PDB is the reference.
            self.dbinsEnd = None

##############################################################################
class Het:    #Each HET within a PDB object will be a Het object
    def __init__(self, TEXT, line_number):
        
        # ('') Het identifier
        if TEXT[7:10].strip() == '':
            self.hetID = None
        else:
            self.hetID = TEXT[7:10].strip()
        
        # ('') Chain identifier
        if TEXT[12:13].strip() == '':
            self.chainID = None
        else:
            self.chainID = TEXT[12:13].strip()
        
        # (#) Sequence number 
        if TEXT[13:17].strip() == '':
            self.seqNum = None
        else:
            self.seqNum = int(TEXT[13:17].strip())
            
        # ('') Insertion code
        if TEXT[17:18].strip() == '':
            self.iCode = None
        else:
            self.iCode = TEXT[17:18].strip()
        
        # (#) Number of HETATM records for the group present in the entry
        if TEXT[20:25].strip() == '':
            self.numHetAtoms = None
        else:
            self.numHetAtoms = int(TEXT[20:25].strip())
        
        # ('') Text describing Het group
        if TEXT[30:70].strip() == '':
            self.text = None
        else:
            self.text = TEXT[30:70].strip()
            
        # (#) Line number
        self.lineNum = line_number

##############################################################################
class Helix:    #Each HELIX within a PDB object will be a Helix object
    def __init__(self, TEXT, line_number):
        
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
        
        # (#) Line number
        self.lineNum = line_number

##############################################################################

class Sheet:    #Each SHEET within a PDB object will be a Sheet object
    def __init__(self, TEXT, line_number):
        
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
        
        # (#) Line number
        self.lineNum = line_number

##############################################################################
class SSbond:    #Each SSBOND within a PDB object will be an SSbond object
    def __init__(self, TEXT, line_number):
        
        # (#) Serial number of the SSBOND.
        if TEXT[7:10].strip() == '':
            self.serNum = None
        else:
            self.serNum = int(TEXT[7:10].strip())
        
        # ('') Residue name.
        if TEXT[11:14].strip() == '':
            self.ResName1 = None
        else:
            self.ResName1 = TEXT[11:14].strip()
        
        # (A-Z) Chain identifier. 
        if TEXT[15:16].strip() == '':
            self.chainID1 = None
        else:
            self.chainID1 = TEXT[15:16].strip()
            
        # (#) Residue sequence number.
        if TEXT[17:21].strip() == '':
            self.seqNum1 = None
        else:
            self.seqNum1 = int(TEXT[17:21].strip())
        
        # Insertion code.
        if TEXT[21:22].strip() == '':
            self.icode1 = None
        else:
            self.icode1 = TEXT[21:22].strip()
        
        # ('') Second residue name.
        if TEXT[25:28].strip() == '':
            self.ResName2 = None
        else:
            self.ResName2 = TEXT[25:28].strip()
        
        # (A-Z) Second chain identifier. 
        if TEXT[29:30].strip() == '':
            self.chainID2 = None
        else:
            self.chainID2 = TEXT[29:30].strip()
            
        # (#) Second residue sequence number.
        if TEXT[31:35].strip() == '':
            self.seqNum2 = None
        else:
            self.seqNum2 = int(TEXT[31:35].strip())
        
        # Second insertion code.
        if TEXT[35:36].strip() == '':
            self.icode2 = None
        else:
            self.icode2 = TEXT[35:36].strip()
            
        # ('') Symmetry operator for residue 1.
        if TEXT[59:65].strip() == '':
            self.symOp1 = None
        else:
            self.symOp1 = TEXT[59:65].strip()
            
        # ('') Symmetry operator for residue 2.
        if TEXT[66:72].strip() == '':
            self.symOp2 = None
        else:
            self.symOp2 = TEXT[66:72].strip()
        
        # (#) Disulfide bond distance
        if TEXT[73:78].strip() == '':
            self.bondLen = None
        else:
            self.bondLen = float(TEXT[73:78].strip())
        
        # (#) Line number
        self.lineNum = line_number

##############################################################################
class Link:    #Each LINK within a PDB object will be an Link object
    def __init__(self, TEXT, line_number):
        
        # ('') Atom name.
        if TEXT[12:16].strip() == '':
            self.atomName1 = None
        else:
            self.atomName1 = TEXT[12:16].strip()
        
        # ('') Alternate location indicator.
        if TEXT[16:17].strip() == '':
            self.altLoc1 = None
        else:
            self.altLoc1 = TEXT[16:17].strip()
        
        # ('') Residue name.
        if TEXT[17:20].strip() == '':
            self.resName1 = None
        else:
            self.resName1 = TEXT[17:20].strip()
        
        # (A-Z) Chain identifier. 
        if TEXT[21:22].strip() == '':
            self.chainID1 = None
        else:
            self.chainID1 = TEXT[21:22].strip()
            
        # (#) Residue sequence number.
        if TEXT[22:26].strip() == '':
            self.seqNum1 = None
        else:
            self.seqNum1 = int(TEXT[22:26].strip())
        
        # Insertion code.
        if TEXT[26:27].strip() == '':
            self.icode1 = None
        else:
            self.icode1 = TEXT[26:27].strip()
        
        # ('') Second atom name.
        if TEXT[42:46].strip() == '':
            self.atomName2 = None
        else:
            self.atomName2 = TEXT[42:46].strip()
        
        # ('') Second alternate location indicator.
        if TEXT[46:47].strip() == '':
            self.altLoc2 = None
        else:
            self.altLoc2 = TEXT[46:47].strip()
        
        # ('') Second residue name.
        if TEXT[47:50].strip() == '':
            self.resName2 = None
        else:
            self.resName2 = TEXT[47:50].strip()
        
        # (A-Z) Second chain identifier. 
        if TEXT[51:52].strip() == '':
            self.chainID2 = None
        else:
            self.chainID2 = TEXT[51:52].strip()
            
        # (#) Second residue sequence number.
        if TEXT[52:56].strip() == '':
            self.seqNum2 = None
        else:
            self.seqNum2 = int(TEXT[52:56].strip())
        
        # Second insertion code.
        if TEXT[56:57].strip() == '':
            self.icode2 = None
        else:
            self.icode2 = TEXT[56:57].strip()
        
        # ('') Symmetry operator for residue 1.
        if TEXT[59:65].strip() == '':
            self.symOp1 = None
        else:
            self.symOp1 = TEXT[59:65].strip()
            
        # ('') Symmetry operator for residue 2.
        if TEXT[66:72].strip() == '':
            self.symOp2 = None
        else:
            self.symOp2 = TEXT[66:72].strip()
        
        # (#) Bond distance
        if TEXT[73:78].strip() == '':
            self.bondLen = None
        else:
            self.bondLen = float(TEXT[73:78].strip())
        
        # (#) Line number
        self.lineNum = line_number

##############################################################################
class Cispep:    #Each CISPEP within a PDB object will be an Cispep object
    def __init__(self, TEXT, line_number):
        
        # (#) Record serial number.
        if TEXT[7:10].strip() == '':
            self.serNum = None
        else:
            self.serNum = int(TEXT[7:10].strip())
        
        # ('') Residue name.
        if TEXT[11:14].strip() == '':
            self.resName1 = None
        else:
            self.resName1 = TEXT[11:14].strip()
        
        # (A-Z) Chain identifier. 
        if TEXT[15:16].strip() == '':
            self.chainID1 = None
        else:
            self.chainID1 = TEXT[15:16].strip()
            
        # (#) Residue sequence number.
        if TEXT[17:21].strip() == '':
            self.seqNum1 = None
        else:
            self.seqNum1 = int(TEXT[17:21].strip())
        
        # Insertion code.
        if TEXT[21:22].strip() == '':
            self.icode1 = None
        else:
            self.icode1 = TEXT[21:22].strip()
        
        # ('') Second residue name.
        if TEXT[25:28].strip() == '':
            self.resName2 = None
        else:
            self.resName2 = TEXT[25:28].strip()
        
        # (A-Z) Second chain identifier. 
        if TEXT[29:30].strip() == '':
            self.chainID2 = None
        else:
            self.chainID2 = TEXT[29:30].strip()
            
        # (#) Second residue sequence number.
        if TEXT[31:35].strip() == '':
            self.seqNum2 = None
        else:
            self.seqNum2 = int(TEXT[31:35].strip())
        
        # Second insertion code.
        if TEXT[35:36].strip() == '':
            self.icode2 = None
        else:
            self.icode2 = TEXT[35:36].strip()

        # ('') Identifies the specific model.
        if TEXT[43:46].strip() == '':
            self.modNum = None
        else:
            self.modNum = int(TEXT[43:46].strip())
        
        # (#) Angle measurement in degrees.
        if TEXT[53:59].strip() == '':
            self.angle = None
        else:
            self.angle = float(TEXT[53:59].strip())
        
        # (#) Line number
        self.lineNum = line_number

##############################################################################

class Atom:    #Each ATOM within a PDB object will be a Atom object
    def __init__(self, TEXT, line_number):
        
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
        
        # (#) Line number
        self.lineNum = line_number

##############################################################################
class Ter:    #Each TER within a PDB object will be an Ter object - Terminal residue
    def __init__(self, TEXT):
        
        # (#) Serial number.
        if TEXT[7:11].strip() == '':
            self.serNum = None
        else:
            self.serNum = int(TEXT[7:11].strip())
        
        # ('') Residue name.
        if TEXT[17:20].strip() == '':
            self.resName = None
        else:
            self.resName = TEXT[17:20].strip()
        
        # (A-Z) Chain identifier. 
        if TEXT[21:22].strip() == '':
            self.chainID = None
        else:
            self.chainID = TEXT[21:22].strip()
            
        # (#) Residue sequence number.
        if TEXT[22:26].strip() == '':
            self.seqNum = None
        else:
            self.seqNum = int(TEXT[22:26].strip())
        
        # Insertion code.
        if TEXT[26:27].strip() == '':
            self.icode = None
        else:
            self.icode = TEXT[26:27].strip()
            
##############################################################################

class Hetatm:    #Each HETATM within a PDB object will be a Hetatm object
    def __init__(self, TEXT, line_number):
        
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
        
        # (#) Line number
        self.lineNum = line_number

##############################################################################
#FUNCTIONS####################################################################
##############################################################################

def readPDBfile(ID=None, Extension=None, Filename=None, Verbose=False,
                IDUpper=False, IDLower=False, ExtUpper=False, ExtLower=False):
    #This code opens the file that has a given PDB ID . Extension, or Filename
    #ID or Extension can be changed to Upper or Lower case by changing IDUpper,
    #IDLower, ExtUpper, ExtLower flags to True
    #Verbose=True will print some more data

    if Filename == '':
        warnings.warn('Empty string was given as filename.')
        sys.exit(0) #Check if string is empty
        
    if ID == '':
        warnings.warn('Empty string was given as PDB ID')
        sys.exit(0) #Check if string is empty

    if IDUpper == True and IDLower == True:
        warnings.warn('Conflict between IDUpper and IDLower')
        sys.exit(0) #Check if both flags are raised
    
    if ExtUpper == True and ExtLower == True:
        warnings.warn('Conflict between ExtUpper and ExtLower')
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
            warnings.warn('There was no PDB ID or Filename given.')
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
        warnings.warn('Opening {} file.'.format(ToOpen))

    with open(ToOpen, 'r') as f:
        File = f.readlines()

    return(File)

##############################################################################

def define_sites(PDB_SITE):
    SITE = {}        
    counter = 0
    for TEXT in PDB_SITE[0]:
        if int(TEXT[7:10].strip()) == 1:
            SITE[counter] = [TEXT[11:14].strip(), int(TEXT[15:17].strip()), []]
            
            residue = []
            if TEXT[18:21].strip() != '':
                residue.append([TEXT[18:21].strip(), TEXT[22:23].strip(),
                                int(TEXT[23:27].strip()), TEXT[27:28].strip()])
                                
            if TEXT[29:32].strip() != '':
                residue.append([TEXT[29:32].strip(), TEXT[33:34].strip(),
                                int(TEXT[34:38].strip()), TEXT[38:39].strip()])
                                
            if TEXT[40:43].strip() != '':
                residue.append([TEXT[40:43].strip(), TEXT[44:45].strip(),
                                int(TEXT[45:49].strip()), TEXT[49:50].strip()])
                                
            if TEXT[51:54].strip() != '':
                residue.append([TEXT[51:54].strip(), TEXT[55:56].strip(),
                                int(TEXT[56:60].strip()), TEXT[60:61].strip()])
            SITE[counter][2].extend(residue)
            counter += 1
        else:
            residue = []
            if TEXT[18:21].strip() != '':
                residue.append([TEXT[18:21].strip(), TEXT[22:23].strip(),
                                int(TEXT[23:27].strip()), TEXT[27:28].strip()])
                                
            if TEXT[29:32].strip() != '':
                residue.append([TEXT[29:32].strip(), TEXT[33:34].strip(),
                                int(TEXT[34:38].strip()), TEXT[38:39].strip()])
                                
            if TEXT[40:43].strip() != '':
                residue.append([TEXT[40:43].strip(), TEXT[44:45].strip(),
                                int(TEXT[45:49].strip()), TEXT[49:50].strip()])
                                
            if TEXT[51:54].strip() != '':
                residue.append([TEXT[51:54].strip(), TEXT[55:56].strip(),
                                int(TEXT[56:60].strip()), TEXT[60:61].strip()])
            
            SITE[counter-1][2].extend(residue)
            
    return(SITE)
    
##############################################################################

def parsePDBfile(File):
    PDB = PDBfile()

    if File == None:
        warnings.warn('Input file variable is == None!')
    elif File == '':
        warnings.warn('Input file variable is empty!')
    
    line_number = 0
    
    for line in File:
        if   len(line) < 7:
            RECORD = line.strip()
        else:
            RECORD = line[0:7].strip()
            
        if RECORD == '':
            warnings.warn('Empty RECORD entry at line {}!!!'.format(line_number))
        
        elif RECORD == 'HEADER':#1x1
            TEXT = line
            if PDB.HEADER == None:
                PDB.HEADER = TEXT[10:50].strip()
                PDB.depDate = datetime.strptime(TEXT[50:59].strip(), '%d-%b-%y')
                PDB.idCode = TEXT[62:66].strip()
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
            TEXT = line
            if PDB.TITLE == None:
                PDB.TITLE = [[[],[]]]
                PDB.TITLE[0][0].append(TEXT[10:80].strip())
                PDB.TITLE[0][1].append(line_number)
                PDB.TITLE[0][1].append(line_number)
            else:
                if PDB.TITLE[-1][1][1] == line_number - 1:
                    PDB.TITLE[-1][0].append(TEXT[10:80].strip())
                    PDB.TITLE[-1][1][1] = line_number
                else:
                    warnings.warn('Another TITLE section found at line {}.'.format(line_number))
                    PDB.TITLE.append([[],[]])
                    PDB.TITLE[-1][0].append(TEXT[10:80].strip())
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
                PDB.CAVEAT[0][0].append(TEXT[19:79].strip())
                PDB.CAVEAT[0][1].append(line_number)
                PDB.CAVEAT[0][1].append(line_number)
            else:
                if PDB.CAVEAT[-1][1][1] == line_number - 1:
                    PDB.CAVEAT[-1][0].append(TEXT[19:79].strip())
                    PDB.CAVEAT[-1][1][1] = line_number
                else:
                    warnings.warn('Another CAVEAT section found at line {}.'.format(line_number))
                    PDB.CAVEAT.append([[],[]])
                    PDB.CAVEAT[-1][0].append(TEXT[19:79].strip())
                    PDB.CAVEAT[-1][1].append(line_number)
                    PDB.CAVEAT[-1][1].append(line_number)
        
        elif RECORD == 'COMPND':#1xN
            TEXT = line
            if PDB.COMPND == None:
                PDB.COMPND = [[[],[]]]
                PDB.COMPND[0][0].append(TEXT[10:80].strip())
                PDB.COMPND[0][1].append(line_number)
                PDB.COMPND[0][1].append(line_number)
            else:
                if PDB.COMPND[-1][1][1] == line_number - 1:
                    PDB.COMPND[-1][0].append(TEXT[10:80].strip())
                    PDB.COMPND[-1][1][1] = line_number
                else:
                    warnings.warn('Another COMPND section found at line {}.'.format(line_number))
                    PDB.COMPND.append([[],[]])
                    PDB.COMPND[-1][0].append(TEXT[10:80].strip())
                    PDB.COMPND[-1][1].append(line_number)
                    PDB.COMPND[-1][1].append(line_number)
        
        elif RECORD == 'SOURCE':#1xN
            TEXT = line
            if PDB.SOURCE == None:
                PDB.SOURCE = [[[],[]]]
                PDB.SOURCE[0][0].append(TEXT[10:79].strip())
                PDB.SOURCE[0][1].append(line_number)
                PDB.SOURCE[0][1].append(line_number)
            else:
                if PDB.SOURCE[-1][1][1] == line_number - 1:
                    PDB.SOURCE[-1][0].append(TEXT[10:79].strip())
                    PDB.SOURCE[-1][1][1] = line_number
                else:
                    warnings.warn('Another SOURCE section found at line {}.'.format(line_number))
                    PDB.SOURCE.append([[],[]])
                    PDB.SOURCE[-1][0].append(TEXT[10:79].strip())
                    PDB.SOURCE[-1][1].append(line_number)
                    PDB.SOURCE[-1][1].append(line_number)
        
        elif RECORD == 'KEYWDS':#1xN
            TEXT = line
            if PDB.KEYWDS == None:
                PDB.KEYWDS = [[[],[]]]
                PDB.KEYWDS[0][0].append(TEXT[10:79].strip())
                PDB.KEYWDS[0][1].append(line_number)
                PDB.KEYWDS[0][1].append(line_number)
            else:
                if PDB.KEYWDS[-1][1][1] == line_number - 1:
                    PDB.KEYWDS[-1][0].append(TEXT[10:79].strip())
                    PDB.KEYWDS[-1][1][1] = line_number
                else:
                    warnings.warn('Another KEYWDS section found at line {}.'.format(line_number))
                    PDB.KEYWDS.append([[],[]])
                    PDB.KEYWDS[-1][0].append(TEXT[10:79].strip())
                    PDB.KEYWDS[-1][1].append(line_number)
                    PDB.KEYWDS[-1][1].append(line_number)
        
        elif RECORD == 'EXPDTA':#1xN
            TEXT = line
            if PDB.EXPDTA == None:
                PDB.EXPDTA = [[[],[]]]
                PDB.EXPDTA[0][0].append(TEXT[10:79].strip())
                PDB.EXPDTA[0][1].append(line_number)
                PDB.EXPDTA[0][1].append(line_number)
            else:
                if PDB.EXPDTA[-1][1][1] == line_number - 1:
                    PDB.EXPDTA[-1][0].append(TEXT[10:79].strip())
                    PDB.EXPDTA[-1][1][1] = line_number
                else:
                    warnings.warn('Another EXPDTA section found at line {}.'.format(line_number))
                    PDB.EXPDTA.append([[],[]])
                    PDB.EXPDTA[-1][0].append(TEXT[10:79].strip())
                    PDB.EXPDTA[-1][1].append(line_number)
                    PDB.EXPDTA[-1][1].append(line_number)
        
        elif RECORD == 'NUMMDL':#1xN
            TEXT = line
            if PDB.NUMMDL == None:
                PDB.NUMMDL = [[[],[]]]
                PDB.NUMMDL[0][0].append(int(TEXT[10:14].strip()))
                PDB.NUMMDL[0][1].append(line_number)
                PDB.NUMMDL[0][1].append(line_number)
            else:
                if PDB.NUMMDL[-1][1][1] == line_number - 1:
                    PDB.NUMMDL[-1][0].append(int(TEXT[10:14].strip()))
                    PDB.NUMMDL[-1][1][1] = line_number
                else:
                    warnings.warn('Another NUMMDL section found at line {}.'.format(line_number))
                    PDB.NUMMDL.append([[],[]])
                    PDB.NUMMDL[-1][0].append(int(TEXT[10:14].strip()))
                    PDB.NUMMDL[-1][1].append(line_number)
                    PDB.NUMMDL[-1][1].append(line_number)
        
        elif RECORD == 'MDLTYP':#1xN
            TEXT = line
            if PDB.MDLTYP == None:
                PDB.MDLTYP = [[[],[]]]
                PDB.MDLTYP[0][0].append(TEXT[10:79].strip())
                PDB.MDLTYP[0][1].append(line_number)
                PDB.MDLTYP[0][1].append(line_number)
            else:
                if PDB.MDLTYP[-1][1][1] == line_number - 1:
                    PDB.MDLTYP[-1][0].append(TEXT[10:79].strip())
                    PDB.MDLTYP[-1][1][1] = line_number
                else:
                    warnings.warn('Another MDLTYP section found at line {}.'.format(line_number))
                    PDB.MDLTYP.append([[],[]])
                    PDB.MDLTYP[-1][0].append(TEXT[10:79].strip())
                    PDB.MDLTYP[-1][1].append(line_number)
                    PDB.MDLTYP[-1][1].append(line_number)
        
        elif RECORD == 'AUTHOR':#1xN
            TEXT = line
            if PDB.AUTHOR == None:
                PDB.AUTHOR = [[[],[]]]
                PDB.AUTHOR[0][0].append(TEXT[10:79].strip())
                PDB.AUTHOR[0][1].append(line_number)
                PDB.AUTHOR[0][1].append(line_number)
            else:
                if PDB.AUTHOR[-1][1][1] == line_number - 1:
                    PDB.AUTHOR[-1][0].append(TEXT[10:79].strip())
                    PDB.AUTHOR[-1][1][1] = line_number
                else:
                    warnings.warn('Another AUTHOR section found at line {}.'.format(line_number))
                    PDB.AUTHOR.append([[],[]])
                    PDB.AUTHOR[-1][0].append(TEXT[10:79].strip())
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
            TEXT = line
            if PDB.JRNL == None:
                PDB.JRNL = [[[],[]]]
                PDB.JRNL[0][0].append(TEXT[12:79].strip())
                PDB.JRNL[0][1].append(line_number)
                PDB.JRNL[0][1].append(line_number)
            else:
                if PDB.JRNL[-1][1][1] == line_number - 1:
                    PDB.JRNL[-1][0].append(TEXT[12:79].strip())
                    PDB.JRNL[-1][1][1] = line_number
                else:
                    warnings.warn('Another JRNL section found at line {}.'.format(line_number))
                    PDB.JRNL.append([[],[]])
                    PDB.JRNL[-1][0].append(TEXT[12:79].strip())
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
            TEXT = line
            if PDB.DBREF == None:
                PDB.DBREF  = [[],[]]
                PDB.DBREF[0].append(TEXT)
                PDB.DBREF[1].append(line_number)
            else:
                PDB.DBREF[0].append(TEXT)
                PDB.DBREF[1].append(line_number)
        
        elif RECORD == 'DBREF1':#Mx1
            TEXT = line
            if PDB.DBREF1 == None:
                PDB.DBREF1  = [[],[]]
                PDB.DBREF1[0].append(TEXT)
                PDB.DBREF1[1].append(line_number)
            else:
                PDB.DBREF1[0].append(TEXT)
                PDB.DBREF1[1].append(line_number)
        
        elif RECORD == 'DBREF2':#Mx1
            TEXT = line
            if PDB.DBREF2 == None:
                PDB.DBREF2  = [[],[]]
                PDB.DBREF2[0].append(TEXT)
                PDB.DBREF2[1].append(line_number)
            else:
                PDB.DBREF2[0].append(TEXT)
                PDB.DBREF2[1].append(line_number)
        
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
            TEXT = line
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
            TEXT = line
            if PDB.HET == None:
                PDB.HET  = [[],[]]
                PDB.HET[0].append(TEXT)
                PDB.HET[1].append(line_number)
            else:
                PDB.HET[0].append(TEXT)
                PDB.HET[1].append(line_number)
        
        elif RECORD == 'HETNAM':#MxN
            TEXT = line
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
            TEXT = line
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
            TEXT = line
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
            TEXT = line
            if PDB.HELIX == None:
                PDB.HELIX  = [[],[]]
                PDB.HELIX[0].append(TEXT)
                PDB.HELIX[1].append(line_number)
            else:
                PDB.HELIX[0].append(TEXT)
                PDB.HELIX[1].append(line_number)
        
        elif RECORD == 'SHEET':#Mx1
            TEXT = line
            if PDB.SHEET == None:
                PDB.SHEET  = [[],[]]
                PDB.SHEET[0].append(TEXT)
                PDB.SHEET[1].append(line_number)
            else:
                PDB.SHEET[0].append(TEXT)
                PDB.SHEET[1].append(line_number)
        
        elif RECORD == 'SSBOND':#Mx1
            TEXT = line
            if PDB.SSBOND == None:
                PDB.SSBOND  = [[],[]]
                PDB.SSBOND[0].append(TEXT)
                PDB.SSBOND[1].append(line_number)
            else:
                PDB.SSBOND[0].append(TEXT)
                PDB.SSBOND[1].append(line_number)
        
        elif RECORD == 'LINK':#Mx1
            TEXT = line
            if PDB.LINK == None:
                PDB.LINK  = [[],[]]
                PDB.LINK[0].append(TEXT)
                PDB.LINK[1].append(line_number)
            else:
                PDB.LINK[0].append(TEXT)
                PDB.LINK[1].append(line_number)
        
        elif RECORD == 'CISPEP':#Mx1
            TEXT = line
            if PDB.CISPEP == None:
                PDB.CISPEP  = [[],[]]
                PDB.CISPEP[0].append(TEXT)
                PDB.CISPEP[1].append(line_number)
            else:
                PDB.CISPEP[0].append(TEXT)
                PDB.CISPEP[1].append(line_number)
        
        elif RECORD == 'SITE':#MxN
            TEXT = line
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
            TEXT = line
            if PDB.TER == None:
                PDB.TER  = [[],[]]
                PDB.TER[0].append(TEXT)
                PDB.TER[1].append(line_number)
            else:
                PDB.TER[0].append(TEXT)
                PDB.TER[1].append(line_number)    
        
        elif RECORD == 'HETATM':#Mx1
            TEXT = line
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
    if PDB.COMPND == None:
        warnings.warn('COMPND entry is missing!')
        sys.exit(0)
    else:
        PDB.COMPND = PDB.COMPND[0]
    
    for line in PDB.COMPND[0]:
        if line.split()[0] == 'CHAIN:':
            if PDB.CHAINS == None:
                chain = []
                for A in line[6:].split():
                    chain.append(A[0])
                PDB.CHAINS = chain
            else:
                chain = []
                for A in line[6:].split():
                    chain.append(A[0])
                PDB.CHAINS.extend(chain)
    
    if PDB.DBREF == None:
        pass
    else:
        DBREF = []
        for line in PDB.DBREF[0]:
            DBREF.append(Dbref([0, line, PDB.DBREF[1][PDB.DBREF[0].index(line)]]))
        PDB.DBREF = DBREF

    if PDB.DBREF1 == None and PDB.DBREF2 == None:
        pass
    else:
        DBREF = []
        for line1, line2 in zip(PDB.DBREF1[0], PDB.DBREF2[0]):
            DBREF.append(Dbref([1, line1, line2, PDB.DBREF1[1][PDB.DBREF1[0].index(line1)], 
                                PDB.DBREF2[1][PDB.DBREF2[0].index(line2)]]))
        PDB.DBREF = DBREF
        PDB.DBREF1= 'Done'
        PDB.DBREF2= 'Done'
    
    if PDB.DBREF == None:
        warnings.warn('DBREF entry is missing!')
    
    if len(PDB.DBREF) == len(PDB.CHAINS):
        pass
    else:
        warnings.warn('Different number of DBREF entries than COMPND CHAIN {} != {}.'.format(len(PDB.DBREF), PDB.NUMMDL))    
    
    if PDB.SEQRES == None:
        warnings.warn('SEQRES entry is missing!')
    else:
        SEQRES = {}
        SEQRESlen = {}
        for seq in PDB.SEQRES[0][0]:
            this_line = []
            
            if seq[19:22].strip() == '':
                pass
            else:
                this_line.append(seq[19:22].strip())
            
            if seq[23:26].strip() == '':
                pass
            else:
                this_line.append(seq[23:26].strip())
            
            if seq[27:30].strip() == '':
                pass
            else:
                this_line.append(seq[27:30].strip())
            
            if seq[31:34].strip() == '':
                pass
            else:
                this_line.append(seq[31:34].strip())
            
            if seq[35:38].strip() == '':
                pass
            else:
                this_line.append(seq[35:38].strip())
            
            if seq[39:42].strip() == '':
                pass
            else:
                this_line.append(seq[39:42].strip())
            
            if seq[43:46].strip() == '':
                pass
            else:
                this_line.append(seq[43:46].strip())
            
            if seq[47:50].strip() == '':
                pass
            else:
                this_line.append(seq[47:50].strip())
            
            if seq[51:54].strip() == '':
                pass
            else:
                this_line.append(seq[51:54].strip())
            
            if seq[55:58].strip() == '':
                pass
            else:
                this_line.append(seq[55:58].strip())
            
            if seq[59:62].strip() == '':
                pass
            else:
                this_line.append(seq[59:62].strip())
            
            if seq[63:66].strip() == '':
                pass
            else:
                this_line.append(seq[63:66].strip())
            
            if seq[67:70].strip() == '':
                pass
            else:
                this_line.append(seq[67:70].strip())
            
            
            if seq[11:12].strip() in SEQRES:
                Value = SEQRES[seq[11:12].strip()]
                Value.extend(this_line)
                SEQRES[seq[11:12].strip()] = Value
            else:
                SEQRES[seq[11:12].strip()] = this_line
                SEQRESlen[seq[11:12].strip()] = int(seq[13:17].strip())
        
        for key, value in SEQRES.items():
            if len(value) == SEQRESlen[key]:
                pass
            else:
                warnings.warn('Sequence length of chain {} is not equal the expected. {} != {}'.format(key, len(value), SEQRESlen[key]))
        
        PDB.SEQRES = SEQRES
        PDB.SEQRESlen = SEQRESlen

    if PDB.HET == None:
        pass
    else:
        HET = []
        for line in PDB.HET[0]:
            HET.append(Het(line, PDB.HET[1][PDB.HET[0].index(line)]))
        PDB.HET = HET
    
    if PDB.HETNAM == None:
        pass
    else:
        PDB.HETNAM = PDB.HETNAM[0]
        HETNAM = {}
        for line in PDB.HETNAM[0]:
            if line[11:14].strip() in HETNAM:
                HETNAM[line[11:14].strip()] = HETNAM[line[11:14].strip()] + line[15:70].strip()
            else:
                HETNAM[line[11:14].strip()] = line[15:70].strip()
        PDB.HETNAM = HETNAM
    
    if PDB.HETSYN == None:
        pass
    else:
        PDB.HETSYN = PDB.HETSYN[0]
        HETSYN = {}
        for line in PDB.HETSYN[0]:
            if line[11:14].strip() in HETSYN:
                HETSYN[line[11:14].strip()] = HETSYN[line[11:14].strip()] + line[15:70].strip()
            else:
                HETSYN[line[11:14].strip()] = line[15:70].strip()
        PDB.HETSYN = HETSYN
    
    if PDB.FORMUL == None:
        pass
    else:
        PDB.FORMUL = PDB.FORMUL[0]
        FORMUL = {}
        for line in PDB.FORMUL[0]:
            if line[12:15].strip() in FORMUL:
                FORMUL[line[12:15].strip()] = FORMUL[line[12:15].strip()] + line[19:70].strip()
            else:
                FORMUL[line[12:15].strip()] = line[19:70].strip()
        PDB.FORMUL = FORMUL
    
    if PDB.HELIX == None:
        pass
    else:
        HELIX = []
        for line in PDB.HELIX[0]:
            HELIX.append(Helix(line, PDB.HELIX[1][PDB.HELIX[0].index(line)]))
        PDB.HELIX = HELIX

    if PDB.SHEET == None:
        pass
    else:
        SHEET = []
        for line in PDB.SHEET[0]:
            SHEET.append(Sheet(line, PDB.SHEET[1][PDB.SHEET[0].index(line)]))
        PDB.SHEET = SHEET
    
    if PDB.SSBOND == None:
        pass
    else:
        SSBOND = []
        for line in PDB.SSBOND[0]:
            SSBOND.append(SSbond(line, PDB.SSBOND[1][PDB.SSBOND[0].index(line)]))
        PDB.SSBOND = SSBOND
    
    if PDB.LINK == None:
        pass
    else:
        LINK = []
        for line in PDB.LINK[0]:
            LINK.append(Link(line, PDB.LINK[1][PDB.LINK[0].index(line)]))
        PDB.LINK = LINK
    
    if PDB.CISPEP == None:
        pass
    else:
        CISPEP = []
        for line in PDB.CISPEP[0]:
            CISPEP.append(Cispep(line, PDB.CISPEP[1][PDB.CISPEP[0].index(line)]))
        PDB.CISPEP = CISPEP
    
    if PDB.SITE == None:
        pass
    else:
        PDB.SITE = define_sites(PDB.SITE[0])
    
    if PDB.ATOM == None:
        warnings.warn('ATOM section is missing!')
        sys.exit(0)
    else:
        ATOM = []
        for line in PDB.ATOM[0]:
            ATOM.append(Atom(line, PDB.ATOM[1][PDB.ATOM[0].index(line)]))
        PDB.ATOM = ATOM

    if PDB.TER == None:
        warnings.warn('TER section is missing!')
        sys.exit(0)
    else:
        TER = []
        for line in PDB.TER[0]:
            TER.append(Ter(line))
        PDB.TER = TER
    
    if PDB.HETATM == None:
        #warnings.warn('HETATM section is missing!')
        pass
    else:
        HETATM = []
        for line in PDB.HETATM[0]:
            HETATM.append(Hetatm(line, PDB.HETATM[1][PDB.HETATM[0].index(line)]))
        PDB.HETATM = HETATM

    return(PDB)



##############################################################################
#MAIN#########################################################################
##############################################################################

if __name__ == '__main__':
    print('PyDB is a package that is designed to read PDB file format.')
    File_1a0k = readPDBfile(ID = '1a0k', Extension = 'pdb')
    PDB_1a0k = parsePDBfile(File_1a0k)
    print(PDB_1a0k.JUNK)
#    '''
    File_1ao6 = readPDBfile(ID = '1ao6', Extension = 'pdb')
    PDB_1ao6 = parsePDBfile(File_1ao6)
    File_5b5t = readPDBfile(ID = '5b5t', Extension = 'pdb')
    PDB_5b5t = parsePDBfile(File_5b5t)
    print(PDB_1a0k.JUNK, PDB_1ao6.JUNK, PDB_5b5t.JUNK)
#    '''
    
    
    