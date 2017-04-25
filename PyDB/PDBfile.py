import sys, warnings, urllib.request, pickle
from datetime import datetime

from PyDB.readPDB       import readPDB
from PyDB.readFromFile  import readFromFile
from PyDB.readFromUrl   import readFromUrl
from PyDB.parsePDB      import parsePDB, parse1xN, parse1x1, parseMx1, parseMxN
from PyDB.processPDB    import processPDB, processPDB_SITE

class PDBfile(): #PDF file class, each PDB file will be a PDB object
    def __init__(self, ID = None, Url = None, File = None):
        self.HEADER = None   #Mandatory 1x1 - First line of the entry, contains PDB ID code,
                                           #classification, and date of deposition.
        self.depDate= None   #Date of deposition - from HEADER line 
        self.idCode = None   #Unique PDB ID code - from HEADER line
        self.OBSLTE = None   #Optional  1xN - Statement that the entry has been removed from
                             #distribution and list of the ID code(s) which replaced it.
        self.TITLE  = None   #Mandatory 1xN - Description of the experiment represented in the entry.
        self.SPLIT  = None   #Optional  1xN - List of PDB entries that compose a larger
                             #macromolecular complexes.
        self.CAVEAT = None   #Optional  1xN - Severe error indicator.    
        self.COMPND = None   #Mandatory 1xN - Description of macromolecular contents of the entry.
        self.CHAINS = None   #                From COMPND, contains list of chains
        self.SOURCE = None   #Mandatory 1xN - Biological source of macromolecules in the entry.
        self.KEYWDS = None   #Mandatory 1xN - List of keywords describing the macromolecule.
        self.EXPDTA = None   #Mandatory 1xN - Experimental technique used for the structure determination.
        self.NUMMDL = None   #Optional  1x1 - Number of models.
        self.MDLTYP = None   #Optional  1xN - Contains additional annotation pertinent
                             #to the coordinates presented in the entry.
        self.AUTHOR = None   #Mandatory 1xN - List of contributors.
        self.REVDAT = None   #Mandatory Mx1 - Revision date and related information.
        self.SPRSDE = None   #Optional  1xN - List of entries obsoleted from public
                             #release and replaced by current entry.
        self.JRNL   = None   #Mandatory 1xN - Literature citation that defines the coordinate set.
        self.REMARK_0 = None #Optional  1xN - Re-refinement notice
        self.REMARK_1 = None #Optional  1xN - Related publications
        self.REMARK_2 = None #Mandatory 1xN - Resolution
        self.REMARK_3 = None #Mandatory 1xN - Final refinement information
        self.REMARK_4 = None #Optional  1xN - Format
        self.REMARK_5 = None #Optional  1xN - Obsolete Statement
        self.REMARK_6 = None #Optional  1xN - 6-99 free text annotation
        self.REMARK_100=None #Optional  1xN - Deposition or Processing Site
        self.REMARK_200=None #Optional  1xN - X-ray Diffraction Experimental Details
        self.REMARK_205=None #Optional  1xN - Fiber Diffraction, Fiber Sample Experiment Details
        self.REMARK_210=None #Optional  1xN - +215/217  NMR Experiment Details
        self.REMARK_230=None #Optional  1xN - Neutron Diffraction Experiment Details
        self.REMARK_240=None #Optional  1xN - Electron Crystallography Experiment Details
        self.REMARK_245=None #Optional  1xN - Electron Microscopy Experiment Details
        self.REMARK_247=None #Optional  1xN - Electron Microscopy details
        self.REMARK_250=None #Optional  1xN - Other Type of Experiment Details
        self.REMARK_265=None #Optional  1xN - Solution Scattering Experiment Details
        self.REMARK_280=None #Optional  1xN - Crystal
        self.REMARK_285=None #Optional  1xN - CRYST1
        self.REMARK_290=None #Optional  1xN - Crystallographic Symmetry
        self.REMARK_300=None #Optional  1xN - Biomolecule
        self.REMARK_350=None #Optional  1xN - Generating the Biomolecule
        self.REMARK_375=None #Optional  1xN - Special Position
        self.REMARK_400=None #Optional  1xN - Compound
        self.REMARK_450=None #Optional  1xN - Source
        self.REMARK_465=None #Optional  1xN - Missing residues
        self.REMARK_470=None #Optional  1xN - Missing Atom(s)
        self.REMARK_475=None #Optional  1xN - Residues modeled with zero occupancy
        self.REMARK_480=None #Optional  1xN - Polymer atoms modeled with zero occupancy
        self.REMARK_500=None #Optional  1xN - Geometry and Stereochemistry
        self.REMARK_525=None #Optional  1xN - Distant Solvent Atoms
        self.REMARK_600=None #Optional  1xN - Heterogen
        self.REMARK_610=None #Optional  1xN - Non-polymer residues with missing atoms
        self.REMARK_615=None #Optional  1xN - Non-polymer residues containing atoms with zero occupancy
        self.REMARK_620=None #Optional  1xN - Metal coordination
        self.REMARK_630=None #Optional  1xN - Inhibitor Description
        self.REMARK_650=None #Optional  1xN - Helix
        self.REMARK_700=None #Optional  1xN - Sheet
        self.REMARK_800=None #Optional  1xN - Important Sites
        self.REMARK_999=None #Optional  1xN - Sequence
        self.DBREF  = None   #Optional  Mx1 - Reference to the entry in the sequence database(s).
        self.DBREF1 = None   #Optional  Mx1 - Same as DBREF but new version, goes with DBREF2
        self.DBREF2 = None   #Optional  Mx1 - goes with DBREF1
        self.SEQADV = None   #Optional  Mx1 - Identification of conflicts between PDB
                             #and the named sequence database.
        self.SEQRES = None   #Mandatory MxN - Primary sequence of backbone residues. As dictionary
        self.SEQRESlen=None  #                Number of amino acids in each chain as dictionary
        self.MODRES = None   #Optional  Mx1 - Identification of modifications to standard residues.
        self.HET    = None   #Optional  Mx1 - Identification of non-standard groups (heterogens).
        self.HETNAM = None   #Optional  MxN - Compound name of the heterogens. (JUPAC?)
        self.HETSYN = None   #Optional  MxN - Synonymous compound names for heterogens. (nonJUPAC?)
        self.FORMUL = None   #Mandatory MxN - Chemical formula of non-standard groups.
        self.HELIX  = None   #Optional  Mx1 - Identification of helical substructures.
        self.SHEET  = None   #Optional  Mx1 - Identification of sheet substructures.
        self.SSBOND = None   #Optional  Mx1 - Identification of disulfide bonds.
        self.LINK   = None   #Optional  Mx1 - Identification of inter-residue bonds.
        self.CISPEP = None   #Optional  Mx1 - Identification of peptide residues in cis conformation.
        self.SITE   = None   #Optional  MxN - Identification of groups comprising important entity sites.
        self.CRYST1 = None   #Mandatory 1x1 - Unit cell parameters, space group, and Z.
        self.ORIGX1 = None   #Mandatory 1x1 - Transformation from orthogonal coordinates
        self.ORIGX2 = None   #Mandatory 1x1 - to the submitted coordinates.
        self.ORIGX3 = None   #Mandatory 1x1
        self.SCALE1 = None   #Mandatory 1x1 - Transformation from orthogonal coordinates
        self.SCALE2 = None   #Mandatory 1x1 - to fractional crystallographic coordinates.
        self.SCALE3 = None   #Mandatory 1x1
        self.MTRIX1 = None   #Optional  Mx1 - Transformations expressing non-crystallographic
        self.MTRIX2 = None   #Optional  Mx1   symmetry.
        self.MTRIX3 = None   #Optional  Mx1
        self.MODEL  = None   #Optional ?MxN - Specification of model number for multiple
                             #structures in a single coordinate entry.
        self.ATOM   = None   #Mandatory Mx1 - Atomic coordinate records for standard groups.
        self.ANISOU = None   #Optional  Mx1 - Anisotropic temperature factors.
        self.TER    = None   #Mandatory Mx1 - Chain terminator.
        self.HETATM = None   #Optional  Mx1 - Atomic coordinate records for heterogens.
        self.ENDMDL = None   #Optional ?MxN - End-of-model record for multiple structures
                             #in a single coordinate entry.
        self.CONECT = None   #Optional  Mx1 - Connectivity records.
        self.MASTER = None   #Mandatory 1x1 - Control record for bookkeeping.
        self.END    = None   #Mandatory 1x1 - Last record in the file.
        self.JUNK   = None   #PyDB property, if it's not None, there were lines starting
                             #with a keyword other than the ones allowed by PDB file format.
        self.LINENUM = []    #PyDB property, it contains a dictionary where the key is the line number
                             #and the value is that lines record name (like HEADER or TITLE)
        self.LINENUMINV = {} #Inverse of LINENUM, dictionary that contains records as key
                             #and array of line numbers as value

        if ID != None:
            FILE = readPDB(ID)
            self.FSOURCE = ID
        elif Url != None:
            FILE = readFromUrl(Url)
            self.FSOURCE = Url
        elif File != None:
            FILE = readFromFile(File)
            self.FSOURCE = File
        else:
            self.FSOURCE = None   #SOURCE will be either the ID code, the Url to the file or the filename of a file on the disk
            warnings.warn("No file source was given. Please use argument: ID='' or Url='' or File=''")
            sys.exit(0)

        parsePDB(self, FILE)
        processPDB(self)
