import sys, warnings

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