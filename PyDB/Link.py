import sys, warnings

class Link:    #Each LINK within a PDB object will be an Link object
    def __init__(self, TEXT):
        
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