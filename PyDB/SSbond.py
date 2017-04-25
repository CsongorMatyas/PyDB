import sys, warnings

class SSbond:    #Each SSBOND within a PDB object will be an SSbond object
    def __init__(self, TEXT):
        
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