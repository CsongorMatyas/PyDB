import sys, warnings

class Cispep:    #Each CISPEP within a PDB object will be an Cispep object
    def __init__(self, TEXT):
        
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