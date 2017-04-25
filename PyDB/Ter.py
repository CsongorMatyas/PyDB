import sys, warnings

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