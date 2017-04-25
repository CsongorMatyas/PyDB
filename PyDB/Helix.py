import sys, warnings

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
