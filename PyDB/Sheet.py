import sys, warnings

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