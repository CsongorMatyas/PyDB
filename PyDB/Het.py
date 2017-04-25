import sys, warnings

class Het:    #Each HET within a PDB object will be a Het object
    def __init__(self, TEXT):
        
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