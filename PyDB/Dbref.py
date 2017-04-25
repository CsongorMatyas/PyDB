import sys, warnings

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