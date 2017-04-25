import sys, warnings

from PyDB.Dbref         import Dbref
from PyDB.Het           import Het
from PyDB.Helix         import Helix
from PyDB.Sheet         import Sheet
from PyDB.SSbond        import SSbond
from PyDB.Link          import Link
from PyDB.Cispep        import Cispep
from PyDB.Atom          import Atom
from PyDB.Ter           import Ter
from PyDB.Hetatm        import Hetatm

def processPDB_SITE(PDB_SITES):
    SITE = {}        
    counter = 0
    for PDB_SITE in PDB_SITES:
        for TEXT in PDB_SITE:
            if int(TEXT[7:10].strip()) == 1:
                SITE[counter] = [TEXT[11:14].strip(), int(TEXT[15:17].strip()), []]
                
                residue = []
                if TEXT[18:21].strip() != '':
                    residue.append([TEXT[18:21].strip(), TEXT[22:23].strip(),
                                    int(TEXT[23:27].strip()), TEXT[27:28].strip()])
                                    
                if TEXT[29:32].strip() != '':
                    residue.append([TEXT[29:32].strip(), TEXT[33:34].strip(),
                                    int(TEXT[34:38].strip()), TEXT[38:39].strip()])
                                    
                if TEXT[40:43].strip() != '':
                    residue.append([TEXT[40:43].strip(), TEXT[44:45].strip(),
                                    int(TEXT[45:49].strip()), TEXT[49:50].strip()])
                                    
                if TEXT[51:54].strip() != '':
                    residue.append([TEXT[51:54].strip(), TEXT[55:56].strip(),
                                    int(TEXT[56:60].strip()), TEXT[60:61].strip()])
                SITE[counter][2].extend(residue)
                counter += 1
            else:
                residue = []
                if TEXT[18:21].strip() != '':
                    residue.append([TEXT[18:21].strip(), TEXT[22:23].strip(),
                                    int(TEXT[23:27].strip()), TEXT[27:28].strip()])
                                    
                if TEXT[29:32].strip() != '':
                    residue.append([TEXT[29:32].strip(), TEXT[33:34].strip(),
                                    int(TEXT[34:38].strip()), TEXT[38:39].strip()])
                                    
                if TEXT[40:43].strip() != '':
                    residue.append([TEXT[40:43].strip(), TEXT[44:45].strip(),
                                    int(TEXT[45:49].strip()), TEXT[49:50].strip()])
                                    
                if TEXT[51:54].strip() != '':
                    residue.append([TEXT[51:54].strip(), TEXT[55:56].strip(),
                                    int(TEXT[56:60].strip()), TEXT[60:61].strip()])
                
                SITE[counter-1][2].extend(residue)
            
    return(SITE)

def processPDB(self):
    if self.COMPND == None:
        warnings.warn('COMPND entry is missing!')
        sys.exit(0)

    for line in self.COMPND:
        if line.split()[0] == 'CHAIN:':
            if self.CHAINS == None:
                chain = []
                for A in line[6:].split():
                    chain.append(A[0])
                self.CHAINS = chain
            else:
                chain = []
                for A in line[6:].split():
                    chain.append(A[0])
                self.CHAINS.extend(chain)

    if self.DBREF != None:
        DBREF = []
        for line in self.DBREF:
            DBREF.append(Dbref([0, line]))
        self.DBREF = DBREF

    if self.DBREF1 == None and self.DBREF2 == None:
        pass
    else:
        DBREF = []
        for line1, line2 in zip(self.DBREF1, self.DBREF2):
            DBREF.append(Dbref([1, line1, line2]))
        self.DBREF = DBREF
        self.DBREF1= 'Done'
        self.DBREF2= 'Done'

    if self.DBREF == None:
        warnings.warn('DBREF entry is missing!')
    
    if len(self.DBREF) != len(self.CHAINS):
        warnings.warn('Different number of DBREF entries than COMPND CHAIN {} != {}.'.format(len(self.DBREF), self.NUMMDL))    

    if self.SEQRES == None:
        warnings.warn('SEQRES entry is missing!')
    else:
        SEQRES = {}
        SEQRESlen = {}
        for seq in self.SEQRES[0]:
            this_line = []
            
            if seq[19:22].strip() == '':
                pass
            else:
                this_line.append(seq[19:22].strip())
            
            if seq[23:26].strip() == '':
                pass
            else:
                this_line.append(seq[23:26].strip())
            
            if seq[27:30].strip() == '':
                pass
            else:
                this_line.append(seq[27:30].strip())
            
            if seq[31:34].strip() == '':
                pass
            else:
                this_line.append(seq[31:34].strip())
            
            if seq[35:38].strip() == '':
                pass
            else:
                this_line.append(seq[35:38].strip())
            
            if seq[39:42].strip() == '':
                pass
            else:
                this_line.append(seq[39:42].strip())
            
            if seq[43:46].strip() == '':
                pass
            else:
                this_line.append(seq[43:46].strip())
            
            if seq[47:50].strip() == '':
                pass
            else:
                this_line.append(seq[47:50].strip())
            
            if seq[51:54].strip() == '':
                pass
            else:
                this_line.append(seq[51:54].strip())
            
            if seq[55:58].strip() == '':
                pass
            else:
                this_line.append(seq[55:58].strip())
            
            if seq[59:62].strip() == '':
                pass
            else:
                this_line.append(seq[59:62].strip())
            
            if seq[63:66].strip() == '':
                pass
            else:
                this_line.append(seq[63:66].strip())
            
            if seq[67:70].strip() == '':
                pass
            else:
                this_line.append(seq[67:70].strip())
            
            
            if seq[11:12].strip() in SEQRES:
                Value = SEQRES[seq[11:12].strip()]
                Value.extend(this_line)
                SEQRES[seq[11:12].strip()] = Value
            else:
                SEQRES[seq[11:12].strip()] = this_line
                SEQRESlen[seq[11:12].strip()] = int(seq[13:17].strip())
        
        for key, value in SEQRES.items():
            if len(value) == SEQRESlen[key]:
                pass
            else:
                warnings.warn('Sequence length of chain {} is not equal the expected. {} != {}'.format(key, len(value), SEQRESlen[key]))
        
        self.SEQRES = SEQRES
        self.SEQRESlen = SEQRESlen

    if self.HET == None:
        pass
    else:
        HET = []
        for line in self.HET:
            HET.append(Het(line))
        self.HET = HET
    
    if self.HETNAM != None:
        HETNAM = {}
        for LIST in self.HETNAM:
            for line in LIST:
                if line[11:14].strip() in HETNAM:
                    HETNAM[line[11:14].strip()] = HETNAM[line[11:14].strip()] + line[15:70].strip()
                else:
                    HETNAM[line[11:14].strip()] = line[15:70].strip()
            self.HETNAM = HETNAM

    if self.HETSYN != None:
        HETSYN = {}
        for LIST in self.HETSYN:
            for line in LIST:
                if line[11:14].strip() in HETSYN:
                    HETSYN[line[11:14].strip()] = HETSYN[line[11:14].strip()] + line[15:70].strip()
                else:
                    HETSYN[line[11:14].strip()] = line[15:70].strip()
            self.HETSYN = HETSYN

    if self.FORMUL != None:
        FORMUL = {}
        for LIST in self.FORMUL:
            for line in LIST:
                if line[12:15].strip() in FORMUL:
                    FORMUL[line[12:15].strip()] = FORMUL[line[12:15].strip()] + line[19:70].strip()
                else:
                    FORMUL[line[12:15].strip()] = line[19:70].strip()
        self.FORMUL = FORMUL
    
    if self.HELIX == None:
        pass
    else:
        HELIX = []
        for line in self.HELIX:
            HELIX.append(Helix(line))
        self.HELIX = HELIX

    if self.SHEET == None:
        pass
    else:
        SHEET = []
        for line in self.SHEET:
            SHEET.append(Sheet(line))
        self.SHEET = SHEET

    if self.SSBOND == None:
        pass
    else:
        SSBOND = []
        for line in self.SSBOND:
            SSBOND.append(SSbond(line))
        self.SSBOND = SSBOND

    if self.LINK == None:
        pass
    else:
        LINK = []
        for line in self.LINK:
            LINK.append(Link(line))
        self.LINK = LINK

    if self.CISPEP == None:
        pass
    else:
        CISPEP = []
        for line in self.CISPEP:
            CISPEP.append(Cispep(line))
        self.CISPEP = CISPEP

    if self.SITE == None:
        pass
    else:
        self.SITE = processPDB_SITE(self.SITE)

    if self.ATOM == None:
        warnings.warn('ATOM section is missing!')
        sys.exit(0)
    else:
        ATOM = []
        for line in self.ATOM:
            ATOM.append(Atom(line))
        self.ATOM = ATOM

    if self.TER == None:
        warnings.warn('TER section is missing!')
        sys.exit(0)
    else:
        TER = []
        for line in self.TER:
            TER.append(Ter(line))
        self.TER = TER

    if self.HETATM == None:
        #warnings.warn('HETATM section is missing!')
        pass
    else:
        HETATM = []
        for line in self.HETATM:
            HETATM.append(Hetatm(line))
        self.HETATM = HETATM