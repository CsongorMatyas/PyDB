import sys, warnings
from datetime import datetime

def parse1x1(obj, text, name, line_number):
    if obj == None:
        obj = text
    else:
        warnings.warn('A second {} section found at line {}!!!'.format(name, line_number))
        sys.exit(0)
    return(obj)

def parse1xN(self, obj, text, name, line_number):
    if obj == None:
        obj = []
    else:
        if self.LINENUM[-1] != name:
            warnings.warn('A second {} section found at line {}.'.format(name, line_number))
    if text != '':
        obj.append(text)
    return(obj)

def parseMx1(obj, text):
    if obj == None:
        obj = []
    if text != '':
        obj.append(text)
    return(obj)

def parseMxN(self, obj, text, name):
    if text != '':
        if obj == None:
            obj = [[]]
        else:
            if self.LINENUM[-1] != name:
                obj[-1].append([])
        obj[-1].append(text)
    return(obj)

def parsePDB(self, FILE):
    line_number = 0
    
    for line in FILE:
        if   len(line) < 7:
            RECORD = line.strip()
        else:
            RECORD = line[0:7].strip()
        
        self.LINENUM.append(RECORD)
        if RECORD in self.LINENUMINV:
            self.LINENUMINV[RECORD].append(line_number)
        else:
            self.LINENUMINV[RECORD] = [line_number]

        if RECORD == '':
            warnings.warn('Empty RECORD entry at line {}!!!'.format(line_number))
        
        elif RECORD == 'HEADER':#1x1
            TEXT = line
            if self.HEADER == None:
                self.HEADER = TEXT[10:50].strip()
                self.depDate = datetime.strptime(TEXT[50:59].strip(), '%d-%b-%y')
                self.idCode = TEXT[62:66].strip()
            else:
                warnings.warn('A second {} section found at line {}!!!'.format(RECORD, line_number))
                sys.exit(0)
        
        elif RECORD == 'OBSLTE':#1xN
            self.OBSLTE = parse1xN(self, self.OBSLTE, line[7:-1].strip(), RECORD, line_number)        

        elif RECORD == 'TITLE':#1xN
            self.TITLE  = parse1xN(self, self.TITLE,  line[10:-1].strip(), RECORD, line_number)

        elif RECORD == 'SPLIT':#1xN
            self.SPLIT  = parse1xN(self, self.SPLIT,  line[7:-1].strip(), RECORD, line_number)
        
        elif RECORD == 'CAVEAT':#1xN
            self.CAVEAT = parse1xN(self, self.CAVEAT, line[7:-1].strip(), RECORD, line_number)
        
        elif RECORD == 'COMPND':#1xN
            self.COMPND = parse1xN(self, self.COMPND, line[10:-1].strip(), RECORD, line_number)
                    
        elif RECORD == 'SOURCE':#1xN
            self.SOURCE = parse1xN(self, self.SOURCE, line[10:-1].strip(), RECORD, line_number)
        
        elif RECORD == 'KEYWDS':#1xN
            self.KEYWDS = parse1xN(self, self.KEYWDS, line[7:-1].strip(), RECORD, line_number)
        
        elif RECORD == 'EXPDTA':#1xN
            self.EXPDTA = parse1xN(self, self.EXPDTA, line[7:-1].strip(), RECORD, line_number)
        
        elif RECORD == 'NUMMDL':#1xN
            self.NUMMDL = parse1xN(self, self.NUMMDL, line[7:-1].strip(), RECORD, line_number)
        
        elif RECORD == 'MDLTYP':#1xN
            self.MDLTYP = parse1xN(self, self.MDLTYP, line[7:-1].strip(), RECORD, line_number)
        
        elif RECORD == 'AUTHOR':#1xN
            self.AUTHOR = parse1xN(self, self.AUTHOR, line[7:-1].strip(), RECORD, line_number)
        
        elif RECORD == 'REVDAT':#Mx1
            self.REVDAT = parseMx1(self.REVDAT, line[10:-1].strip())
                
        elif RECORD == 'SPRSDE':#1xN
            self.SPRSDE = parse1xN(self, self.SPRSDE, line[7:-1].strip(), RECORD, line_number)
        
        elif RECORD == 'JRNL':#1xN
            self.JRNL   = parse1xN(self, self.JRNL,   line[7:-1].strip(), RECORD, line_number)
        
        elif RECORD == 'REMARK':
            RNUM = int(line[7:10].strip())
            if RNUM == 0:
                self.REMARK_0 = parse1xN(self, self.REMARK_0, line[10:-1].strip(), RECORD, line_number)                        
            elif RNUM == 1:
                self.REMARK_1 = parse1xN(self, self.REMARK_1, line[10:-1].strip(), RECORD, line_number)
            elif RNUM == 2:
                self.REMARK_2 = parse1xN(self, self.REMARK_2, line[10:-1].strip(), RECORD, line_number)
            elif RNUM == 3:
                self.REMARK_3 = parse1xN(self, self.REMARK_3, line[10:-1].strip(), RECORD, line_number)
            elif RNUM == 4:
                self.REMARK_4 = parse1xN(self, self.REMARK_4, line[10:-1].strip(), RECORD, line_number)
            elif RNUM == 5:
                self.REMARK_5 = parse1xN(self, self.REMARK_5, line[10:-1].strip(), RECORD, line_number)
            elif RNUM < 100:
                self.REMARK_6 = parse1xN(self, self.REMARK_6, line[10:-1].strip(), RECORD, line_number)
            elif RNUM == 100:
                self.REMARK_100 = parse1xN(self, self.REMARK_100, line[10:-1].strip(), RECORD, line_number)
            elif RNUM == 200:
                self.REMARK_200 = parse1xN(self, self.REMARK_200, line[10:-1].strip(), RECORD, line_number)
            elif RNUM == 205:
                self.REMARK_205 = parse1xN(self, self.REMARK_205, line[10:-1].strip(), RECORD, line_number)
            elif RNUM < 218:
                self.REMARK_210 = parse1xN(self, self.REMARK_210, line[10:-1].strip(), RECORD, line_number)
            elif RNUM == 230:
                self.REMARK_230 = parse1xN(self, self.REMARK_230, line[10:-1].strip(), RECORD, line_number)
            elif RNUM == 240:
                self.REMARK_240 = parse1xN(self, self.REMARK_240, line[10:-1].strip(), RECORD, line_number)
            elif RNUM == 245:
                self.REMARK_245 = parse1xN(self, self.REMARK_245, line[10:-1].strip(), RECORD, line_number)
            elif RNUM == 247:
                self.REMARK_247 = parse1xN(self, self.REMARK_247, line[10:-1].strip(), RECORD, line_number)
            elif RNUM == 250:
                self.REMARK_250 = parse1xN(self, self.REMARK_250, line[10:-1].strip(), RECORD, line_number)
            elif RNUM == 265:
                self.REMARK_265 = parse1xN(self, self.REMARK_265, line[10:-1].strip(), RECORD, line_number)
            elif RNUM == 280:
                self.REMARK_280 = parse1xN(self, self.REMARK_280, line[10:-1].strip(), RECORD, line_number)
            elif RNUM == 285:
                self.REMARK_285 = parse1xN(self, self.REMARK_285, line[10:-1].strip(), RECORD, line_number)
            elif RNUM == 290:
                self.REMARK_290 = parse1xN(self, self.REMARK_290, line[10:-1].strip(), RECORD, line_number)
            elif RNUM == 300:
                self.REMARK_300 = parse1xN(self, self.REMARK_300, line[10:-1].strip(), RECORD, line_number)
            elif RNUM == 350:
                self.REMARK_350 = parse1xN(self, self.REMARK_350, line[10:-1].strip(), RECORD, line_number)
            elif RNUM == 375:
                self.REMARK_375 = parse1xN(self, self.REMARK_375, line[10:-1].strip(), RECORD, line_number)
            elif RNUM == 400:
                self.REMARK_400 = parse1xN(self, self.REMARK_400, line[10:-1].strip(), RECORD, line_number)
            elif RNUM == 450:
                self.REMARK_450 = parse1xN(self, self.REMARK_450, line[10:-1].strip(), RECORD, line_number)
            elif RNUM == 465:
                self.REMARK_465 = parse1xN(self, self.REMARK_465, line[10:-1].strip(), RECORD, line_number)
            elif RNUM == 470:
                self.REMARK_470 = parse1xN(self, self.REMARK_470, line[10:-1].strip(), RECORD, line_number)
            elif RNUM == 475:
                self.REMARK_475 = parse1xN(self, self.REMARK_475, line[10:-1].strip(), RECORD, line_number)
            elif RNUM == 480:
                self.REMARK_480 = parse1xN(self, self.REMARK_480, line[10:-1].strip(), RECORD, line_number)
            elif RNUM == 500:
                self.REMARK_500 = parse1xN(self, self.REMARK_500, line[10:-1].strip(), RECORD, line_number)
            elif RNUM == 525:
                self.REMARK_525 = parse1xN(self, self.REMARK_525, line[10:-1].strip(), RECORD, line_number)
            elif RNUM == 600:
                self.REMARK_600 = parse1xN(self, self.REMARK_600, line[10:-1].strip(), RECORD, line_number)
            elif RNUM == 610:
                self.REMARK_610 = parse1xN(self, self.REMARK_610, line[10:-1].strip(), RECORD, line_number)
            elif RNUM == 615:
                self.REMARK_615 = parse1xN(self, self.REMARK_615, line[10:-1].strip(), RECORD, line_number)
            elif RNUM == 620:
                self.REMARK_620 = parse1xN(self, self.REMARK_620, line[10:-1].strip(), RECORD, line_number)
            elif RNUM == 630:
                self.REMARK_630 = parse1xN(self, self.REMARK_630, line[10:-1].strip(), RECORD, line_number)
            elif RNUM == 650:
                self.REMARK_650 = parse1xN(self, self.REMARK_650, line[10:-1].strip(), RECORD, line_number)
            elif RNUM == 700:
                self.REMARK_700 = parse1xN(self, self.REMARK_700, line[10:-1].strip(), RECORD, line_number)
            elif RNUM == 800:
                self.REMARK_800 = parse1xN(self, self.REMARK_800, line[10:-1].strip(), RECORD, line_number)
            elif RNUM == 999:
                self.REMARK_999 = parse1xN(self, self.REMARK_999, line[10:-1].strip(), RECORD, line_number)
                        
        elif RECORD == 'DBREF':#Mx1
            self.DBREF  = parseMx1(self.DBREF,  line)
        elif RECORD == 'DBREF1':#Mx1
            self.DBREF1 = parseMx1(self.DBREF1, line)
        elif RECORD == 'DBREF2':#Mx1
            self.DBREF2 = parseMx1(self.DBREF2, line)
        elif RECORD == 'SEQADV':#Mx1
            self.SEQADV = parseMx1(self.SEQADV, line[7:-1].strip())

        elif RECORD == 'SEQRES':#MxN
            self.SEQRES = parseMxN(self, self.SEQRES, line, RECORD)

        elif RECORD == 'MODRES':#Mx1
            self.MODRES = parseMx1(self.MODRES, line[7:-1].strip())
        elif RECORD == 'HET':#Mx1
            self.HET    = parseMx1(self.HET,    line[7:-1].strip())

        elif RECORD == 'HETNAM':#MxN
            self.HETNAM = parseMxN(self, self.HETNAM, line, RECORD)
        elif RECORD == 'HETSYN':#MxN
            self.HETSYN = parseMxN(self, self.HETSYN, line, RECORD)
        elif RECORD == 'FORMUL':#MxN
            self.FORMUL = parseMxN(self, self.FORMUL, line, RECORD)
        
        elif RECORD == 'HELIX':#Mx1
            self.HELIX  = parseMx1(self.HELIX,  line)
        elif RECORD == 'SHEET':#Mx1
            self.SHEET  = parseMx1(self.SHEET,  line)
        elif RECORD == 'SSBOND':#Mx1
            self.SSBOND = parseMx1(self.SSBOND, line)
        elif RECORD == 'LINK':#Mx1
            self.LINK   = parseMx1(self.LINK,   line[7:-1].strip())
        elif RECORD == 'CISPEP':#Mx1
            self.CISPEP = parseMx1(self.CISPEP, line)

        elif RECORD == 'SITE':#MxN
            self.SITE   = parseMxN(self, self.SITE, line, RECORD)
        
        elif RECORD == 'CRYST1':#1x1
            self.CRYST1 = parse1x1(self.CRYST1, line[7:-1].strip(), RECORD, line_number)
        elif RECORD == 'ORIGX1':#1x1
            self.ORIGX1 = parse1x1(self.ORIGX1, line[7:-1].strip(), RECORD, line_number)
        elif RECORD == 'ORIGX2':#1x1
            self.ORIGX2 = parse1x1(self.ORIGX2, line[7:-1].strip(), RECORD, line_number)
        elif RECORD == 'ORIGX3':#1x1
            self.ORIGX3 = parse1x1(self.ORIGX3, line[7:-1].strip(), RECORD, line_number)
        elif RECORD == 'SCALE1':#1x1
            self.SCALE1 = parse1x1(self.SCALE1, line[7:-1].strip(), RECORD, line_number)
        elif RECORD == 'SCALE2':#1x1
            self.SCALE2 = parse1x1(self.SCALE2, line[7:-1].strip(), RECORD, line_number)
        elif RECORD == 'SCALE3':#1x1
            self.SCALE3 = parse1x1(self.SCALE3, line[7:-1].strip(), RECORD, line_number)

        elif RECORD == 'MTRIX1':#Mx1
            self.MTRIX1 = parseMx1(self.MTRIX1, line[7:-1].strip())
        elif RECORD == 'MTRIX2':#Mx1
            self.MTRIX2 = parseMx1(self.MTRIX2, line[7:-1].strip())
        elif RECORD == 'MTRIX3':#Mx1
            self.MTRIX3 = parseMx1(self.MTRIX3, line[7:-1].strip())
        elif RECORD == 'MODEL':#MxN
            self.MODEL  = parseMxN(self, self.MODEL,  line[7:-1].strip(), RECORD)
        
        elif RECORD == 'ATOM':#Mx1
            self.ATOM   = parseMx1(self.ATOM,   line)
        elif RECORD == 'ANISOU':#Mx1
            self.ANISOU = parseMx1(self.ANISOU, line[7:-1].strip())
        elif RECORD == 'TER':#Mx1
            self.TER    = parseMx1(self.TER,    line)
        elif RECORD == 'HETATM':#Mx1
            self.HETATM = parseMx1(self.HETATM, line)

        elif RECORD == 'ENDMDL':#MxN
            self.ENDMDL = parseMxN(self, self.ENDMDL, line[7:-1].strip(), RECORD)
        
        elif RECORD == 'CONECT':#Mx1
            self.CONECT = parseMx1(self.CONECT, line[7:-1].strip())
        
        elif RECORD == 'MASTER':#1x1
            self.MASTER = parse1x1(self.MASTER, line[7:-1].strip(), RECORD, line_number)

        elif RECORD == 'END':#1x1
            if self.END == None:
                self.END = 'Done'
            else:
                warnings.warn('A second {} section found at line {}!!!'.format(RECORD, line_number))
                sys.exit(0)
        
        else:
            TEXT = line[:-1]
            if self.JUNK == None:
                self.JUNK  = [[],[]]
                self.JUNK[0].append(TEXT)
                self.JUNK[1].append(line_number)
                warnings.warn('Unrecognized section found at line {}.'.format(line_number))
            else:
                self.JUNK[0].append(TEXT)
                self.JUNK[1].append(line_number)
                warnings.warn('Unrecognized section found at line {}.'.format(line_number))
        
        line_number += 1