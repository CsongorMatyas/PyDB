def translate3to1(LIST):
	TEXT = ''
	DICT = {'ALA' : 'A', 'CYS' : 'C', 'ASP' : 'D', 'GLU' : 'E', 'PHE' : 'F', 'GLY' : 'G',
	        'HIS' : 'H', 'ILE' : 'I', 'LYS' : 'K', 'LEU' : 'L', 'MET' : 'M', 'ASN' : 'N',
	        'PRO' : 'P', 'GLN' : 'Q', 'ARG' : 'R', 'SER' : 'S', 'THR' : 'T', 'SEC' : 'U',
	        'VAL' : 'V', 'TRP' : 'W', 'TYR' : 'Y'}
	for aa in LIST:
		if aa in DICT:
			TEXT += DICT[aa]
		else:
			warnings.warn('Unrecognized amino acid: {} was replaced with X'.format(aa))
			TEXT += 'X'
	return(TEXT)
