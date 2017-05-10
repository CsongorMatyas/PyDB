def translate1to3(TEXT):
	LIST = []
	DICT = {'A' : 'ALA', 'C' : 'CYS', 'D' : 'ASP', 'E' : 'GLU', 'F' : 'PHE', 'G' : 'GLY',
	        'H' : 'HIS', 'I' : 'ILE', 'K' : 'LYS', 'L' : 'LEU', 'M' : 'MET', 'N' : 'ASN',
	        'P' : 'PRO', 'Q' : 'GLN', 'R' : 'ARG', 'S' : 'SER', 'T' : 'THR', 'U' : 'SEC',
	        'V' : 'VAL', 'W' : 'TRP', 'Y' : 'TYR'}
	for aa in TEXT:
		if aa in DICT:
			LIST.append(DICT[aa])
		else:
			warnings.warn('Unrecognized amino acid: "{}" was replaced with "UNK"'.format(aa))
			LIST.append('UNK')
	return(LIST)
