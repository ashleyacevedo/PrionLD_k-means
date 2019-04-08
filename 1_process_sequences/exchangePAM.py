"""
This module replaces an amino acid sequence with a sequence of exchange groups.

Author: Ashley Acevedo
"""
def exchangePAM(seq):

	"""
	The exchange PAM function requires an input string composed of single
	letter amino acids and returns a sequence of U, O, J, Z, X and B signifying
	distinct exchange groups. Exchange groups are based on similarity in the
	PAM substitution matrix.
	"""
	
	seq = str(seq)

	seq = seq.replace("H", "U")
	seq = seq.replace("R", "U")
	seq = seq.replace("K", "U")
	seq = seq.replace("D", "O")
	seq = seq.replace("E", "O")
	seq = seq.replace("N", "O")
	seq = seq.replace("Q", "O")
	seq = seq.replace("C", "J")
	seq = seq.replace("S", "Z")
	seq = seq.replace("T", "Z")
	seq = seq.replace("P", "Z")
	seq = seq.replace("A", "Z")
	seq = seq.replace("G", "Z")
	seq = seq.replace("M", "X")
	seq = seq.replace("I", "X")
	seq = seq.replace("L", "X")
	seq = seq.replace("V", "X")
	seq = seq.replace("F", "B")
	seq = seq.replace("Y", "B")
	seq = seq.replace("W", "B")

	return(seq)
