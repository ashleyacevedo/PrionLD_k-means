"""
generate_sequence_windows_and_features.py reads fasta formatted sequences and
breaks them into windows with a given size and step length. From these windows
the composition of ngrams of varying length are computed and written to file.
Ngrams can be composed of amino acids alone or amino acids and exchange groups.

Author: Ashley Acevedo
"""

# Import required packages
from Bio import SeqIO
from entropy import shannon_entropy
from ngram import ngram
from exchangePAM import exchangePAM
import gzip

#-------------------------------------------------------------------------------
# Inputs:
window_size      = 50
step_size        = 10
max_ngram_length = 2
exchange         = False
fasta_file       = "../input_data/FUS_family.fasta"
out_file_prefix  = "../processed_sequences/FUS_family_window_" +\
				   str(window_size) + "_step_" + str(step_size) + "_ngram_" +\
				   str(max_ngram_length) + "_"
#-------------------------------------------------------------------------------

amino_acids = ["A", "C", "D", "E", "F",
               "G", "H", "I", "K", "L",
               "M", "N", "P", "Q", "R",
               "S", "T", "V", "W", "Y"]

# Generate lists of all possible amino acid combinations of length n
ngram_lists = []
for n in range(1, max_ngram_length + 1):
	ngram_lists.append(ngram(amino_acids, amino_acids, n))

# Generate lists of all possible exchange group combinations of length n
if exchange:
	exchange_groups = ["U", "O", "J", "Z", "X", "B"]

	ex_ngram_lists = []
	for n in range(1, max_ngram_length + 1):
		ex_ngram_lists.append(ngram(exchange_groups, exchange_groups, n))

out_num = 1
iter = 1
for record in SeqIO.parse(fasta_file, "fasta"):
	
	# Open an indexed outfile for 200 genes to be written and add a header
	if iter == 1:
		# Open file and add header to outfile
		out = gzip.open(out_file_prefix + str(out_num).zfill(3) + ".txt.gz", "wb")
		out.write("ProteinName\tStartPosition\tSequence\t")
		for n in range(1, max_ngram_length + 1):
			out.write("\t".join(ngram_lists[n - 1]) + "\tEntropy_" + str(n) + "\t")
			if exchange:
				out.write("\t".join(ex_ngram_lists[n - 1]) + "\texEntropy_" + str(n) + "\t")
		out.write("\n")
	
	seq_len = len(record.seq)
	i = 0
	while i < seq_len - window_size + 1:
		window = record.seq[i:(i + window_size)]
		out.write(record.id + "\t" + str(i + 1) + "\t" + str(window))
		for n in range(1, max_ngram_length + 1):
			
			# Compute ngram probability and entropy for amino acids
			probs = []
			for ngram in ngram_lists[n - 1]:
				p = window.count(ngram)/float(window_size - (n - 1))
				probs.append(p)
				out.write("\t" + str(p))
			out.write("\t" + str(shannon_entropy(probs)))
			
			if exchange:
				# Compute ngram probability and entropy for exchange groups
				ex_window = exchangePAM(window)
				probs = []
				for ngram in ex_ngram_lists[n - 1]:
					p = ex_window.count(ngram)/float(window_size - (n - 1))
					probs.append(p)
					out.write("\t" + str(p))
				out.write("\t" + str(shannon_entropy(probs)))
		
		out.write("\n")
		i += step_size

	# Reset iteration count, increase file index and close outfile so the next
	# file can be opened and written to
	if iter == 200:
		out_num += 1
		iter = 0
		out.close()

	iter += 1
