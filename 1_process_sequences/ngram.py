"""
This module contains the function ngram which recursively generates a list of
all possible ngrams.

Author: Ashley Acevedo
"""

def ngram(current_list, types, n):

	"""
	The ngram function recursively generates a list of all possible ngrams of
	length n. The ngrams are composed of components in the list, types. In each
	round of ngram, each 'type' is added to each component in current_list.
	
	The ngram function should be called as: ngram(types, types, n)
	"""
	
	if n == 1:
		return current_list
	else:
		new_list = []
		for c in current_list:
			for t in types:
				new_list.append(c + t)
		return(ngram(new_list, types, n - 1))
