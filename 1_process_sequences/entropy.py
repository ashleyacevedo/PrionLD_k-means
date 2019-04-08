"""
This module contains the function shannon_entropy which calculates the shannon
entropy for sequences.

Author: Ashley Acevedo
"""

import numpy as np

def shannon_entropy(p):

    """
    The shannon_entropy function takes as input a list of probabilities
    (fraction) of each type (i.e. nucleic acid base or amino acid) and returns
    the shannon entropy.
	
    Shannon entropy = - sum(p * log(p)).
    """

    # Remove zeros entries from the list of probabilities
    p = [i for i in p if i != 0]
	
    p = np.asarray(p)
    logp = np.log(p)

    return -1 * np.sum(p * logp)
