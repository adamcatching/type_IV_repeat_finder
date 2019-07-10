# Benjamin Adam Catching
# 2018-02-05
# Bondy-Denomy Lab
# Type IV CRISPR Project

"""
The goal of this code is to make a simple python code that any person
with an installation of python 2.7 or 3 on their computer can run. This code
does not have any package dependencies besides the basics so shouldn't break 
that often.

This code has been tested on a Mac running OS Mojave 10.14.5 without python
packages installed
 """

def hamming(seq1, seq2):
	"""
	From two seqs of a given length, find the hamming distance.

	:param seq1: first string of letters to compare against
	:type seq1: any subscriptable type; string, list...

	:param seq2: second string of letters to compare against
	:type seq2: any subscriptable type

	return: distance between the two sequences
	rtype: integer
	"""

	# Initiate hamming distance value
	ham = 0
	# Ensure that both sequences are of equal length
	if len(seq1) != len(seq2):
		print('Sequences are not of equal length')
	else:
		# Iterate over the values in seq1
		for l1, letter1 in enumerate(seq1):
			# Get the compareable 
			letter2 = seq2[l1]
			# If the two values are not equal, add to the hamming distance
			if letter1 != letter2:
				ham += 1

		return ham

def rev_comp(seq):
	"""
	From a string of DNA letters, return the reverse complement

	param seq: a string of letters only including A, T, C, G
	type seq: string

	return: the reverse complement sequence of the input
	rtype: string
	"""

	# List to store nucleotides in
	result = []

	# Work from last letter to first letter finding complement
	for letter in seq[::-1]:
		if letter == 'C':
			result.append('G')
		elif letter == 'G':
			result.append('C')
		elif letter == 'T':
			result.append('A')
		elif letter == 'A':
			result.append('T')
	return ''.join(result)

def palindromic_finder(target_seq, 
	template_palin='CCCCGC',
	template_mismatch=2,
	palindrome_mismatch=2,
	up_down_save=10):
    """
    A target sequence thought to have palindromic repeats is searched 
    through for repets that match the template sequence and has a 
    palindrome four base-pairs downstream that matches the palindrome of the
    target sequence. All potential palindromes are returned.

    param seq:
    """
    
    # Length of the sequence that will be searched for palindromes
    target_len = len(target_seq)
    # Number of upstream and downstream basepairs to keep for analysis
    repeat_end_len = 10
    # Length of the first half of the palindrome
    template_len = len(template_palin)
    
    # Dictionary to store location and sequence of palindromic repeats  
    potential_repeats = {}
    
    # Sliding window search for k-mers matching the template with a hamming
    # distance equal to or less than template_mismatch
    for i in range(up_down_save,(target_len-23)):
    	# Get k-mer at ith position to test against
        temp_pal_1 = target_seq[i:(i+template_len)]
        if hamming(temp_pal_1, template_palin) <= template_mismatch:
            # New position to test against original k-mer
            temp_pal_2_loc = i+repeat_end_len
            temp_pal_2_loc_end = temp_pal_2_loc + template_len
            # New k-mer to test against original k-mer
            temp_pal_2 = target_seq[(temp_pal_2_loc):(temp_pal_2_loc_end)]
            # Reverse-complement of k-mer which will match the original k-mer
            # if the total sequence is a palindrome
            true_pal_2 = rev_comp(temp_pal_1)
            # Test if the hamming distance between the original k-mer and the
            # reverse-complement of the downstream k-mer are below the given
            # threshold 
            if hamming(true_pal_2, temp_pal_2) <= palindrome_mismatch:
                # Region begin and ending position
                temp_reg_begin = i-repeat_end_len
                temp_reg_end = i + 4 + 2 * template_len + repeat_end_len
                # Slice out palindromic region
                temp_repeat = target_seq[temp_reg_begin:temp_reg_end]
                # Add position and sequence to dictionary
                potential_repeats[i] = temp_repeat

    return potential_repeats

"""This is a test of palindrome recognition below"""

# Palindromic template used for paper's repeats, can be changed
template_pal_1 = 'CCCCGC'

# Test sequence containing palindromes
test_seq = 'TCTCATGCTCGATCCATGGGCCGCTCACCCCCGCACATGCGGGGAACCCATCGGCTGACTCA\
TATGACGTTCGTCAGAATCCCGCTCACCCCCGCGCACGCGGGAAATACCAACGCAGAGACGACTGGGACGACCT\
GCCAGACCGCTCACCCCCGCACACGCGGGAAATACCCTTATCCGCCAAATGCGGCCTCAGCATGATGCCGCTCA\
CCCCCGCGCACGCGGGAAATACTTCGCGCTTGCTCATTTCTCACCCCCACAGCCCCGCTCACCCCCGCGCACGC\
GGGAAATACAACCGTCAATACCCCTTCAAATCCCCAAAATACCGTTCACCCCCGCTTGCGCAGGGAATACGTTG\
CTGCGGTTTCGCAACGTGCGCTACAGGACCGCTTACCCCTACGGCGCGGGGAACACGGGTTCATCAACGGCGAT\
GACGGGCGCAGCTACCGAATACCTCCACGCAGGTGGACAACACCAGAACCCTGATTTTCTCCCTGGATATATAA\
GCCGCTTACCCCCGCGCACGCGGGAAATACCCGAAAACAGCGCCACCGCTCTCATCGTATGCCCGCTCACCCCC\
GCGCACGCGGGAAATACTGGCACCGCATGCTGTTCAGTTCGCAAACATGCCGCTCACCCCCGCGCACGCGGGAA\
ATACGGCATCGTCGCACTGCTGATCGCCTACTTCGCCCGCTCACCCCCGCGCACGCGGGAAATACAAATGATGG\
GAGCTGAGAAGATCGCCACAAACCCGCTCACCCCCGCGCACGCGGGAAATACGGAGCTGATCGAGGGCTGGTTT\
ACCCAAGAGAACCGCTCACCCCCGCGCACGCGGGAAATACGTAAGTGGTTTGCCGCCAATCTCAGGCATCTGCC\
GCTCACCCCCGCGCACGCGGGAAATACGTTGTGGTGCTTGGGTGGGACAGAGTTTGGGCCCGCTCACCCCCGCG\
CACGCGGGAAATACACATAGATATTGATCTTCCGCATAGTGCGAACACTACATCTGGTAACTTGCATGGTTAAG'

test_palindromes = palindromic_finder(test_seq, template_mismatch=2)
print(test_palindromes)
