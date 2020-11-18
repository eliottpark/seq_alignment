""" 
    RNA Alignment Assignment
    
    Implement each of the functions below using the algorithms covered in class.
    You can construct additional functions and data structures but you should not
    change the functions' APIs.

    You will be graded on the helper function implementations as well as the RNA alignment, although
    you do not have to use your helper function.
    
    *** Make sure to comment out any print statement so as not to interfere with the grading script
"""

import sys # DO NOT EDIT THIS
import numpy as np
import time

from pprint import pprint
from shared import *

ALPHABET = [TERMINATOR] + BASES

def get_suffix_array(s):
    """
    Naive implementation of suffix array generation (0-indexed). You do not have to implement the
    KS Algorithm. Make this code fast enough so you have enough time in Aligner.__init__ (see bottom).

    Input:
        s: a string of the alphabet ['A', 'C', 'G', 'T'] already terminated by a unique delimiter '$'
    
    Output: list of indices representing the suffix array

    >>> get_suffix_array('GATAGACA$')
    [8, 7, 5, 3, 1, 6, 4, 0, 2]
    """
    
    # suffixes = {}
    # for i in range(len(s)):
    #     suffixes[s[i:]] = i

    # sorted_indices = [suffixes[x] for x in sorted(suffixes)]
    
    # K = number of groups to split into
    # k = 10
    # chunk_size = len(s)/k
    # groups = []
    # # split, sort, and merge
    # for j in range(k):
    #     groups.append[np.sort(suffixes[j*chunk_size:(j + 1)*chunk_size])]

    # sorted_indices = []

    # # Try radix sort

    # while groups:
    #     c = np.concatenate(groups.pop(), groups.pop())
    #     c = np.sort(c, kind="mergesort")

    # Naive approach ##########################################
    # suffixes = []
    # for i in range(len(s)):
    #     suffixes = np.append(suffixes, s[i:])
    
    # sorted_indices = np.argsort(suffixes, kind='mergesort') # supposedly uses radix sort under the hood

    # Radix sort with buckets of len = k #######################
    k = 49
    r = 20
    suffixes = {}
    str_len = len(s)
    # Construct dictionary that represents buckets defined by kth prefix of each suffix
    print("str_len: ",str_len)
    for i in range(str_len):
        if s[i:i + k] in suffixes.keys():
            suffixes[s[i:i + k]] += [(i, s[i + k:])]
        else:
            suffixes[s[i:i + k]] = [(i, s[i + k:])]
    
    print("suffix len: ", len(suffixes))
    sorted_suffixes = k_radix_sort(suffixes, k, r)

    # pprint(sorted_suffixes)
    sout = [0 for _ in range(str_len)]
    print("str_len: ", len(sout))
    for i in range(str_len):
        sout[i] = sorted_suffixes[i] 

    print("sa len: ", len(sout))
    return sout

def k_radix_sort(suffixes, k, r):
    """
    Radix sort according to the kth prefix of each suffix.

    suffixes: a dictionary of suffixes {[:k] prefix: (index in original string, [k:] suffix)}
    k: size of prefix
    r: recursion depth

    Returns a sorted dictionary of suffixes {relative order of suffix: index in original string}
    """
    # Sort dictionary by keys and replace by order 
    # sorted_keys = sorted(suffixes.keys())
    for i, key in enumerate(sorted(suffixes.keys())):
        suffixes[i] = suffixes.pop(key)


    outer_suffixes = {}
    total_suffixes = 0
    # Recurse if there are more than one strings in a bucket
    for i in range(len(suffixes)):
        # if r == 20 and i < 20:
            # print("outer suffixes len  at i = {}: ".format(i), len(outer_suffixes))
        # print("\ni: ", i) 
        if len(suffixes[i]) > 1 and r > 0:
            inner_count = 0
            recurse_suffix = {}
            # Create suffix buckets from inner pairs
            for pair in suffixes[i]:
                total_suffixes += 1
                # print("\npair: ", pair)
                inner_count += 1
                s = pair[1]
                if s[i:i + k] in recurse_suffix.keys():
                    recurse_suffix[s[i:i + k]] += [(pair[0], s[i + k:])]
                else:
                    recurse_suffix[s[i:i + k]] = [(pair[0], s[i + k:])]
            # if r == 20 and i < 20:
                # print("    inner count: ", inner_count)
                # print("    recurse suf len: ", len(recurse_suffix))
            # Return sorted suffix dictionary
            sorted_suffixes = k_radix_sort(recurse_suffix, k, r-1)
            # if r == 20 and i < 20:
                # print("    sorted_suf len: ", len(sorted_suffixes))
            # Correct dictionary keys after index i to account for order of new suffixes
            delta = len(sorted_suffixes)
            # if r == 20 and i < 20:
            #     pprint(outer_suffixes)
            # print("sorted suffixes:", suffixes)
            outer_suffixes = {(key + delta - 1 if key > i else key): v for key, v in outer_suffixes.items()}
            # Adjust key pairs in newly sorted dictionary
            sorted_suffixes = { key + i: v for key, v in sorted_suffixes.items()}
            outer_suffixes.update(sorted_suffixes)
            # if r == 20 and i < 20:
            #     pprint(sorted_suffixes)
            #     pprint(outer_suffixes)
            
        else:
            outer_suffixes[i] = suffixes[i][0][0]
            total_suffixes += 1
    
    # if r == 20:
        # print("total suffixes: ", total_suffixes)
    return outer_suffixes
    
    # # Recurse if there are more than one strings in a bucket
    # for i, key in enumerate(sorted(suffixes.keys())):
    #     if len(suffixes[key]) > 1:
    #         recurse_suffix = {}
    #         # Create suffix buckets from inner pairs
    #         for pair in suffixes[key]:
    #             s = pair[1]
    #             if s[i:i + k] in recurse_suffix.keys():
    #                 recurse_suffix[s[i:i + k]] += [(i, s[i + k:])]
    #             else:
    #                 recurse_suffix[s[i:i + k]] = [(i, s[i + k:])]
    #         # Return sorted suffix dictionary
    #         sorted_suffixes = k_radix_sort(recurse_suffix, k)
    #         # Correct dictionary keys after index i to account for order of new suffixes
    #         delta = len(sorted_suffixes)
    #         print("\n\ni: ", i)
    #         print(suffixes)
    #         suffixes = { (j + delta if j > i else j): v for j, v in suffixes.items()}
    #         # Adjust key pairs in newly sorted dictionary
    #         sorted_suffixes = { j + i: v for j, v in sorted_suffixes.items()}
    #         suffixes.update(sorted_suffixes)
    #     else:
    #         suffixes[key] = suffixes[key][0]
    
    # # Sort dictionary by keys and replace by order 
    # for i, key in enumerate(sorted(suffixes.keys())):
    #     suffixes[i] = suffixes.pop(key)

def get_bwt(s, sa):
    """
    Input:
        s: a string terminated by a unique delimiter '$'
        sa: the suffix array of s

    Output:
        L: BWT of s as a string

    >>> s = 'GATAGACA$'
    >>> sa = get_suffix_array(s)
    >>> get_bwt(s, sa)
    'ACGTGAA$A'
    """
    bwt = ""
    for i in sa:
        bwt += s[i-1]
    
    return bwt
    

def get_F(L):
    """
    Input: L = get_bwt(s)

    Output: F, first column in Pi_sorted

    >>> s = 'BARBARA$'
    >>> sa = get_suffix_array(s)
    >>> get_F(get_bwt(s, sa))
    '$AAABBRR'
    """
    F = ''
    return F.join(sorted(list(L)))


def get_M(F):
    """
    Returns the helper data structure M (using the notation from class). M is a dictionary that maps character
    strings to start indices. i.e. M[c] is the first occurrence of "c" in F.

    If a character "c" does not exist in F, you may set M[c] = -1
    
    >>> s = 'BARBARA$'
    >>> sa = get_suffix_array(s)
    >>> f = get_F(get_bwt(s, sa))
    >>> get_M(f)
    {'$': 0, 'A': 1, 'B': 4, 'R': 6}
    """
    M = {}
    for i in range(len(F)):
        if F[i] not in M.keys():
            M[F[i]] = i
    
    return M

def get_occ(L):
    """
    Returns the helper data structure OCC (using the notation from class). OCC should be a dictionary that maps 
    string character to a list of integers. If c is a string character and i is an integer, then OCC[c][i] gives
    the number of occurrences of character "c" in the bwt string up to and including index i

    >>> s = 'BARBARA$'
    >>> sa = get_suffix_array(s)
    >>> l = get_bwt(s, sa)
    >>> get_occ(l)
    {'$': [0, 0, 0, 0, 0, 1, 1, 1], 'A': [1, 1, 1, 1, 1, 1, 2, 3], 'B': [0, 0, 1, 2, 2, 2, 2, 2], 'R': [0, 1, 1, 1, 2, 2, 2, 2]}
    """
    occ = {}

    for c in sorted(set(L)):
        occ[c] = [0]*len(L)

    for i in range(len(L)):

        for c in occ.keys():
            if i == 0:
                if L[i] == c:
                    occ[c][i] = 1
            else:
                if L[i] == c:
                    occ[c][i] = occ[c][i-1] + 1
                else:
                    occ[c][i] = occ[c][i-1]

    return occ
    

def exact_suffix_matches(p, M, occ):
    """
    Find the positions within the suffix array sa of the longest possible suffix of p 
    that is a substring of s (the original string).
    
    Note that such positions must be consecutive, so we want the range of positions.

    Input:
        p: the pattern string
        M, occ: buckets and repeats information used by sp, ep

    Output: a tuple (range, length)
        range: a tuple (start inclusive, end exclusive) of the indices in sa that contains
            the longest suffix of p as a prefix. range=None if no indices matches any suffix of p
        length: length of the longest suffix of p found in s. length=0 if no indices matches any suffix of p

        An example return value would be ((2, 5), 7). This means that p[len(p) - 7 : len(p)] is
        found in s and matches positions 2, 3, and 4 in the suffix array.

    >>> s = 'ACGT' * 10 + '$'
    >>> sa = get_suffix_array(s)
    >>> sa
    [40, 36, 32, 28, 24, 20, 16, 12, 8, 4, 0, 37, 33, 29, 25, 21, 17, 13, 9, 5, 1, 38, 34, 30, 26, 22, 18, 14, 10, 6, 2, 39, 35, 31, 27, 23, 19, 15, 11, 7, 3]
    >>> L = get_bwt(s, sa)
    >>> L
    'TTTTTTTTTT$AAAAAAAAAACCCCCCCCCCGGGGGGGGGG'
    >>> F = get_F(L)
    >>> F
    '$AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT'
    >>> M = get_M(F)
    >>> sorted(M.items())
    [('$', 0), ('A', 1), ('C', 11), ('G', 21), ('T', 31)]
    >>> occ = get_occ(L)
    >>> type(occ) == dict, type(occ['$']) == list, type(occ['$'][0]) == int
    (True, True, True)
    >>> occ['$']
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    >>> exact_suffix_matches('ACTGA', M, occ)
    ((1, 11), 1)
    >>> exact_suffix_matches('$', M, occ)
    ((0, 1), 1)
    >>> exact_suffix_matches('AA', M, occ)
    ((1, 11), 1)
    >>> exact_suffix_matches('AC', M, occ)
    ((1, 11), 2)
    >>> s = 'BARBARA$'
    >>> sa = get_suffix_array(s)
    >>> sa
    [7, 6, 4, 1, 3, 0, 5, 2]
    >>> l = get_bwt(s, sa)
    >>> f = get_F(l)
    >>> m = get_M(f)
    >>> occ = get_occ(l)
    >>> exact_suffix_matches('AR', m, occ)
    ((2, 4), 2)
    """
    
    iterator = iter(M.keys())
    
    while True:
        char = next(iterator)
        if char == p[-1]:
            sp = M[char]
            try:
                ep = M[next(iterator)] - 1
            except StopIteration:
                ep = len(occ[char]) - 1
            break

    if ep < sp:
        return ((None), 0)

    i = len(p) - 2
    while i >= 0:
        temp_sp = M[p[i]] + occ[p[i]][sp - 1]
        temp_ep = M[p[i]] + occ[p[i]][ep] - 1
        if temp_sp >= temp_ep:
            break
        sp, ep = temp_sp, temp_ep
        i -= 1

    return ((sp, ep + 1), len(p) - i - 1)

MIN_INTRON_SIZE = 20
MAX_INTRON_SIZE = 10000

class Aligner:
    def __init__(self, genome_sequence, known_genes):
        """
        Initializes the aligner. Do all time intensive set up here. i.e. build suffix array.

        genome_sequence: a string (NOT TERMINATED BY '$') representing the bases of the of the genome
        known_genes: a python set of Gene objects (see shared.py) that represent known genes. You can get the isoforms 
                     and exons from a Gene object

        Time limit: 500 seconds maximum on the provided data. Note that our server is probably faster than your machine, 
                    so don't stress if you are close. Server is 1.25 times faster than the i7 CPU on my computer

        """
        # Making a set of isoforms to represent transcriptome. TODO: could be a hash table
        self._isoforms = {}
        # Iterate through all genes
        print("iterating thru genes")
        for gene in known_genes:
            # Iterate through all isoforms for each gene
            print("Gene: ", gene.id)
            for isoform in gene.isoforms:
                # Build isoform from individual exons
                full_isoform = ''
                for exon in isoform.exons: # If ordered, ok. if not, need to make sure
                    full_isoform += genome_sequence[exon.start:exon.end]
                full_isoform += '$'
                print("full isoform len: ", len(full_isoform))
                print("isoform id: ", isoform.id)
                if isoform.id == 'ENST00000475864':
                    print(full_isoform)
                
                # Find SA, M, and OCC of the isoform
                print("    getting suffix array")
                start_time = time.time()
                sa = get_suffix_array(full_isoform)
                print("    --- %s seconds ---" % (time.time() - start_time))
                # print("    suffix array: ", sa)
                self._isoforms[isoform.id] = [get_bwt(full_isoform, sa), sa,
                                                get_M(full_isoform), get_occ(full_isoform), 
                                                full_isoform, isoform]
                    
        # Build SA, M, OCC for whole genome
        print("building sa, m, occ, bwt for whole genomes")
        self._sa = get_suffix_array(genome_sequence)
        self._bwt = get_bwt(genome_sequence, self._sa)
        self._m = get_M(genome_sequence)
        self._occ = get_occ(genome_sequence)
        self._genome_seq = genome_sequence
            
        # To make transcriptome: for each known_gene, find isoform; for each isoform, find exons; concatenate together.
        # Build suffix array, M array, Occ array.

    
    def greedy_inexact_alignment(self, p, M, occ, isoform_id=None):
        """
        Run a greedy inexact alignment algorithm to determine an alignment to the reference text
        with less than or equal to MAX_NUM_MISMATCHES. 

        p: Query string
        M: M array of the reference text
        occ: occ array of the reference text
        isoform_id: Isoform id to refer

        Returns an alignment to reference text if one can be found with less than MAX_NUM_MISMATCHES
        in the form ((<range of matches in suffix array>, <length of longest match>), <number of mismatches>)
        """
        num_mismatches = 0
        longest_match_length = 0
        query = p
        query_len = len(query)

        # Iterate through until there are too many mismatches or we get an alignment
        while num_mismatches < MAX_NUM_MISMATCHES:
            locations = exact_suffix_matches(query, M, occ)
            location_range, match_len = locations[0], locations[1]

            # If the match is the length of the query, return
            if match_len == query_len:
                return locations, num_mismatches
            # Else, mutate the query at the point of mismatch to match the reference text
            else:
                # Location of mismatch in the query
                i = query_len - match_len - 1

                # Align to genome or isoform TODO: Figure out if this is an accurate way to choose characters
                # ex when we have multiple alignments up to this point. Maybe we need to keep track of the longest alignment?
                # --> This way we can select the last mismatched character for the most promising/longest alignment instead of the
                # first in the range.
                if isoform_id:
                    j = self._isoforms[isoform_id][1][location_range[0]] # location of first aligned char in suffix array
                    query[i] = self._isoforms[isoform_id][4][j] # Set query mismatch to reference text value at corresponding index
                else:
                    j = self._sa[location_range[0]] # location of first aligned char in suffix array
                    query[i] = self._genome_seq[j] # Set query mismatch to reference text value at corresponding index

                num_mismatches += 1

        return None, None


    def align_to_isoforms(self, read_sequence):
        """
        Returns the best alignment of the read sequence to the isoform database. 
        Prioritize matching to known isoforms. Minimize number of mismatches

        read_sequence: input read text to be aligned

        Return value must be in the form (also see the project pdf):
        [(<read_start_1>, <reference_start_1, length_1), (<read_start_2>, <reference_start_2, length_2), ...]

        If no good matches are found: return the best match you can find or return []
        """

        # Initialize return value and matches
        isoform_matches = [] # [(isoform.id, locations, num_mismatches)...]
        exon_matches = [] # [(<read_start_1>, <reference_start_1, length_1), (<read_start_2>, <reference_start_2, length_2), ...]

        # Matches are in the form (read start index, genome start index, length)
        for isoform_id in self._isoforms.keys():
            isoform = self._isoforms[isoform_id]
            m, occ =  isoform[2], isoform[3]
            locations, num_mismatches = self.greedy_inexact_alignment(read_sequence, m, occ, isoform_id)
            if locations:
                # Want to add to list of matches and then find the best match.
                isoform_matches.append((isoform_id, locations, num_mismatches))

        
        # Check to see if there are any returned matches.
        if isoform_matches:
            # Finding the match with the least number of mismatches and converting it to the correct format.
            best_match = (0, 0, 7)
            for match in isoform_matches:
                if match[2] < best_match[2]:
                    best_match = match
            # Now need to convert into correct format.
            # Format: [(index of read, index of genome, length of match)]
            isoform_arr = self._isoforms[best_match[0]]

            # Start and end indices for the isoform of interest
            sa = isoform_arr[1]
            ifs_index = sa[best_match[1][0][0]] # index of first match in fully concatenated isoform string
            end = sa[best_match[1][0][1]] # index of last match in fully concatenated isoform string
        
            # check to see if start = end (only a single match)
            if ifs_index != end:
                pass
                # TODO figure out what happens in this situation

            # Iterate through exons to find indices
            read_index = 0
            genome_index = 0
            read_len = len(read_sequence)
            ifs_offset = 0 # length of all exons preceding concatenated isoform string index

            # iterate through all exons until we finish the read
            for exon in isoform_arr[5].exons: # This should be ordered
                genome_index = exon.start + ifs_index - ifs_offset

                # If the genome index is in the current exon
                if genome_index < exon.end:
                    # Check to see if the read ends in this exon window
                    if exon.end > (exon.start + read_len - read_index):
                        length = read_len - read_index
                        
                        # Append to list
                        exon_matches.append((read_index, genome_index, length))

                        return exon_matches

                    # If not, extend match to end of exon window
                    else:
                        length = exon.end - genome_index

                        # Append to list
                        exon_matches.append((read_index, genome_index, length))

                    read_index += length
                    ifs_index += length
                # Else, progress to the next exon
                else:
                    length = exon.end - exon.start
                ifs_offset += length
        
        return exon_matches
                                

    def align_to_genome(self, read_sequence):
        """
        Returns the best alignment of the read sequence to the genome. This will cover reads
        that originate from new exons.

        read_sequence: input read text to be aligned

        Return value must be in the form (also see the project pdf):
        [(<read_start_1>, <reference_start_1, length_1), (<read_start_2>, <reference_start_2, length_2), ...]

        If no good matches are found: return the best match you can find or return []
        """
        # Return value
        matches = []
        # Retrieve the best inexact alignment to the genome
        # locations = (<range of matches in suffix array>, <length of longest match>)
        locations, num_mismatches = self.greedy_inexact_alignment(read_sequence, self._m, self._occ)

        if locations:
            matches.append((0, self._sa[locations[0][0]]))
        
        return matches

    def align_seeds(self, read_sequence):
        """
        Split into seeds and align.

        read_sequence: input read text to be aligned

        Return value must be in the form (also see the project pdf):
        [(<read_start_1>, <reference_start_1, length_1), (<read_start_2>, <reference_start_2, length_2), ...]

        If no good matches are found: return the best match you can find or return []
        """
        pass
    

    def align(self, read_sequence):
        """
        Returns an alignment to the genome sequence. An alignment is a list of pieces. 
        Each piece consists of a start index in the read, a start index in the genome, and a length 
        indicating how many bases are aligned in this piece. Note that mismatches are count as "aligned".

        Note that <read_start_2> >= <read_start_1> + <length_1>. If your algorithm produces an alignment that 
        violates this, we will remove pieces from your alignment arbitrarily until consecutive pieces 
        satisfy <read_start_2> >= <read_start_1> + <length_1>

        Return value must be in the form (also see the project pdf):
        [(<read_start_1>, <reference_start_1, length_1), (<read_start_2>, <reference_start_2, length_2), ...]

        If no good matches are found: return the best match you can find or return []

        Time limit: 0.5 seconds per read on average on the provided data.
        """

        reads = []

        # Try aligning to the isoform database (case 1 in tophat2)
        reads = self.align_to_isoforms(read_sequence)

        # Align to genome if no isoform matches (case 2 in tophat2)
        if not reads:
            reads = self.align_to_genome(read_sequence)

        # If no isoform matches or direct genome matches, split into seeds and
        # align seed to try to find new splice junctions (case 3 in tophat2)
        if not reads:
            reads = self.align_seeds(read_sequence)
                
        return reads


            
            
            ## OLD BOWTIE1 CODE #################################################################
                # if index_ifs >= exon.start and index_ifs < exon.end:
                #     exon_len = exon.end - start
                #     exon_matches.append((read_index, )
            # num_backtracks = 0
            # while num_backtracks < k:
            #     locations = exact_suffix_matches(query, M, occ)

            #     # If we improve our match length, then we have successfully introduced a mismatch
            #     if locations[1] > longest_match_length:
            #         longest_match_length = locations[1]
            #         best_match = locations
            #         num_mismatches += 1
                    
            #     if locations[1] == 0: # might fold into else 
            #         return [], 0
            #     elif locations[1] == query_len and num_mismatches <= MAX_NUM_MISMATCHES:
            #         return locations, num_mismatches
            #     else:
            #         i = query_len - locations[1] - 1
            #         mismatched_char = query[i]
            #         # leftmost just-visited position with minimal quality
            #         j = i + 1
            #         if j not in R.keys():
            #             c = random.choice(set(nucleotides) - set(query[j]))
            #             R[j] = set(c)
                        
            #         else:
            #             c = random.choice(set(nucleotides) - set(R[j]))
            #             R[j] += c
            #         query[j] = c # Introducing a mismatch

            #         num_backtracks += 1                   
            # return [], 0
    