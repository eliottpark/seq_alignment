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
import gc
import zlib
import numpy as np
import heapq
import time
import random
# import sufarray

# from pysuffixarray.core import SuffixArray
# from pprint import pprint
# from memory_profiler import profile

from shared import *

ALPHABET = [TERMINATOR] + BASES

# @profile
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
    k = 280
    r = 2
    suffixes = {}
    str_len = len(s)
    # Construct dictionary that represents buckets defined by kth prefix of each suffix
    for i in range(str_len):
        # # Only split if the bucket will store data
        prefix = s[i:i + k]
        if prefix in suffixes.keys():
            # suffixes[prefix] += [(i, zlib.compress(bytes(suffix, encoding='utf-8')))]
            suffixes[prefix] += [i]
        else:
            # suffixes[prefix] = [(i, zlib.compress(bytes(suffix, encoding='utf-8')))]
             suffixes[prefix] = [i]
        if i % 500000 == 0:
            print("slice ", i)

    # Sort long suffixes
    sorted_suffixes = k_radix_sort(s, suffixes, k, r)

    suffixes.clear()
    # gc.collect()

    return sorted_suffixes

# @profile
def k_radix_sort(s, suffixes, k, r):
    """
    Radix sort according to the kth prefix of each suffix.

    suffixes: a dictionary of suffixes {kth prefix: [(index in original string, [k:] suffix), ...]}
    k: size of prefix
    r: recursion depth

    Returns a sorted dictionary of suffixes {relative order of suffix: index in original string}
    """
    # print("r = ", r)
    # Sort dictionary by keys and replace by order 
    s_len = len(suffixes)
    for i, key in enumerate(sorted(suffixes.keys())):
        suffixes[i] = suffixes.pop(key) # switch keys from kth prefix to index/order
    # print(suffixes)

    outer_suffixes = []
    # offset = 0
    # Recurse if there are more than one strings in a bucket
    for i in range(len(suffixes)):
        if i % 10000 == 0 and i != 0:
            print("bucket ", i)
        if len(suffixes[i]) > 1 and r > 0:
            inner_count = 0
            recurse_suffix = {}
            # short_suffix = {}
            j = 0
            # Create suffix buckets from inner pairs
            for index in suffixes.pop(i):
                inner_count += 1
                # s = str(zlib.decompress(pair[1])).lstrip('b\'').rstrip('\'') # Rest of the string available.
                suf = s[index + k*r:]
                # if pair[0] in missing:
                #     print("missing")
                # If length is greater than k, we want to recurse 
                # if len(suf) > k:
                pre = s[:k]
                # if pair[0] in missing:
                #     print("missing")
                if pre in recurse_suffix.keys():
                    # recurse_suffix[prefix] += [(pair[0], zlib.compress(bytes(suffix, encoding='utf-8')))]
                    recurse_suffix[pre] += [index]
                else:
                    # recurse_suffix[prefix] = [(pair[0], zlib.compress(bytes(suffix, encoding='utf-8')))]
                    recurse_suffix[pre] = [index]
                # Otherwise, just sort directly on what we have
                # else:
                #     short_suffix[pair[1]] = pair[0]
                    # if pair[0] in missing:
                    #     print("missing")
            
            # if short_suffix:
            #     ##### Sort and add in short suffixes #####
            #     for j, key in enumerate(sorted(short_suffix.keys())):
            #         short_suffix[i + j] = short_suffix.pop(key)
            #     # Add in short suffixes to outer suffixes
            #     # for key in short_suffix.keys():
            #         # if key in outer_suffixes.keys():
            #         #     raise KeyError("must add unique dict elements")
            #     outer_suffixes += list(short_suffix.values())
            
            ##### Sort and add in recursive suffixes #####
            sorted_suffixes = k_radix_sort(s, recurse_suffix, k, r-1)
            outer_suffixes += sorted_suffixes
            recurse_suffix.clear()
            # short_suffix.clear()

        elif len(suffixes[i]) > 1:
            # Sort the bucket directly at maximum recursion depth
            max_r_suffixes = {}
            # Turn list of tuples (index in original string, remainder of string) into sorted dict
            for index in suffixes.pop(i):
                # s = str(zlib.decompress(pair[1])).lstrip('b\'').rstrip('\'') 
                suf = s[index + k*r:]
                max_r_suffixes[suf] = index
                # if pair[0] in missing:
                #     print("missing")
                # offset += 1
            # Sort by keys, and reassign keys to be overall index
            for j, key in enumerate(sorted(max_r_suffixes.keys())):
                max_r_suffixes[i + j] = max_r_suffixes.pop(key)
            # Combine with outer dictionary
            outer_suffixes += list(max_r_suffixes.values())
            max_r_suffixes.clear()

        else:
            # Append suffix index to output array
            outer_suffixes.append(suffixes.pop(i)[0])
    
    suffixes.clear()

    return outer_suffixes
    

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
        if temp_sp > temp_ep:
            break
        sp, ep = temp_sp, temp_ep
        i -= 1

    return ((sp, ep + 1), len(p) - i - 1)

MIN_INTRON_SIZE = 20
MAX_INTRON_SIZE = 10000

class Aligner:
    # @profile
    def __init__(self, genome_sequence, known_genes):
        """
        Initializes the aligner. Do all time intensive set up here. i.e. build suffix array.

        genome_sequence: a string (NOT TERMINATED BY '$') representing the bases of the of the genome
        known_genes: a python set of Gene objects (see shared.py) that represent known genes. You can get the isoforms 
                     and exons from a Gene object

        Time limit: 500 seconds maximum on the provided data. Note that our server is probably faster than your machine, 
                    so don't stress if you are close. Server is 1.25 times faster than the i7 CPU on my computer

        """
        overall_start_time = time.time()
        # Making a set of isoforms to represent transcriptome. TODO: could be a hash table
        self._isoforms = {}
        self._transcriptome = ''
        # Iterate through all genes
        print("iterating thru genes")
        transcriptome_index = 0
        for gene in known_genes:
            # Iterate through all isoforms for each gene
            # print("Gene: ", gene.id)
            for isoform in gene.isoforms:
                start = transcriptome_index
                # # Build isoform from individual exons
                # full_isoform = ''
                # for exon in isoform.exons: # If ordered, ok. if not, need to make sure
                #     full_isoform += genome_sequence[exon.start:exon.end]
                # full_isoform += '$'
                # Add to full transcriptome with delimiter
                for exon in isoform.exons: # If ordered, ok. if not, need to make sure
                    self._transcriptome += genome_sequence[exon.start:exon.end]
                self._transcriptome += '$'
                transcriptome_index = len(self._transcriptome) - transcriptome_index
                # print("isoform id: ", isoform.id)
                # print("    full isoform len: ", len(full_isoform))
                # if isoform.id == "ENST00000475864":
                #     print(full_isoform)
                
                # Find SA, M, and OCC of the isoform
                # print("    getting suffix array")
                # start_time = time.time()
                # sa = get_suffix_array(full_isoform)
                # print("ours:    --- %s seconds ---" % (time.time() - start_time))
                # start_time = time.time()
                # sa_theirs = sufarray.SufArray(full_isoform)
                # print("thiers:  --- %s seconds ---" % (time.time() - start_time))
                # if (sa != sa_theirs.get_array()):
                #     raise KeyError("Wrong suffix array")
                # print("Match? ", sa == sa_theirs.get_array())
                # self._isoforms[isoform.id] = [get_bwt(full_isoform, sa), # 0
                #                                 sa,                      # 1
                #                                 get_M(full_isoform),     # 2
                #                                 get_occ(full_isoform),   # 3
                #                                 full_isoform,            # 4
                #                                 isoform]                 # 5

                self._isoforms[isoform.id] = [isoform,                     # 0 - full isoform object
                                              start,                       # 1 - start index
                                              transcriptome_index - 1      # 2 - end index
                                             ]                 

        # FM Index for Transcriptome
        print("getting sa for full transcriptome")
        start_time = time.time()
        self._t_sa = get_suffix_array(self._transcriptome)
        print("    --- %s seconds ---" % (time.time() - start_time))
        print("getting bwt for full transcriptome")
        self._t_bwt = get_bwt(self._transcriptome, self._t_sa)
        print("    --- %s seconds ---" % (time.time() - start_time))
        print("getting m array for full transcriptome")
        self._t_m = get_M(get_F(self._t_bwt))
        print("    --- %s seconds ---" % (time.time() - start_time))
        print("getting occ array for full transcriptome")
        self._t_occ = get_occ(self._t_bwt)
        print("    --- %s seconds ---" % (time.time() - start_time))
                    
        print("isoforms:    --- %s seconds ---" % (time.time() - overall_start_time))
        # Build SA, M, OCC for whole genome
        print("building sa, m, occ, bwt for whole genomes")
        start_time = time.time()
        self._sa = get_suffix_array(genome_sequence)
        print("sa:    --- %s seconds ---" % (time.time() - start_time))
        start_time = time.time()
        print("bwt:    --- %s seconds ---" % (time.time() - start_time))
        self._bwt = get_bwt(genome_sequence, self._sa)
        self._m = get_M(get_F(self._bwt))
        self._occ = get_occ(genome_sequence)
        print("bwt, m, occ:    --- %s seconds ---" % (time.time() - start_time))
        self._genome_seq = genome_sequence
        print("full init:    --- %s seconds ---" % (time.time() - overall_start_time))
            
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
        in the form (<index of match in reference text>, <number of mismatches>)

        Make entire transcriptome --> concat w delimiters --> keep track of limits within the isoform 
        array w/isoform id of every index in transcriptome
        """
        num_mismatches = 0
        longest_match_length = 0
        query = p
        query_len = len(query)
        choices = [] # Substitution priority queue <index pattern: substitutions we know of but haven't tried> 
                     # based on the length of the match up to that point
        subbed = {} # Already substituted values 

        # Need to check both the forward and the reverse of the isoforms.
        # rev = self._isoforms[isoform_id][4][::-1] # sa
        # rev_m = get_M(rev)
        # rev_occ = get_occ(rev)
        # need_to_check = [[M, occ, 0], [rev_m, rev_occ, 0]] # [M array, occ array, number of mismatches]

        return_val = (0, 7) # Format: (index of match, number of mismatches)

        # Need to compare the forward and reverse matches. --> change the 
        fewest_mismatches = 0
        # for isoform in need_to_check:
        #     M = isoform[0]
        #     occ = isoform[1]
        #     num_mismatches = 0

        # Iterate through until there are too many mismatches or we get an alignment
        while num_mismatches < MAX_NUM_MISMATCHES:
            locations = exact_suffix_matches(query, M, occ)
            location_range, longest_match_len = locations[0], locations[1]
            print("longest match len: ", locations)

            # If the match is the length of the query, return
            if longest_match_len == query_len:
                if num_mismatches < return_val[1]:
                    return (locations, num_mismatches) # Still need to convert this to index in reference text.
                #return locations, num_mismatches # incorrect --> going to return a range of locations
            # Else, mutate the query, at the point of mismatch with longest match length, to match the reference text
            else:
                # Location of longest mismatch in the query. This is where we will start from to check for len
                q_match_index = query_len - longest_match_len - 1
                # Iterate through all matches and add to a priority queue based off of length of match
                for k in location_range: # location of suffix that matches
                    
                    # location and string of match in reference text
                    if isoform_id:
                        j = self._isoforms[isoform_id][1][k] # location of aligned suffix in reference text
                        ref_suffix = self._isoforms[isoform_id][4][j:] # suffix at j in reference text

                    else:
                        j = self._sa[k] # location of aligned suffix in reference text
                        ref_suffix = self._genome_seq[j:] # suffix at j in reference text

                    curr_match_len = 0 
                    curr_query_index = 0
                    n, m = q_match_index, j  # location in query, location in reference 
                    # Start at latest possible location in query, then do direct char comparison 
                    while (n < query_len):
                        if isoform_id:
                            if (self._isoforms[isoform_id][4][m] == query[n]):
                                if curr_match_len == 0:
                                    curr_query_index = n
                                curr_match_len += 1
                                n += 1
                                m += 1
                            else:
                                n += 1
                        else:
                            if (self._genome_seq[m] == query[n]):
                                if curr_match_len == 0:
                                    curr_query_index = n
                                curr_match_len += 1
                                n += 1
                                m += 1
                            else:
                                n += 1

                    # Prioritize by the inverse of the length. Store (<location in query>, <location in reference>)
                    heapq.heappush(choices, (-1*curr_match_len, (curr_query_index, j)))

                # Pop best option off stack and edit query
                curr_best = heapq.heappop(choices)
                print(curr_best)
                query_index, ref_index = curr_best[1][0], curr_best[1][1]

                # Set query mismatch for longest match in reference text
                if isoform_id:
                    sub = self._isoforms[isoform_id][4][ref_index - 1] # Set query mismatch to reference text value at corresponding index
                    query = query[:query_index - 1] + sub + query[query_index:]
                else:
                    sub = self._genome_seq[ref_index - 1] # Set query mismatch to reference text value at corresponding index
                    query = query[:query_index - 1] + sub + query[query_index:]
                
                print("\nnew query: ", query)

                # Save character substitution TODO may need to save immediately previous subbed index as well
                subbed[num_mismatches] = sub

                num_mismatches += 1

                # # Align to genome or isoform TODO: Figure out if this is an accurate way to choose characters
                # # ex when we have multiple alignments up to this point. Maybe we need to keep track of the longest alignment?
                # # --> This way we can select the last mismatched character for the most promising/longest alignment instead of the
                # # first in the range.
                # if isoform_id:
                #     j = self._isoforms[isoform_id][1][location_range[0]] # location of first aligned char in suffix array
                #     query[i] = self._isoforms[isoform_id][4][j] # Set query mismatch to reference text value at corresponding index
                # else:
                #     j = self._sa[location_range[0]] # location of first aligned char in suffix array
                #     query[i] = self._genome_seq[j] # Set query mismatch to reference text value at corresponding index



        return None, None

    def recursive_inexact_alignment(self, p, M, occ, num_mismatches, isoform_id=None):
        """
        Run a greedy inexact alignment algorithm to determine an alignment to the reference text
        with less than or equal to MAX_NUM_MISMATCHES. 

        p: Query string
        M: M array of the reference text
        occ: occ array of the reference text
        isoform_id: Isoform id to refer

        Returns an alignment to reference text if one can be found with less than MAX_NUM_MISMATCHES
        in the form (<index of match in reference text>, <number of mismatches>)

        Make entire transcriptome --> concat w delimiters --> keep track of limits within the isoform 
        array w/isoform id of every index in transcriptome
        """

        iterator = iter(M.keys())
        r = {} # Dictionary (index, [bases tried at index])
        nucleotides = ['A', 'C', 'T', 'G']
        while True:
            char = next(iterator)
            if char == p[-1]:
                sp = M[char]
                try:
                    ep = M[next(iterator)] - 1
                except StopIteration:
                    ep = len(occ[char]) - 1
                break
        if ep < sp: # Need to fix this.
            return ((None), 0)
        i = len(p) - 2
        while i >= 0:
            print(p[i])
            temp_sp = M[p[i]] + occ[p[i]][sp - 1]
            temp_ep = M[p[i]] + occ[p[i]][ep] - 1
            while temp_sp > temp_ep:
                if i in r.keys():
                    r[i].append(p[i])
                else:
                    r[i] = [p[i]]
                    num_mismatches += 1
                if len(r[i]) == 4:
                    return (None, None)
                else:
                    p = p[:i] + str(random.choice(list(set(nucleotides) - set(r[i])))) + p[i+1:]
                    temp_sp = M[p[i]] + occ[p[i]][sp - 1]
                    temp_ep = M[p[i]] + occ[p[i]][ep] - 1
            sp, ep = temp_sp, temp_ep
            i -= 1
        return ((sp, ep + 1), num_mismatches)

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
            ifs_index = sa[best_match[1]] # index of first match in fully concatenated isoform string
        
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
        locations, num_mismatches = self.recursive_inexact_alignment(read_sequence, self._m, self._occ)

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

        # Making the seeds.
        l = len(read_sequence)
        n = l / 4
        i = 0        
        seeds = []
        while i < l:
            if i + n < l:
                seeds.append(read_sequence[i:i + n])
            else:
                seeds.append(read_sequence[i:l])

        # Aligning the seeds against the genome.
        alignments = []
        for seed in seeds:
            alignments.append(self.greedy_inexact_alignment(seed, self._m, self._occ))

        # Making sure that they are in the right order, need intron size range.
        MIN_INTRON_SIZE = 20
        MAX_INTRON_SIZE = 10000


        for i in range(0, len(alignments)):
            for j in range(i + 1, len(alignments)):
                # Making sure it's the difference between the end of the first and beginning of second.
                diff = alignments[j][0] - alignments[i][0] + len(seeds[i])
                # Doing this accounts for the seeds being in the necessary order.
                if diff > MIN_INTRON_SIZE and diff < MAX_INTRON_SIZE:
                    # Need to figure out the ID stuff --> do we need to make new isoform?
                    exon1 = genome_sequence[alignments[i][0] - (l - len(seeds[i])):alignments[i][0] + len(seeds[i])]
                    exon2 = genome_sequence[alignments[j[0]]:alignments[j[0]] + l]
                    isoform_string = exon1 + exon2
                    # Now trying to align read to this new isoform.
                    val = self.greedy_inexact_alignment(read_sequence, get_m(isoform_string), get_occ(isoform_string))
                    if val is not None:
                        return val

        

        # Making offset seeds in order to find all unknown exons.
        i = 6
        reads.append(read_sequence[0:i])
        while i < l:
            if i + n < l:
                reads.append(read_sequence[i:i + n])
            else:
                reads.append(read_sequence[i:l])
        
        # Finding the reads.
        for read in reads:
            greedy_inexact_alignment

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
    