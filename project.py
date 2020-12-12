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
import heapq
import time
import random
import difflib
import math

from shared import *
from test import get_diff

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
    k = 500
    r = 2
    suffixes = {}
    str_len = len(s)
    # Construct dictionary that represents buckets defined by kth prefix of each suffix
    for i in range(str_len):
        # Only split if the bucket will store data
        prefix = s[i:i + k]
        if prefix in suffixes.keys():
            suffixes[prefix] += [i]
        else:
             suffixes[prefix] = [i]

    # Sort long suffixes
    sorted_suffixes = k_radix_sort(s, suffixes, k, r)

    suffixes.clear()

    return sorted_suffixes

def k_radix_sort(s, suffixes, k, r):
    """
    Radix sort according to the kth prefix of each suffix.

    suffixes: a dictionary of suffixes {kth prefix: [(index in original string, [k:] suffix), ...]}
    k: size of prefix
    r: recursion depth

    Returns a sorted dictionary of suffixes {relative order of suffix: index in original string}
    """
    # Sort dictionary by keys and replace by order 
    s_len = len(suffixes)
    for i, key in enumerate(sorted(suffixes.keys())):
        suffixes[i] = suffixes.pop(key) # switch keys from kth prefix to index/order

    outer_suffixes = []
    # Recurse if there are more than one strings in a bucket
    for i in range(len(suffixes)):
        if len(suffixes[i]) > 1 and r > 0:
            inner_count = 0
            recurse_suffix = {}
            j = 0
            # Create suffix buckets from inner pairs
            for index in suffixes.pop(i):
                inner_count += 1
                suf = s[index + k*r:]
                pre = s[:k]
                if pre in recurse_suffix.keys():
                    recurse_suffix[pre] += [index]
                else:
                    recurse_suffix[pre] = [index]
            
            ##### Sort and add in recursive suffixes #####
            sorted_suffixes = k_radix_sort(s, recurse_suffix, k, r-1)
            outer_suffixes += sorted_suffixes
            recurse_suffix.clear()

        elif len(suffixes[i]) > 1:
            # Sort the bucket directly at maximum recursion depth
            max_r_suffixes = {}
            # Turn list of tuples (index in original string, remainder of string) into sorted dict
            for index in suffixes.pop(i):
                suf = s[index + k*r:]
                max_r_suffixes[suf] = index
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
MAX_INTRON_SIZE = 45000

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
        self._transcriptome = ''
        # Iterate through all genes
        # print("iterating thru genes")
        transcriptome_index = 0
        for gene in known_genes:
            # Iterate through all isoforms for each gene
            # print("Gene: ", gene.id)
            for isoform in gene.isoforms:
                start = transcriptome_index # start in transcripome string
                # Add to full transcriptome with delimiter
                q = len(self._transcriptome)
                temp = ''
                for exon in isoform.exons: 
                    temp += genome_sequence[exon.start:exon.end]
                temp += '$'
                transcriptome_index += len(temp)    
                self._transcriptome += temp          

                self._isoforms[isoform.id] = [isoform,                       # 0 - full isoform object
                                              start,                         # 1 - start (inclusive) index in stranscriptome string
                                              transcriptome_index            # 2 - end (exclusive) index in transcriptome string (at $)
                                             ]     

        # FM Index for Transcriptome
        self._t_sa = get_suffix_array(self._transcriptome)
        self._t_sa_r = get_suffix_array(self._transcriptome[::-1])
        self._t_bwt = get_bwt(self._transcriptome, self._t_sa)
        self._t_bwt_r = get_bwt(self._transcriptome[::-1], self._t_sa_r)
        self._t_m = get_M(get_F(self._t_bwt))
        self._t_m_r = get_M(get_F(self._t_bwt_r))
        self._t_occ = get_occ(self._t_bwt)
        self._t_occ_r = get_occ(self._t_bwt_r)
        self._transcriptome_len = len(self._transcriptome)
                    
        # Build SA, M, OCC for whole genome
        self._genome_seq = genome_sequence + '$'
        self._sa = get_suffix_array(self._genome_seq)
        self._sa_r = get_suffix_array(self._genome_seq[::-1])
        self._bwt = get_bwt(self._genome_seq, self._sa)
        self._bwt_r = get_bwt(self._genome_seq[::-1], self._sa_r)
        self._m = get_M(get_F(self._bwt))
        self._m_r = get_M(get_F(self._bwt_r))
        self._occ = get_occ(self._bwt)
        self._occ_r = get_occ(self._bwt_r)
        self._genome_seq_len = len(self._genome_seq)
            
       
    def recursive_inexact_alignment(self, query, M, occ, isoform_id=None, num_mismatches=0):
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
        r = {} # Dictionary (index, [bases tried at index])
        p = query

        # Count algorithm
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
        match = [self._genome_seq[self._sa[x]:self._sa[x]+len(p)-1-i] for x in range(sp, ep+1)]
        c = p[i] 
        while i >= 0 and num_mismatches < MAX_NUM_MISMATCHES:
            match = [self._genome_seq[self._sa[x]:self._sa[x]+len(p)-1-i] for x in range(sp, ep+1)]
            c = p[i]
            temp_sp = M[c] + occ[c][sp - 1] # update rule
            temp_ep = M[c] + occ[c][ep] - 1 # update rule
            # if there are no matches, try to correct. 
            # If no correction can be made, return nothing
            while (temp_sp > temp_ep) and num_mismatches < MAX_NUM_MISMATCHES: 
                # Check to see if character already subbed
                if i in r.keys():
                    r[i].append(c)
                else:
                    r[i] = [c]
                    num_mismatches += 1
                # If we run out of characters and we cannot add more mismatches, exit
                if len(r[i]) == 4 or num_mismatches == MAX_NUM_MISMATCHES:
                    return ((None), 0)
                # Else, substitute in a random character into p and continue
                else:
                    c = str(random.choice(list(set(BASES) - set(r[i]))))
                    p = p[:i] + c + p[i+1:]
                    temp_sp = M[c] + occ[c][sp - 1]
                    temp_ep = M[c] + occ[c][ep] - 1
            # finalize update of pointers
            sp, ep = temp_sp, temp_ep
            i -= 1
        return (sp, ep + 1), num_mismatches

    def align_to_transcriptome(self, read_sequence):
        """
        Returns the best alignment of the read sequence to the isoform database. 
        Prioritize matching to known isoforms. Minimize number of mismatches

        read_sequence: input read text to be aligned

        Return value must be in the form (also see the project pdf):
        [(<read_start_1>, <reference_start_1, length_1), (<read_start_2>, <reference_start_2, length_2), ...]

        If no good matches are found: return the best match you can find or return []
        """
        exon_matches = []

        # Align forward and reverse to transcriptome
        range_for, num_mismatches_for = self.recursive_inexact_alignment(read_sequence, self._t_m, self._t_occ ) # range in transcriptome
        range_rev, num_mismatches_rev = self.recursive_inexact_alignment(read_sequence[::-1], self._t_m_r, self._t_occ_r) # range in r transcriptome

        # Check for any alignments
        if num_mismatches_for <= num_mismatches_rev and range_for:
            # Read sequence indices in the transcriptome
            start_r_t = self._t_sa[range_for[0]] 
            end_r_t = start_r_t + len(read_sequence)
        elif num_mismatches_for >= num_mismatches_rev and range_rev:
            # Read sequence indices in the transcriptome
            start_r_t = self._t_sa_r[range_for[0]]
            end_r_t = start_r_t + len(read_sequence)
        else:
            return []

        # Read sequence indices in original read
        start_r_o, end_r_o = 0, len(read_sequence) - 1

        # Find corresponding isoform and location in isoform
        isoform_id, start_i_t, end_i_t = '', 0, 0 # isoform id, start of alignment in transcriptome, end of alignment in transcriptome
        for _id in self._isoforms.keys():
            s = self._isoforms[_id][1] # start of isoform in transcriptome string
            e = self._isoforms[_id][2] # end of isoform in transcriptome string
            # check to see if range in isoform
            if s <= start_r_t and end_r_t <= e:
                start_i_t = s # start of isoform of interest in the transcriptome
                end_i_t = e # end of isoform of interest in the transcriptome
                isoform_id = _id
                break

        if isoform_id:

            # Now need to convert into correct format.
            # Format: [(index of read, index of genome, length of match)]

            read_index = 0 # Current index in read
            genome_index = 0 # Current index in genome
            t_index = start_i_t # Current index in transcriptome 
            t_offset = start_r_t - start_i_t # offset left in isoform

            # Iterate thru exons
            for exon in self._isoforms[isoform_id][0].exons:
                remaining_read = end_r_o - read_index # amount remaining in original read
                if remaining_read > 0:
                    # Progress to exon in isoform where read starts
                    # if length of offset greater than that of exon, go to next exon
                    if (exon.end - exon.start - 1) < t_offset:
                        # subtract exon length from offset
                        t_offset -= (exon.end - exon.start - 1)
                        continue
                    
                    # Progress to location in exon if any t_offset left
                    genome_index = exon.start + t_offset

                    # Check to get amount of read left
                    if remaining_read > (exon.end - genome_index - 1):
                        length = exon.end - genome_index - 1# length of this read
                    else:
                        length = remaining_read

                    # Add to matches
                    if length > 0:
                        exon_matches.append((read_index, genome_index, length))

                    # Progress read_index
                    read_index += length + 1
    
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

        range_for, num_mismatches_for = self.recursive_inexact_alignment(read_sequence, self._m, self._occ)
        if not range_for or num_mismatches_for > 0:
            range_rev, num_mismatches_rev = self.recursive_inexact_alignment(read_sequence[::-1], self._m_r, self._occ_r)
            if num_mismatches_for <= num_mismatches_rev and range_for:
                return [(0, self._sa[range_for[0]], len(read_sequence))]
            elif num_mismatches_for >= num_mismatches_rev and range_rev:
                return [(0, self._genome_seq_len - self._sa_r[range_rev[0]], len(read_sequence))]
            else:
                return []
        return [(0, self._sa[range_for[0]], len(read_sequence))]
        
    def align_seeds(self, read_sequence):
        """
        Split into seeds and align.

        read_sequence: input read text to be aligned

        Return value must be in the form (also see the project pdf):
        [(<read_start_1>, <reference_start_1, length_1), (<read_start_2>, <reference_start_2, length_2), ...]

        If no good matches are found: return the best match you can find or return []
        """

        def exon_creator(read_sequence, genome, m, occ, sa):
            # Making the seeds.
            l = len(read_sequence)
            n = math.ceil(l /8)
            i = 0        
            seeds = []
            while i < l:
                if i + n < l:
                    seeds.append(read_sequence[i:i + n])
                    i += n
                else:
                    seeds.append(read_sequence[i:l])
                    i = l

            # Aligning the seeds against the genome. May need to keep track of overall # of mismatches
            alignments = []
            total_mismatches = 0
            for seed in seeds:
                seed_aligned = self.recursive_inexact_alignment(seed, m, occ, num_mismatches=total_mismatches)
                if seed_aligned[0]:
                    alignments.append(seed_aligned)
                    total_mismatches += seed_aligned[1]
                else:
                    return [], 0
                
            isoform_string = '' # overall isoform string for all seed matches
            read_index = 0
            matches = []
            next_index = 0
            curr_match_index, curr_match_len = -1, -1

            # Need to have all seeds aligned sequentially for a true alignment
            for i in range(0, len(alignments)):
                range_i = list(range(alignments[i][0][0], alignments[i][0][1]))
                # If we are at the final seed, append (check between seeds occurs at the prior seed index)
                if i == len(alignments) - 1:
                    isoform_string += genome[sa[range_i[j]]:(sa[range_i[j]] + len(seeds[i]))]
                    matches.append((read_index, sa[range_i[j]], len(seeds[i])))
                else:
                    # Test if this seed is close enough to the next
                    # Test for all alignments in ranges j, k
                    range_i_1 = list(range(alignments[i+1][0][0], alignments[i+1][0][1]))
                    offset1, offset2 = alignments[i+1][0][0], alignments[i][0][0]
                    if i == 0:
                        this_sa_start, this_sa_end = alignments[i][0][0], alignments[i][0][1]
                    else:
                        # if using previous index
                        this_sa_start, this_sa_end = next_index, next_index + 1
                    for k in range(this_sa_start, this_sa_end):
                        for j in range(alignments[i+1][0][0], alignments[i+1][0][1]):
                            diff = sa[range_i_1[j-offset1]] - (sa[range_i[k-offset2]] + len(seeds[i]))
                            if diff > MIN_INTRON_SIZE and diff < MAX_INTRON_SIZE:
                                # If the next is within the allowed distance, append current to list
                                start_index = sa[range_i[k-offset2]]
                                isoform_string += genome[sa[start_index]:(sa[start_index] + len(seeds[i]))]
                                # Check if next val should be appended together (dont add to matches yet)
                                if curr_match_len > 0:
                                    # Add to runing match
                                    curr_match_len += len(seeds[i])
                                else:
                                    # initialize running match
                                    curr_match_index = sa[start_index]
                                    curr_match_len = len(seeds[i])
                                # if we are not at an extension
                                if diff > 0:
                                    matches.append((read_index, curr_match_index, curr_match_len))
                                    read_index += len(seeds[i])
                                    next_index = j
                                    # Close running match
                                    curr_match_index, curr_match_len = 0, 0
                                break
                        else:
                            continue
                        break
                    else:
                        # we must align all seeds - if one is not close enought then we must exit
                        return [], 0
            
            return matches, self.recursive_inexact_alignment(read_sequence, get_M(isoform_string), get_occ(isoform_string))[1]

        matches_for, num_mismatches_for = exon_creator(read_sequence, self._genome_seq, self._m, self._occ, self._sa)
        matches_rev, num_mismatches_rev = exon_creator(read_sequence[::-1], self._genome_seq[::-1], self._m_r, self._occ_r, self._sa_r)
         
        if num_mismatches_for <= num_mismatches_rev and matches_for:
            return matches_for
        elif num_mismatches_for >= num_mismatches_rev and matches_rev:
            return matches_rev
        else:
            return []

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
        reads = self.align_to_transcriptome(read_sequence)

        # Align to genome if no isoform matches (case 2 in tophat2)
        if not reads:
            reads = self.align_to_genome(read_sequence)

        # If no isoform matches or direct genome matches, split into seeds and
        # align seed to try to find new splice junctions (case 3 in tophat2)
        # if not reads:
        #     reads = self.align_seeds(read_sequence)
                
        return reads