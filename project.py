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

# @profile
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
                # If length is greater than k, we want to recurse 
                # if len(suf) > k:
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
                for exon in isoform.exons: # If ordered, ok. if not, need to make sure
                    self._transcriptome += genome_sequence[exon.start:exon.end]
                self._transcriptome += '$'
                diff = len(self._transcriptome) - q
                transcriptome_index += diff               

                self._isoforms[isoform.id] = [isoform,                       # 0 - full isoform object
                                              start,                          # 1 - start index in stranscriptome string
                                              transcriptome_index - 1       # 2 - end index in transcriptome string (at $)
                                             ]     

        # FM Index for Transcriptome
        self._t_sa_r = get_suffix_array(self._transcriptome[::-1])
        self._t_sa = get_suffix_array(self._transcriptome)
        self._t_bwt = get_bwt(self._transcriptome, self._t_sa)
        self._t_bwt_r = get_bwt(self._transcriptome[::-1], self._t_sa_r)
        self._t_m = get_M(get_F(self._t_bwt))
        self._t_m_r = get_M(get_F(self._t_bwt_r))
        self._t_occ = get_occ(self._t_bwt)
        self._t_occ_r = get_occ(self._t_bwt_r)
                    
        # Build SA, M, OCC for whole genome
        self._sa = get_suffix_array(genome_sequence)
        self._sa_r = get_suffix_array(genome_sequence[::-1])
        self._bwt = get_bwt(genome_sequence, self._sa)
        self._bwt_r = get_bwt(genome_sequence[::-1], self._sa_r)
        self._m = get_M(get_F(self._bwt))
        self._m_r = get_M(get_F(self._bwt_r))
        self._occ = get_occ(self._bwt)
        self._occ_r = get_occ(self._bwt_r)
        self._genome_seq = genome_sequence
            
       
    def recursive_inexact_alignment(self, query, M, occ, isoform_id=None):
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
        num_mismatches=0
        r = {} # Dictionary (index, [bases tried at index])
        string_matched = ''
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

        string_matched = char
        i = len(p) - 2
        c = p[i]
        while i >= 0 and num_mismatches < MAX_NUM_MISMATCHES:
            match = [self._transcriptome[self._t_sa[x]:self._t_sa[x]+len(p) - i] for x in range(sp, ep+1)]
            c = p[i]
            temp_sp = M[c] + occ[c][sp - 1] # update rule
            temp_ep = M[c] + occ[c][ep] - 1 # update rule
            # if there are no matches, try to correct
            while (temp_sp > temp_ep): 
                if i in r.keys():
                    r[i].append(c)
                else:
                    r[i] = [c]
                    num_mismatches += 1
                # If we run out of characters and we cannot add more mismatches, exit
                if len(r[i]) == 4 and num_mismatches == MAX_NUM_MISMATCHES:
                    return ((None), 0)
                # Else, substitute in a random character into p and continue
                else:
                    c_prime = str(random.choice(list(set(BASES) - set(r[i]))))
                    p = p[:i] + c_prime + p[i+1:]
                    temp_sp = M[c_prime] + occ[c_prime][sp - 1]
                    temp_ep = M[c_prime] + occ[c_prime][ep] - 1
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
        range_for, num_mismatches_for = self.recursive_inexact_alignment(read_sequence, self._t_m, self._t_occ, ) # range in transcriptome
        range_rev, num_mismatches_rev = self.recursive_inexact_alignment(read_sequence[::-1], self._t_m_r, self._t_occ_r) # range in r transcriptome

        # Check for any alignments
        if range_for[0] or range_rev[0]:        
            # Select best alignment (f or r)
            if num_mismatches_for <= num_mismatches_rev and range_for[0]:
                # Read sequence indices in the transcriptome
                start_r_t = self._t_sa[range_for[0] - 1] 
                end_r_t = start_r_t + len(read_sequence)
            elif range_rev[0]:
                # Read sequence indices in the transcriptome
                start_r_t = self._t_sa_r[range_for[0] - 1]
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
                    # Progress to exon in isoform where read starts
                    # if length of offset greater than that of exon, go to next exon
                    if (exon.end - exon.start) < t_offset:
                        # subtract exon length from offset
                        t_offset -= (exon.end - exon.start)
                        continue
                    
                    # Progress to location in exon if any t_offset left
                    genome_index = exon.start + t_offset

                    # Check to get amount of read left
                    remaining_read = end_r_o - read_index # in original read
                    if remaining_read > (exon.end - genome_index):
                        length = exon.end - genome_index # length of this read
                    else:
                        length = remaining_read

                    # Add to matches
                    exon_matches.append((read_index, genome_index, length))

                    # Progress read_index
                    read_index += length
        
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
        range_rev, num_mismatches_rev = self.recursive_inexact_alignment(read_sequence[::-1], self._m_r, self._occ_r)
        print("forward start: ", self._sa[range_for[0]])
        print("reverse start: ", self._sa_r[range_rev[0]])
        if range_for[0] or range_rev[0]:     
            if num_mismatches_for <= num_mismatches_rev:
                # range_for[0] - 1
                location = self._sa[range_for[0]]
            else:
                location = self._sa_r[range_rev[0]]

            match = self._genome_seq[location:location+len(read_sequence)]
            
            diff = 0
            for i,s in enumerate(difflib.ndiff(match, read_sequence)):
                if s[0]==' ': continue
                elif s[0]=='-':
                    diff += 1
                    print(u'Delete "{}" from position {}'.format(s[-1],i))
                elif s[0]=='+':
                    diff += 1
                    print(u'Add "{}" to position {}'.format(s[-1],i))
            print("diff: ",  diff)
            
            return (0, location, len(read_sequence))

        return []
        
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
            n = int(l / 4)
            i = 0        
            seeds = []
            while i < l:
                if i + n < l:
                    seeds.append(read_sequence[i:i + n])
                    i += n
                else:
                    seeds.append(read_sequence[i:l])
                    i = l

            # Aligning the seeds against the genome.
            alignments = []
            for seed in seeds:
                alignments.append(self.recursive_inexact_alignment(seed, m, occ, 0))

            # Making sure that they are in the right order, need intron size range.
            for i in range(0, len(alignments)):
                for j in range(i + 1, len(alignments)):
                    # Making sure it's the difference between the end of the first and beginning of second.
                    diff = sa[alignments[j][0][0] - 1] - sa[alignments[i][0][0] - 1] + len(seeds[i])
                    # Doing this accounts for the seeds being in the necessary order.
                    if diff > MIN_INTRON_SIZE and diff < MAX_INTRON_SIZE:
                        # Need to figure out the ID stuff --> do we need to make new isoform?
                        exon1 = genome[sa[alignments[i][0][0] - 1] - (l - len(seeds[i])):sa[alignments[i][0][0] - 1] + len(seeds[i])]
                        exon2 = genome[sa[alignments[j][0][0] - 1]:sa[alignments[j][0][0] - 1] + l]
                        isoform_string = exon1 + exon2
                        # Now trying to align read to this new isoform.
                        return self.recursive_inexact_alignment(read_sequence, get_M(isoform_string), get_occ(isoform_string), 0)

        range_for, num_mismatches_for = exon_creator(read_sequence, self._genome_seq, self._m, self._occ, self._sa)
        range_rev, num_mismatches_rev = exon_creator(read_sequence[::-1], self._genome_seq[::-1], self._m_r, self._occ_r, self._sa_r)
        
        if range_for[0] or range_rev[0]:  
            if num_mismatches_for <= num_mismatches_rev:
                location = self._sa[range_for[0] - 1]
            else:
                location = self._sa_r[range_rev[0] - 1]
            
            return (0, location, len(read_sequence))
        
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
    