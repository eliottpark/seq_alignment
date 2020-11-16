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
    
    suffixes = {}
    for i in range(len(s)):
        suffixes[s[i:]] = i
    
    sorted_indices = [suffixes[x] for x in sorted(suffixes)]

    return sorted_indices

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
        for gene in known_genes:
            # Iterate through all isoforms for each gene
            for isoform in gene.isoforms:
                # Build isoform from individual exons
                full_isoform = ''
                for exon in isoform.exons: # If ordered, ok. if not, need to make sure
                    full_isoform += genome_sequence[exon.start:exon.end]
                # Find SA, M, and OCC of the isoform
                sa = get_suffix_array(full_isoform)
                self._isoforms[isoform.id] = [get_bwt(full_isoform, sa), sa,
                                                get_M(full_isoform), get_occ(full_isoform), 
                                                full_isoform, isoform]
                    
        # Build SA, M, OCC for whole genome
        self._sa = get_suffix_array(genome_sequence)
        self._bwt = get_bwt(genome_sequence, self._sa)
        self._m = get_M(genome_sequence)
        self._occ = get_occ(genome_sequence)
            
        # To make transcriptome: for each known_gene, find isoform; for each isoform, find exons; concatenate together.
        # Build suffix array, M array, Occ array.

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

        def simplified_bowtie(p, M, occ, isoform_id):
            """
            p: Query string
            T: BWT of the text (for isoforms, pass in BWT)
            K: max # of backtracks

            Iterate from last to 
            """
            nucleotides = ['A', 'C', 'T', 'G']
            R = {} # keep track of all the substitutions so far; R{i} is the set of characters already considered for position i
            num_mismatches = 0
            longest_match_length = 0
            query = p
            query_len = len(query)

            while num_mismatches < MAX_NUM_MISMATCHES:
                locations = exact_suffix_matches(query, M, occ)

                if locations[1] == query_len:
                    return locations, num_mismatches
                else:
                    i = query_len - locations[1] - 1
                    # Is there a way to access the index in the reference string just to set it equal?
                    # Line below is supposed to get the index number in the original isoform that our query doesn't align to
                    # so we can assign it directly to the query.
                    j = self._isoforms[isoform_id][1][locations[0][0]]
                    
                    query[i] = self._isoforms[isoform_id][4][j]

                    num_mismatches += 1
            return None, None

        # Prioritize matching to known isoforms. Minimize number of mismatches
        matches = [] # (isoform.id, locations, num_mismatches)
        # Matches are in the form (read start index, genome start index, length)
        for isoform in self._isoforms:
            locations, num_mismatches = simplified_bowtie(read_sequence, isoform[2], isoform[3])
            if locations not None and num_mismatches not None:
                # Want to add to list of matches and then find the best match.
                matches.append((isoform.id, locations, num_mismatches))

        
        # Check to see if there are any returned matches.
        if matches:
            # Finding the match with the least number of mismatches and converting it to the correct format.
            best_match = (0, 0, 6)
            for match in matches:
                if match[2] < best_match[2]:
                    best_match = match
            # Now need to convert into correct format.
            # Format: [(index of read, index of genome, length of match)]
            isoform = self._isoforms[best_match[0]]

            # Start and end indices for the isoform of interest
            sa = isoform[1]
            index_ifs = sa[best_match[1][0][0]] # index of first match in fully concatenated isoform string
            end = sa[best_match[1][0][1]] # index of last match in fully concatenated isoform string
        
            # check to see if start = end (only a single match)
            if index_ifs != end:
                # TODO figure out what happens in this situation
            # Iterate through exons to find indices
            exon_matches = []
            read_index = 0
            genome_index = 0
            offset = 0
            for exon in isoform[5].exons: # This should be ordered
                genome_index = index_ifs + exon.start
                if genome_index < exon.end:
                    length = exon.end - genome_index
                    exon_matches.append(read_index, genome_index, length)
                    read_index += length
                else:
                    length = exon.end - exon.start
                index_ifs += length
                                



                if index_ifs >= exon.start and index_ifs < exon.end:
                    exon_len = exon.end - start
                    exon_matches.append((read_index, ))

            


            

                



            
            
            
            num_backtracks = 0
            while num_backtracks < k:
                locations = exact_suffix_matches(query, M, occ)

                # If we improve our match length, then we have successfully introduced a mismatch
                if locations[1] > longest_match_length:
                    longest_match_length = locations[1]
                    best_match = locations
                    num_mismatches += 1
                    
                if locations[1] == 0: # might fold into else 
                    return [], 0
                elif locations[1] == query_len and num_mismatches <= MAX_NUM_MISMATCHES:
                    return locations, num_mismatches
                else:
                    i = query_len - locations[1] - 1
                    mismatched_char = query[i]
                    # leftmost just-visited position with minimal quality
                    j = i + 1
                    if j not in R.keys():
                        c = random.choice(set(nucleotides) - set(query[j]))
                        R[j] = set(c)
                        
                    else:
                        c = random.choice(set(nucleotides) - set(R[j]))
                        R[j] += c
                    query[j] = c # Introducing a mismatch

                    num_backtracks += 1                   
            return [], 0
                        

        




        pass
