"""
Testing file. Parse test files into Gene, Isoform, and Exon objects and test 
implementation.
"""

import sys
import re
import time
# import unittest

from pprint import pprint
from itertools import zip_longest, islice

from shared import *
from project import *

def parse_reads():
    # Dict of reads
    reads = {}
    filepath = 'reads.fa'
    with open(filepath) as fp:
        line = fp.readline().rstrip()
        key = ""
        while line:
            if re.search("^>", line):
                key = line[1:]
            else:
                reads[key] = line
            line = fp.readline().rstrip()
    return reads

def parse_genes():
    # list of genes
    genes = []
    filepath = 'genes.tab'
    with open(filepath) as fp:
        line = fp.readline().rstrip()
        gene_id = ""
        curr_isoforms = []
        isoform_id = ""
        curr_exons = []
        # Iterate through all lines
        while line:
            # Split into array by whitespace
            line = line.split()

            if line[0] == "gene" or line[0] == "unknown_gene":
                if gene_id:
                    curr_isoforms.append(Isoform(isoform_id, curr_exons))
                    genes.append(Gene(gene_id, curr_isoforms))
                    isoform_id = ""
                    curr_isoforms = []
                gene_id = line[1]
                isoforms = line[2].split(";")
            elif line[0] == "isoform" or line[0] == "unknown_isoform":
                if isoform_id:
                    curr_isoforms.append(Isoform(isoform_id, curr_exons))
                    curr_exons = []
                isoform_id = line[1]
                exons = line[2].split(";")
            elif line[0] == "exon" or line[0] == "unknown_exon":
                curr_exons.append(Exon(line[1], int(line[2]), int(line[3])))

            line = fp.readline().rstrip()

        if gene_id:
            curr_isoforms.append(Isoform(isoform_id, curr_exons))
            genes.append(Gene(gene_id, curr_isoforms))
            isoform_id = ""
            curr_isoforms = []

    return genes

def parse_genome():
    # Dict of reads
    genome = ""
    filepath = 'genome.fa'
    with open(filepath) as fp:
        index = fp.readline().rstrip()
        genome = fp.readline().rstrip()
    return genome

def test_init(sequence, genes):
    # Test initialization time
    start_time = time.time()
    aligner = Aligner(sequence, genes)
    print("--- %s seconds ---" % (time.time() - start_time))
    
    return aligner

def test_greedy(aligner, read_sequence):
    # Test align time
    start_time = time.time()

    isoform_id = list(aligner._isoforms.keys())[0]
    isoform = aligner._isoforms[isoform_id]
    m, occ =  isoform[2], isoform[3]
    locations, num_mismatches = aligner.greedy_inexact_alignment(read_sequence, m, occ, isoform_id)
    print("\nlocations: ", locations, "\nnum mismatches: ", num_mismatches)
    print("--- %s seconds ---" % (time.time() - start_time))


def test_exact_matches(p, s):
    sa = get_suffix_array(s)
    L = get_bwt(s, sa)
    F = get_F(L)
    M = get_M(F)
    occ = get_occ(L)
    print("\ntesting exact suffix matches")
    locations = exact_suffix_matches(p, M, occ)
    print("locations: ", locations)


def main():
    print("Testing Sequence alignment")

    print("\nParsing reads")
    reads = parse_reads()

    print("\nParsing genes")
    genes = parse_genes()
    # pprint(genes)

    print("\nParsing genome")
    genome = parse_genome()
    

    print("\nTesting suffix array")
    start_time = time.time()
    suffix_array = get_suffix_array("GAAAATGACATTTGTAATAGGAATTTTAAGATTTTATTACTTCTTCTTTCATCTCCTCAAAGGAATGAAGAAATGTTGGTGGGCAGATCTTCCTTTCCCATGCAGTATTTTCCAATTAATAATACTTCTTTACCTCAAGAATATTTTCTCTTTCAGAAGATGAGTGGCATGTGTATAATCCAGTGCTGCAGCTTCCTTACTATGAAATGACAGCTCCACTTCCTAATAGTGCATCCGTGTCTTCCTCACTGAATCATGTTCCAGATCTTGAGGCTGGACCCAGCTCATATAAATGGACTCACCAACAACCAAGTGACTCTGACCTTTATCAGATGACAGCTCCACTTCCTAATAGTGCATCCGTGTCTTCCTCACTGAATCATGTTCCAGATCTTGAGGCTGGACCCAGCTCATATAAATGGACTCACCAACAACCAAGTGACTCTGACCTTTATCAGATGAATAAACGAAAGAGACAAAAGCAAACAAGTGATAGTGATAGTAGCACAGACAACAACAGAGGCAACGAATGTAGCCAAAAGTTCCGAAAGTCTAAGAAGAAGAAAAGATACTAGTATTACCTACAAATGAAACTTACCTACACTGATCTTAGTTCTCTTATGAAAAAAATAAGATGTTATCCCATCAAATAAACAATGTCATGGC$")
    print("--- %s seconds ---" % (time.time() - start_time))
    # print("suffix array len, max: ", len(suffix_array), max(suffix_array))

    dummy_genome = "AAAATGACATT$"

    print("\nTesting init function")
    aligner = test_init(genome, genes)

    print("\nTesting greedy: ")
    # test_greedy(aligner, "GCAA")

    test_exact_matches("GGGATG","AAAATGACATTAAAATGACATTAAAATGACATTAAAATGACATTAAAATGACATTAAAAATTGGGAC$")
    

if __name__ == "__main__":
    main()