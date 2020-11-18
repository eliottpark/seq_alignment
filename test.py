"""
Testing file. Parse test files into Gene, Isoform, and Exon objects and test 
implementation.
"""

import sys
import re
import time
# import unittest

from pprint import pprint

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

    # Test isoforms, etc.


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

    print("\nTesting init function")
    # test_init(genome, genes)


    

if __name__ == "__main__":
    main()