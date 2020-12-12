"""
Testing file. Parse test files into Gene, Isoform, and Exon objects and test 
implementation.
"""

import sys
import re
import time
import zlib
import signal
import pickle
# import unittest

import sufarray 
from multiprocessing import Process
from pprint import pprint
# from itertools import zip_longest, islice

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
    unknown = []
    filepath = 'genes.tab'
    with open(filepath) as fp:
        line = fp.readline().rstrip()
        gene_id = ""
        u_gene_id = ""
        curr_isoforms = []
        u_curr_isoforms = []
        isoform_id = ""
        u_isoform_id = ""
        curr_exons = []
        u_curr_exons = []
        # Iterate through all lines
        while line:
            # Split into array by whitespace
            line = line.split()

            # Check for known genes
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

            # Check for unknown genes
            if line[0] == "unknown_gene":
                if gene_id:
                    curr_isoforms.append(Isoform(isoform_id, curr_exons))
                    genes.append(Gene(gene_id, curr_isoforms))
                    isoform_id = ""
                    curr_isoforms = []
                if u_gene_id:
                    u_curr_isoforms.append(Isoform(u_isoform_id, u_curr_exons))
                    unknown.append(Gene(u_gene_id, u_curr_isoforms))
                    u_isoform_id = ""
                    u_curr_isoforms = []
                gene_id = ""
                u_gene_id = line[1]
                u_isoforms = line[2].split(";")
            elif line[0] == "unknown_isoform":
                if u_isoform_id:
                    u_curr_isoforms.append(Isoform(u_isoform_id, u_curr_exons))
                    u_curr_exons = []
                u_isoform_id = line[1]
                u_exons = line[2].split(";")
            elif line[0] == "unknown_exon":
                u_curr_exons.append(Exon(line[1], int(line[2]), int(line[3])))

            line = fp.readline().rstrip()

        if u_gene_id:
            u_curr_isoforms.append(Isoform(u_isoform_id, u_curr_exons))
            unknown.append(Gene(u_gene_id, u_curr_isoforms))
            u_isoform_id = ""
            u_curr_isoforms = []

    return genes, unknown

def parse_genome():
    # Dict of reads
    genome = ""
    filepath = 'genome.fa'
    with open(filepath) as fp:
        index = fp.readline().rstrip()
        genome = fp.readline().rstrip()
    return genome

def timeout_handler(signum, frame):
    raise Exception("function timeout")

def test_init(sequence, genes):
    # Test initialization time
    # signal.signal(signal.SIGALRM, timeout_handler)
    # signal.alarm(500)
    # try:
    start_time = time.time()
    aligner = Aligner(sequence, genes)
    print("--- %s seconds ---" % (time.time() - start_time))
    # except Exception as exc:
    #     print(exc)
    
    return aligner

def test_greedy(aligner, read_sequence):
    # Test align time
    start_time = time.time()

    # isoform_id = list(aligner._isoforms.keys())[0]
    # isoform = aligner._isoforms[isoform_id]
    m, occ, sa =  aligner._m, aligner._occ, aligner._sa
    locations, num_mismatches = aligner.recursive_inexact_alignment(read_sequence, m, occ, 0, None)
    print("\nlocations: ", locations, "\nnum mismatches: ", num_mismatches)
    index = sa[locations[0]-1]
    print(aligner._genome_seq[index:index+len(read_sequence)])
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

def test_exact_matches_two(p, genome_sequence):
    # Build SA, M, OCC for whole genome
    _sa = get_suffix_array(genome_sequence)
    _sa_r = get_suffix_array(genome_sequence[::-1])
    _bwt = get_bwt(genome_sequence, _sa)
    _bwt_r = get_bwt(genome_sequence[::-1], _sa_r)
    _m = get_M(get_F(_bwt))
    _m_r = get_M(get_F(_bwt_r))
    _occ = get_occ(_bwt)
    _occ_r = get_occ(_bwt_r)
    _genome_seq = genome_sequence
    print("\ntesting exact suffix matches")
    locations = exact_suffix_matches(p, _m, _occ)
    print("locations: ", locations)

def test_align_to_transcriptome(aligner, p):
    return aligner.align_to_transcriptome(p)

def test_align_seeds(aligner, p):
    m, occ =  aligner._m, aligner._occ
    values = aligner.align_seeds(p)

def test_align(aligner, p):
    return aligner.align(p)

def test_align_to_genome(aligner, p):
    return aligner.align_to_genome(p)

def get_diff(s1, s2):
    diff = 0
    for i,s in enumerate(difflib.ndiff(s1, s2)):
        if s[0]==' ': continue
        elif s[0]=='-':
            diff += 1
            print(u'Delete "{}" from position {}'.format(s[-1],i))
        elif s[0]=='+':
            diff += 1
            print(u'Add "{}" to position {}'.format(s[-1],i))
    return diff

example = "GAAAAGGACATTTGTAATAGGAATTTTAAGATTTTATTACTTCTTCTTTCATCTCCTCAAAGGAATGAAGAAATGTTGGTGGGCAGATCTTCCTTTCCCATGCAGTATTTTCCAATTAATAATACTTCTTTACCTCAAGAATATTTTCTCTTTCAGAAGATGAGTGGCATGTGTATAATCCAGTGCTGCAGCTTCCTTACTATGAAATGACAGCTCCACTTCCTAATAGTGCATCCGTGTCTTCCTCACTGAATCATGTTCCAGATCTTGAGGCTGGACCCAGCTCATATAAATGGACTCACCAACAACCAAGTGACTCTGACCTTTATCAGATGACAGCTCCACTTCCTAATAGTGCATCCGTGTCTTCCTCACTGAATCATGTTCCAGATCTTGAGGCTGGACCCAGCTCATATAAATGGACTCACCAACAACCAAGTGACTCTGACCTTTATCAGATGAATAAACGAAAGAGACAAAAGCAAACAAGTGATAGTGATAGTAGCACAGACAACAACAGAGGCAACGAATGTAGCCAAAAGTTCCGAAAGTCTAAGAAGAAGAAAAGATACTAGTATTACCTACAAATGAAACTTACCTACACTGATCTTAGTTCTCTTATGAAAAAAATAAGATGTTATCCCATCAAATAAACAATGTCATGGCGAAAATGACATTTGTAATAGGAATTTTAAGATTTTATTACTTCTTCTTTCATCTCCTCAAAGGAATGAAGAAATGTTGGTGGGCAGATCTTCCTTTCCCATGCAGTATTTTCCAATTAATAATACTTCTTTACCTCAAGAATATTTTCTCTTTCAGAAGATGAGTGGCATGTGTATAATCCAGTGCTGCAGCTTCCTTACTATGAAATGACAGCTCCACTTCCTAATAGTGCATCCGTGTCTTCCTCACTGAATCATGTTCCAGATCTTGAGGCTGGACCCAGCTCATATAAATGGACTCACCAACAACCAAGTGACTCTGACCTTTATCAGATGACAGCTCCACTTCCTAATAGTGCATCCGTGTCTTCCTCACTGAATCATGTTCCAGATCTTGAGGCTGGACCCAGCTCATATAAATGGACTCACCAACAACCAAGTGACTCTGACCTTTATCAGATGAATAAACGAAAGAGACAAAAGCAAACAAGTGATAGTGATAGTAGCACAGACAACAACAGAGGCAACGAATGTAGCCAAAAGTTCCGAAAGTCTAAGAAGAAGAAAAGATACTAGTATTACCTACAAATGAAACTTACCTACACTGATCTTAGTTCTCTTATGAAAAAAATAAGATGTTATCCCATCAAATAAACAATGTCATGGCGAAAATGACATTTGTAATAGGAATTTTAAGATTTTATTACTTCTTCTTTCATCTCCTCAAAGGAATGAAGAAATGTTGGTGGGCAGATCTTCCTTTCCCATGCAGTATTTTCCAATTAATAATACTTCTTTACCTCAAGAATATTTTCTCTTTCAGAAGATGAGTGGCATGTGTATAATCCAGTGCTGCAGCTTCCTTACTATGAAATGACAGCTCCACTTCCTAATAGTGCATCCGTGTCTTCCTCACTGAATCATGTTCCAGATCTTGAGGCTGGACCCAGCTCATATAAATGGACTCACCAACAACCAAGTGACTCTGACCTTTATCAGATGACAGCTCCACTTCCTAATAGTGCATCCGTGTCTTCCTCACTGAATCATGTTCCAGATCTTGAGGCTGGACCCAGCTCATATAAATGGACTCACCAACAACCAAGTGACTCTGACCTTTATCAGATGAATAAACGAAAGAGACAAAAGCAAACAAGTGATAGTGATAGTAGCACAGACAACAACAGAGGCAACGAATGTAGCCAAAAGTTCCGAAAGTCTAAGAAGAAGAAAAGATACTAGTATTACCTACAAATGAAACTTACCTACACTGATCTTAGTTCTCTTATGAAAAAAATAAGATGTTATCCCATCAAATAAACAATGTCATGGCGAAAATGACATTTGTAATAGGAATTTTAAGATTTTATTACTTCTTCTTTCATCTCCTCAAAGGAATGAAGAAATGTTGGTGGGCAGATCTTCCTTTCCCATGCAGTATTTTCCAATTAATAATACTTCTTTACCTCAAGAATATTTTCTCTTTCAGAAGATGAGTGGCATGTGTATAATCCAGTGCTGCAGCTTCCTTACTATGAAATGACAGCTCCACTTCCTAATAGTGCATCCGTGTCTTCCTCACTGAATCATGTTCCAGATCTTGAGGCTGGACCCAGCTCATATAAATGGACTCACCAACAACCAAGTGACTCTGACCTTTATCAGATGACAGCTCCACTTCCTAATAGTGCATCCGTGTCTTCCTCACTGAATCATGTTCCAGATCTTGAGGCTGGACCCAGCTCATATAAATGGACTCACCAACAACCAAGTGACTCTGACCTTTATCAGATGAATAAACGAAAGAGACAAAAGCAAACAAGTGATAGTGATAGTAGCACAGACAACAACAGAGGCAACGAATGTAGCCAAAAGTTCCGAAAGTCTAAGAAGAAGAAAAGATACTAGTATTACCTACAAATGAAACTTACCTACACTGATCTTAGTTCTCTTATGAAAAAAATAAGATGTTATCCCATCAAATAAACAATGTCATGGCGAAAATGACATTTGTAATAGGAATTTTAAGATTTTATTACTTCTTCTTTCATCTCCTCAAAGGAATGAAGAAATGTTGGTGGGCAGATCTTCCTTTCCCATGCAGTATTTTCCAATTAATAATACTTCTTTACCTCAAGAATATTTTCTCTTTCAGAAGATGAGTGGCATGTGTATAATCCAGTGCTGCAGCTTCCTTACTATGAAATGACAGCTCCACTTCCTAATAGTGCATCCGTGTCTTCCTCACTGAATCATGTTCCAGATCTTGAGGCTGGACCCAGCTCATATAAATGGACTCACCAACAACCAAGTGACTCTGACCTTTATCAGATGACAGCTCCACTTCCTAATAGTGCATCCGTGTCTTCCTCACTGAATCATGTTCCAGATCTTGAGGCTGGACCCAGCTCATATAAATGGACTCACCAACAACCAAGTGACTCTGACCTTTATCAGATGAATAAACGAAAGAGACAAAAGCAAACAAGTGATAGTGATAGTAGCACAGACAACAACAGAGGCAACGAATGTAGCCAAAAGTTCCGAAAGTCTAAGAAGAAGAAAAGATACTAGTATTACCTACAAATGAAACTTACCTACACTGATCTTAGTTCTCTTATGAAAAAAATAAGATGTTATCCCATCAAATAAACAATGTCATGGCGAAAATGACATTTGTAATAGGAATTTTAAGATTTTATTACTTCTTCTTTCATCTCCTCAAAGGAATGAAGAAATGTTGGTGGGCAGATCTTCCTTTCCCATGCAGTATTTTCCAATTAATAATACTTCTTTACCTCAAGAATATTTTCTCTTTCAGAAGATGAGTGGCATGTGTATAATCCAGTGCTGCAGCTTCCTTACTATGAAATGACAGCTCCACTTCCTAATAGTGCATCCGTGTCTTCCTCACTGAATCATGTTCCAGATCTTGAGGCTGGACCCAGCTCATATAAATGGACTCACCAACAACCAAGTGACTCTGACCTTTATCAGATGACAGCTCCACTTCCTAATAGTGCATCCGTGTCTTCCTCACTGAATCATGTTCCAGATCTTGAGGCTGGACCCAGCTCATATAAATGGACTCACCAACAACCAAGTGACTCTGACCTTTATCAGATGAATAAACGAAAGAGACAAAAGCAAACAAGTGATAGTGATAGTAGCACAGACAACAACAGAGGCAACGAATGTAGCCAAAAGTTCCGAAAGTCTAAGAAGAAGAAAAGATACTAGTATTACCTACAAATGAAACTTACCTACACTGATCTTAGTTCTCTTATGAAAAAAATAAGATGTTATCCCATCAAATAAACAATGTCATGGCGAAAATGACATTTGTAATAGGAATTTTAAGATTTTATTACTTCTTCTTTCATCTCCTCAAAGGAATGAAGAAATGTTGGTGGGCAGATCTTCCTTTCCCATGCAGTATTTTCCAATTAATAATACTTCTTTACCTCAAGAATATTTTCTCTTTCAGAAGATGAGTGGCATGTGTATAATCCAGTGCTGCAGCTTCCTTACTATGAAATGACAGCTCCACTTCCTAATAGTGCATCCGTGTCTTCCTCACTGAATCATGTTCCAGATCTTGAGGCTGGACCCAGCTCATATAAATGGACTCACCAACAACCAAGTGACTCTGACCTTTATCAGATGACAGCTCCACTTCCTAATAGTGCATCCGTGTCTTCCTCACTGAATCATGTTCCAGATCTTGAGGCTGGACCCAGCTCATATAAATGGACTCACCAACAACCAAGTGACTCTGACCTTTATCAGATGAATAAACGAAAGAGACAAAAGCAAACAAGTGATAGTGATAGTAGCACAGACAACAACAGAGGCAACGAATGTAGCCAAAAGTTCCGAAAGTCTAAGAAGAAGAAAAGATACTAGTATTACCTACAAATGAAACTTACCTACACTGATCTTAGTTCTCTTATGAAAAAAATAAGATGTTATCCCATCAAATAAACAATGTCATGGCGAAAATGACATTTGTAATAGGAATTTTAAGATTTTATTACTTCTTCTTTCATCTCCTCAAAGGAATGAAGAAATGTTGGTGGGCAGATCTTCCTTTCCCATGCAGTATTTTCCAATTAATAATACTTCTTTACCTCAAGAATATTTTCTCTTTCAGAAGATGAGTGGCATGTGTATAATCCAGTGCTGCAGCTTCCTTACTATGAAATGACAGCTCCACTTCCTAATAGTGCATCCGTGTCTTCCTCACTGAATCATGTTCCAGATCTTGAGGCTGGACCCAGCTCATATAAATGGACTCACCAACAACCAAGTGACTCTGACCTTTATCAGATGACAGCTCCACTTCCTAATAGTGCATCCGTGTCTTCCTCACTGAATCATGTTCCAGATCTTGAGGCTGGACCCAGCTCATATAAATGGACTCACCAACAACCAAGTGACTCTGACCTTTATCAGATGAATAAACGAAAGAGACAAAAGCAAACAAGTGATAGTGATAGTAGCACAGACAACAACAGAGGCAACGAATGTAGCCAAAAGTTCCGAAAGTCTAAGAAGAAGAAAAGATACTAGTATTACCTACAAATGAAACTTACCTACACTGATCTTAGTTCTCTTATGAAAAAAATAAGATGTTATCCCATCAAATAAACAATGTCATGGCGAAAATGACATTTGTAATAGGAATTTTAAGATTTTATTACTTCTTCTTTCATCTCCTCAAAGGAATGAAGAAATGTTGGTGGGCAGATCTTCCTTTCCCATGCAGTATTTTCCAATTAATAATACTTCTTTACCTCAAGAATATTTTCTCTTTCAGAAGATGAGTGGCATGTGTATAATCCAGTGCTGCAGCTTCCTTACTATGAAATGACAGCTCCACTTCCTAATAGTGCATCCGTGTCTTCCTCACTGAATCATGTTCCAGATCTTGAGGCTGGACCCAGCTCATATAAATGGACTCACCAACAACCAAGTGACTCTGACCTTTATCAGATGACAGCTCCACTTCCTAATAGTGCATCCGTGTCTTCCTCACTGAATCATGTTCCAGATCTTGAGGCTGGACCCAGCTCATATAAATGGACTCACCAACAACCAAGTGACTCTGACCTTTATCAGATGAATAAACGAAAGAGACAAAAGCAAACAAGTGATAGTGATAGTAGCACAGACAACAACAGAGGCAACGAATGTAGCCAAAAGTTCCGAAAGTCTAAGAAGAAGAAAAGATACTAGTATTACCTACAAATGAAACTTACCTACACTGATCTTAGTTCTCTTATGAAAAAAATAAGATGTTATCCCATCAAATAAACAATGTCATGGCGAAAATGACATTTGTAATAGGAATTTTAAGATTTTATTACTTCTTCTTTCATCTCCTCAAAGGAATGAAGAAATGTTGGTGGGCAGATCTTCCTTTCCCATGCAGTATTTTCCAATTAATAATACTTCTTTACCTCAAGAATATTTTCTCTTTCAGAAGATGAGTGGCATGTGTATAATCCAGTGCTGCAGCTTCCTTACTATGAAATGACAGCTCCACTTCCTAATAGTGCATCCGTGTCTTCCTCACTGAATCATGTTCCAGATCTTGAGGCTGGACCCAGCTCATATAAATGGACTCACCAACAACCAAGTGACTCTGACCTTTATCAGATGACAGCTCCACTTCCTAATAGTGCATCCGTGTCTTCCTCACTGAATCATGTTCCAGATCTTGAGGCTGGACCCAGCTCATATAAATGGACTCACCAACAACCAAGTGACTCTGACCTTTATCAGATGAATAAACGAAAGAGACAAAAGCAAACAAGTGATAGTGATAGTAGCACAGACAACAACAGAGGCAACGAATGTAGCCAAAAGTTCCGAAAGTCTAAGAAGAAGAAAAGATACTAGTATTACCTACAAATGAAACTTACCTACACTGATCTTAGTTCTCTTATGAAAAAAATAAGATGTTATCCCATCAAATAAACAATGTCATGGCGAAAATGACATTTGTAATAGGAATTTTAAGATTTTATTACTTCTTCTTTCATCTCCTCAAAGGAATGAAGAAATGTTGGTGGGCAGATCTTCCTTTCCCATGCAGTATTTTCCAATTAATAATACTTCTTTACCTCAAGAATATTTTCTCTTTCAGAAGATGAGTGGCATGTGTATAATCCAGTGCTGCAGCTTCCTTACTATGAAATGACAGCTCCACTTCCTAATAGTGCATCCGTGTCTTCCTCACTGAATCATGTTCCAGATCTTGAGGCTGGACCCAGCTCATATAAATGGACTCACCAACAACCAAGTGACTCTGACCTTTATCAGATGACAGCTCCACTTCCTAATAGTGCATCCGTGTCTTCCTCACTGAATCATGTTCCAGATCTTGAGGCTGGACCCAGCTCATATAAATGGACTCACCAACAACCAAGTGACTCTGACCTTTATCAGATGAATAAACGAAAGAGACAAAAGCAAACAAGTGATAGTGATAGTAGCACAGACAACAACAGAGGCAACGAATGTAGCCAAAAGTTCCGAAAGTCTAAGAAGAAGAAAAGATACTAGTATTACCTACAAATGAAACTTACCTACACTGATCTTAGTTCTCTTATGAAAAAAATAAGATGTTATCCCATCAAATAAACAATGTCATGGCGAAAATGACATTTGTAATAGGAATTTTAAGATTTTATTACTTCTTCTTTCATCTCCTCAAAGGAATGAAGAAATGTTGGTGGGCAGATCTTCCTTTCCCATGCAGTATTTTCCAATTAATAATACTTCTTTACCTCAAGAATATTTTCTCTTTCAGAAGATGAGTGGCATGTGTATAATCCAGTGCTGCAGCTTCCTTACTATGAAATGACAGCTCCACTTCCTAATAGTGCATCCGTGTCTTCCTCACTGAATCATGTTCCAGATCTTGAGGCTGGACCCAGCTCATATAAATGGACTCACCAACAACCAAGTGACTCTGACCTTTATCAGATGACAGCTCCACTTCCTAATAGTGCATCCGTGTCTTCCTCACTGAATCATGTTCCAGATCTTGAGGCTGGACCCAGCTCATATAAATGGACTCACCAACAACCAAGTGACTCTGACCTTTATCAGATGAATAAACGAAAGAGACAAAAGCAAACAAGTGATAGTGATAGTAGCACAGACAACAACAGAGGCAACGAATGTAGCCAAAAGTTCCGAAAGTCTAAGAAGAAGAAAAGATACTAGTATTACCTACAAATGAAACTTACCTACACTGATCTTAGTTCTCTTATGAAAAAAATAAGATGTTATCCCATCAAATAAACAATGTCATGGC$"
example2 = "GAAAATGACATTTGTAATAGGAATTTTAAGATTTTATTACTTCTTCTTTCATCTCCTCAAAGGAATGAAGAAATGTTGGTGGGCAGATCTTCCTTTCCCATGCAGTATTTTCCAATTAATAATACTTCTTTACCTCAAGAATA$"
test_iso = "GAAAATGACATTTGTAATAGGAATTTTAAGATTTTATTACTTCTTCTTTCATCTCCTCAAAGGAATGAAGAAATGTTGGTGGGCAGATCTTCCTTTCCCATGCAGTATTTTCCAATTAATAATACTTCTTTACCTCAAGAATATTTTCTCTTTCAGAAGATGAGTGGCATGTGTATAATCCAGTGCTGCAGCTTCCTTACTATGAAATGACAGCTCCACTTCCTAATAGTGCATCCGTGTCTTCCTCACTGAATCATGTTCCAGATCTTGAGGCTGGACCCAGCTCATATAAATGGACTCACCAACAACCAAGTGACTCTGACCTTTATCAGATGACAGCTCCACTTCCTAATAGTGCATCCGTGTCTTCCTCACTGAATCATGTTCCAGATCTTGAGGCTGGACCCAGCTCATATAAATGGACTCACCAACAACCAAGTGACTCTGACCTTTATCAGATGAATAAACGAAAGAGACAAAAGCAAACAAGTGATAGTGATAGTAGCACAGACAACAACAGAGGCAACGAATGTAGCCAAAAGTTCCGAAAGTCTAAGAAGAAGAAAAGATACTAGTATTACCTACAAATGAAACTTACCTACACTGATCTTAGTTCTCTTATGAAAAAAATAAGATGTTATCCCATCAAATAAACAATGTCATGGC$"

# @profile
def main():
    print("Testing Sequence alignment")

    print("\nParsing reads")
    reads = parse_reads()

    print("\nParsing genome")
    genome = parse_genome()

    print("\nParsing genes")
    genes, unknown = parse_genes()

    # First known gene
    first = genes[0]
    iso = ''
    indices = []
    i = 0
    for exon in first.isoforms[0].exons:
        if i < 400:
        # print((exon.start, exon.end-exon.start)) in genome
            indices.append((i, exon.start,  exon.end-exon.start - 1))
            iso += genome[exon.start:exon.end]
            i += exon.end-exon.start
    # print(first.id, ": ", iso)
    # print(indices)

    # First unknown gene
    u_first = unknown[0]
    u_iso = ''
    u_indices = []
    i = 0
    for exon in u_first.isoforms[0].exons:
        if i < 600:
        # print((exon.start, exon.end-exon.start)) in genome
            u_indices.append((i, exon.start,  exon.end-exon.start - 1))
            u_iso += genome[exon.start:exon.end]
            i += exon.end-exon.start
    # print(u_first.id, ": ", u_iso)
    # print(u_indices)

    # print("\nTesting suffix array")
    # start_time = time.time()
    # suffix_array = get_suffix_array(example)
    # print("--- %s seconds ---" % (time.time() - start_time))
    # print("\nTesting thiers")
    # start_time = time.time()
    # sa_theirs = sufarray.SufArray(example).get_array()
    # print(len(suffix_array))
    # print(len(sa_theirs))
    # # diff =  set(sa_theirs) - set(suffix_array)
    # # print(list(diff))
    # # print(suffix_array)
    # print("--- %s seconds ---" % (time.time() - start_time))
    # print("accurate? ", suffix_array == sa_theirs)


    # print("Pickle")
    # aligner = test_init(genome, genes)
    # pickle.dump(aligner, open('pickle/aligner_onject', 'wb'))
    print("Unpickle")
    aligner = pickle.load(open('pickle/aligner_onject', 'rb'))

    print("\nTesting align to transcriptome: ")
    start_time = time.time()
    out = test_align(aligner,iso[:10]+'C'+iso[11:40]+'A'+iso[41:80]+'A'+iso[81:])
    # out = test_align(aligner,iso)    
    print("--- %s seconds ---" % (time.time() - start_time))
    print("true: ", indices)
    print("out: ", out)
    print(f'accurate? { indices == out}')

    print("\nTesting align to genome: ")
    ex = [(0, 50, 100)]
    start_time = time.time()
    inp = "AATTTCTTGGTGATGTGCACGTTTGTCACACGGAATTGAACCCTTCTTCTGATTGAGCAGTTTGGAATC"
    op = [(0, 0, len(inp))]
    # out = test_align(aligner, inp)
    # out = test_align(aligner, genome[ex[0][1]:ex[0][1]+ex[0][2]])
    out = test_align_to_genome(aligner, genome[ex[0][1]:ex[0][1]+ex[0][2]])
    print("--- %s seconds ---" % (time.time() - start_time))
    print("true: ", ex)
    print("out: ", out)
    print(f'accurate? { ex == out}')

    print("\nTesting align to unknown exons: ")
    start_time = time.time()
    out = test_align(aligner,u_iso[:10]+'C'+u_iso[11:40]+'A'+u_iso[41:80]+'A'+u_iso[81:])
    # out = test_align(aligner,u_iso)
    print("--- %s seconds ---" % (time.time() - start_time))
    print("true: ", u_indices)
    print("out: ", out)
    print(f'accurate? { u_indices == out}')

if __name__ == "__main__":
    main()