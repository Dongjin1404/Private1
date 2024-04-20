import unittest

from ch03_assignments import *

s1 = 'AGGCTAGTACGGACTTACGCACAACGCTTGCGGA'
s2 = 'AGGCTAGTACGGACTTGCGCACAAGGCTTGCGGA'
s3 = 'ACGTACGTACGTACGT'
 
class TestAssignment1(unittest.TestCase):
    def test_sorted_kmers(self):
        x = sorted_kmers(s1, 3)
        assert x == ['AAC', 'ACA', 'ACG', 'ACG', 
         'ACG', 'ACT', 'AGG', 'AGT',
         'CAA', 'CAC', 'CGC', 'CGC', 
         'CGG', 'CGG', 'CTA', 'CTT',
         'CTT', 'GAC', 'GCA', 'GCG',
         'GCT', 'GCT', 'GGA', 'GGA', 
         'GGC', 'GTA', 'TAC', 'TAC', 
         'TAG', 'TGC', 'TTA', 'TTG']

    def test_all_genomes(self):
        x = set(all_genomes(sorted_kmers(s1, 5)))
        assert len(x) == 2
        assert s1 in x
        x = set(all_genomes(sorted_kmers(s1, 6)))
        assert len(x) == 1
        assert s1 in x

class TestAssignment2(unittest.TestCase):
    def test_de_bruijn_graph(self):
        kmers = sorted_kmers(s3, 15)
        nodes, edges = get_deBruijn_graph(kmers)
        assert sorted(nodes) == ['ACGTACGTACGTAC', 'CGTACGTACGTACG', 'GTACGTACGTACGT']
        assert sorted(edges) == [('ACGTACGTACGTAC', 'CGTACGTACGTACG'), ('CGTACGTACGTACG', 'GTACGTACGTACGT')]

    def test_get_genome_deBrujin(self):
        kmers = sorted_kmers(s1, 6)
        path = get_genome_deBruijn(kmers)
        assert path == s1

if __name__ == '__main__':
    unittest.main()

