import unittest

from ch02_assignments import *

implanted_3mers = ['TTACCTTAAC',
                   'GATATCTGTC', 
                   'ACGGCGTTCG',
                   'CCCTAAAGAG',
                   'CGTCAGAGGT']

implanted_8mers = ['TCAGATGGGCTGCTTGCAGGTTTCTTTTGACCCGGGCCGCGCGGTAGCACCCCGAGGACGCTATCTGAGGGATAC',
                   'CAGTTAAGCAGTTTCCTTGCTCGCCGGAACCCGACTCGCAAGCCAACCGTTTAGTGGAAGGAACCCACCGCGGGG', 
                   'CCCCGGAGGACCTACGACTTGGGAGGGTAATGCACTTTTTCCACACAGTCGGTCCAAAGTCGGGAAGAACTTACC',
                   'CTGCCAGTCCAAGGTATCGTATAGACCGTCAGTAATTGTATACCACGGGGGGTGCAGCTCTGTGCCGGTCGGTGC',
                   'GTTCAGAACAATAGCCCCGGTGCATACGCGTGAAATAAATCCGCTAGCTTCGGTTTTTGCCGAGCAGCTCTTAAT']
 
class TestAssignment1(unittest.TestCase):
    def test_kd_neighbors(self): 
        assert type(kd_neighbors('A', 0)) == list
        assert len(kd_neighbors('ACG', 0)) == 1

    def test_motifenumeration(self):
        assert len(motif_enumeration(implanted_3mers, 3, 1)) == 18
        assert len(motif_enumeration(implanted_8mers, 8, 1)) == 1
        assert motif_enumeration(implanted_8mers, 8, 1) == set(['CCCCGGGG'])

motifs = ['AAAACC', 
          'GACAAA', 
          'AGACAA', 
          'CAGGAA',
          'ACAAGG']

class TestAssignment2(unittest.TestCase):
    def test_score(self):
        assert motif_score(motifs) == 12

    def test_profile(self):
        assert round(motif_profile(motifs)[0][0], 2) == 0.44
        assert round(motif_profile(motifs)[0][-1], 2) == 0.44
        assert round(motif_profile(motifs)[-1][-1], 2) == 0.11
                                       
    def test_entropy(self):
        p = motif_profile(motifs)
        assert round(profile_entropy(p), 5) == 11.01955

    def test_get_motifs(self):
        p = motif_profile(motifs)
        m = get_motifs_from_profile(p, implanted_8mers)
        assert m == ['AGCACC', 'AAGGAA', 'GAAGAA', 'CCAAGG', 'AGAACA']

if __name__ == '__main__':
    unittest.main()

