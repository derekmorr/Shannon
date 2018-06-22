import unittest
from dna import DNA

class DNATests(unittest.TestCase):

  def test_reverse_complement(self):
    dna = DNA()
    output = dna.reverse_complement("ATGC")
    self.assertEquals("GCAT", output)
  
  def test_reverse_complement_no_n(self):
    dna = DNA()
    output = dna.reverse_complement_no_n("ATGC")
    self.assertEquals("GCAT", output)

  def test_reverse_complement_n(self):
    dna = DNA()
    output = dna.reverse_complement("GNA")
    self.assertEquals("TNC", output)
