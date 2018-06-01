import unittest
from multiprocessing import Queue
from extension_correction import rc, rc_mate_ds

class ExtensionCorrectionTests(unittest.TestCase):

  def test_rc(self):
    q = Queue()
    lines = ['AGCT', 'TCGA']
    rc(lines, q)
    output = q.get()
    self.assertEquals(output, ['AGCT\n', 'TCGA\n'])

  def test_rc_mate_ds(self):
    q = Queue()
    reads1 = ['AGCT', 'TCGA']
    reads2 = ['AGCT', 'TCGA']
    rc_mate_ds(reads1, reads2, False, q)
    output = q.get()
    self.assertEquals(output, [['AGCT', 'AGCT\n'], ['TCGA', 'TCGA\n']])

  def test_rc_mate_ds_double_stranded(self):
    q = Queue()
    reads1 = ['AGCT', 'TCGA']
    reads2 = ['AGCT', 'TCGA']
    rc_mate_ds(reads1, reads2, True, q)
    output = q.get()
    self.assertEquals(output, [['AGCT', 'AGCT\n'], ['TCGA', 'TCGA\n'], ['AGCT', 'AGCT\n'], ['TCGA', 'TCGA\n']])

  