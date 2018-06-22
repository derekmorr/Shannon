import unittest
from junkdrawer import JunkDrawer

class JunkDrawerTests(unittest.TestCase):

  def incr(self, x):
    return x + 1

  def test_argmax(self):
    output = JunkDrawer.argmax([1,2,3], self.incr)
    self.assertEqual(output, 3)