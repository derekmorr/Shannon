import unittest
from cli import ArgumentParser

class ArgParserTests(unittest.TestCase):
    """Unit tests for argument parsing."""

    def test(self):
      self.assertEqual(True, True)

    def test_version(self):
      args = ["--version"]
      namespace = ArgumentParser.parse(args)
      self.assertTrue(namespace.version)
    
    def test_indisk(self):
        args = ["--inDisk"]
        namespace = ArgumentParser.parse(args)
        self.assertTrue(namespace.inDisk)
        self.assertFalse(namespace.inMem)  # if inDisk, the inMem must be False
    
    def test_inmemory(self):
        args = ["--inMem"]
        namespace = ArgumentParser.parse(args)
        self.assertTrue(namespace.inMem)
        self.assertFalse(namespace.inDisk)  # if inMem, then inDisk must be False

    # def test_parallel(self):
    #     args = ["-p 64"]
    #     namespace = ArgumentParser.parse(args)
    #     # self.assertTrue(namespace.run_parallel)
    #     self.assertEquals(namespace.nJobs, 64)
    
    # def test_parallel_rejects_nonnumeric_input(self):
    #     args = ["-p two"]
    #     # self.assertRaises
    #     namespace = ArgumentParser.parse(args)
    #     self.assertEquals(namespace.nJobs, 64)

    # def test_left_and_right(self):
    #     args = ["--left", "left_file", "--right", "right_file"]
    #     namespace = ArgumentParser.parse(args)
    #     self.assertEquals(namespace.left, "left_file")
    #     self.assertEquals(namespace.right, "right_file")
    
    # def test_single(self):
    #     args = ["--single", "single_file"]
    #     namespace = ArgumentParser.parse(args)
    #     self.assertEquals(namespace.single, "single_file")

    # # def test_both_left_and_right_must_be_present(self):
    # #     """if one of '--left' or '--right' are present, both must be present."""
    # #     args = ["--left", "left_file"]
    # #     namespace = ArgumentParser.parse(args)
    # #     # XXX how to test SystemExit?

    def test_strand_specific(self):
        args = ["--strand_specific"]
        namespace = ArgumentParser.parse(args)
        self.assertTrue(namespace.strand_specific)

    def test_strand_specific_ouput(self):
        args = ["--strand_specific"]
        namespace = ArgumentParser.parse(args)
        output = ArgumentParser.generate_output(namespace)
        self.assertIn("OPTIONS --strand_specific: Single-stranded mode enabled", output)

    def test_output_dir(self):
        args = ["-o", "/blah"]
        namespace = ArgumentParser.parse(args)
        self.assertEquals(namespace.output, "/blah")