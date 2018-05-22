import unittest
from cli import ArgumentParser

class ArgParserTests(unittest.TestCase):
    """Unit tests for argument parsing."""

    def test(self):
      self.assertEqual(True, True)

    def test_version(self):
      args = ["--version"]
      parser = ArgumentParser(args)
      self.assertTrue(parser.namespace.version)
    
    def test_indisk(self):
        args = ["--inDisk"]
        parser = ArgumentParser(args)
        self.assertTrue(parser.namespace.inDisk)
        self.assertFalse(parser.namespace.inMem)  # if inDisk, the inMem must be False
    
    def test_inmemory(self):
        args = ["--inMem"]
        parser = ArgumentParser(args)
        self.assertTrue(parser.namespace.inMem)
        self.assertFalse(parser.namespace.inDisk)  # if inMem, then inDisk must be False

    def test_parallel(self):
        args = ["-p", "64"]
        parser = ArgumentParser(args)
        self.assertTrue(parser.namespace.run_parallel)
        self.assertEquals(parser.namespace.nJobs, 64)

    def test_njobs_default(self):
        """if -p is not passed, the parser defaults to 1 job."""
        args = []
        parser = ArgumentParser(args)
        self.assertEquals(1, parser.namespace.nJobs)
    
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
        parser = ArgumentParser(args)
        self.assertTrue(parser.namespace.strand_specific)

    def test_strand_specific_ouput(self):
        args = ["--strand_specific"]
        parser = ArgumentParser(args)
        self.assertIn("OPTIONS --strand_specific: Single-stranded mode enabled", parser.output)

    def test_output_dir(self):
        args = ["-o", "/blah"]
        parser = ArgumentParser(args)
        self.assertEquals(parser.namespace.output_dir, "/blah")

    def test_exit_now_if_output_dir_exists(self):
        """because / exists, exit_now should be True"""
        args = ["-o", "/"]
        parser = ArgumentParser(args)
        self.assertTrue(parser.exit_now)

    def test_not_exit_now_if_output_dir_doesnt_exists(self):
        """because (presumably) /foobar does not exist, exit_not should be False"""
        args = ["-o", "/foobar"]
        parser = ArgumentParser(args)
        self.assertFalse(parser.exit_now)
    
    def test_both_inmem_and_indisk(self):
        """if both inMem and inDisk are set, inDisk trumps."""
        args = ['--inMem', '--inDisk']
        parser = ArgumentParser(args)
        self.assertTrue(parser.namespace.inDisk)
        self.assertFalse(parser.namespace.inMem)

    def test_fastq(self):
        args = ["--fastq"]
        parser = ArgumentParser(args)
        self.assertTrue(parser.run_quorum)
        self.assertTrue(parser.fastq)

    def test_fasta(self):
        args = ["--fasta"]
        parser = ArgumentParser(args)
        self.assertFalse(parser.run_quorum)
        self.assertFalse(parser.fastq)

    def test_both_fasta_and_fastq(self):
        """If both --fasta and --fastq are specified, fasta trumps."""
        args = ["--fasta", "--fastq"]
        parser = ArgumentParser(args)
        self.assertFalse(parser.run_quorum)
        self.assertFalse(parser.fastq)

    def test_kmer_size(self):
        args = ["-K", "12"]
        parser = ArgumentParser(args)
        self.assertEquals(parser.namespace.kmer_size, 12)

        args = ["--kmer-size", "16"]
        parser = ArgumentParser(args)
        self.assertEquals(parser.namespace.kmer_size, 16)
    
    def test_filter_false_positives(self):
        args = ["--filter_FP"]
        parser = ArgumentParser(args)
        self.assertTrue(parser.namespace.filter_FP_flag)
    
    def test_filter_false_positives_inMem(self):
        args = ["--filter_FP", "--inMem"]
        parser = ArgumentParser(args)
        self.assertIn('OPTIONS --filter_FP: INCOMPATIBLE with inMem mode', parser.output)
    
    def test_filter_false_positives_default(self):
        args = []
        parser = ArgumentParser(args)
        self.assertFalse(parser.namespace.filter_FP_flag)

    def test_partition(self):
        args = ['--partition', '8']
        parser = ArgumentParser(args)
        self.assertEquals(parser.namespace.partition_size, 8)
    
    def test_partition_default_value(self):
        """if --partition is not specified, a default must be used."""
        args = []
        parser = ArgumentParser(args)
        self.assertTrue(parser.namespace.partition_size)  # check that a default is present, not the specific value.

    def test_single_stranded(self):
        args = ["-s"]
        parser = ArgumentParser(args)
        self.assertTrue(parser.namespace.single_stranded)

        args = ["--ss"]
        parser = ArgumentParser(args)
        self.assertTrue(parser.namespace.single_stranded)

    def test_defaults_to_double_stranded(self):
        """if neither -s or --ss are passed, the parser defaults to double_stranded mode."""
        args = []
        parser = ArgumentParser(args)
        self.assertFalse(parser.namespace.single_stranded)
    
    def test_only_reads(self):
        args = ['--only_reads']
        parser = ArgumentParser(args)
        self.assertTrue(parser.namespace.only_reads)
    
    def test_only_reads_default(self):
        """if --only_reads isn't set, it should default to False."""
        args = []
        parser = ArgumentParser(args)
        self.assertFalse(parser.namespace.only_reads)
