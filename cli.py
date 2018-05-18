import argparse
import os

class ArgumentParser(object):


    def __init__(self, args):
        self.args = args
        self.output = []
        self.exit_now = False
        self.parse()
        self.generate_output()


    def parse(self):
        """returns an argument parser."""
        parser = argparse.ArgumentParser()

        # parser.add_argument("--help", dest="help", required=False)
        parser.add_argument("--version", action="store_true", help="display program version.")
        parser.add_argument("--inDisk", action="store_true", help="use disk to transfer data between components.")
        parser.add_argument("--inMem", action="store_true", help="use memory to transfer data between components.")

        parser.add_argument("-p", type=int, dest="nJobs", help="number of parallel worker processes.")
        parser.add_argument("--strand_specific", action="store_true", help="Enable single-stranded mode")

        parser.add_argument("--fastq", action="store_true", help="Use FASTQ format")
        parser.add_argument("--fasta", action="store_true", help="Use FASTA format")

        parser.add_argument("-K", "--kmer-size", type=int, help="Kmer size")
        parser.add_argument("--partition", type=int, dest="partition_size", default=500, help="Partition size.")

        # XXX: think about using required=True
        parser.add_argument("-o", dest="output_dir", help="output directory. Must not exist.")

        parser.add_argument("--filter_FP", action="store_true", dest="filter_FP_flag",
                            help="Enable false-positive filtering. This option is incompatible with --inMem.")

        self.namespace = parser.parse_args(self.args)
        return self.namespace


    def generate_output(self):
        self.namespace.run_parallel = bool(self.namespace.nJobs)
        if self.namespace.nJobs:
            self.output.append('OPTIONS -p: Running parallel with ' + str(self.namespace.nJobs) + ' jobs.')

        if self.namespace.partition_size:
            self.output.append('OPTIONS --partition: Partition size set to ' + str(self.namespace.partition_size))

        if self.namespace.strand_specific:
            self.output.append('OPTIONS --strand_specific: Single-stranded mode enabled')
        
        if self.namespace.inMem:
            self.output.append('OPTIONS --inMem: In Memory mode enabled')
        
        if self.namespace.inDisk:
            self.output.append('OPTIONS --inDisk: In Memory mode disabled')

        if self.namespace.inMem and self.namespace.inDisk:
            self.namespace.inMem = False
            self.output.append('OPTIONS: Both --inMem and --inDisk specified. Using inDisk mode.')

        if self.namespace.output_dir:
            self.comp_directory_name = os.path.abspath(self.namespace.output_dir)
            if os.path.isdir(self.comp_directory_name) and os.listdir(self.comp_directory_name):
                self.output.append('ERROR: Output directory specified with -o needs to be an empty or non-existent directory')
                self.exit_now = True
        else:
            self.output.append('ERROR: Output directory needed. Use -o flag, which is mandatory.')
            self.exit_now = True

        if self.namespace.fastq:
            self.run_quorum = True
            self.fastq = True
            self.output.append('OPTIONS --fastq: Input is fastq format')

        if self.namespace.fasta:
            self.run_quorum = False
            self.fastq = False
            self.output.append('OPTIONS --fasta: Input is fasta format')

        if self.namespace.fasta and self.namespace.fastq:
            self.output.append("OPTIONS both --fasta and --fastq specified. Using FASTA.")

        if self.namespace.kmer_size:
            self.output.append('OPTIONS -K: Kmer size set to ' + str(self.namespace.kmer_size))

        if self.namespace.filter_FP_flag:
            self.output.append('OPTIONS --filter_FP: False-positive filtering enabled')
        
        if self.namespace.filter_FP_flag and self.namespace.inMem:
            self.output.append('OPTIONS --filter_FP: INCOMPATIBLE with inMem mode')

        return self.output