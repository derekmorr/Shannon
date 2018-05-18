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

        parser.add_argument("-p", type=int, help="number of parallel worker processes.")
        parser.add_argument("--strand_specific", action="store_true", help="Enable single-stranded mode")

        # XXX: think about using required=True
        parser.add_argument("-o", dest="output_dir", help="output directory. Must not exist.")

        self.namespace = parser.parse_args(self.args)
        return self.namespace


    def generate_output(self):
        if self.namespace.strand_specific:
            self.output.append('OPTIONS --strand_specific: Single-stranded mode enabled')
        
        if self.namespace.output_dir:
            self.comp_directory_name = os.path.abspath(self.namespace.output_dir)
            if os.path.isdir(self.comp_directory_name) and os.listdir(self.comp_directory_name):
                self.output.append('ERROR: Output directory specified with -o needs to be an empty or non-existent directory')
                self.exit_now = True
        else:
            self.output.append('ERROR: Output directory needed. Use -o flag, which is mandatory.')
            self.exit_now = True

        return self.output