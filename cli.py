import argparse

class ArgumentParser(object):

    


    @staticmethod
    def parse(args):
        """returns an argument parser."""
        parser = argparse.ArgumentParser()

        # parser.add_argument("--help", dest="help", required=False)
        parser.add_argument("--version", action="store_true", help="display program version.")
        parser.add_argument("--inDisk", action="store_true", help="use disk to transfer data between components.")
        parser.add_argument("--inMem", action="store_true", help="use memory to transfer data between components.")

        parser.add_argument("-p", type=int, help="number of parallel worker processes.")

        return parser.parse_args(args)