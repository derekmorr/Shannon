import time

class Counter(object):
    def __init__(self, name, report_length, noun=''):
        self.name = name
        self.count = 0
        self.report_length = report_length

    def increment(self):
        self.count += 1
        if self.count % self.report_length == 0:
            print("{:s}: {:s}, processed {:d} kmers".format(time.asctime(), self.name, self.count))
