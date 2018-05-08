import copy
import logging
import math
import multiprocessing
import numpy
import sys
import time
from operator import itemgetter
from collections import defaultdict

BASES = ['A', 'G', 'C', 'T']
correct_errors = True


class Counter(object):
    def __init__(self, name, report_length):
        self.name = name
        self.count = 0
        self.report_length = report_length

    def increment(self):
        self.count += 1
        if self.count % self.report_length == 0:
            print "{:s}: {:s}, processed {:d} kmers".format(time.asctime(), self.name, self.count)


c1 = Counter("Loading", 10**6)
c2 = Counter("Correction", 10**6)

reverse_complement = \
    lambda x: ''.join([{'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}[B] for B in x][::-1])


def rc(lines, out_q):
    nl = copy.deepcopy(lines)
    for (i, line) in enumerate(lines):
        nl[i] = (reverse_complement(line.strip())+'\n')
    out_q.put(nl)


def rc_mate_ds(reads_1, reads_2, double_stranded, out_q):
    nr1 = copy.deepcopy(reads_1)
    if double_stranded:
        nr2 = copy.deepcopy(reads_2)
    for (i, read_1) in enumerate(reads_1):
        nr1[i] = [reads_1[i], reverse_complement(reads_2[i].strip())+'\n']
        if double_stranded:
            nr2[i] = [reads_2[i], reverse_complement(reads_1[i].strip())+'\n']
    if double_stranded:
        nr1.extend(nr2)
    out_q.put(nr1)


def par_read_ns(reads_files, double_stranded, nJobs, ns):
    if len(reads_files) == 1:
        with open(reads_files[0]) as f:
            lines = f.readlines()
        reads = lines[1::2]
        if double_stranded:
            chunk = int(math.ceil(len(reads)/float(nJobs)))
            temp_q = multiprocessing.Queue()
            procs = [multiprocessing.Process(target=rc, args=(reads[x*chunk:(x+1)*chunk], temp_q)) for x in range(nJobs)]
            for p in procs:
                p.start()
            for i in range(nJobs):
                reads.extend(temp_q.get())
            for p in procs:
                p.join()
        ns.x = reads
    elif len(reads_files) == 2:
        with open(reads_files[0]) as f:
            lines_1 = f.readlines()
        with open(reads_files[1]) as f:
            lines_2 = f.readlines()
        assert len(lines_1) == len(lines_2)
        reads_1 = lines_1[1::2]
        reads_2 = lines_2[1::2]

        # double_stranded
        chunk = int(math.ceil(len(reads_1)/float(nJobs)))
        temp_q = multiprocessing.Queue()
        procs = [multiprocessing.Process(target=rc_mate_ds, args=(reads_1[x*chunk:(x+1)*chunk], reads_2[x*chunk:(x+1)*chunk], double_stranded, temp_q)) for x in range(nJobs)]
        for p in procs:
            p.start()
        reads = []
        for i in range(nJobs):
            reads.extend(temp_q.get())
        for p in procs:
            p.join()
        r1 = [r[0] for r in reads]
        r2 = [r[1] for r in reads]
        reads = [r1, r2]
        ns.x = reads
    else:
        ns.x = []


def par_read(reads_files, double_stranded, nJobs, out_q):
    if len(reads_files) == 1:
        with open(reads_files[0]) as f:
            lines = f.readlines()
        reads = lines[1::2]
        if double_stranded:
            chunk = int(math.ceil(len(reads)/float(nJobs)))
            temp_q = multiprocessing.Queue()
            procs = [multiprocessing.Process(target=rc, args=(reads[x*chunk:(x+1)*chunk], temp_q)) for x in range(nJobs)]
            for p in procs:
                p.start()
            for i in range(nJobs):
                reads.extend(temp_q.get())
            for p in procs:
                p.join()
        out_q.put(reads)
    elif len(reads_files) == 2:
        with open(reads_files[0]) as f:
            lines_1 = f.readlines()
        with open(reads_files[1]) as f:
            lines_2 = f.readlines()
        assert len(lines_1) == len(lines_2)
        reads_1 = lines_1[1::2]
        reads_2 = lines_2[1::2]

        # double_stranded
        chunk = int(math.ceil(len(reads_1)/float(nJobs)))
        temp_q = multiprocessing.Queue()
        procs = [multiprocessing.Process(target=rc_mate_ds, args=(reads_1[x*chunk:(x+1)*chunk], reads_2[x*chunk:(x+1)*chunk], double_stranded, temp_q)) for x in range(nJobs)]
        for p in procs:
            p.start()
        reads = []
        for i in range(nJobs):
            reads.extend(temp_q.get())
        for p in procs:
            p.join()
        r1 = [r[0] for r in reads]
        r2 = [r[1] for r in reads]
        reads = [r1, r2]
        out_q.put(reads)
    else:
        out_q.put([])


def lowComplexity(kmer):
    """Input: kmer k, Output: whether the kmer is within hamming distance of 2 from all A or all T"""
    nA = sum(c == 'A' for c in kmer)
    nT = sum(c == 'T' for c in kmer)
    nC = sum(c == 'C' for c in kmer)
    nG = sum(c == 'G' for c in kmer)
    return max(nA, nC, nG, nT) >= len(kmer) - 2


def argmax(lst, key):
    """Returns the element x in LST that maximizes KEY(x)."""
    best = lst[0]
    for x in lst[1:]:
        if key(x) > key(best):
            best = x
    return best


def par_load(lines, double_stranded, polyA_del, out_q):
    """Loads the list of K-mers and copycounts and determines K.
    Returns (kmers, K) in out_q.
    """
    reverse_complement = \
        lambda x: ''.join([{'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}[B] for B in x][::-1])
    d = defaultdict(lambda: 0)
    for (i, line) in enumerate(lines):
        kmer, weight = line.strip().split()
        if polyA_del and lowComplexity(kmer):
            continue
        d[kmer] += float(weight)
        if double_stranded:
            rc = reverse_complement(kmer)
            d[rc] += float(weight)
    out_q.put(d)


def load_kmers_parallel(infile, double_stranded, polyA_del=True, nJobs=1):
    kmers = defaultdict(lambda: 0)
    with open(infile) as f:
        lines = f.readlines()
    chunk = int(math.ceil(len(lines) / float(nJobs)))
    out_q = multiprocessing.Queue()
    procs = [multiprocessing.Process(target=par_load, args=(lines[x*chunk:(x+1)*chunk], double_stranded, polyA_del, out_q)) for x in range(nJobs)]
    for p in procs:
        p.start()
    K_set = False
    for i in range(nJobs):
        par_dict = out_q.get()
        for key in par_dict:
            if not K_set:
                K = len(key)
                K_set = True
            kmers[key] += par_dict[key]

    for p in procs:
        p.join()
    print(len(kmers))
    return kmers, K


def load_kmers(infile, double_stranded, polyA_del=True):
    """Loads the list of K-mers and copycounts and determines K.
    Returns (kmers, K).
    """

    kmers = defaultdict(lambda: 0)
    with open(infile) as f:
        for line in f:
            if not line:
                continue
            c1.increment()
            kmer, weight = line.split()
            kmer = kmer.upper()
            if polyA_del and lowComplexity(kmer):
                continue
            weight = float(weight)
            kmers[kmer] += weight
            if double_stranded:
                kmers[reverse_complement(kmer)] += weight
    K = len(kmers.keys()[0])
    return kmers, K


def extend(start, extended, unextended, traversed, kmers):
    last = unextended(start)
    tot_weight = 0
    tot_kmer = 0
    extension = []

    while True:
        next_bases = [b for b in BASES if extended(last, b) in kmers and extended(last, b) not in traversed]
        if not next_bases:
            return [extension, tot_weight, tot_kmer]

        next_base = argmax(next_bases, lambda b: kmers[extended(last, b)])
        extension.append(next_base)
        tot_weight += kmers[extended(last, next_base)]
        tot_kmer += 1
        last = extended(last, next_base)
        traversed.add(last)
        last = unextended(last)
        c2.increment()


def extend_right(start, traversed, kmers, K):
    return extend(start, lambda last, b: last + b, lambda kmer: kmer[-(K - 1):], traversed, kmers)


def extend_left(start, traversed, kmers, K):
    return extend(start, lambda last, b: b + last, lambda kmer: kmer[:K - 1], traversed, kmers)


def duplicate_check(contig, rmer_to_contig, r=15, f=0.5):
    # To add: if rmer in the contig multiple times,
    # only increment the dup-contig once for each time its in dup-contig
    dup_count = defaultdict(lambda: 1)
    max_till_now = 0
    max_contig_index = -1
    for i in range(0, len(contig)-r+1):
        if contig[i:i+r] in rmer_to_contig:
            for dup in rmer_to_contig[contig[i:i+r]]:
                dup_count[dup] += 1
                if dup_count[dup] >= max_till_now:
                    max_till_now = dup_count[dup]
                    max_contig_index = dup
    a = numpy.zeros(len(contig))
    for i in range(0, len(contig)-r+1):
        if contig[i:i+r] in rmer_to_contig:
            if max_contig_index in rmer_to_contig[contig[i:i+r]]:
                a[i:i+r] = 1

    return sum(a) > f * float(len(contig))


def duplicate_check2(contig, rmer_to_contig, r=15, f=0.5):
    # To add: if rmer in the contig multiple times,
    # only increment the dup-contig once for each time its in dup-contig
    dup_count = defaultdict(lambda: 1)
    max_till_now = 0
    max_contig_index = -1
    for i in range(0, len(contig)-r+1):
        # parallel map over this:
        contig_slice = contig[i:i+r]
        if contig_slice in rmer_to_contig:
            for dup in rmer_to_contig[contig_slice]:
                # return a struct with:
                # - dup_count[dup]'s new increment (0 or 1)
                # - max_till_now
                # - max_contig_index
                # then in the driver, after the processes have joined,
                # merge dup_count[dup] values
                # # how to process max_till_now and max_config_index ?
                dup_count[dup] += 1
                if dup_count[dup] >= max_till_now:
                    max_till_now = dup_count[dup]
                    max_contig_index = dup
    a = 0
    for i in range(0, len(contig)-r+1):
        if contig_slice in rmer_to_contig:
            if max_contig_index in rmer_to_contig[contig_slice]:
                a += r

    return a > f * float(len(contig))


def trim_polyA(contig):
    """Trim polyA tails if at least last minLen bases are polyA or first minLen bases are polyT"""
    startLen = 0
    endLen = 0
    startPt = 0

    maxMiss = 1
    startMiss = 0
    endMiss = 0
    endPt = 0
    for i in range(len(contig)):
        if contig[i] != 'T':
            startMiss += 1
        else:
            startPt = startLen + 1
        if startMiss > maxMiss:
            break
        startLen += 1

    totLen = len(contig)

    for j in range(len(contig)):
        if contig[totLen-1-j] != 'A':
            endMiss += 1
        else:
            endPt = endLen + 1
        if endMiss > maxMiss:
            break
        endLen += 1

    return contig[startPt:totLen-endPt]


def run_correction(infile, outfile, min_weight, min_length, double_stranded, comp_directory_name,
                   comp_size_threshold, polyA_del, inMem, nJobs, reads_files):

    # setup logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    log_file_name = comp_directory_name + "/before_sp_log.txt"
    file_handler = logging.FileHandler(log_file_name)
    file_handler.setLevel(logging.INFO)

    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    stdout_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(stdout_handler)

    logger.info('nJobs: %d', nJobs)
    logger.info('reads_files: %s', ' '.join(reads_files))
    logger.info("Starting Kmer error correction.")

    if nJobs == 1:
        kmers, K = load_kmers(infile, double_stranded, polyA_del)
    elif nJobs > 1:
        kmers, K = load_kmers_parallel(infile, double_stranded, polyA_del, nJobs)

    logger.info("%d K-mers loaded.", len(kmers))

    logger.info("Reads loading in background process.")

    heaviest = sorted(kmers.items(), key=itemgetter(1))
    traversed = set()
    allowed = set()
    f1 = open(outfile+'_contig', 'w')
    contig_index = 0
    contig_connections = {}
    contigs = ["buffer"]
    rmer_to_contig = {}
    cmer_to_contig = {}
    while len(heaviest) > 0:
        start_kmer, w = heaviest.pop()
        if w < min_weight:
            break  # min_weight is used here
        if start_kmer in traversed:
            continue
        traversed.add(start_kmer)

        [right_extension, right_kmer_wt, right_kmer_no] = \
            extend_right(start_kmer, traversed, kmers, K)
        [left_extension, left_kmer_wt, left_kmer_no] = \
            extend_left(start_kmer, traversed, kmers, K)
        tot_wt = right_kmer_wt + left_kmer_wt + kmers[start_kmer]
        tot_kmer = right_kmer_no + left_kmer_no + 1
        avg_wt = tot_wt / max(1, tot_kmer)
        contig = ''.join(reversed(left_extension)) + start_kmer + ''.join(right_extension)

        r = 15
        # Check whether this contig is significantly represented in a earlier contig.
        duplicate_suspect = duplicate_check(contig, rmer_to_contig, r)

        # The line below is the hyperbola error correction.
        if (not correct_errors) or \
                (len(contig) >= min_length and
                 len(contig)*math.pow(avg_wt, 1/4.0) >= 2*min_length*math.pow(min_weight, 1/4.0)
                 and not duplicate_suspect):
            f1.write("{:s}\n".format(contig))
            contig_index += 1
            contigs.append(contig)
            if contig_index not in contig_connections:
                contig_connections[contig_index] = {}
            # For adding new kmers
            for i in range(len(contig) - K + 1):
                allowed.add(contig[i:i+K])

            # For graph partitioning
            C = K-1

            for i in range(len(contig) - C + 1):
                if contig[i:i+C] in cmer_to_contig:
                    for contig2_index in cmer_to_contig[contig[i:i+C]]:
                        if contig2_index in contig_connections[contig_index] and \
                                contig2_index != contig_index:
                            contig_connections[contig_index][contig2_index] += 1
                        elif contig2_index != contig_index:
                            contig_connections[contig_index][contig2_index] = 1
                        if contig_index in contig_connections[contig2_index] and \
                                contig2_index != contig_index:
                            contig_connections[contig2_index][contig_index] += 1
                        elif contig2_index != contig_index:
                            contig_connections[contig2_index][contig_index] = 1
                else:
                    cmer_to_contig[contig[i:i+C]] = []
                cmer_to_contig[contig[i:i+C]].append(contig_index)

            # For error correction
            for i in range(len(contig) - r + 1):
                if contig[i:i+r] in rmer_to_contig:
                    rmer_to_contig[contig[i:i+r]].append(contig_index)
                else:
                    rmer_to_contig[contig[i:i+r]] = [contig_index]
    f1.close()
    logger.info("%d K-mers remaining after error correction.", len(allowed))

    # Writes out kmers from all allowed contigs
    allowed_kmer_dict = {}
    with open(outfile, 'w') as f:
        for kmer in allowed:
            if not inMem:
                f.write("{:s}\t{:d}\n".format(kmer, int(kmers[kmer])))
            allowed_kmer_dict[kmer] = int(kmers[kmer])
    del kmers
    logger.info("%d K-mers written to file.", len(allowed))

    # Depth First Search to find components of contig graph.

    logger.debug("Before dfs")

    contig2component = {}
    component2contig = {}
    seen_before = {}
    # Some
    for contig_index in contig_connections:
        if contig_index not in contig2component:
            component2contig[contig_index] = []
            stack1 = [contig_index]
            seen_before[contig_index] = True
            while stack1:
                curr = stack1.pop()
                contig2component[curr] = contig_index
                component2contig[contig_index].append(curr)
                for connected in contig_connections[curr]:
                    if connected not in seen_before:
                        stack1.append(connected)
                        seen_before[connected] = True

    logger.debug("After dfs")
    # Finds all connections for Metis graph file.

    connections_drawn = {}
    for component in component2contig:
        connections_drawn[component] = {}
    for contig in contig_connections:
        for contig2 in contig_connections[contig]:
            key1 = [contig, contig2]
            key1.sort()
            key1 = tuple(key1)
            component = contig2component[contig]
            if key1 not in connections_drawn[component]:
                # write to file here
                connections_drawn[component][key1] = True
    logger.debug("After Edges Loaded")

    # Build Metis Graph file.
    new_comp_num = 1

    remaining_file_curr_size = 0
    remaining_file_num = 1
    single_contig_index = 0
    single_contigs = \
        open(comp_directory_name+"/reconstructed_single_contigs.fasta", 'w')
    non_comp_contigs = \
        open(comp_directory_name+"/remaining_contigs" + str(remaining_file_num) + ".txt", 'w')
    for component in component2contig:
        if len(component2contig[component]) == 1:
            contig_ind = component2contig[component][0]
            contig = contigs[contig_ind]
            single_contigs.write('>Single_' + str(single_contig_index)+'\n')
            single_contigs.write(contig + '\n')
            single_contig_index += 1
            continue

        if len(component2contig[component]) > comp_size_threshold:
            with open(comp_directory_name+"/component" + str(new_comp_num) + ".txt", 'w') as f:
                num_nodes = len(component2contig[component])
                num_edges = len(connections_drawn[component])
                f.write(str(num_nodes) + "\t" + str(num_edges) + "\t" + "001" + "\n")
                code = {}
                i = 1
                for contig in component2contig[component]:
                    code[contig] = i
                    i += 1
                for contig in component2contig[component]:
                    for contig2 in contig_connections[contig]:
                        f.write(str(code[contig2]) + "\t" +
                                str(contig_connections[contig][contig2]) + "\t")
                    f.write("\n")

            with open(comp_directory_name+"/component" +
                      str(new_comp_num) + "contigs" + ".txt", 'w') as f:
                for contig in component2contig[component]:
                    f.write(contigs[contig] + "\n")

            new_comp_num += 1
        else:
            for contig_ind in component2contig[component]:
                contig = contigs[contig_ind]
                non_comp_contigs.write(contig + "\n")

            remaining_file_curr_size += len(component2contig[component])
            if remaining_file_curr_size > comp_size_threshold:
                remaining_file_num += 1
                non_comp_contigs.close()
                non_comp_contigs = open(comp_directory_name + "/remaining_contigs" +
                                        str(remaining_file_num) + ".txt", 'w')
                remaining_file_curr_size = 0

    logger.info("Metis Input File Created")
    logger.info("Read-loader in background process joining back.")

    reads = []
    logger.info("%d Reads loaded in background process.", len(reads))
    return allowed_kmer_dict, reads


def extension_correction(arguments, inMem=False):
    double_stranded = '-d' in arguments
    arguments = [a for a in arguments if len(a) > 0 and a[0] != '-']

    infile, outfile = arguments[:2]
    min_weight = int(arguments[2])
    min_length = int(arguments[3])
    comp_directory_name = arguments[4]
    comp_size_threshold = int(arguments[5])
    if len(arguments) > 6:
        nJobs = int(arguments[6])
    else:
        nJobs = 1
    if len(arguments) > 7:
        reads_files = [arguments[7]]
        if len(arguments) > 8:
            reads_files.append(arguments[8])
    else:
        reads_files = []

    allowed_kmer_dict, reads = run_correction(infile, outfile, min_weight, min_length,
                                              double_stranded, comp_directory_name,
                                              comp_size_threshold, True, inMem, nJobs, reads_files)
    return allowed_kmer_dict, reads


if __name__ == '__main__':
    if len(sys.argv) == 1:
        arguments = ['kmers.dict', 'allowed_kmers.dict', '1', '1', '-d']
    else:
        arguments = sys.argv[1:]

    extension_correction(arguments)
