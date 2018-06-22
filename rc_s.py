import sys,pdb,time, math, multiprocessing, copy
D={'A':'T','C':'G','G':'C','T':'A','N':'N'}
reverse_complement = lambda x: ''.join([D[B] for B in x][::-1])
nJobs = 20


def reverse_complement_serial(infile,outfile):
    reads = []
    if 1: #with open(outfile,'w') as of:
        for line in open(infile):
            if not line.strip(): continue
            fields = line.strip().split()
            if fields[0][0]=='>': 
                reads.append(line.strip())
            else:
                reads.append(reverse_complement(fields[0]))
    with open(outfile,'w') as of:
        of.write('\n'.join(reads) + '\n')


def main():
    if len(sys.argv) == 1:
        arguments = ['asd', 'in_fasta', 'out_fasta', '-d']
    else:
        arguments = sys.argv[1:]

    infile, outfile = arguments[:2]
    reverse_complement_serial(infile, outfile)

if __name__ == '__main__':
    main()





