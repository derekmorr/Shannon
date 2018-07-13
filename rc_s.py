import sys

D = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
reverse_complement = lambda x: ''.join([D[B] for B in x][::-1])

def reverse_complement_serial(infile, outfile):
    with open(outfile, 'w') as of:
        for line in open(infile):
            if not line.strip():
                continue
            fields = line.strip().split()
            if fields[0][0] == '>':
                of.write(line.strip())
                of.write('\n')
            else:
                of.write(reverse_complement(fields[0]))
                of.write('\n')

def main():
    if len(sys.argv) == 1:
        arguments = ['asd', 'in_fasta', 'out_fasta', '-d']
    else:
        arguments = sys.argv[1:]

    infile, outfile = arguments[:2]
    reverse_complement_serial(infile, outfile)

if __name__ == '__main__':
    main()





