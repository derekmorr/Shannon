import os, pdb
import sys, time

def run_cmd(s1):
	print(s1)
	os.system(s1)

def find_representatives(file_in,file_out):
	'''Find representative transcripts from file_in and return them in file_out.'''
	blat_file = file_out+ '_temp_blat.psl'
	thresh_factor = 0.99 
	os.system('python parallel_blat.py '+ file_in + ' '+ file_in + ' ' + blat_file + ' 100')
	allTr ={}
	with open(blat_file,'r') as fp:
		i = 0
		for line in fp:
			if i <0:
				i+=1; continue;
			tokens = line.split();qName = tokens[9]; tName = tokens[13]; qLen = int(tokens[10]); tLen = int(tokens[14]); score = int(tokens[0]); 
			if qName == tName:
				continue
			if qLen < 200:
				allTr[qName] = -1
			if score > thresh_factor*qLen:
				if (tLen > qLen) or ((tLen == qLen) and (tName > qName)):
					allTr[qName] = -1
			i+=1
	
	with open(file_out,'w') as fp_out:
		write_now = False
		with open(file_in,'r') as fp_in:
			for line in fp_in:
				tokens=line.strip().split();
				if tokens[0][0]=='>':
					tr_name=tokens[0][1:]
					if allTr.get(tr_name,0) < 0:
						write_now = False
					else:
						write_now = True
				if write_now:
					fp_out.write(line)


def main():
    if len(sys.argv) == 1:
        arguments = ['asd', 'in_fasta', 'out_fasta', '-d']
    else:
        arguments = sys.argv
    ds = '-d' in arguments
    arguments = [a for a in arguments if len(a) > 0 and a[0] != '-'][1:]

    infile, outfile = arguments[:2]
    find_representatives(infile, outfile)


if __name__ == '__main__':
    main()
