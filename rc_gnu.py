import os, subprocess, sys,pdb,time, math, multiprocessing, copy
import run_parallel_cmds
D={'A':'T','C':'G','G':'C','T':'A','N':'N'}
reverse_complement = lambda x: ''.join([D[B] for B in x][::-1])

def run_cmd(cmd):
    os.system(cmd)


def cut_file(in_name,out_name,line_start,line_end):
    os.system('awk \'NR > ' + str(line_end) + ' { exit } NR >= ' + str(line_start) +  '\' '+ in_name + ' > ' + out_name )


def find_L(readfile):
        N = 0; L = 0; even_line = False
        for line in open(readfile):
                if even_line: N+=1; L+=len(line)-1
                even_line = not even_line
        return (N,float(L)/N)


def rc_gnu(infile,tempfile,outfile,nCPU,python_path='python ',shannon_dir=''):
    chunks = nCPU
    if chunks ==1: 
        run_cmd(python_path + ' ' + shannon_dir + 'rc_s.py ' + infile + ' ' + outfile); 
        return find_L(infile)


    file_length = float(subprocess.check_output('grep -c \'>\' ' + infile,shell=True))
    split_size = int(math.ceil(float(file_length)/chunks))
    infile_piece = open(tempfile+'_1','w'); piece_no = 1;  curr_seqs = []
    read_tot = 0; no_reads = 0
    for line in open(infile):
	curr_seqs.append(line);
        fields = line.strip().split();
        if fields and fields[0][0]!='>': read_tot+=len(fields[0]); no_reads+=1
	if len(curr_seqs)==split_size*2:
                infile_piece = open(tempfile+'_'+str(piece_no),'w') 
		infile_piece.write(''.join(curr_seqs))
		infile_piece.close()
		piece_no +=1
	        curr_seqs = [];  

    if curr_seqs:
        infile_piece = open(tempfile+'_'+str(piece_no),'w') 
        infile_piece.write(''.join(curr_seqs))
        infile_piece.close()
    else:
        piece_no-=1


    N=no_reads; L= (read_tot)/ max(N,1)
    '''for n in range(chunks):
        if n==chunks-1:
            cut_file(infile,infile+'_'+str(n+1),2*(n)*split_size+1,2*file_length)
        else:
            cut_file(infile,infile+'_'+str(n+1),2*(n)*split_size+1,2*(n+1)*split_size)'''
    chunks = piece_no; c_range = range(chunks)
    x = [int(i)+1 for i in c_range]
    c_str = " ".join(map(str,x))
    cmds = []
    for i in range(chunks):
        cmds.append(python_path + ' ' + shannon_dir + 'rc_s.py ' + tempfile + '_' + str(i+1) + ' ' + outfile + '_' + str(i+1))
    run_parallel_cmds.run_cmds(cmds,chunks)

    #os.system('parallel --bibtex ' + python_path + ' ' + shannon_dir + 'rc_s.py ' + tempfile + '_{} ' + outfile + '_{} ' + ' ::: ' + c_str  )

    file_list = ' '.join([outfile+'_'+str(i+1) for i in range(chunks)])

    run_cmd('cat ' + file_list+' > ' + outfile)
    run_cmd('rm ' + outfile+'_* ')
    run_cmd('rm ' + tempfile+'_*  ')
    print(N); print(L)
    return (N,L)

def main():
    if len(sys.argv) == 1:
        arguments = ['asd', 'in_fasta', 'out_fasta', '-d']
    else:
        arguments = sys.argv[1:]

    infile, tempfile, outfile = arguments[:3]
    nCPU = int(arguments[3])
    rc_gnu(infile,tempfile,outfile,nCPU)

if __name__ == '__main__':
    main()


