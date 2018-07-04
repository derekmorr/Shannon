import subprocess
import os
import math

def cut_file(in_name, out_name, line_start, line_end):
    os.system('awk \'NR > ' + str(line_end) + ' { exit } NR >= ' + str(line_start) \
        +  '\' '+ in_name + ' > ' + out_name)

def parallel_blat(target_fasta, query_fasta, out_file, QUERY_SPLIT):
    '''Function takes in target,query and output file. parallelizes blat by running GNU parallel
    - Currently only parallelizes on query space
    - Also assumes that query fasta file takes two lines per sequence (not wrapped)'''
    target_length = float(subprocess.check_output('grep -c \'>\' ' + target_fasta, shell=True))
    query_length = float(subprocess.check_output('grep -c \'>\' ' + query_fasta, shell=True))
    os.system('awk \'/^>/{print s? s"\\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}\' ' \
        + query_fasta + ' > ' + query_fasta +'_nospace')
    query_fasta = query_fasta + '_nospace'
    print('Query length: ' + str(query_length) \
        + ' Target length: ' + str(target_length) \
        + ' Query Split: ' + str(QUERY_SPLIT))
    split_size = int(math.floor(float(query_length)/QUERY_SPLIT))

    for n in range(QUERY_SPLIT):
        if n == QUERY_SPLIT-1:
            line_end = 2 * query_length
        else:
            line_end = 2*(n+1)*split_size
        cut_file(query_fasta, query_fasta + '_' + str(n+1), 2*(n)*split_size+1, line_end)

    q_range = range(QUERY_SPLIT)
    x = [int(i)+1 for i in q_range]
    q_str = " ".join(map(str, x))
    os.system('rm ' + out_file + '_*')
    blat_command = 'parallel blat -noHead ' + target_fasta + ' ' + query_fasta + '_{} ' \
        + out_file + '_{} ::: ' + q_str
    print(blat_command)
    os.system('time ' + blat_command)
    os.system('cat ' + out_file + '_* > ' + out_file)
    os.system('rm ' + out_file + '_*')
    os.system('rm ' + query_fasta + '_*')

def main():
    import sys
    if len(sys.argv) == 5:
        Query_split = int(sys.argv[4])
    else:
        Query_split = 100
    parallel_blat(sys.argv[1], sys.argv[2], sys.argv[3], Query_split)


if __name__ == '__main__':
    main()
