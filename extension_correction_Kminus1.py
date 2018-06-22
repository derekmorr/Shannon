import time
import sys
import re
import pdb,math
import numpy
from counter import Counter

BASES = ['A', 'G', 'C', 'T']
dont_correct_errors = False
rmer_to_contig = {}
contig_to_rmer = {}
cmer_to_contig = {}

c1 = Counter("Loading", 10**6)
c2 = Counter("Correction", 10**6)


def reverse_complement(bases):
    """Return the reverse complement of BASES. Assumes BASES is
    all uppercase.
    """
    replacements = [('A', 't'), ('T', 'a'), ('C', 'g'), ('G', 'c')]
    for ch1, ch2 in replacements:
        bases = re.sub(ch1, ch2, bases)
    return bases[::-1].upper()

def argmax(lst, key):
    """Returns the element x in LST that maximizes KEY(x).
    """
    best = lst[0]
    for x in lst[1:]:
        if key(x) > key(best):
            best = x
    return best

def load_kmers(infile, double_stranded):
    """Loads the list of K-mers and copycounts and determines K.
    Returns (kmers, K).
    """

    kmers = {}
    with open(infile) as f:
        for line in f:
            if len(line) == 0: continue
            c1.increment()
            kmer, weight = line.split()
            kmer = kmer.upper()
            weight = (float(weight))
            kmers[kmer] = kmers.get(kmer,0)+weight
            if double_stranded:
                kmers[reverse_complement(kmer)] = kmers.get(reverse_complement(kmer),0)+weight

    K = len(kmers.keys()[0])
    return kmers, K

def extend(start, extended, unextended, traversed, kmers, K):
    last = unextended(start)
    tot_weight = 0; tot_kmer=0
    extension = []
    
    while True:
        next_bases = [b for b in BASES if extended(last, b) in kmers and extended(last, b) not in traversed]
        if len(next_bases) == 0: return [extension,tot_weight,tot_kmer]
        
        next_base = argmax(next_bases, lambda b : kmers[extended(last, b)])
        extension.append(next_base); tot_weight += kmers[extended(last,next_base)]; tot_kmer +=1
        last = extended(last, next_base)
        traversed.add(last)
        last = unextended(last)
        c2.increment()

def extend_right(start, traversed, kmers, K):
    return extend(start, lambda last, b: last + b, lambda kmer: kmer[-(K - 1):],
        traversed, kmers, K)

def extend_left(start, traversed, kmers, K):
    return extend(start, lambda last, b: b + last, lambda kmer: kmer[:K - 1],
        traversed, kmers, K)

def duplicate_check(contig, r = 12, f = 0.5):
    # To add: if rmer in the contig multiple times, only increment the dup-contig once for each time its in dup-contig
    dup_count = {}
    max_till_now = 0; max_contig_index = -1
    for i in range(0, len(contig)-r+1):
        if contig[i:i+r] in rmer_to_contig:
            for dup in rmer_to_contig[contig[i:i+r]]:
                if dup not in dup_count:
                    dup_count[dup] = 1
                else:
                    dup_count[dup] += 1
                if dup_count[dup] >= max_till_now:
                    max_till_now = dup_count[dup]; max_contig_index = dup
    a = numpy.zeros(len(contig))
    for i in range(0, len(contig)-r+1):
        if contig[i:i+r] in rmer_to_contig:
            if max_contig_index in rmer_to_contig[contig[i:i+r]]:
                a[i:i+r] = 1

    if 1: #for dup in dup_count:
        if sum(a) > f * float(len(contig)):
            return True
        else:
            return False
           


def trim_polyA(contig):
    '''Trim polyA tails if atleast last minLen bases are polyA or first minLen bases are polyT'''
    minLen = 13; startPt = 0; endPt = 0
    maxMiss = 0; startMiss = 0; endMiss = 0; 
    for i in range(len(contig)):
        startPt +=1
        if contig[i] != 'A':
            startMiss+=1
        if startMiss > maxMiss:
            break
    
    totLen = len(contig)
    

    for j in range(len(contig)):
        endPt +=1
        if contig[totLen-1-j] != 'T':
            endMiss+=1
        if endMiss > maxMiss:
            break;
        
    if startPt < minLen:
        startPt = 0
    else:
        startPt -=2
 
    if endPt < minLen:
        endPt = 0
    else:
        endPt -=2
 
    return contig[startPt:totLen-endPt]


def duplicate_check_ends(contig,C):
        qLen = len(contig);
        first_cmer = contig[:C]; last_cmer = contig[-C:]
        if first_cmer in cmer_to_contig and last_cmer in cmer_to_contig:
                contig_dict = {}
                for c,p in cmer_to_contig[first_cmer]:
                        if c in contig_dict:
                                contig_dict[c][0] = p
                        else:
                                contig_dict[c] = [p,-1]  
                for c,p in cmer_to_contig[last_cmer]:
                        if c in contig_dict:
                                contig_dict[c][1] = p
                        else:
                                contig_dict[c] = [-1,p]
                for c in contig_dict:
                        if contig_dict[c][0] >= 0 and contig_dict[c][1]>=0:
                                diff = contig_dict[c][1] - contig_dict[c][0]
                                if abs(diff - (len(contig)-C)) < 3: #we are assuming that if two kmers match and if the length is the same, then the sequences must be the same
                                                #if hamming_dist(contig,)
                                                #print(contig_name + ' is already in '+ c)
                                                #print(contig)
                                                #print(contigs[c][contig_dict[c][0]: contig_dict[c][1]+r])
                                                return True
        return False

        
def run_correction(infile, outfile, min_weight, min_length,double_stranded, comp_directory_name, comp_size_threshold):
    f_log = open(comp_directory_name+"/before_sp_log.txt", 'w')
    #pdb.set_trace()
    print "{:s}: Starting Kmer error correction..".format(time.asctime())
    f_log.write("{:s}: Starting..".format(time.asctime()) + "\n")
    kmers, K = load_kmers(infile, double_stranded)
    print "{:s}: {:d} K-mers loaded.".format(time.asctime(), len(kmers))
    f_log.write("{:s}: {:d} K-mers loaded.".format(time.asctime(), len(kmers)) + "\n")
    heaviest = sorted(kmers.items(), key=lambda kv: kv[1])
    heaviest = [(k, w) for k, w in heaviest if w >= min_weight]
    traversed, allowed = set(), set()
    f1 = open(outfile+'_contig','w')
    contig_index = 0;
    contig_connections = {}
    components_numbers = {}
    contigs = ["buffer"]
    #pdb.set_trace()
    while len(heaviest) > 0:
        start_kmer, _ = heaviest.pop()
        if start_kmer in traversed: continue
        traversed.add(start_kmer)

        [right_extension,right_kmer_wt,right_kmer_no] = extend_right(start_kmer, traversed, kmers, K)
        [left_extension,left_kmer_wt,left_kmer_no] = extend_left(start_kmer, traversed, kmers, K)
        tot_wt = right_kmer_wt + left_kmer_wt + kmers[start_kmer];
        tot_kmer = right_kmer_no + left_kmer_no + 1;
        avg_wt = tot_wt/max(1,tot_kmer);
        contig = ''.join(reversed(left_extension)) + start_kmer + ''.join(right_extension)
        #contig = trim_polyA(contig)

        r = K-1
        duplicate_suspect = duplicate_check_ends(contig, r) #Check whether this contig is significantly represented in a earlier contig. 
        #pdb.set_trace()
        # The line below is the hyperbola error correction.
        if dont_correct_errors or (len(contig) >= min_length and len(contig)*math.pow(avg_wt,1/4.0) >= 2*min_length*math.pow(min_weight,1/4.0) and not duplicate_suspect):
            f1.write("{:s}\n".format(contig));  contig_index+=1
            contigs.append(contig)
            if contig_index not in contig_connections:
                contig_connections[contig_index] = {}
            # For adding new kmers
            for i in range(len(contig) - K +1):
                allowed.add(contig[i:i+K])
                
            # For graph partitioning
            C = K-1

            #pdb.set_trace()
            for i in range(len(contig) - C +1):
                if contig[i:i+C] in cmer_to_contig:
                    for contig2_index,pos in cmer_to_contig[contig[i:i+C]]:

                        #pdb.set_trace()
                        if contig2_index in contig_connections[contig_index] and contig2_index != contig_index:
                            contig_connections[contig_index][contig2_index] += 1
                        elif contig2_index != contig_index:
                            contig_connections[contig_index][contig2_index] = 1
                        if contig_index in contig_connections[contig2_index] and contig2_index != contig_index:
                            contig_connections[contig2_index][contig_index] += 1
                        elif contig2_index != contig_index:
                            contig_connections[contig2_index][contig_index] = 1
                else:
                    cmer_to_contig[contig[i:i+C]] = []
                cmer_to_contig[contig[i:i+C]].append([contig_index,i])     
                #pdb.set_trace()
                
            # For error correction
            # for i in range(len(contig) -r + 1):
            #     rmer = contig[i:i+r]
            #     if rmer in rmer_to_contig:
            #         rmer_to_contig[rmer].append([contig_index,i])
            #     else:
            #         rmer_to_contig[rmer]=[[contig_index,i]]
    f1.close()
    print "{:s}: {:d} K-mers remaining after error correction.".format(time.asctime(), len(allowed))
    f_log.write("{:s}: {:d} K-mers remaining after error correction.".format(time.asctime(), len(allowed)) + " \n")
    
    # Writes out kmers from all allowed contigs
    allowed_kmer_dict = {}
    with open(outfile, 'w') as f:
        for kmer in allowed:
            f.write("{:s}\t{:d}\n".format(kmer, int(kmers[kmer])))
            allowed_kmer_dict[kmer] = int(kmers[kmer])
    del kmers
    f_log.write("{:s}: {:d} K-mers written to file.".format(time.asctime(), len(allowed)) + " \n") 
    #pdb.set_trace()

    # Depth First Search to find components of contig graph.

    f_log.write(str(time.asctime()) + ": " +"Before dfs " + "\n")
    
    contig2component = {}
    component2contig = {}
    seen_before = {}
    # Some 
    for contig_index in contig_connections:
        if contig_index not in contig2component:
            component2contig[contig_index] = []
            stack1 = [contig_index]
            seen_before[contig_index] = True
            #pdb.set_trace()
            while len(stack1) != 0:
                curr = stack1.pop()
                contig2component[curr] = contig_index
                component2contig[contig_index].append(curr)
                for connected in contig_connections[curr]:
                    if connected not in seen_before:
                        stack1.append(connected)
                        seen_before[connected] = True

    f_log.write(str(time.asctime()) + ": " + "After dfs " + "\n")                       
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
    f_log.write(str(time.asctime()) + ": " + "After Edges Loaded " + "\n")
    
    #comp_size_threshold = 1000
    #comp_directory_name = 'S15_SE_contig_partition2' #'NM_metis_test'
    
    # Build Metis Graph file.
    #pdb.set_trace()
    new_comp_num = 1
    
    remaining_file_curr_size = 0
    remaining_file_num = 1
    single_contig_index = 0
    single_contigs = open(comp_directory_name+"/reconstructed_single_contigs.fasta", 'w')
    non_comp_contigs = open(comp_directory_name+"/remaining_contigs"+str(remaining_file_num)+".txt", 'w')
    if 1:
        for component in component2contig:
            if len(component2contig[component]) == 1:
                contig_ind = component2contig[component][0]
                contig = contigs[contig_ind]
                single_contigs.write('>Single_'+ str(single_contig_index)+'\n')
                single_contigs.write(contig+'\n')
                single_contig_index+=1
                continue

            #pdb.set_trace()
            if len(component2contig[component]) > comp_size_threshold:
                with open(comp_directory_name+"/component" + str(new_comp_num) + ".txt" , 'w') as f:
                    num_nodes = len(component2contig[component])
                    #pdb.set_trace()
                    num_edges = len(connections_drawn[component])
                    f.write(str(num_nodes) + "\t" + str(num_edges) + "\t" + "001" + "\n")
                    code = {}
                    i = 1
                    for contig in component2contig[component]:
                        code[contig] = i
                        i += 1
                    for contig in component2contig[component]:
                        #f.write(str(contig) + "\t")
                        for contig2 in contig_connections[contig]:
                            f.write(str(code[contig2]) + "\t" + str(contig_connections[contig][contig2]) + "\t")
                        f.write("\n")
                        
                with open(comp_directory_name+"/component" + str(new_comp_num) + "contigs" + ".txt" , 'w') as f:
                    for contig in component2contig[component]:
                        f.write(contigs[contig] + "\n")
                        
                new_comp_num += 1
            else:
                for contig_ind in component2contig[component]:
                    #pdb.set_trace()
                    contig = contigs[contig_ind]
                    non_comp_contigs.write(contig + "\n")
                    
                remaining_file_curr_size += len(component2contig[component])
                if remaining_file_curr_size > comp_size_threshold:
                    remaining_file_num += 1
                    non_comp_contigs.close()
                    non_comp_contigs = open(comp_directory_name+"/remaining_contigs"+str(remaining_file_num)+".txt", 'w')
                    remaining_file_curr_size = 0

                    '''for i in range(len(contig) - K + 1):
                        kmer = contig[i:i+K]
                        non_comp_kmers.write(kmer + "\t" + str(int(kmers[kmer])) + "\n")'''
    
    f_log.write(str(time.asctime()) + ": " + "Metis Input File Created " + "\n")
    f_log.close()
    return allowed_kmer_dict
    #pdb.set_trace()
            
                
                

def extension_correction(arguments):
    #pdb.set_trace()
    double_stranded = '-d' in arguments
    arguments = [a for a in arguments if len(a) > 0 and a[0] != '-']

    infile, outfile = arguments[:2]
    min_weight, min_length = int(arguments[2]), int(arguments[3])
    comp_directory_name, comp_size_threshold = arguments[4], int(arguments[5])
    #pdb.set_trace()
    allowed_kmer_dict = run_correction(infile, outfile, min_weight, min_length, double_stranded, comp_directory_name, comp_size_threshold)
    return allowed_kmer_dict

if __name__ == '__main__':
    if len(sys.argv) == 1:
        arguments = ['kmers.dict', 'allowed_kmers.dict', '1', '1', '-d']
    else:
        arguments = sys.argv[1:]

    extension_correction(arguments)
