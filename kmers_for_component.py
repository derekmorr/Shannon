import time
import pdb,math
import os
import os.path
import copy
import multiprocessing
from weight_updated_graph import weight_updated_graph

def run_cmd(s1):
    print(s1)
    os.system(s1)

class Counter():
    '''Used for printing the number of kmers processed'''
    def __init__(self, name, report_length):
        self.name = name
        self.count = 0
        self.report_length = report_length

    def increment(self):
        self.count += 1
        if self.count % self.report_length == 0:
            print "{:s}: {:s}, processed {:d}".format(time.asctime(), self.name, self.count)

D = {'A':'T','C':'G','G':'C','T':'A', 'N':'N'}
reverse_complement = lambda x: ''.join([D[B] for B in x[::-1]])
'''def reverse_complement(bases):
    replacements = [('A', 't'), ('T', 'a'), ('C', 'g'), ('G', 'c')]
    for ch1, ch2 in replacements:
        bases = re.sub(ch1, ch2, bases)
    return bases[::-1].upper()'''

def rc(lines,out_q):
        nl = copy.deepcopy(lines)
        for (i,line) in enumerate(lines):
                nl[i]=(reverse_complement(line.strip()))
        out_q.put(nl)

def rc_mate_ds(reads_1,reads_2,ds, out_q):
        nr1 = copy.deepcopy(reads_1)
        if ds: nr2 = copy.deepcopy(reads_2)
        for (i,read_1) in enumerate(reads_1):
                nr1[i]=[reads_1[i],reverse_complement(reads_2[i].strip())]
                if ds: nr2[i]=[reads_2[i],reverse_complement(reads_1[i].strip())]
        if ds: nr1.extend(nr2)
        out_q.put(nr1)

def par_read(reads_files,double_stranded, nJobs):
    reverse_complement = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A'}[B] for B in x][::-1])
    if len(reads_files)==1:
        with open(reads_files[0]) as f:
            lines = f.readlines()
        reads = lines[1::2]
        if double_stranded:
            chunk = min(1000000, len(reads))
            nProcs = int(math.ceil(len(reads)/float(chunk)))
            temp_q = multiprocessing.Queue()
            procs = [multiprocessing.Process(target=rc,args=(reads[x*chunk:(x+1)*chunk],temp_q)) for x in range(nProcs)]
            split_procs = []

            for i in range(int(math.ceil(float(nProcs)/nJobs))):
                split_procs.append(procs[(i)*nJobs:(i+1)*nJobs])


            for curr_procs in split_procs:
                for p in curr_procs:
                    p.start()
                for i in range(len(curr_procs)):
                    reads.extend(temp_q.get())
                for p in curr_procs:
                    p.join()
        print(str(len(reads)) + ' Reads loaded')
        return reads
    elif len(reads_files)==2:
        with open(reads_files[0]) as f:
            lines_1 = f.readlines()
        with open(reads_files[1]) as f:
            lines_2 = f.readlines()
        assert len(lines_1)==len(lines_2)
        reads_1 = lines_1[1::2]
        reads_2 = lines_2[1::2]
        if 1: #double_stranded:
            chunk = min(1000000, len(reads_1))
            nProcs = int(math.ceil(len(reads_1)/float(chunk)))

            temp_q = multiprocessing.Queue()
            procs = [multiprocessing.Process(target=rc_mate_ds,args=(reads_1[x*chunk:(x+1)*chunk],reads_2[x*chunk:(x+1)*chunk],double_stranded,temp_q)) for x in range(nProcs)]
            split_procs = []
            for i in range(int(math.ceil(float(nProcs)/nJobs))):
                split_procs.append(procs[(i)*nJobs:(i+1)*nJobs])
            reads = []
            for curr_procs in split_procs:
                for p in curr_procs:
                    p.start()
                for i in range(len(curr_procs)):
                    reads.extend(temp_q.get())
                for p in curr_procs:
                    p.join()
            print(str(len(reads)) + ' Reads loaded')
            r1 = [r[0] for r in reads]
            r2 = [r[1] for r in reads]
            reads = [r1,r2]    
        return reads
    else:
        return []

def par_SE_rc(reads,nJobs):
    chunk = int(math.ceil(len(reads)/float(nJobs)))
    temp_q = multiprocessing.Queue()
    procs = [multiprocessing.Process(target=rc,args=(reads[x*chunk:(x+1)*chunk],temp_q)) for x in range(nJobs)]
    for proc in procs:
        proc.start()
    for proc in procs:
        reads.extend(temp_q.get())
    for p in procs:
        p.join()
    return reads

def par_PE_rc(reads_1,reads_2,double_stranded,nJobs):
    chunk = int(math.ceil(len(reads_1)/float(nJobs)))
    temp_q = multiprocessing.Queue()
    procs = [multiprocessing.Process(target=rc_mate_ds,args=(reads_1[x*chunk:(x+1)*chunk],reads_2[x*chunk:(x+1)*chunk],double_stranded,temp_q)) for x in range(nJobs)]
    reads = []
    for proc in procs:
        proc.start()
    for proc in procs:
        reads.extend(temp_q.get())
    for p in procs:
        p.join()

    return reads


def kmers_for_component(k1mer_dictionary, kmer_directory, reads, reads_files, directory_name,
                        contig_file_extension, get_partition_k1mers, double_stranded=True,
                        paired_end=False, repartition=False, partition_size=500, overload=1.5,
                        K=24, gpmetis_path='gpmetis', penalty=5, only_reads=False, inMem=False, nJobs=1):
    """This fuction runs gpmetis on the components above a threshold size.  It then creates a dictionary called 
    k1mers2component {k1mer : component}.  It then sends any reads that share a kmer with a component to that component.
    It then cycles through the k1mer file and outputs the k1mers along with their weights to a file for each component.
    It then creates a kmers2component dictionary.  It then outputs a kmers file for each component. 
    Inputs: 
    k1mer_dictionary
    kmer_directory: used instead of k1mer_dictionary, currently unused
    reads_files: where reads are
    directory_name: where files should be stored
    contig_file_extension: the contig files are in this extension
    """
    if os.path.exists(directory_name+"/before_sp_log.txt"):
        f_log = open(directory_name+"/before_sp_log.txt", 'a')
    else:
        f_log = open(directory_name+"/before_sp_log.txt", 'w')
    
    def write_log(s):
        f_log.write(s + "\n")
        print(s)
    
    # Counts number of components above size threshold
    Num_Components = 0
    i = 1
    while os.path.exists(directory_name+"/component" + str(i) + "contigs.txt"):
        i += 1
        Num_Components += 1  

    # Counts number of components below size threshold
    Num_Remaining_Components = 0
    i = 1
    while os.path.exists(directory_name+"/remaining_contigs" + str(i) + ".txt"):
        i += 1
        Num_Remaining_Components += 1   
    
    # def find_comps(comp_dict, tries):
    #     #Only the exact 
    #     (comp,hits) = max(comp_dict.iteritems(), key=operator.itemgetter(1))
    #     if hits>=tries:
    #         return set(comp)
    #     else:
    #         return set()

    def get_rmers(read,R):
        i = 0
        k1mers = []
        while i < len(read) - R:
            k1mers.append(read[i:i+R])
            i+=R
        k1mers.append(read[-R:])
        return k1mers

    def get_comps(read,k1mers2component):
        k1mers = get_rmers(read,K+1)
        comps = [k1mers2component.get(k1mer, [-1,-1]) for k1mer in k1mers]
        comp_list = [set(comp[0]) for comp in comps if comp[0] != -1]
        if comp_list:
                assigned_comp = set.union(*comp_list)
        else:
                assigned_comp = set()
        return assigned_comp

    def get_comps_paired(read1,read2,k1mers2component):
        return set.union(get_comps(read1,k1mers2component),get_comps(read2,k1mers2component))

    if get_partition_k1mers:
        ufactor = int(1000.0*overload - 1000.0)
        components_broken = {}
        temp_string = ""
        
        # Runs gpmetis
        for i in range(Num_Components):
            with open(directory_name+"/component" + str(i+1) + contig_file_extension , 'r') as f:
                lines = f.readlines()
                num_contigs = len(lines)
                Partitions = min(int(math.ceil(float(num_contigs)/float(partition_size))),100)
                components_broken[i] = Partitions          
                temp_string += "Component " + str(i) + ": " + str(Partitions) + " partitions, "  
                if len(lines) >= 2:
                    run_cmd(gpmetis_path + " -ufactor="+str(ufactor) + " "+directory_name+"/component" + str(i+1) + ".txt" + " " +str(Partitions))
                    if repartition:
                        #Create GPMETIS file for repartition
                        partition_file = "/component" + str(i+1) + ".txt" + ".part." + str(components_broken[i])
                        og_graph_file = "/component" + str(i+1) + ".txt"
                        new_graph_file = "/component" + str(i+1) + "r2.txt"
                        contig_file = "/component" + str(i+1) + contig_file_extension
                        new_contig_file = contig_file #Currently unused option
                        randomize = False
                        write_log(str(time.asctime()) + ": " + "Creating graph for repartition ")
                        weight_updated_graph(directory_name, partition_file, og_graph_file, new_graph_file, contig_file, new_contig_file, penalty, randomize)
                        write_log(str(time.asctime()) + ": " + "Created graph for repartition ")
                        #Rerun GPMETIS file for repartition
                        run_cmd(gpmetis_path + " -ufactor="+str(ufactor) + " "+directory_name+"/component" + str(i+1)  + "r2.txt " +str(Partitions))
                    
        
        write_log(str(time.asctime()) + ": " + "gpmetis for partitioning is complete \n " + temp_string)
        
        new_components = {}
        k1mers2component = {}
        kmers2component = {}
                
        # Builds k1mer2component dictionary
        for i in components_broken:
            with open(directory_name+"/component" + str(i+1) + contig_file_extension, 'r') as f_contigs:
                contig_lines = f_contigs.readlines()
                with open(directory_name+"/component" + str(i+1)  + ".txt.part." + str(components_broken[i]) , 'r') as f_component:
                    j = 0
                    for line in f_component:
                        tokens = line.split()
                        comp = 'c' + str(i+1) + "_" + tokens[0]
                        contig = contig_lines[j].split()[0]
                        if comp not in new_components:
                            new_components[comp] = [contig]
                        else:
                            new_components[comp].append(contig)
                        for each in range(len(contig)-(K+1) + 1):
                            k1mer = contig[each:each+(K+1)]

                            if k1mer not in k1mers2component:
                                k1mers2component[k1mer] = [set([comp]), k1mer_dictionary.get(k1mer,0)]
                            else:
                                k1mers2component[k1mer][0].add(comp)
                        j += 1
                if repartition:
                    with open(directory_name+"/component" + str(i+1)  + "r2.txt.part." + str(components_broken[i]) , 'r') as f_component:
                        j = 0
                        for line in f_component:
                            tokens = line.split()
                            comp = 'r2_c' + str(i+1) + "_" + tokens[0]
                            contig = contig_lines[j].split()[0]
                            if comp not in new_components:
                                new_components[comp] = [contig]
                            else:
                                new_components[comp].append(contig)
                            for each in range(len(contig)-(K+1) + 1):
                                k1mer = contig[each:each+(K+1)]
                                if k1mer not in k1mers2component:
                                    k1mers2component[k1mer] = [set([comp]), k1mer_dictionary.get(k1mer,0)]
                                else:
                                    k1mers2component[k1mer][0].add(comp)
                            j += 1                        
                         
        # Adds remaining components to k1mers2component
        if 1:
            for i in range(Num_Remaining_Components):
                with open(directory_name+"/remaining_contigs"+str(i+1)+".txt", 'r') as remaining_contigs_file:
                    if 1:
                        lines = remaining_contigs_file.readlines()
                        comp = "cremaining"+str(i+1)
                        j = 0
                        for line in lines:
                            tokens = line.split()
                            contig = tokens[0]
                            if comp not in new_components:
                                new_components[comp] = [contig]
                            else:
                                new_components[comp].append(contig)
                            for each in range(len(contig)-(K+1) + 1):
                                k1mer = contig[each:each+(K+1)]
                                if k1mer not in k1mers2component:
                                    k1mers2component[k1mer] = [set([comp]), k1mer_dictionary.get(k1mer,0)]
                                else:
                                    k1mers2component[k1mer][0].add(comp)
                            j += 1                             

        write_log(str(time.asctime()) + ": " + "k1mers2component dictionary created ")
        
        '''no_kmers_in_comp = {}
        for comp in new_components:
            temp = 0
            for contig in new_components[comp]:
                temp += len(contig)
            no_kmers_in_comp[comp] = temp'''

        '''iter_tag = "_c"
        if second_iteration:
            iter_tag = "_r2_c"'''
        read_line =''   
        
        # Assigns reads to components in the non paired end case
        NR=10000000
        if not paired_end:
            read_ctr = 0
            offset = {}
            read_part_seq = {}   
            for comp in new_components:
                read_part_seq[comp] = []
                offset[comp] = 0
            f = open(reads_files[0],'r')
            while 1:
                reads = []
                while 1:
                    readname = f.readline()[:-1]
                    if not readname:
                        last_read = readname
                        break
                    read = f.readline()[:-1]
                    if read.strip('ACTG'):
                        continue
                    read_ctr +=1
                    reads.append(read)
                    last_read = read
                    if (read_ctr % NR) == (0) or (not read):
                        break
                
                if double_stranded: reads = par_SE_rc(reads,nJobs)
                for read in reads:
                        assigned_comp = get_comps(read,k1mers2component)
                        for each_comp in assigned_comp:
                            read_part_seq[each_comp].append(read)
                if not inMem:
                    for comp in new_components:
                        read_part_file = open(directory_name+"/reads"+str(comp)+".fasta", 'a')
                        read_part_file.write("".join(['>' + str(e+offset.get(comp,0)) + '\n' + read + '\n' for (e,read) in enumerate(read_part_seq[comp])]))
                        read_part_file.close()  
                        offset[comp] = offset.get(comp,0) + len(read_part_seq[comp])
                        read_part_seq[comp][:] = []
                if not last_read: break

        # Assigns reads to components in the paired end case
        elif paired_end:
            read_ctr = 0
            offset = {}
            read1_part_seq = {}
            read2_part_seq = {}
            for comp in new_components:
                read1_part_seq[comp] = []
                read2_part_seq[comp] = []
                offset[comp] = 0
            f1= open(reads_files[0],'r')
            f2 = open(reads_files[1],'r')
            while 1:
                reads_1 = []
                reads_2 = []
                while 1:
                    readname_1 = f1.readline()[:-1]
                    readname_2 = f2.readline()[:-1]
                    if not readname_1:
                        last_read = readname_1
                        break
                    read_1 = f1.readline()[:-1]
                    read_2 = f2.readline()[:-1]
                    if read_1.strip('ACTG') or read_2.strip('ACTG'):
                        continue
                    read_ctr +=1
                    reads_1.append(read_1)
                    reads_2.append(read_2)
                    last_read = read_1
                    if (read_ctr % NR) == (0) or (not read_1):
                        break
                if double_stranded:
                    reads = par_PE_rc(reads_1,reads_2,double_stranded,nJobs)
                if not double_stranded:
                    reads = [[reads_1[i],reads_2[i]] for i in range(len(reads_1))]
                del reads_1, reads_2

                for read in reads:
                    assigned_comp = get_comps_paired(read[0],read[1],k1mers2component)

                    for each_comp in assigned_comp:
                        read1_part_seq[each_comp].append(read[0])
                        read2_part_seq[each_comp].append(read[1])

                if not inMem: 
                    for comp in new_components:
                        read1_part_file = open(directory_name+"/reads"+str(comp)+"_1.fasta", 'a')
                        read2_part_file = open(directory_name+"/reads"+str(comp)+"_2.fasta", 'a')
                        read1_part_file.write("".join(['>' + str(e+offset.get(comp,0)) + '_1\n' + read + '\n' for (e,read) in enumerate(read1_part_seq[comp])]))
                        read2_part_file.write("".join(['>' + str(e+offset.get(comp,0)) + '_2\n' + read + '\n' for (e,read) in enumerate(read2_part_seq[comp])]))
                        read1_part_file.close()
                        read2_part_file.close()
                        offset[comp] = offset.get(comp,0) + len(read1_part_seq[comp])
                        read1_part_seq[comp][:] = []
                        read2_part_seq[comp][:] = []
                if not last_read:
                    break


            if not inMem: 
                for comp in new_components:
                    read1_part_file = open(directory_name+"/reads"+str(comp)+"_1.fasta", 'a')
                    read2_part_file = open(directory_name+"/reads"+str(comp)+"_2.fasta", 'a')
                    read1_part_file.write("".join(['>' + str(e+offset.get(comp,0)) + '_1\n' + read + '\n' for (e,read) in enumerate(read1_part_seq[comp])]))
                    read2_part_file.write("".join(['>' + str(e+offset.get(comp,0)) + '_2\n' + read + '\n' for (e,read) in enumerate(read2_part_seq[comp])]))
                    read1_part_file.close()
                    read2_part_file.close()
                    offset[comp] = offset.get(comp,0) + len(read1_part_seq[comp])
            elif inMem:
                for comp in new_components:
                    rps1 = [read1_part_seq[comp]]
                    rps2 = [read2_part_seq[comp]]
                    read1_part_seq[comp] = rps1
                    read2_part_seq[comp] = rps2
        
        write_log(str(time.asctime()) + ": " + "reads partititoned ")        

        contig_weights = {}
        if not only_reads: #If only_reads, no need to write k1mers
            write_log(str(time.asctime()) + ": Writing k1mers to file")
            # Writes out k1mers with weights for each partition
            for comp in new_components:
                with open(directory_name+"/component" + comp + "k1mers_allowed.dict" , 'w') as k1mer_file:
                    k1mer_file_data = []
                    contig_weights[comp] = []
                    for contig in new_components[comp]:
                        weight_list = []
                        for i in range(len(contig)-(K+1) + 1):
                            k1mer = contig[i:i+(K+1)]
                            k1mer_wt = k1mers2component[k1mer][1]
                            weight_list.append(k1mer_wt)
                            if not inMem: 
                                k1mer_file_data.append(k1mer + "\t" + str(k1mer_wt) + "\n")
                        if inMem:
                            contig_weights[comp].append(weight_list)

                    k1mer_file.writelines(k1mer_file_data)
            write_log(str(time.asctime()) + ": k1mers written to file ")

        write_log(str(time.asctime()) + ": kmers written to file " + "\n")
        
        if inMem:
            new_comps = new_components
        else:
            new_comps = [c for c in new_components]
            contig_weights = []
        rps = {}
        if inMem:
            if paired_end:
                for comp in new_components:
                    rps[comp] = [read1_part_seq[comp], read2_part_seq[comp]]
            else:
                for comp in new_components:
                    rps[comp] = [read_part_seq[comp]]
        return [components_broken, new_comps, contig_weights, rps]
