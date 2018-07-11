import os
import sys
import shutil
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import warnings
warnings.simplefilter("ignore")

import matplotlib as mplt
mplt.use('Agg')
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

min_k = 11

def ReadFasta(fasta_fname):
    records = []
    for r in SeqIO.parse(fasta_fname, 'fasta'):
        r.seq = str(r.seq).upper()
        records.append(r)
    return records

def GetDDict(d_records):
    d_dict = dict()
    for r in d_records:
        basename = r.id.split('*')[0]
        if basename not in d_dict:
            d_dict[basename] = []
        d_dict[basename].append(r.seq)
    return d_dict

def GetOriginalCDR3sLength(cdr3_header, cdr3_seq):
    prefix = 'ORIG_LEN:'
    if cdr3_header.find(prefix) == -1:
        return len(cdr3_seq)
    splits = cdr3_header.split('|')
    return int(splits[len(splits) - 1][len(prefix): ])

############################################
def FindLongestDMatches(cdr3, d_genes):
    d_dict = dict()
    for d in d_genes:
        forbidden_pos = set()
        for i in range(0, len(d.seq) - min_k + 1):
            if i in forbidden_pos:
                continue
            d_kmer = d.seq[i : i + min_k]
            if cdr3.find(d_kmer) == -1:
                continue
            d_pos = i + min_k
            cdr3_pos = cdr3.index(d_kmer) + min_k
            while d_pos < len(d.seq) and cdr3_pos < len(cdr3):
                if cdr3[cdr3_pos] != d.seq[d_pos]:
                    break
                d_kmer += cdr3[cdr3_pos]
                cdr3_pos += 1
                d_pos += 1
            if d_kmer not in d_dict:
                d_dict[d_kmer] = []
            d_dict[d_kmer].append(d.id)
            for j in range(i + 1, i + len(d_kmer)):
                forbidden_pos.add(j)
#    print d_dict
#    if cdr3 == 'TGACGTATTACTATGGTTCAGGGAGTTATGCGGGCAACC':
#        print d_dict
    return d_dict

def CreateLenDict(seq_dict):
    len_dict = dict()
    for d in seq_dict:
        if len(d) not in len_dict:
            len_dict[len(d)] = []
        len_dict[len(d)].append(d)
    return len_dict

def CreateNonRedundantDDict(d_dict):
    len_dict = CreateLenDict(d_dict)
    sorted_lens = sorted(len_dict.keys(), reverse = True)
    non_red_dict = dict()
    for l in sorted_lens:
        ds = len_dict[l]
        for d in ds:
            is_subseq = False
            for large_d in non_red_dict:
                if large_d.find(d) != -1:
                    is_subseq = True
                    break
            if not is_subseq:
                non_red_dict[d] = d_dict[d]
    return non_red_dict

########################################
class SimpleTandem:
    def __init__(self, cdr3 = "", d1_seq = "", d1_name = "", d2_seq = "", d2_name = "", insertion = ""):
        self.cdr3 = cdr3
        self.d1_seq = d1_seq
        self.d1_name = d1_name
        self.d2_seq = d2_seq
        self.d2_name = d2_name
        self.ins = insertion

    def Empty(self):
        return self.d1_seq == ''

    def DD(self):
#        return (self.d1_name, self.d2_name)
        return (self.d1_name.split('*')[0], self.d2_name.split('*')[0])

    def Seq(self):
        return self.d1_seq + self.ins + self.d2_seq

def DListIsAmbiguous(d_name_list):
    d_set = set([d.split('*')[0] for d in d_name_list])
    return len(d_set) == 1

def AnalyzeSimpleTandems(dd_dict, cdr3):
    if len(dd_dict) != 2:  
        return SimpleTandem()
    d_seqs = list(dd_dict.keys())
    d1 = d_seqs[0]  
    d2 = d_seqs[1]
    index1 = cdr3.index(d1)
    index2 = cdr3.index(d2)
    if len(dd_dict[d1]) != 1 or len(dd_dict[d2]) != 1: # D segments cannot inambiguously identifies
        return SimpleTandem()
    if index1 > index2:
        d1 = d_seqs[1]
        tmp = index1
        index1 = index2
        index2 = tmp
        d2 = d_seqs[0]
#    print cdr3
#    print index1, d1, index2, d2
    if index1 + len(d1) > index2: # sequences overlap
        return SimpleTandem()
    insertion = cdr3[index1 + len(d1) : index2]
    tandem = SimpleTandem(cdr3, d1, dd_dict[d1][0], d2, dd_dict[d2][0], insertion)
#    if tandem.DD()[0] == tandem.DD()[1]:
#        print tandem.DD(), d1, d2, cdr3
    return tandem

def AnalyzeTripleTandems(dd_dict, cdr3):
    return
    if len(dd_dict) != 3:
        return 
    d_seqs = list(dd_dict.keys())
    d1 = d_seqs[0]
    d2 = d_seqs[1]
    d3 = d_seqs[2]
    index1 = cdr3.index(d1)
    index2 = cdr3.index(d2)
    index3 = cdr3.index(d3)
    if len(dd_dict[d1]) != 1 or len(dd_dict[d2]) != 1 or len(dd_dict[d3]) != 1: # D segments cannot inambiguously identifies
        print "Triple tandem: ambiguous Ds"
        return 
    sorted_indices = dict()
    sorted_indices[index1] = d1
    sorted_indices[index2] = d2
    sorted_indices[index3] = d3
    sorted_list = sorted(sorted_indices)
    print sorted_list[0], len(sorted_indices[sorted_list[0]])
    print sorted_list[1], len(sorted_indices[sorted_list[1]])
    print sorted_list[2], len(sorted_indices[sorted_list[2]])
    max_overlap = 3
    for i in range(len(sorted_list) - 1):
        cur_index = sorted_list[i]
        next_index = sorted_list[i + 1]
        if cur_index + len(sorted_indices[cur_index]) - next_index > max_overlap:
            print "Triple tandem: overlap!"
            return 
    print "!!! Correct triple tandem was found"
    d1 = sorted_indices[sorted_list[0]]
    d2 = sorted_indices[sorted_list[1]]
    d3 = sorted_indices[sorted_list[2]]
    print dd_dict[d1][0], dd_dict[d2][0], dd_dict[d3][0]

############################################
def OutputDDMatrix(d_genes, all_ds, tandem_ds, output_fname, label):
    matrix = []
    annot_matrix = []
    for d in all_ds:
        matrix.append([0] * len(all_ds))
        annot_matrix.append([''] * len(all_ds))
    total_sign = 0
    for dd in tandem_ds:
        num_pairs = 0
        for p in tandem_ds[dd]:
            dist, d_gene = HammingDistanceOverAllDs(p.Seq(), d_genes)
            if dist <= 3:
                continue
            num_pairs += 1
        matrix[all_ds.index(dd[0])][all_ds.index(dd[1])] += num_pairs
    num_dd_pairs = 0
    for i in range(0, len(matrix)):
        for j in range(0, len(matrix[i])):
            if matrix[i][j] > 0:
                num_dd_pairs += 1
                annot_matrix[i][j] = str(matrix[i][j])
    print "# DD pairs: " + str(num_dd_pairs)
    sns.heatmap(matrix, cmap = 'jet', xticklabels = all_ds, yticklabels = all_ds, annot = np.array(annot_matrix), fmt = '', cbar = False, square = True, linewidth = .1, linecolor = 'grey', annot_kws = {'size' : 10})
    plt.yticks(rotation = 0, fontsize = 10)
    plt.ylabel('Start D gene', fontsize = 12)
    plt.xticks(rotation = 60, fontsize = 10)
    plt.xlabel('End D gene', fontsize = 12)    
    plt.title(label, fontsize = 14)
    pp = PdfPages(output_fname)
    pp.savefig()
    pp.close()
    plt.clf()    

def ComputeTandemSignificancy(tandem_ds):
    ins_dict = dict()
    for t in tandem_ds:
        if len(t.ins) not in ins_dict:
            ins_dict[len(t.ins)] = 0
        ins_dict[len(t.ins)] += 1
#    return len(ins_dict) 
    num_unique_ins = 0
    for ins_l in ins_dict:
        if ins_dict[ins_l] == 1:
            num_unique_ins += 1
    return num_unique_ins

def OutputDDMatrixInTxt(tandem_ds, d_genes, output_base):
    output_txt = output_base + ".txt"
    fh = open(output_txt, 'w')
    for dd in tandem_ds:
        cur_tandems = tandem_ds[dd]
        significancy = ComputeTandemSignificancy(tandem_ds[dd])
        num_good_cdr3s = 0
        good_ins = []
        for p in cur_tandems:
            dist, d_gene = HammingDistanceOverAllDs(p.Seq(), d_genes)
            if dist > 3:
                num_good_cdr3s += 1
                good_ins.append(p.ins)
        if num_good_cdr3s == 0:
            continue
        fh.write(dd[0] + "\t" + dd[1] + '\t' + str(len(tandem_ds[dd])) + "\t" + str(num_good_cdr3s) + '\t' + ','.join(good_ins) + "\n")
    fh.close()
#    os.system("Rscript visualize_circle_tandem_plot.r " + output_txt + " " + output_base + "_circle.pdf")

def OutputDDMatrixInDot(tandem_ds, output_base):
    output_dot = output_base + '.dot'
    fh = open(output_dot, 'w')
    fh.write('digraph {\n')
    for dd in tandem_ds:
        num_pairs = 0
        for p in tandem_ds[dd]:
            dist, d_gene = HammingDistanceOverAllDs(p.Seq(), d_genes)
            if dist <= 3:
                continue
            num_pairs += 1
        if num_pairs == 0:
            continue
        color = 'black'
        fh.write(dd[0] + ' -> ' + dd[1] + ' [color = ' + str(color) + "]\n")
    fh.write('}')
    fh.close()
    output_pdf = output_base + ".pdf"
    os.system('dot -Tpdf ' + output_dot + " -o " + output_pdf)

def OutputTandemCDR3s(tandem_ds, output_fname):
    fh = open(output_fname, 'w')
    fh.write('CDR3\tD1_name\tD1_seq\tD2_name\tD2_seq\n')
    for dd in tandem_ds:
        cur_tandems = tandem_ds[dd]
        for p in cur_tandems:
            dist, d_gene = HammingDistanceOverAllDs(p.Seq(), d_genes)
            if dist < 3:
                continue
            fh.write(p.cdr3 + '\t' + dd[0] + '\t' + p.d1_seq + '\t' + dd[1] + '\t' + p.d2_seq + '\n')
    fh.close()

############################################
def UpdateSeqCoverage(cov_list, seqs, subseq):
    index = -1
    for s in seqs:
        if s.find(subseq) != -1:
            index = s.index(subseq)
            break
    for i in range(0, len(subseq)):
        cov_list[index + i] += 1
    return cov_list

def OutputDDPairStats(d1_seqs, d2_seqs, tandem_list, output_fname):
    ins_lens = [len(t.ins) for t in tandem_list]
    d1_lens = [len(s) for s in d1_seqs]
    d2_lens = [len(s) for s in d2_seqs]
    d1_cov = [0] * max(d1_lens)
    d2_cov = [0] * max(d2_lens)
    for t in tandem_list:
        d1_cov = UpdateSeqCoverage(d1_cov, d1_seqs, t.d1_seq)
        d2_cov = UpdateSeqCoverage(d2_cov, d2_seqs, t.d2_seq)
    fig, (ax0, ax1, ax2) = plt.subplots(nrows = 3, figsize=(7, 10))
    # plot 1
    ax0.set_title('Start D coverage')
    ax0.bar(range(0, len(d1_cov)), d1_cov, color = 'green')
    # plot 2
    ax1.set_title('End D coverage')
    ax1.bar(range(0, len(d2_cov)), d2_cov, color = 'blue')    
    # plot 3
    ax2.set_title('Insertion lengths')
    ax2.hist(ins_lens, bins = 100)
#    plt.title(str(tandem_list) + " supporting CDR3s")
    pp = PdfPages(output_fname)
    pp.savefig() 
    pp.close()
    plt.clf()    

############################################
class SingleDMatch:
    def __init__(self, d_name = "", d_seq = ""):
        self.d_name = d_name
        self.d_seq = d_seq

    def D(self):
        return self.d_name.split('*')[0]

    def Seq(self):
        return self.d_seq

    def Empty(self):
        return self.d_name == '' or self.d_seq == ''

def AnalyzeSingleDMatch(d_dict):
    if len(d_dict) > 1:
        return SingleDMatch() # match is not single
    d_seq = ''
    for d in d_dict:
        d_seq = d
    if len(d_dict[d_seq]) != 1:
        return SingleDMatch() # match is ambiguous
    return SingleDMatch(d_dict[d_seq][0], d_seq)

def OutputSingleDUsage(single_d_usage, output_fname):
    sum = 0
    for d in single_d_usage:
        sum += len(single_d_usage[d])
    fh = open(output_fname, 'w')
    for d in single_d_usage:
        fh.write(d + '\t' + str(len(single_d_usage[d])) + '\t' + str(float(len(single_d_usage[d])) / float(sum) * 100) + '\n')
    fh.close()

def OutputSingleDUsageInPdf(all_ds, single_d_usage, num_cdr3s, output_fname, label):
#    all_ds = ['D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7', 'D8', 'D9', 'D10', 'D11', 'D12', 'D13', 'D14', 'D15', 'D16', 'D17', 'D18', 'D19', 'D20', 'D21', 'D22','D23', 'D24', 'D25']
    usage = []
    for d in all_ds:
        cur_usage = 0
        if d in single_d_usage:
            cur_usage = float(len(single_d_usage[d])) / float(num_cdr3s) * 100
        usage.append(cur_usage)
    plt.figure()
    plt.bar(range(0, len(all_ds)), usage)
    plt.xticks(range(0, len(all_ds)), all_ds, fontsize = 10, rotation = 60)
    plt.ylabel('Usage of D gene (%)', fontsize = 12)
    plt.title(label, fontsize = 14)
    pp = PdfPages(output_fname)
    pp.savefig()
    pp.close()
    plt.clf()

def OutputDCoverage(d_name, d_seqs, d_subseqs, output_fname):
    d_lens = [len(d) for d in d_seqs]
    d_cov = [0] * max(d_lens)
    for s in d_subseqs:
        d_cov = UpdateSeqCoverage(d_cov, d_seqs, s)
    x = range(0, len(d_cov))
    plt.bar(x, d_cov, color = 'green')
    plt.xticks(x, [str(i + 1) for i in x], fontsize = 12)
    plt.xlabel('D gene position', fontsize = 12)
    plt.ylabel('# CDR3s', fontsize = 12)
    plt.title(d_name, fontsize = 16)
    pp = PdfPages(output_fname)
    pp.savefig()
    pp.close()
    plt.clf()

############################################
def HammingDistance(seq1, seq2):
    dist = 0
    for i in range(0, len(seq1)):
        if seq1[i] != seq2[i]:
            dist += 1
    return dist

def ExtHammingDistance(short_seq, long_seq, k = 30):
    k = len(short_seq) - 5
    min_len = min(k, len(short_seq))
#    min_len = min(k, len(long_seq))
    if min_len > len(long_seq):
        return len(long_seq)
    min_dist = min_len
    for i in range(0, len(long_seq) - min_len + 1):
        long_substr = long_seq[i : i + min_len]
        for j in range(0, len(short_seq) - min_len + 1):
            short_substr = short_seq[j : j + min_len]
            cur_dist = HammingDistance(short_substr, long_substr)
            min_dist = min(min_dist, cur_dist)
    return min_dist

def HammingDistanceOverAllDs(superstr, d_genes, min_allowed_dist = 0):
    min_dist = len(superstr)
    closest_d = []
    for d in d_genes:
#        if len(superstr) > len(d.seq):
#            continue
        cur_dist = ExtHammingDistance(superstr, d.seq)
        if min_dist > cur_dist and cur_dist >= min_allowed_dist:
            closest_d = [d.id]
            min_dist = cur_dist
        elif min_dist == cur_dist:
            closest_d.append(d.id)
    return min_dist, closest_d

############################################
def GetAllDs(d_genes):
    all_ds = []
    for d in d_genes:
        d_base = d.id.split('*')[0]
        if len(all_ds) == 0 or all_ds[len(all_ds) - 1] != d_base:
            all_ds.append(d_base)
    return all_ds

############################################
def CDR3sIsGood(cdr3, d_genes, k):
    for d in d_genes:
        for i in range(0, len(d.seq) - k + 1):
            kmer = d.seq[i : i + k]
            if cdr3.find(kmer) != -1:
                return True
    return False

def GetNumGoodCDR3s(cdr3s, d_genes):
    num_good_cdr3s = 0
    for cdr3 in cdr3s:
        if CDR3sIsGood(cdr3.seq, d_genes, min_k):
            num_good_cdr3s += 1
    return num_good_cdr3s

############################################
def DDPairIsUpper(d1, d2):
    d1 = d1.split('*')[0]
    d2 = d2.split('*')[0]
    return int(d1[1 : ]) < int(d2[1 : ])

def ClassifyPairs(dd, start_ds, start_dist, end_ds, end_dist):
    pair_types = set()
    if start_dist == 1:
        for sd in start_ds:
            type = 'lower'
            if DDPairIsUpper(sd, dd[1]):
                type = 'upper'
            pair_types.add(type)
    if end_dist == 1:
        for ed in end_ds:
            type = 'lower'
            if DDPairIsUpper(dd[0], ed):
                type = 'upper'
            pair_types.add(type)
    return pair_types

############################################

d_fasta = sys.argv[1]
cdr3_fasta = sys.argv[2]
d_genes = ReadFasta(d_fasta)
d_dict = GetDDict(d_genes)
cdr3s = ReadFasta(cdr3_fasta)
output_dir = sys.argv[3]
if os.path.exists(output_dir):
    shutil.rmtree(output_dir)
os.mkdir(output_dir)
label = sys.argv[4]
all_counts = []

all_ds = GetAllDs(d_genes)
tandem_ds = dict()
num_tandem_cdr3s = 0
tandem_lens = []
tandem_counts = []

tandem_d_lens = dict()

single_ds = dict()
num_single_cdr3s = 0
traceable_lens = []

simple_long_cdr3s = []
min_cdr3_length = 70

non_traceable_lens = []

for cdr3 in cdr3s:
#    all_counts.append(int(cdr3.id.split('|')[1][len('COUNT:') : ]))
    d_seg_dict = FindLongestDMatches(cdr3.seq, d_genes)
    original_length = GetOriginalCDR3sLength(cdr3.id, cdr3.seq)
    if len(d_seg_dict) == 0:
        non_traceable_lens.append(original_length)
        continue
    non_red_dict = CreateNonRedundantDDict(d_seg_dict)
    if len(non_red_dict) == 1:
        single_d = AnalyzeSingleDMatch(non_red_dict)
        if single_d.Empty():
            non_traceable_lens.append(original_length)
            continue
        num_single_cdr3s += 1
        traceable_lens.append(original_length)   
        d = single_d.D()
        if d not in single_ds:
            single_ds[d] = []
        single_ds[d].append(single_d.Seq())
        if original_length >= min_cdr3_length:
            simple_long_cdr3s.append(cdr3.seq)
        continue
    tandem_d = AnalyzeSimpleTandems(non_red_dict, cdr3.seq)
    if tandem_d.Empty():
        non_traceable_lens.append(original_length)
        AnalyzeTripleTandems(non_red_dict, cdr3.seq)
        continue
    num_tandem_cdr3s += 1
    tandem_lens.append(original_length)
#    tandem_counts.append(int(cdr3.id.split('|')[1][len('COUNT:') : ]))
    dd = tandem_d.DD()
    if dd not in tandem_ds:
        tandem_ds[dd] = []
        tandem_d_lens[dd] = []
    tandem_ds[dd].append(tandem_d)
    tandem_d_lens[dd].append(original_length)

#print simple_long_cdr3s

single_cdr3_perc = float(len(traceable_lens)) / float(len(cdr3s)) * 100
print str(len(traceable_lens)) + ' (' + str(single_cdr3_perc) + '%) CDR3s with single known D genes were found'
print str(len(single_ds)) + ' D genes were used in single CDR3s'
print 'Average single CDR3 length: ' + str(np.mean(traceable_lens)) + ' nt'
  
tandem_cdr3_perc = float(num_tandem_cdr3s) / float(len(cdr3s)) * 100
print str(num_tandem_cdr3s) + " (" +  str(tandem_cdr3_perc) + "%) tandem CDR3s were found in " + cdr3_fasta
print str(len(tandem_ds)) + " pairs of D segments participate in tandem CDR3s"
print "Average tandem CDR3 length: " + str(np.mean(tandem_lens)) + " nt"

print "# trivial CDR3s: " + str(len(non_traceable_lens)) + " (" + str(float(len(non_traceable_lens)) / float(len(cdr3s)) * 100) + "%)"
print "Average length of trivial CDR3s: " + str(np.mean(non_traceable_lens))

print str(len(traceable_lens)) + ' (' + str(single_cdr3_perc) + '%),' + str(np.mean(traceable_lens)) + ',' + str(num_tandem_cdr3s) + " (" +  str(tandem_cdr3_perc) + "%)," + str(np.mean(tandem_lens)) + ',' + str(len(non_traceable_lens)) + " (" + str(float(len(non_traceable_lens)) / float(len(cdr3s)) * 100) + "%)," + str(np.mean(non_traceable_lens))

OutputSingleDUsage(single_ds, os.path.join(output_dir, 'single_d_usage.txt'))
OutputSingleDUsageInPdf(all_ds, single_ds, len(cdr3s), os.path.join(output_dir, 'single_d_usage.pdf'), label)
single_output_dir = os.path.join(output_dir, "single_usage")
os.mkdir(single_output_dir)
for d in single_ds:
    OutputDCoverage(d, d_dict[d], single_ds[d], os.path.join(single_output_dir, 'single_' + d + '.pdf'))

OutputDDMatrix(d_genes, all_ds, tandem_ds, os.path.join(output_dir, "dd_heatmap.pdf"), label)
OutputDDMatrixInTxt(tandem_ds, d_genes, os.path.join(output_dir, "dd_usage"))
OutputDDMatrixInDot(tandem_ds, os.path.join(output_dir, "dd_graph"))
OutputTandemCDR3s(tandem_ds, os.path.join(output_dir, "tandem_cdr3s.txt"))

d1_lens = []
d2_lens = []
distances = []
min_dist = 3
num_single_pairs = 0
num_low_distance_pairs = 0
num_good_dist_pairs = 0
num_upper_pairs = 0
num_all_tandems = 0
double_output_dir = os.path.join(output_dir, 'double_usage')
os.mkdir(double_output_dir)

good_tandems = []
for dd in tandem_ds:
    if len(tandem_ds[dd]) < 1:
        continue
    if len(tandem_ds[dd]) == 1:
        num_single_pairs += 1
#    print '====', dd[0], dd[1]
    insertions = []
    dds = tandem_ds[dd]
    for p in dds:
        dist, closest_d = HammingDistanceOverAllDs(p.Seq(), d_genes)
        if dist <= min_dist:
            num_low_distance_pairs += 1
        else:
            num_good_dist_pairs += 1
            good_tandems.append(p)
            if int(dd[0][1:]) < int(dd[1][1:]):
                num_upper_pairs += 1
#        print p.d1_seq, p.ins, p.d2_seq, dist
        distances.append(dist)
        insertions.append(p.ins)
        d1_lens.append(len(p.d1_seq))
        d2_lens.append(len(p.d2_seq))
        num_all_tandems += 1
#    print set(insertions)
#    OutputDDPairStats(d_dict[dd[0]], d_dict[dd[1]], tandem_ds[dd], os.path.join(double_output_dir, dd[0] + "_" + dd[1] + "_mult_" + str(len(tandem_ds[dd])) + ".pdf"))

print 'Total # tandem CDR3s: ' + str(num_good_dist_pairs) + ", # upper: " + str(num_upper_pairs) + ", # lower: " + str(num_good_dist_pairs - num_upper_pairs)
print '% of tandem CDR3s: ' + str(float(num_good_dist_pairs) / float(len(cdr3s)) * 100) + ", % wrt to good CDR3s: " + str(float(num_good_dist_pairs) / float(GetNumGoodCDR3s(cdr3s, d_genes)) * 100)
print "Average abundance of tandem CDR3s: " + str(np.mean(tandem_counts)) + " (all CDR3s: " + str(np.mean(all_counts)) + ")"
print 'Tandem bias: ' + str(float(num_good_dist_pairs - num_upper_pairs) / float(num_good_dist_pairs))

#print tandem_counts

#print min(d1_lens), np.mean(d1_lens)
#print min(d2_lens), np.mean(d2_lens)

plt.figure()
small_dist = [d for d in distances if d <= 3]
large_dist = [d for d in distances if d > 3]
#plt.hist([small_dist, large_dist], bins = 100, color = ['red', 'blue'])
plt.hist(distances, bins = 100)
plt.xticks(fontsize = 12)
plt.xlabel('$\Delta$-distance to closest D gene', fontsize = 14)
plt.ylabel('# spans', fontsize = 14)
plt.xticks(range(min(distances), max(distances) + 1), fontsize = 12)
plt.title(label, fontsize = 14)
pp = PdfPages(os.path.join(output_dir, "distance.pdf"))
pp.savefig()
pp.close()
plt.clf()

################## Analysis of mutation tandems ###############
sys.exit(1)
print "Analysis of " + str(len(good_tandems)) + ' good tandem CDR3s starts'
all_lower = 0
lower_ambigious = 0
all_upper = 0
upper_ambigious = 0
for cdr3 in good_tandems:
    start_dist, start_d = HammingDistanceOverAllDs(cdr3.d1_seq, d_genes, 1)
    end_dist, end_d = HammingDistanceOverAllDs(cdr3.d2_seq, d_genes, 1)    
    print str(cdr3.DD()) + ": " + str(start_d) + " (" + str(start_dist) + ") & " + str(end_d) + ' (' + str(end_dist) + ')'
    if start_dist != 1 and end_dist != 1:
        continue
    pair_types = ClassifyPairs(cdr3.DD(), start_d, start_dist, end_d, end_dist)
    if DDPairIsUpper(cdr3.DD()[0], cdr3.DD()[1]):
        all_upper += 1
        if 'lower' in pair_types:
            upper_ambigious += 1
    if not DDPairIsUpper(cdr3.DD()[0], cdr3.DD()[1]):
        all_lower += 1
        if 'upper' in pair_types:
            lower_ambigious += 1

print "Upper: " + str(upper_ambigious) + " / " + str(all_upper)
print "Lower: " + str(lower_ambigious) + ' / ' + str(all_lower)
