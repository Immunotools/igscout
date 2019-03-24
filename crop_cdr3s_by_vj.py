import os
import sys
from Bio import SeqIO

def SeqBelongsToRecords(seq, records):
    for r in records:
        if r.seq.find(seq) != -1:
            return r.seq.index(seq)
    return -1

def ReadFasta(fasta_fname):
    records = []
    for r in SeqIO.parse(fasta_fname, 'fasta'):
        r.seq = str(r.seq).upper()
        records.append(r)
    return records

cdr3_fasta = sys.argv[1]
v_fasta = sys.argv[2]
j_fasta = sys.argv[3]
output_fasta = sys.argv[4]

cdr3_records = ReadFasta(cdr3_fasta)
v_genes = ReadFasta(v_fasta)
j_genes = ReadFasta(j_fasta)

k_min = 9
left_trim = 10
right_trim = 15
fh = open(output_fasta, "w")
num_processed = 0
num_trimmed_cdr3s = 0
min_length = 11
for cdr3 in cdr3_records:
    seq = cdr3.seq
    start_pos = 0
    for i in range(k_min, len(seq)):
        prefix = seq[ : i]
        v_pos = SeqBelongsToRecords(prefix, v_genes)
        if v_pos > 280: # v_pos != -1:
#            print "Pos of V k-mer: " + str(v_pos) + ", CDR3 prefix: " + prefix
            start_pos = i
        else:
            break
    if start_pos == 0:
        start_pos = left_trim
    end_pos = len(cdr3)
    for i in range(start_pos, len(seq)):
        suffix = seq[len(seq) - i : ]
        j_pos = SeqBelongsToRecords(suffix, j_genes)
        if j_pos != -1:
#            print "Pos of J k-mer: " + str(j_pos) + ", CDR3 suffix: " + suffix
            end_pos = len(seq) - i
        else:
            break
    if end_pos == len(cdr3):
        end_pos = len(cdr3) - right_trim
#    print start_pos, end_pos, seq[start_pos : end_pos]
    trimmed_cdr3 = seq[start_pos : end_pos]
    if len(trimmed_cdr3) >= min_length:
        num_trimmed_cdr3s += 1
        fh.write(">" + cdr3.id + "|ORF:" + str(start_pos % 3) + "|ORIG_LEN:" + str(len(seq)) + "\n")
        fh.write(seq[start_pos : end_pos] + "\n")
    num_processed += 1
fh.close()
print str(num_trimmed_cdr3s) + " out of " + str(num_processed) + " trimmed CDR3s have length >= 11 nt"
