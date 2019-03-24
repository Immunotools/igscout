import os
import sys
import shutil

input_cdr3_fname = sys.argv[1]
output_dir = sys.argv[2]

if os.path.exists(output_dir):
    shutil.rmtree(output_dir)
os.mkdir(output_dir)

graph_fname = os.path.join(output_dir, 'compressed_cdr3s.fasta')
print "Constructing Hamming graph on CDR3s..."
os.system('/Poppy/ysafonova/ig_repertoire_constructor/build/release/bin/./ig_swgraph_construct -i ' + input_cdr3_fname + ' --tau 3 -T 0 -t 20 -k 5 -o ' + graph_fname)
consensus_fname = os.path.join(output_dir, 'consensus_cdr3s.fasta')
print "Computing consensus for connected components of Hamming graph..."
os.system('python compute_consensus_using_hg.py ' + input_cdr3_fname + ' ' + graph_fname + ' ' + consensus_fname)
trimmed_fname = os.path.join(output_dir, 'trimmed_cdr3s.fasta')
print "Trimming CDR3s by V suffixes and J prefixes..."
os.system('python crop_cdr3s_by_vj.py ' + consensus_fname + ' germline_data/human/IGHV.fa germline_data/human/IGHJ.fa ' + trimmed_fname) 
