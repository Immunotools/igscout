# igscout
IgScout is a tool for de novo inference of D and J genes using Rep-seq data; analysis of CDR3 diversity; finding VDDJ recombination. 

Please note that this is a very first version of the tool and its usability will be improved in the future. Any suggestions for improvements are welcomed!

Usage:

python igscout.py -i input_cdr3s.fasta -k kmer_size -o output_dir

python tandem_cdr3_finder.py d_genes.fasta input_cdr3s.fasta output_dir


Dependencies:

Biopython

Matplotlib

seaborn

NumPy
