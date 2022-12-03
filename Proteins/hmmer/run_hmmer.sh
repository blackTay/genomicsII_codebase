#!/bin/bash

# simulans grimshawi sechellia willistoni
species="willistoni"

hmmerscan=/local/data/public/eadc2/Genomics_1/programs/hmmer-3.1b2-linux-intel-x86_64/src/hmmscan
pfam=/local/data/public/eadc2/Genomics_1/assignments/assignment_2/pfam_database/Pfam-A.hmm

in="/local/data/mphilcompbio/2022/mw894/genomics/a2/truth_data/uniprot_$species.fasta"
out="/local/data/mphilcompbio/2022/mw894/genomics/a2/hmmer/hmmer_out_$species.tsv"

$hmmerscan --tblout $out -E 0.001 $pfam $in 

#hmmer_domain_out=/home/mw894/comp_bio/genomics/a2/hmmer/hmmer_domain_output.txt
# --domtblout $hmmer_domain_out