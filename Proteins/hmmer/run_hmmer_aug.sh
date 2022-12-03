#!/bin/bash

# simulans grimshawi sechellia willistoni
species="willistoni"

hmmerscan=/local/data/public/eadc2/Genomics_1/programs/hmmer-3.1b2-linux-intel-x86_64/src/hmmscan
pfam=/local/data/public/eadc2/Genomics_1/assignments/assignment_2/pfam_database/Pfam-A.hmm


blastp=/local/data/mphilcompbio/2022/mw894/genomics/ncbi-blast-2.13.0+/bin/blastp

in="/local/data/mphilcompbio/2022/mw894/genomics/a2/input_prots/$species\_augustus_v2.fasta"
out="/local/data/mphilcompbio/2022/mw894/genomics/a2/hmmer/hmmer_out_aug_$species.tsv"

#rm $out
$hmmerscan --tblout $out -E 0.001 $pfam $in 
