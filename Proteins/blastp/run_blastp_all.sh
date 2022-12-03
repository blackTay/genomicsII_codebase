#!/bin/bash

# Align to database
# simulans grimshawi sechellia willistoni
species="willistoni"

export BLASTDB=/local/data/mphilcompbio/2022/mw894/genomics/swissprot/

blastp=/local/data/mphilcompbio/2022/mw894/genomics/ncbi-blast-2.13.0+/bin/blastp

in="/local/data/mphilcompbio/2022/mw894/genomics/a2/truth_data/swissprot_$species.fasta"
out="/local/data/mphilcompbio/2022/mw894/genomics/a2/blastp/blastp_out_all_$species.tsv"

rm $out
$blastp -query $in -db swissprot -out $out -outfmt 6
