#!/bin/bash

# Align to database

species="willistoni"

simulansid=7240
grimshawiid=7222
sechelliaid=7238
willistoniid=7260

export BLASTDB=/local/data/mphilcompbio/2022/mw894/genomics/swissprot/

blastp=/local/data/mphilcompbio/2022/mw894/genomics/ncbi-blast-2.13.0+/bin/blastp

in="/local/data/mphilcompbio/2022/mw894/genomics/a2/truth_data/swissprot_$species.fasta"
out="/local/data/mphilcompbio/2022/mw894/genomics/a2/blastp/blastp_out_$species.tsv"

rm $out
$blastp -query $in -db swissprot -out $out -outfmt 6 -taxids $willistoniid
