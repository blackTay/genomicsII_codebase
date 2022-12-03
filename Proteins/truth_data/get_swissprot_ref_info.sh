#!/bin/bash

# get pfam to go annotation as per the database
simulans_out=/local/data/mphilcompbio/2022/mw894/genomics/a2/truth_data/swissprot_simulans.tsv
grimshawi_out=/local/data/mphilcompbio/2022/mw894/genomics/a2/truth_data/swissprot_grimshawi.tsv
sechellia_out=/local/data/mphilcompbio/2022/mw894/genomics/a2/truth_data/swissprot_sechellia.tsv
willistoni_out=/local/data/mphilcompbio/2022/mw894/genomics/a2/truth_data/swissprot_willistoni.tsv

simulans_id=7240
grimshawi_id=7222
sechellia_id=7238
willistoni_id=7260

# simulans
wget -O $simulans_out   "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Clength%2Cgo_id%2Cxref_pfam&format=tsv&query=%28%28taxonomy_id%3A$simulans_id%29%20AND%20%28reviewed%3Atrue%29%29"
wget -O $grimshawi_out  "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Clength%2Cgo_id%2Cxref_pfam&format=tsv&query=%28%28taxonomy_id%3A$grimshawi_id%29%20AND%20%28reviewed%3Atrue%29%29"
wget -O $sechellia_out  "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Clength%2Cgo_id%2Cxref_pfam&format=tsv&query=%28%28taxonomy_id%3A$sechellia_id%29%20AND%20%28reviewed%3Atrue%29%29"
wget -O $willistoni_out "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Clength%2Cgo_id%2Cxref_pfam&format=tsv&query=%28%28taxonomy_id%3A$willistoni_id%29%20AND%20%28reviewed%3Atrue%29%29"
            