using CSV, DataFrames

function annotate_blastp(s, b_p, pi_cutoff=30)
    blastp_output_p = b_p * s * ".tsv"
    blastp_annoated_p = b_p * s * "_annotated.csv"

    # get headers
    headers = [:qseqid, :sseqid, :pident, :length, :mismatch,
        :gapopen, :qstart, :qend, :sstart, :send, :evalue, :bitscore]
    df = DataFrame(CSV.File(blastp_output_p, header=false))
    rename!(df, headers)

    # set p_ident cutoff
    df = filter(:pident => p -> (p >= pi_cutoff), df)

    # only retain alignment with max p ident for each query
    df = combine(g -> g[argmax(g.pident), :], groupby(df, :qseqid),)

    # extract swissprot ids from matches
    get(sseqid) = match(r"(.*)(?=\.)", sseqid).captures[1]
    df[!, :blastp_swissprot_id] = get.(df[:, :sseqid])

    CSV.write(blastp_annoated_p, df)
end

species = ["grimshawi", "simulans", "sechellia", "willistoni"]

b_p = "/local/data/mphilcompbio/2022/mw894/genomics/a2/blastp/blastp_out_"
map(s -> annotate_blastp(s, b_p,), species)

# 1]
b_p = "/local/data/mphilcompbio/2022/mw894/genomics/a2/blastp/blastp_out_all_"
map(s -> annotate_blastp(s, b_p,), species)

# 2]
b_p = "/local/data/mphilcompbio/2022/mw894/genomics/a2/blastp/blastp_out_all_uni_"
map(s -> annotate_blastp(s, b_p,), species)