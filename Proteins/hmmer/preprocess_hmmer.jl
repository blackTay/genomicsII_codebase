using DataFrames, CSV

function get_hmmer_output(s, b_p)
    hmmer_out_p = b_p * s * ".tsv"
    hmmer_annoated_p = b_p * s * "_annotated.csv"

    # reading in the hmmer output file, requires manual because of delim
    df = DataFrame(
        target_name=String[], accession=String[], query_name=String[],
        accession_bdom=String[], evalue=Float64[], score=Float64[], bias=Float64[],
        Evalue=Float64[], score_bdom=Float64[], bias_bdom=Float64[], exp=Float64[],
        reg=Float64[], clu=Integer[], ov=Integer[], env=Integer[], dom=Integer[],
        rep=Integer[], inc=Integer[], description_of_target=String[]
    )

    function parse_row(r)
        r = split.(r, r"\s+")
        r = push!(r[1:(19-1)], join(r[19:end], " "))
        r = convert(Array{Any,1}, r)
        r[5:12] = parse.(Float64, r[5:12])
        r[13:18] = parse.(Int64, r[13:18])
        push!(df, r)
    end

    # read lines, but first and 10 are info
    map(parse_row, readlines(hmmer_out_p)[4:(end-10)])

    # extract pfam id from target profile
    get_pfam_id(s) = (match(r"(.*)(?=\.)", s).captures[1])
    df[!, :hmmer_pfam_id] = map(get_pfam_id, df[!, :accession])

    CSV.write(hmmer_annoated_p, df)
end

species = ["grimshawi", "simulans", "sechellia", "willistoni"]

b_p = "/local/data/mphilcompbio/2022/mw894/genomics/a2/hmmer/hmmer_out_"
map(s -> get_hmmer_output(s, b_p), species)

b_p_aug = "/local/data/mphilcompbio/2022/mw894/genomics/a2/hmmer/hmmer_out_aug_"
map(s -> get_hmmer_output(s, b_p_aug), species)