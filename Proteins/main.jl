using CSV, DataFrames, JLD2

include("./helpers/swiss_to_goid.jl")
include("./helpers/pfam_to_goid.jl")

species = ["simulans", "sechellia", "grimshawi", "willistoni"]
s = "simulans"

function main(s)
    println("...processing $s")
    results = Dict{Symbol,Any}(:species => s)

    p_uniprot = "/local/data/mphilcompbio/2022/mw894/genomics/a2/truth_data/uniprot_" * s * ".tsv"
    df = DataFrame(CSV.File(p_uniprot, select=["Entry", "Pfam", "Gene Ontology IDs", "Reviewed"]))
    rename!(df, ["Entry" => :gt_uni_id, "Gene Ontology IDs" => :gt_go_ids, "Pfam" => :gt_pfam_ids, "Reviewed" => :reviewed])

    p_blastp_quni = "/local/data/mphilcompbio/2022/mw894/genomics/a2/blastp/blastp_out_all_uni_" * s * "_annotated.csv"
    df_blastp_quni = DataFrame(CSV.File(p_blastp_quni, select=["qseqid", "blastp_swissprot_id"]))
    rename!(df_blastp_quni, ["qseqid" => :blastp_quni_qid, "blastp_swissprot_id" => :blastp_quni_tid])

    p_hmmer_quni = "/local/data/mphilcompbio/2022/mw894/genomics/a2/hmmer/hmmer_out_" * s * "_annotated.csv"
    df_hmmer_quni = DataFrame(CSV.File(p_hmmer_quni, select=["query_name", "hmmer_pfam_id"]))
    rename!(df_hmmer_quni, ["query_name" => :hmmer_quni_qid, "hmmer_pfam_id" => :hmmer_quni_tpfamid])

    # tp
    df.gt_go_ids = map(s -> ismissing(s) ? String[] : split(s, "; "), df.gt_go_ids)
    df.gt_pfam_ids = map(s -> ismissing(s) ? String[] : split(rstrip(s, ';'), ";"), df.gt_pfam_ids)

    # 2]
    df_blastp_quni.blastp_quni_qid = map(s -> match(r"(?<=\|).*?(?=\|)", s).match, df_blastp_quni.blastp_quni_qid)
    df = leftjoin(df, df_blastp_quni, on=[:gt_uni_id => :blastp_quni_qid])
    df.blastp_quni_tp = (df.gt_uni_id .== df.blastp_quni_tid)
    df.blastp_quni_goids = get_goids_from_swissprotid.(df.blastp_quni_tid)

    # prot pred sucess
    df_swiss = filter(:reviewed => r -> (r == "reviewed"), df)
    prots_tps = count(==(true), df_swiss[!, :blastp_quni_tp])
    prots_fps = count(==(false), df_swiss[!, :blastp_quni_tp])
    results[:prots_recall] = prots_tps / size(df_swiss)[1]
    results[:prots_precision] = prots_tps / (prots_tps + prots_fps)

    # prot annotation success
    df.blastp_quni_goids_tps = map(
        r -> ismissing(r.blastp_quni_goids) ? missing : intersect(r.gt_go_ids, r.blastp_quni_goids),
        eachrow(df)
    )
    df.blastp_quni_goids_fps = map(
        r -> ismissing(r.blastp_quni_goids) ? missing : setdiff(r.blastp_quni_goids, r.gt_go_ids),
        eachrow(df)
    )

    prot_goids_tps = sum(length.(filter(inters -> !ismissing(inters), df.blastp_quni_goids_tps)))
    prot_goids_fps = sum(length.(filter(setdif -> !ismissing(setdif), df.blastp_quni_goids_fps)))
    results[:prot_goids_recall] = prot_goids_tps / sum(length.(df.gt_go_ids))
    results[:prot_goids_precision] = prot_goids_tps / (prot_goids_tps + prot_goids_fps)

    # 3]
    df_hmmer_quni.hmmer_quni_qid = map(s -> match(r"(?<=\|).*?(?=\|)", s).match, df_hmmer_quni.hmmer_quni_qid)
    df_hmmer_quni = combine(g -> [g[!, :hmmer_quni_tpfamid]], groupby(df_hmmer_quni, :hmmer_quni_qid))
    rename!(df_hmmer_quni, ["x1" => :hmmer_quni_tpfamids])
    df = leftjoin(df, df_hmmer_quni, on=[:gt_uni_id => :hmmer_quni_qid])
    df.hmmer_quni_goids = get_goids_from_pfamids.(df.hmmer_quni_tpfamids)

    # prot fam pred sucess
    df.hmmer_quni_tps = map(
        r -> ismissing(r.hmmer_quni_tpfamids) ? missing : intersect(r.gt_pfam_ids, r.hmmer_quni_tpfamids),
        eachrow(df)
    )
    df.hmmer_quni_fps = map(
        r -> ismissing(r.hmmer_quni_tpfamids) ? missing : setdiff(r.hmmer_quni_tpfamids, r.gt_pfam_ids),
        eachrow(df)
    )

    protfams_tps = sum(length.(filter(inters -> !ismissing(inters), df.hmmer_quni_tps)))
    protfams_fps = sum(length.(filter(setdif -> !ismissing(setdif), df.hmmer_quni_fps)))
    results[:protfams_recall] = protfams_tps / sum(length.(df.gt_pfam_ids))
    results[:protfams_precision] = protfams_tps / (protfams_tps + protfams_fps)

    # prot fam annotation success
    df.hmmer_quni_goids_tps = map(
        r -> ismissing(r.hmmer_quni_goids) ? missing : intersect(r.gt_go_ids, r.hmmer_quni_goids),
        eachrow(df)
    )
    df.hmmer_quni_goids_fps = map(
        r -> ismissing(r.hmmer_quni_goids) ? missing : setdiff(r.hmmer_quni_goids, r.gt_go_ids),
        eachrow(df)
    )

    protfams_goids_tps = sum(length.(filter(inters -> !ismissing(inters), df.hmmer_quni_goids_tps)))
    protfams_goids_fps = sum(length.(filter(setdif -> !ismissing(setdif), df.hmmer_quni_goids_fps)))
    results[:protfams_goids_recall] = protfams_goids_tps / sum(length.(df.gt_go_ids))
    results[:protfams_goids_precision] = protfams_goids_tps / (protfams_goids_tps + protfams_goids_fps)

    save("/local/data/mphilcompbio/2022/mw894/genomics/a2/" * s * "_df.jld2", Dict("df" => df))
    results
end

results_vec = map(s -> main(s), species)
results_df = DataFrame([NamedTuple{Tuple(keys(d))}(values(d)) for d in results_vec])
save("./results.jld2", Dict("results_df" => results_df))
