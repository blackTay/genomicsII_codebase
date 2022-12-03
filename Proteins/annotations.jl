using JLD2, DataFrames, InlineStrings, CSV

include("./helpers/pfam_to_goid.jl")

s = "simulans"

df = load("/local/data/mphilcompbio/2022/mw894/genomics/a2/" * s * "_df.jld2", "df")

# blastp
open("./goids/blastpquni_anno_$s.txt", "w") do f
    n = 0
    for goids in df.blastp_quni_goids
        if !ismissing(goids)
            n += 1
            write(f, ("N_$n" * "_protein\t") * join(goids, "\t") * "\n")
        end
    end
end

# hmmer
open("./goids/hmmerquni_anno_$s.txt", "w") do f
    n = 0
    for goids in df.hmmer_quni_goids
        if !ismissing(goids)
            n += 1
            write(f, ("N_$n" * "_proteinfamily\t") * join(goids, "\t") * "\n")
        end
    end
end

# hmmer aug
p_hmmer_qaug = "/local/data/mphilcompbio/2022/mw894/genomics/a2/hmmer/hmmer_out_aug_" * s * "_annotated.csv"
df_hmmer_qaug = DataFrame(CSV.File(p_hmmer_qaug, select=["query_name", "hmmer_pfam_id"]))
rename!(df_hmmer_qaug, ["query_name" => :hmmer_qaug_qid, "hmmer_pfam_id" => :hmmer_qaug_tpfamid])

df_hmmer_qaug = combine(g -> [g[!, :hmmer_qaug_tpfamid]], groupby(df_hmmer_qaug, :hmmer_qaug_qid))
rename!(df_hmmer_qaug, ["x1" => :hmmer_qaug_tpfamids])
df_hmmer_qaug.hmmer_qaug_goids = get_goids_from_pfamids.(df_hmmer_qaug.hmmer_qaug_tpfamids)

open("./goids/hmmerqaug_anno_$s.txt", "w") do f
    n = 0
    for goids in df_hmmer_qaug.hmmer_qaug_goids
        if !ismissing(goids)
            n += 1
            write(f, ("N_$n" * "_proteinfamily\t") * join(goids, "\t") * "\n")
        end
    end
end