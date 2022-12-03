using CSV, DataFrames, StatsBase

s = "simulans"

# get go ids from uniprot
p_uniprot = "/local/data/mphilcompbio/2022/mw894/genomics/a2/truth_data/uniprot_" * s * ".tsv"
df = DataFrame(CSV.File(p_uniprot, select=["Entry", "Pfam", "Gene Ontology IDs", "Reviewed"]))
rename!(df, ["Entry" => :gt_uni_id, "Gene Ontology IDs" => :gt_go_ids, "Pfam" => :gt_pfam_ids, "Reviewed" => :reviewed])

df.gt_go_ids = map(s -> ismissing(s) ? String[] : split(s, "; "), df.gt_go_ids)
df.gt_pfam_ids = map(s -> ismissing(s) ? String[] : split(rstrip(s, ';'), ";"), df.gt_pfam_ids)

p_pfam_to_goid = "/local/data/public/eadc2/Genomics_1/assignments/assignment_2/pfam_database/pfam2go.txt"
mapping_pfam_to_go_id = DataFrame(CSV.File(p_pfam_to_goid, delim=" > ", comment="!", header=[:pfam_id, :go_id]))
mapping_pfam_to_go_id.pfam_id = map(s -> match(r"(?<=Pfam:)(.*)(?= )", s).match, mapping_pfam_to_go_id.pfam_id)
mapping_pfam_to_go_id.go_id = map(s -> match(r"(?<= ; )(.*)", s).match, mapping_pfam_to_go_id.go_id)

function get_goids_from_pfamids(pfam_ids)
    if ismissing(pfam_ids)
        return missing
    end

    collect(Iterators.flatten(map(
        pfam_id -> filter(:pfam_id => pfam_id_mapping -> (pfam_id_mapping == pfam_id),
            mapping_pfam_to_go_id
        ).go_id,
        pfam_ids
    )))
end

df.hmmer_quni_goids = get_goids_from_pfamids.(df.gt_pfam_ids)

df.check = map(r -> intersect(r.hmmer_quni_goids, r.gt_go_ids),
    eachrow(df)
)

println(length(collect(Iterators.flatten(df.check))))