# SWISSPROT -> GO ID
p_swissprot_to_goid = "/local/data/public/eadc2/Genomics_1/assignments/assignment_2/swissprot_database/swissprot_database-2022.11.07.tsv"
mapping_swissprot_to_goid = DataFrame(CSV.File(p_swissprot_to_goid, select=["Entry", "Gene Ontology IDs"]))
rename!(mapping_swissprot_to_goid, ["Entry" => :swissprot_id, "Gene Ontology IDs" => :go_ids])
mapping_swissprot_to_goid.go_ids = map(
    s -> ismissing(s) ? String[] : split(s, "; "),
    mapping_swissprot_to_goid.go_ids
)

function get_goids_from_swissprotid(swissprot_id)
    if ismissing(swissprot_id)
        return missing
    elseif !(swissprot_id in mapping_swissprot_to_goid.swissprot_id)
        return String[]
    end

    filter(
        :swissprot_id => mapping_id -> mapping_id == swissprot_id,
        mapping_swissprot_to_goid
    ).go_ids[1]
end