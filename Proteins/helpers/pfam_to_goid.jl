# PFAM -> GO ID
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