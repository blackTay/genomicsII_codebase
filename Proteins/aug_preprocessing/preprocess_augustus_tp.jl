using JLD2, FASTX

# paths
bp_augustus_out = "/local/data/mphilcompbio/2022/mw894/genomics/a2/input_prots/"
bp_goodmatches = "/home/mw894/comp_bio/gi/genomicsII/Genes/output/augustus/"
species = ["grimshawi", "simulans", "sechellia", "willistoni"]

function load_true_cds_match_locs(p)
    # loads the CDS locations in the genome identified as significant
    good_matches = load_object(p)

    true_cds_match_locs = String[]
    map(m -> push!(true_cds_match_locs, join([
                m[2].name_of_species_contig,
                m[2].start_pos,
                m[2].stop_pos], "->")),
        good_matches)

    true_cds_match_locs
end

function prep_augustus(p, true_cds_match_locs)
    reg_CDS_pos = r"(?<=CDS\t)([0-9]+)(\t)([0-9]+)"
    reg_CDS_chrom = r"^([\S\-]+)"
    reg_protein_seq = r"(?<=protein sequence = \[)(.*?)(?=])"
    prots = FASTA.Record[]

    open(p) do f
        raw_CDS = ""
        raw_protein_seq = ""
        prev_CDS_line = true
        line_count = 0

        for l in readlines(f)
            # breaking out when all the transcripts are repeated
            line_count += 1
            (line_count < 14) && (continue)
            (occursin(">", l)) && (break)

            # if CDS line
            if !(startswith(l, '#'))
                # if previous was !CDS line => check CDS & prot 
                if !prev_CDS_line
                    cds_chrom_match = match(reg_CDS_chrom, raw_CDS).match
                    cds_matches = [cds_chrom_match * "->" * replace(m.match, "\t" => "->") for m in collect(
                        eachmatch(reg_CDS_pos, raw_CDS))]

                    # if true positive save the seq
                    if !isdisjoint(cds_matches, true_cds_match_locs)
                        id = string(length(prots) + 1) * "_" * join(cds_matches, ";")
                        protein_seq = match(reg_protein_seq, raw_protein_seq).match
                        push!(prots, FASTA.Record(id, protein_seq))
                    end

                    # reset
                    raw_protein_seq = ""
                    raw_CDS = ""
                end

                # else just collect CDS
                raw_CDS = raw_CDS * l
                prev_CDS_line = true
            else
                # collect the protein seq
                raw_protein_seq = raw_protein_seq * l[3:end]
                prev_CDS_line = false
            end
        end
    end

    prots
end

function augustus_out_to_fasta(s)
    s_bp_goodmatches = bp_goodmatches * s * "_augustus_to_groundtruth_good_matches.jld2"
    s_bp_augustus_out = bp_augustus_out * s * "_augustus_v2.txt"
    s_p_fasta = bp_augustus_out * s * "_augustus_v2.fasta"

    true_cds_match_locs = load_true_cds_match_locs(s_bp_goodmatches)

    prots = prep_augustus(s_bp_augustus_out, true_cds_match_locs)

    FASTX.FASTA.Writer(open(s_p_fasta, "w")) do writer
        for r in prots
            write(writer, r)
        end
    end
end

map(augustus_out_to_fasta, species)
