using FASTX

bp_augustus_out = "/local/data/mphilcompbio/2022/mw894/genomics/a2/input_prots/"
species = ["grimshawi", "simulans", "sechellia", "willistoni"]

function prep_augustus(p)
    reg_CDS_score = r"(?<=CDS)(\t)(\d+)(\t)(\d+)(\t)(\d\.\d\d)"
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
                    cds_scores = [parse(Float64, m.captures[end]) for m in collect(eachmatch(reg_CDS_score, raw_CDS))]

                    # if true positive save the seq
                    if any(cds_scores .> 0.8) #score high enough
                        protein_seq = match(reg_protein_seq, raw_protein_seq).match
                        id = "aug_seq_" * string(length(prots) + 1)
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
    s_bp_augustus_out = bp_augustus_out * s * "_augustus_v2.txt"
    s_p_fasta = bp_augustus_out * s * "_augustus_v2.fasta"

    prots = prep_augustus(s_bp_augustus_out)

    FASTX.FASTA.Writer(open(s_p_fasta, "w")) do writer
        for r in prots
            write(writer, r)
        end
    end
end

map(augustus_out_to_fasta, species)
