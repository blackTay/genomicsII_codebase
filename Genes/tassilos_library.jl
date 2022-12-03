# parse file "simulans_prot_excerpt.txt" and return a dictionary of the form
# {gene_id: [hit_id, hit_id, ...]}

# define own type for a match containing 
# - the name where the matching contig is from in the species (eg species simulans, etc.)
# - the score of the match
# - start_position 
# - end_position

using Plots

struct Match
    name_of_species_contig::String
    score_of_match::Float64
    start_pos::Int64
    stop_pos::Int64
end

struct GroundTruth
    chromosome_id::String
    start_pos::Int64
    stop_pos::Int64
end

# overwrite the == operator for GroundTruth
Base.:(==)(a::GroundTruth, b::GroundTruth) = a.chromosome_id == b.chromosome_id && a.start_pos == b.start_pos && a.stop_pos == b.stop_pos


    

function parse_blast(filename)
    file = open(filename)
    dict = Dict{String, Vector{Match}}()
    
    currently_in_a_block = false
    current_protein_id_mel = ""

    # the start/ stop position in our species (simulans, ...)
    current_nuc_seq_start_species = -1
    current_nuc_seq_stop_species = -1
    current_name_of_species_contig = ""
    current_score = -1
    update_end_pos=false

    for line in eachline(file)
        # check if line contains "Query= "
        if occursin( "***** No hits found *****",line)
            # reset variables
            current_nuc_seq_start_species = -1
            current_nuc_seq_stop_species = -1
            current_name_of_species_contig = ""
            current_score = -1
        end

        if occursin("Query= ", line)
            # extract gene_id, keep old to store
            old_protein_id_mel = current_protein_id_mel
            current_protein_id_mel = split(line, " ")[2]

            # currently_in_a_block should be false
            if currently_in_a_block
                # safe the old block
                if haskey(dict, old_protein_id_mel)
                    push!(dict[old_protein_id_mel], Match(current_name_of_species_contig, current_score, current_nuc_seq_start_species, current_nuc_seq_stop_species))
                else
                    dict[old_protein_id_mel] = [Match(current_name_of_species_contig, current_score, current_nuc_seq_start_species, current_nuc_seq_stop_species)]
                end

                # reset variables
                current_nuc_seq_start_species = -1
                current_nuc_seq_stop_species = -1
                current_name_of_species_contig = ""
                current_score = -1
            end
            currently_in_a_block = true
            continue
        end


        if currently_in_a_block
            # get the name of species
            if current_name_of_species_contig == "" && occursin(">", line)
                current_name_of_species_contig = line[2:end]
                continue
            end

            # get the current_score (E value)
            if  startswith(line, " Score = ")
                if current_score == -1
                    # extract score
                    current_score = parse(Float64, replace(split(line, " ")[10],","=>""))
                    update_end_pos = true
                else
                    update_end_pos = false # don't update. we are now outside of the first match
                end
                continue
            end

            # get the start position in species
            if current_nuc_seq_start_species == -1 && occursin("Sbjct  ", line)
                # extract start pos
                current_nuc_seq_start_species = parse(Int64, split(line, "  ")[2])
            end
            
            # always update stop pos if update_end_pos is true. The last one when a new block starts is the correct one
            if update_end_pos && occursin("Sbjct  ", line)
                # extract stop pos
                current_nuc_seq_stop_species = parse(Int64, split(line, "  ")[4])
            end
        end

    end

    # safe the last block
    if(current_nuc_seq_stop_species != -1)
        if haskey(dict, current_protein_id_mel)
            push!(dict[current_protein_id_mel], Match(current_name_of_species_contig, current_score, current_nuc_seq_start_species, current_nuc_seq_stop_species))
        else
            dict[current_protein_id_mel] = [Match(current_name_of_species_contig, current_score, current_nuc_seq_start_species, current_nuc_seq_stop_species)]
        end
    end

    # close file
    close(file)
    # return dictionary
    return dict
end


function parse_blast_fmt6(filename)
    file = open(filename)
    matches = Vector{Match}()

    for line in eachline(file)
        split_line = split(line, "\t")
        start_pos = parse(Int64, split_line[9])
        stop_pos = parse(Int64, split_line[10])
        score = parse(Float64, split_line[11])
        species_gene_id = String(split_line[2])

        push!(matches,  Match(species_gene_id,score,min(start_pos,stop_pos),max(start_pos,stop_pos)))
    end

    return matches
end

"""
Parse augustus file

Input: 
- filename
- threshold_augustus_score: if set to -1, take all values. Else, take only those augustus results with a score (the score from augustus) above threshold_augustus_score
"""
function parse_augustus(filename,threshold_augustus_score = -1)
    file = open(filename)
    matches = Vector{Match}()

    for line in eachline(file)
        # check if line contains "AUGUSTUS	gene"
        if occursin("AUGUSTUS	CDS", line)
            # extract start_pos, stop_pos, score, augustus_gene_id
            split_line = split(line, "\t")
            start_pos = parse(Int64, split_line[4])
            stop_pos = parse(Int64, split_line[5])
            score = parse(Float64, split_line[6])
            contig_id_in_species = split_line[1]

            if threshold_augustus_score == -1 || score > threshold_augustus_score
                push!(matches,  Match(contig_id_in_species,score,min(start_pos,stop_pos),max(start_pos,stop_pos)))
            end
        end
    end
    return matches
end


function parse_groundtrouth(filename)
    file = open(filename)
    matches = Vector{GroundTruth}()

    for line in eachline(file)
        # check if line contains "AUGUSTUS	gene"
        if occursin("Gnomon	exon", line)
            split_line = split(line, "\t")
            name_of_gene = split_line[1]
            start_pos = parse(Int64, split_line[4])
            stop_pos = parse(Int64, split_line[5])
            unique_id = split(split_line[9],";")[1]

            push!(matches, GroundTruth(name_of_gene, min(start_pos, stop_pos), max(start_pos,stop_pos)))
        end
    end
    return matches
end



function parse_groundtrouth_CDS(filename)
    file = open(filename)
    matches = Vector{GroundTruth}()

    for line in eachline(file)
        if occursin("Gnomon	CDS", line)
            split_line = split(line, "\t")
            name_of_gene = split_line[1]
            start_pos = parse(Int64, split_line[4])
            stop_pos = parse(Int64, split_line[5])
            unique_id = split(split_line[9],";")[1]
            push!(matches, GroundTruth(name_of_gene, min(start_pos, stop_pos),max(start_pos, stop_pos)))
        end
    end
    return matches
end