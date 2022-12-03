using BioSequences
using FASTX
using JLD2
# using Debugger

include("tassilos_library.jl")

threshold_augustus_score = 0.8 # THIS NEEDS TO BE SET IN count_FP_FN.jl as well!!!!

# species is either of simulans willistoni sechellia grimshawi
# tool is either of blastn tblastn augustus
function compare_tool_to_groundtruth(species, tool)
    println("species: ",species)
    println("tool: ",tool)
    println("parsing...")
    # TODO remove "_mRNA" from command when comparing augustus
    if tool == "augustus"
        groundtruth = parse_groundtrouth(species*".groundtruuuth")
    else
        groundtruth = parse_groundtrouth_CDS(species*".groundtruuuth")
    end

    if tool == "augustus"
        candidate_matches = parse_augustus("/Users/black/Documents/programming/julia/genI_ass2/input/augustus/"*species*"/augustus.gtf",threshold_augustus_score) 
        println("aug length ",length(candidate_matches))
    elseif tool == "blastn"
        # blastn
        candidate_matches = parse_blast_fmt6("/Users/black/Documents/programming/julia/genI_ass2/input/blastn/"*species*"_out.txt")
    elseif tool == "tblastn"
        # tblastn
        candidate_matches = parse_blast_fmt6("/Users/black/Documents/programming/julia/genI_ass2/input/tblastn/"*species*"_out.txt")
    else
        println("tool not supported")
        exit(-1)
        return
    end



    # define list of pairs
    groundtruth_to_candidate_pairs = Vector{Tuple{GroundTruth,Match,Float64}}()

    println("comparing to groundtruth...")

    # create a counter to monitor progress of for loop
    counter = 0
    counter_max = length(candidate_matches)

    # for sensitivity first groundtruth loop, then candidate loop
    # TODO REVERSE ORDER OF TWO FOR LOOPS HERE AGAIN
    for groundtruth_element in groundtruth
        # store the best overlap as type Tuple{GroundTruth,Match,Float64}
        best_overlap = (nothing,nothing,0.0)
            for candidate_match in candidate_matches
            # check that chromosome matches
            if groundtruth_element.chromosome_id != candidate_match.name_of_species_contig
                    continue
            end

            # check how much aug_match and groundtruth_element overlap
            overlap = max(0,min(candidate_match.stop_pos, groundtruth_element.stop_pos) - max(candidate_match.start_pos, groundtruth_element.start_pos))

            # get the proportional overlap to both aug_match and groundtruth_element
            overlap_prop_aug_match = overlap / (candidate_match.stop_pos - candidate_match.start_pos)
            overlap_prop_groundtruth = overlap / (groundtruth_element.stop_pos - groundtruth_element.start_pos)

            # the final overlap proportion is the minimum of the overlaps
            overlap_prop = min(overlap_prop_aug_match, overlap_prop_groundtruth)
            # if(overlap_prop > 0)
            #     print(overlap_prop)
            # end

            # if overlap is better than the best overlap so far, update best_overlap
            if overlap_prop > best_overlap[3]
                best_overlap = (groundtruth_element,candidate_match,overlap_prop)
            end

        end

        # if the overlap is at least 10% add the pair to the list
        if best_overlap[3] >= 0.1
            push!(groundtruth_to_candidate_pairs, best_overlap)
        end

        # update counter
        counter += 1
        if counter % 1000 == 0
            println("progress: ",counter/counter_max)
        end
    
    end

    save_object("output/"*tool*"/"*species*"_"*tool*"_to_groundtruth_for_sensitivity.jld2",groundtruth_to_candidate_pairs)
end


compare_tool_to_groundtruth("grimshawi","tblastn") # one of simulans willistoni sechellia grimshawi
#compare_tool_to_groundtruth("simulans","tblastn") 
# compare_tool_to_groundtruth("willistoni","tblastn") 
#compare_tool_to_groundtruth("sechellia","tblastn") 