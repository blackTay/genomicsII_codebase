using JLD2

# THIS FILE OUTPUT S THE SENSITIVTY CORRECTLY 

include("tassilos_library.jl")

println("for calculating sensitivity.")

species = "grimshawi" # one of simulans willistoni sechellia grimshawi
threshold_score = 0.97 # the thresold for our own percentage overlap score when a hit should be considered a match

threshold_augustus_score = 0.8


# TODO CHANGE WHEN BLAST / AUG CHANGE
# groundtruth_data = parse_groundtrouth_CDS(species*".groundtruuuth")
groundtruth_data = parse_groundtrouth(species*".groundtruuuth")



# TODO CHANGE THIS FOR AUG/ BLASTN
# check for each element in groundtruth_data if it has a match
# tblastn_species_comparision_to_groundtruth = load_object("output/blastn/"*species*"_blastn_to_groundtruth.jld2")
# tblastn_species_comparision_to_groundtruth_filtered = [elem for elem in tblastn_species_comparision_to_groundtruth if elem[3]>=threshold_score]
aug_species_comparision_to_groundtruth_filtered = load_object("output/augustus/"*species*"_augustus_to_groundtruth_good_matches.jld2")

println("species: ", species)

let prog_cnt = 0
    let counter = 0
        @inbounds for elem in groundtruth_data
            prog_cnt += 1
            if prog_cnt % 10 == 0
                println(prog_cnt/length(groundtruth_data))
            end
            # todo change this for augustus / BLASTN
            @inbounds for comp in aug_species_comparision_to_groundtruth_filtered
            # @inbounds for comp in tblastn_species_comparision_to_groundtruth_filtered
                if elem.chromosome_id == comp[1].chromosome_id && elem.start_pos == comp[1].start_pos && elem.stop_pos == comp[1].stop_pos
                    counter += 1
                    break
                end
            end
        end
        println("augustus")
        println("counter: ", counter)
        println("sensitiviy",round(counter/length(groundtruth_data), digits=3))
        println("species: ", species)
        println("threshold_score: ", threshold_score)
    end
end
