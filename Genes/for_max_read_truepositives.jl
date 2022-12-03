using JLD2

# include("tassilos_library.jl")

# config begin
path = "output/augustus/"
species = "simulans"
# config end

# this file contains only true postives according to our score measure.
simulans_good_matches = load_object("./Genes/" * path * species * "_augustus_to_groundtruth_good_matches.jld2")

cds_match_locations = Tuple[]

# prints all information on each match
for match in simulans_good_matches
    # println("in annotated file (groundtruth):")
    # println(match[1])
    # println("in nuceltide file of " * species * ":")
    push!(cds_match_locations, (match[2].start_pos, match[2].stop_pos,))
    # println("our overlap percentage (score):")
    # println(match[3])
    # println("")
end

