using JLD2

include("tassilos_library.jl")

# config begin
path = "output/augustus/"
species = "simulans"
# config end

simulans_good_matches = load_object(path*species*"_augustus_to_groundtruth_good_matches.jld2") # this file contains only true postives according to our score measure.

# prints all information on each match
for match in simulans_good_matches
    println("in annotated file (groundtruth):")
    println(match[1])
    println("in nuceltide file of "*species*":")
    println(match[2])
    println("our overlap percentage (score):")
    println(match[3])
    println("")
end

