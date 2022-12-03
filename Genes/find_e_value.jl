
using BioSequences
using FASTX
using JLD2

include("tassilos_library.jl")


blastn_dict = parse_blast("/Users/black/Documents/programming/julia/genI_ass2/simulans_blasted_all.txt")
candidate_matches = [y[1] for (x,y) in collect(blastn_dict) if y[1].score_of_match != -1]

candidate_scores = [y.score_of_match for y in candidate_matches]

candidate_scores_nonzero = [x for x in candidate_scores if x != 0]




# plot histogram of candidate_scores with 10 bins
scores_hist = histogram(candidate_scores_nonzero, bins=100)


plot!(scores_hist)

