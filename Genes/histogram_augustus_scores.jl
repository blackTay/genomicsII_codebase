using Plots

species = "sechellia" # simulans willistoni sechellia grimshawi
candidate_matches = parse_augustus("/Users/black/Documents/programming/julia/genI_ass2/input/augustus/"*species*"/augustus.gtf")

# plot histogram of scores
scores = [elem.score_of_match for elem in candidate_matches]
histogram(scores, bins=100, xlabel="score", ylabel="count", title="Histogram of Augustus scores for "*species, legend=false)

