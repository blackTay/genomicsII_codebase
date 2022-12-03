using JLD2
using Plots
# augustus = load_object("comparison_simulans/candidate_augustus_matches_to_groundtruth.jld2")

simulans = load_object("output/blastn/simulans_blastn_to_groundtruth.jld2")
willistoni = load_object("output/blastn/willistoni_blastn_to_groundtruth.jld2")
sechellia = load_object("output/blastn/sechellia_blastn_to_groundtruth.jld2")
grimshawi = load_object("output/blastn/grimshawi_blastn_to_groundtruth.jld2")


xs = range(0, 1, length=50)
ys_simulans = [length([elem for elem in simulans if elem[3]>=x]) for x in xs]
ys_willistoni = [length([elem for elem in willistoni if elem[3]>=x]) for x in xs]
ys_sechellia = [length([elem for elem in sechellia if elem[3]>=x]) for x in xs]
ys_grimshawi = [length([elem for elem in grimshawi if elem[3]>=x]) for x in xs]


plot(xs,[ys_simulans,ys_willistoni,ys_sechellia,ys_grimshawi],label=["simulans" "willistoni" "sechellia" "grimshawi"])
title!("Gene finding using blastn")
xlabel!("score: overlap between predicition and ground truth")
ylabel!("# sequences with score >= x")