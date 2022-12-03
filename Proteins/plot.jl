using JLD2, Plots

results_df = load("./results.jld2", "results_df")

# prediction success
plot(
    size=(800, 600),
    title="Prediction success: Swiss-Prot proteins and Pfam profiles identified",
    xlabel="Recall",
    ylabel="Precision"
)

scatter!(results_df.prots_recall, results_df.prots_precision,
    group=["BLASTp_qUNI_$s" for s in results_df.species],
    legend=:topleft,
    shape=[:circle],
    color=[:green, :red, :blue, :orange],
    markersize=10,
    dpi=300
)

scatter!(results_df.protfams_recall, results_df.protfams_precision,
    group=["HMMERp_qUNI_$s" for s in results_df.species],
    legend=:topleft,
    shape=[:diamond],
    color=[:green, :red, :blue, :orange],
    markersize=10,
    dpi=300
)

savefig("/home/mw894/comp_bio/gi/a2/plots/pred_sucess.png")

# annotation success
plot(
    size=(800, 600),
    title="Annotation success: Common Benchmark UniProtKB annotation",
    xlabel="Recall",
    ylabel="Precision"
)

scatter!(results_df.prot_goids_recall, results_df.prot_goids_precision,
    group=["BLASTp_qUNI_$s" for s in results_df.species],
    legend=:bottomleft,
    shape=[:circle],
    color=[:green, :red, :blue, :orange],
    markersize=10,
    dpi=300
)

scatter!(results_df.protfams_goids_recall, results_df.protfams_goids_precision,
    group=["HMMERp_qUNI_$s" for s in results_df.species],
    legend=:bottomleft,
    shape=[:diamond],
    color=[:green, :red, :blue, :orange],
    markersize=10,
    dpi=300
)

savefig("/home/mw894/comp_bio/gi/a2/plots/anno_sucess.png")