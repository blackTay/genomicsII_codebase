using JLD2

# config begin
threshold_score = 0.97 # 0.7 # the thresold for our own percentage overlap score when a hit should be considered a match

threshold_augustus_score = 0.8
# THIS NEEDS TO BE SET IN compare_outputs.jl as well!!!!
 # the threshold form the  augustus score when we should take it as an augustus proposed entry
do_blastn = true
do_tblastn = true
do_augustus = true
# confing end

include("tassilos_library.jl")

println("for calculating precision.")

function calculate_sensitivity_and_precision(species)

    println("\nSpecies: ",species)

    if do_blastn || do_tblastn
        blast_groundtruth = length(parse_groundtrouth_CDS(species*".groundtruuuth"))
    end

    if do_blastn
        # blastn
        blastn_numbermatches = length(parse_blast_fmt6("/Users/black/Documents/programming/julia/genI_ass2/input/blastn/"*species*"_out.txt"))

        blastn_species_comparision_to_groundtruth = load_object("output/blastn/"*species*"_blastn_to_groundtruth.jld2")

        blastn_truematches = length([elem for elem in blastn_species_comparision_to_groundtruth if elem[3]>=threshold_score])

        # aka true positive rate
        blastn_sensitivity = blastn_truematches / blast_groundtruth

        # TP / (TP + FP)
        blastn_precision = blastn_truematches / blastn_numbermatches

        println("blastn:")
        # println("sensitiviy:\t",round(blastn_sensitivity,digits=3))
        println("precision:\t", round(blastn_precision,digits=3))
    end


    if do_tblastn
        # tblastn
        tblastn_numbermatches = length(parse_blast_fmt6("/Users/black/Documents/programming/julia/genI_ass2/input/tblastn/"*species*"_out.txt"))

        tblastn_species_comparision_to_groundtruth = load_object("output/tblastn/"*species*"_tblastn_to_groundtruth.jld2")

        tblastn_truematches = length([elem for elem in tblastn_species_comparision_to_groundtruth if elem[3]>=threshold_score])

        # aka true positive rate
        tblastn_sensitivity = tblastn_truematches/ blast_groundtruth

        # TP / (TP + FP)
        tblastn_precision = tblastn_truematches/tblastn_numbermatches

        println("tblastn:")
        # println("sensitiviy:\t",round(tblastn_sensitivity, digits=3))
        println("precision:\t", round(tblastn_precision, digits=3))
    end



    if do_augustus
        # augustus
        filename = "/Users/black/Documents/programming/julia/genI_ass2/input/augustus/"*species*"/augustus.gtf"

        augustus_numbermatches = length(parse_augustus(filename, threshold_augustus_score))

        augustus_species_comparision_to_groundtruth = load_object("output/augustus/"*species*"_augustus_to_groundtruth.jld2")
        augustus_species_comparision_to_groundtruth_good_matches = [elem for elem in augustus_species_comparision_to_groundtruth if elem[3]>=threshold_score]



        save_object("output/augustus/"*species*"_augustus_to_groundtruth_good_matches.jld2",augustus_species_comparision_to_groundtruth_good_matches)

        augustus_truematches = length(augustus_species_comparision_to_groundtruth_good_matches)

        augustus_groundtruth = length(parse_groundtrouth(species*".groundtruuuth"))

        # aka true positive rate
        augustus_sensitivity = augustus_truematches/ augustus_groundtruth

        # TP / (TP + FP)
        augustus_precision = augustus_truematches / augustus_numbermatches

        println("augustus:")
        # print two decimals of augustus_sensitivity
        # println("sensitiviy:\t",round(augustus_sensitivity, digits=3))
        println("precision:\t",round(augustus_precision, digits=3))
    end
end

calculate_sensitivity_and_precision("simulans")
calculate_sensitivity_and_precision("sechellia")
calculate_sensitivity_and_precision("willistoni")
calculate_sensitivity_and_precision("grimshawi") # simulans willistoni sechellia grimshawi)





