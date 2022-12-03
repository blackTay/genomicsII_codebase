using FASTX

in_path = "/home/mw894/comp_bio/genomics/a2/data/drosophila_simulans.fasta"
out_path = "/home/mw894/comp_bio/genomics/a2/data/drosophila_simulans_sample.fasta"

sample_first_n = function (in_path, out_path, n)
    i = 0
    reader = open(FASTA.Reader, in_path)
    writer = open(FASTA.Writer, out_path)
    for rec in reader
        if i == n
            break
        end
        write(writer, rec)
        i += 1
    end

    close(writer)
end

sample_first_n(in_path, out_path, 10)