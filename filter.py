
def functie(snakemakeIN, snakemakeOUT):
    # Read file and make list for every line
    file = open(snakemakeIN.data, 'r')
    content = file.readlines()
    list = []
    for line in content:
        list.append(line.strip("\n").split("\t"))
    #file.close()

    # Check anti-codons for chars other than A, T, C or G and add Genome name to list
    delete = []
    for item in list:
        try:
            anticodon = item[5]
            check = 0
            for char in anticodon:
                if char in ['A', 'T', 'C', 'G']:
                    check += 1
            if check != 3 and item[0] not in delete:
                delete.append(item[0])
        except IndexError:
            pass

    # Clean data and write to new file
    file_clean = open(snakemakeOUT.filtered_data, "w")
    file_clean.write("\t".join(list[0]))
    file_clean.write("\n")
    for item in list:
        if item[0] not in delete:
            line = "\t".join(item)
            file_clean.write(line)
            file_clean.write("\n")
    file_clean.close()


if __name__ == "__main__":
    functie(snakemake.input, snakemake.output)
