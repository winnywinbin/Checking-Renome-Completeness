from itertools import islice
import numpy as np
import xlsxwriter
import pandas as pd
import os


def filesMaken(files):
    """
    This function adds the files that are being used in this script. It
    receives a list with names of files that are needed.
    :param files:
    :return:
    """
    for i in files:
        file = open(i, "w")
        file.close()


def opschonen(inputfile, outputfile):
    """
    In this function the output from the R-script is cleaned. This needs
    to be done in order to get a set amount of lines per test instead
    of an irregular amount so it is more easy to work with. The
    filtering is done by looking what is in the line, if it's usefull
    the line gets appended to a list. When the loop reaches the next
    test the list is send to a function for formatting and writes it to
    a file.
    :param outputfile:
    :return:
    """
    test = []
    outfile = open(outputfile, "a")
    with open(inputfile, "r") as file:
        while True:
            line = file.readline()
            if not line:
                break
            if "alternative" in line or len(line) == 1 or "Fisher" in line \
                    or "NULL" in line or "data" in line:
                pass
            else:
                test.append(line.strip("\n"))
            if "Fisher" in line:
                if len(test) > 0:
                    outfile.write(testFormat(test))
                test = []


def testFormat(testList):
    """
    Here the list from the function opschonen is send to get formatted
    and is returned in the correct format.
    :param testList:
    :return formatted list:
    """
    teststr = ""
    for i in testList:
        teststr += i + "\n"
    return teststr


def fillDataFrame(presentcodons, filenaam, sheets):
    """
    In this function the data from the cleaned file is put in seperate
    dataframes. One for the p-value and four others for the values from
    the fisher exact test. This is done with a while loop which grabs
    a seperate test from the cleaned file an processes it.
    :param presentcodons:
    :param filenaam:
    :param sheets:
    :return a list of dataframes:
    """
    dataFrames = []
    for sheet in sheets:
        corelations = {}
        for i in presentcodons:
            corelations[i] = [None] * len(presentcodons)
        with open(filenaam, "r") as file:
            while True:
                test = list(islice(file, 10))
                if not test:
                    break
                codons = get_codons(test)
                if float(test[0].split(" ")[2]) < 1.1 and "N" not in codons[
                    0] and "N" not in codons[1]:
                    if sheet == "P-value":
                        negPos = calcCor(testAssesment(test))
                        pvalue = str(float(test[0].split(" ")[2]))
                        if "e" in pvalue:
                            x = pvalue.split("e")
                            if len(x[0]) > 4:
                                pvalue = x[0][0:4] + "E" + x[1]
                        corelations[codons[0]][presentcodons.index(codons[1])] \
                            = negPos + pvalue.upper()
                    else:
                        values = testAssesment(test)
                        corelations[codons[0]][presentcodons.index(codons[1])] \
                            = int(values[(sheets.index(sheet)-1)])
        dataFrames.append(pd.DataFrame(corelations, index=presentcodons))
    return dataFrames


def calcCor(values):
    """
    Here is calculated if the corelation is either positive of negative. This
    is done by using the Fi-coefficient and is indicated by + or -.
    :param values:
    :return + or -:
    """
    if int(values[0])+int(values[3]) == 0:
        return "-"
    elif int(values[1])+int(values[2]) == 0:
        return ""
    else:
        var = ((int(values[0])+int(values[3]))/(int(values[1])+int(values[2])))
        if var > 1:
            return ""
        else:
            return "-"


def get_codons(test):
    """
    Here the present codons in the fisher exact test are filtered out
    from the other character in the line.
    :param test:
    :return a list of the two codons:
    """
    codon1 = test[6].strip(" ").strip("\n")
    codon2 = test[7].split(" ")[0]
    return [codon1, codon2]


def testAssesment(test):
    """
    In this function the table which is created with the fisher exact
    test is filtered. This way only the values are left.
    :param test:
    :return list with the values from the table:
    """
    kruis = []
    for i in [test[8], test[9]]:
        temp_list = []
        regel = i.strip("\n").split(" ")
        for y in regel:
            if y != "":
                temp_list.append(y)
        kruis.append(list(temp_list[1:]))
    kruislist = []
    for i in kruis:
        for y in i:
            kruislist.append(y)
    return kruislist


def naarExcel(dataFrames, file, sheets):
    """
    In this function the list with dataframes is writen to a excel file.
    Each dataframe gets a seperate sheet in the excel file.
    :param dataFrames:
    :param file:
    :param sheets:
    :return:
    """
    writer = pd.ExcelWriter(file, engine='xlsxwriter')
    for df, sheet in zip(dataFrames, sheets):
        df.to_excel(writer, sheet_name=sheet)
        workbook = writer.book
        worksheet = writer.sheets[sheet]
        if sheet != "P-value":
            worksheet.conditional_format('C3:BN66', {'type': '2_color_scale'})
    writer.save()


def delTempFiles(file):
    """
    Here the temporary file with the cleaned data is deleted after it
    has been used.
    :param file:
    :return:
    """
    if os.path.exists(file):
        os.remove(file)


def main():
    """
    This is the main function. Here al the functions above are called in
    the correct order.
    :return:
    """

    alle_codons = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT',
                       'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT',
                       'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT',
                       'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT',
                       'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT',
                       'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT',
                       'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT',
                       'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG',
                       'TTT']
    files = ["tempfile", outputfile]
    sheets = ["P-value", "Lb", "Rb", "Lo", "Ro"]
    filesMaken(files)
    opschonen(inputfile, files[0])
    dataFrames = fillDataFrame(alle_codons, files[0], sheets)
    naarExcel(dataFrames, files[1], sheets)
    delTempFiles(files[0])


"""
The if-statement below is used so this script can be called from a 
snakefile.
"""

if __name__ == "__main__":
    inputfile = snakemake.input.results
    outputfile = snakemake.output.excel
    main()
