#csv file aanmaken voor aanwezigheid en afwezigheid anticodon.
anticodons="AAA AAC AAG AAT ACA ACC ACG ACT AGA AGC AGG AGT ATA ATC ATG ATT CAA CAC CAG CAT CCA CCC CCG CCT CGA CGC CGG CGT CTA CTC CTG CTT GAA GAC GAG GAT GCA GCC GCG GCT GGA GGC GGG GGT GTA GTC GTG GTT TAA TAC TAG TAT TCA TCC TCG TCT TGA TGC TGG TGT TTA TTC TTG TTT"

cat accessions_genomes.txt | while read line
do {
  data=""
  data=${data}${line};
  for i in ${anticodons}
  do {
    i+="_accessions"
    if grep -Rq $line tree_accesions/${i};
    then {
      nummer=",1.0"
      data=${data}${nummer};
    } else {
      nummer=",0.0"
      data=${data}${nummer};
    }
    fi;
  }
  done;
  echo ${data} >> tdata.csv
}
done;

#header aanmaken op basis aantal anticodons
nummer=$(echo $anticodons | wc -w)
header=""
for ((x=1; x<=${nummer}; x++ ))
do {
  header+=",V${x}";
}
done;

#data opschonen en header toevoegen aan de csv file
echo $header > tdata_clean.csv
cat tdata.csv | egrep -v "^,1.0" >> tdata_clean.csv

