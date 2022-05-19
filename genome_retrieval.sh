#domain="-A"
#checking which domain is given by the user.
echo $1
if [ $1 = "-A" ]
then {
  name="archaea"
} elif [ $1 = "-B"]
then {
  name="bacteria"
} else {
  name="eukaryote"
}
fi

#download assembly summary of the domain given by the user.
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt 

# Filter for complete genomes.
awk -F "\t" '$12=="Complete Genome" && $11=="latest"{print $20}' assembly_summary.txt > ftpdirpaths

# Identify the FASTA files (.fna.) other files may also be downloaded here.
awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' ftpdirpaths > ftpfilepaths

# Download all genomes in parallel in the directory all. The directory can be changed.
mkdir -p all
while read l; do
	#wget -P all/ $l
	echo " "
done <ftpfilepaths

for i in all/*.gz
do
	gzip -dc $i > $i
done
#set +u; cat ftpfilepaths | parallel -j 20 --verbose --progress "cd all && wget -c --continue {}";
#cat ftpfilepaths | cd all && curl -O {}
#gunzip all/*.fna.gz
