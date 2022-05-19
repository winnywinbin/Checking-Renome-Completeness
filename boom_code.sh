#the files with accesions numbers of a domain are made here.
anticodons="AAA AAC AAG AAT ACA ACC ACG ACT AGA AGC AGG AGT ATA ATC ATG ATT CAA CAC CAG CAT CCA CCC CCG CCT CGA CGC CGG CGT CTA CTC CTG CTT GAA GAC GAG GAT GCA GCC GCG GCT GGA GGC GGG GGT GTA GTC GTG GTT TAA TAC TAG TAT TCA TCC TCG TCT TGA TGC TGG TGT TTA TTC TTG TTT"

mkdir tree_accesions/
for x in ${anticodons}
do {
  cat output_genomes_tRNA.txt | egrep "${x}" | awk -F '\t' '{print $1}' | sort | uniq > tree_accesions/"${x}_accessions"
}
done;



#16S gene sequence is retrieved and put in a large fasta file
for i in all/*.fna
do {
  barrnap $i > test_bar.gff
  cat test_bar.gff | egrep "16S" | head -n 1 > barr_results.gff
  bedtools getfasta -fi $i -bed barr_results.gff > fasta.fa
  code=$(cat fasta.fa | egrep '>' | awk -F ':' '{print $1}')
  echo $code | cut -c2- >> accessions_genomes.txt
  sed -i "1s/.*/$code/" fasta.fa
  cat fasta.fa >> fasta_compleet.fa
}
done;

# multiple sequence aligment and phylogenetic tree
clustalw -INFILE=fasta_compleet.fa -ALIGN -OUTFILE=domain.aln
clustalw -INFILE=domain.aln -TREE -OUTFILE=domain.ph

