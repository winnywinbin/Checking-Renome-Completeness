#!/bin/bash


touch rapport_gc_perc.txt
printf "%s\t%s\n" "Genome" "GC%" >> rapport_gc_perc.txt;
for file in GCF_*
do {
	seq=$(cat ${file} | egrep '[ATCG]{10,}' | tr -d '\n');
	gc=$(cat ${file} | egrep '[ATCG]{10,}' | tr -d '\n' | egrep -o '[GC]' | wc -w);
	perc=$( echo ${gc} ${#seq} | awk '{print $1/$2*100}');
	/usr/bin/printf "%s\t%.1f\n" ${file} ${perc} >> rapport_gc_perc.txt;
}
done;
