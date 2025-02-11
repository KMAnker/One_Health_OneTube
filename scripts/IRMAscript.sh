#!/bin/bash

read -p "Are you running this script from one of the following folders: run/human or run/avian or run/swine? Or where you have the IRMA output folders (yes/no)"$'\n' yn
    case $yn in
        yes ) echo "ok, we proceed";;
        no ) echo "exiting..."$'\n'"You better go to the right directory :)"; 
	    return ;;
        * ) echo "Please answer yes or no.";;
    esac


echo "sample: "$'\t'"Min_coverage"$'\t'"Max_coverage"$'\t'"Average_coverage"  > cov.txt
echo  "Min_quality"$'\t'"Max_quality"$'\t'"Average_quality" > qual.txt
echo "%_Missing_sites(N)" > Ns.txt

for dir in */; do
cd $dir
for file in tables/*coverage.txt; do 
SAMP="$(basename "`pwd`")"
SEG=$(basename $file -coverage.txt )
awk -v var=${SAMP}_${SEG} ' NR==2 {min=$3; max=$3} $3 > max {max=$3} $3 < min {min=$3} {s+=$3} END {print var"\t"min"\t"max"\t"s/(NR)}' $file >> ../cov.txt
awk ' NR==2 {min=$8; max=$8} $8 > max {max=$8} $8 < min {min=$8} {s+=$8} END {print min"\t"max"\t"s/(NR)}' $file >> ../qual.txt
awk 'BEGIN {counter = 0}{ if($4 == "N") {counter++}}END{print "\t"(counter/NR) *100}' $file >> ../Ns.txt
done; cd ../ ; done 

paste cov.txt qual.txt Ns.txt > stats.txt
 
echo "Sample "$'\t'"Min_coverage"$'\t'"Max_coverage"$'\t'"Average_coverage"$'\t'"Min_quality"$'\t'"Max_quality"$'\t'"Average_quality"$'\t'"%_Missing_sites(N)" > mapping_stats.xls

awk 'NR>1 {printf "%s\t %s\t %s\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\n",$1,$2,$3,$4,$5,$6,$7,$8}' stats.txt >> mapping_stats.xls

rm stats.txt
rm cov.txt 
rm qual.txt
rm Ns.txt

awk 'NR>1 {if ($8>10 || $7< 30) print $1"\t""=>""\t""Fail"; else print $1"\t""=>""\t""Pass"}' mapping_stats.xls > segments_validation.txt
