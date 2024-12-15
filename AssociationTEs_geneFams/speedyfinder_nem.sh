#!/bin/bash

for i in *; do
  cd $i;
  for file in *.EDTA.TEanno.gff3.mod ; do
    echo $file;
    sort -k1,1 -k4,4n -k5,5n "$file" > "$file.sorted";
    bedtools closest -k 10 -D a -a gene.bed.sorted -b "$file.sorted" | \
    awk '$16 >= - 500 && $16 <= 500' | sort -k4,4 | awk '$8 != - 1' > "${file}.sorted.closest";
    rm *.mod.sorted ;
    join -1 2 -2 4 -t$'\t' -o2.9,2.12,2.15,2.16,1.2,1.3,1.4 pfam.tsv <(sort -k4 ${file}.sorted.closest) | sort -u > ${file}.sorted.closest.join ;
    rm *.closest ;
    sort -k6V ${file}.sorted.closest.join | awk 'function p(){print l,c; delete a; b=c=0} NR!=1&&l!=$6{p()} ++a[$3]==1{c++} {l=$6} END{p()}' | grep -v "^ " > ${file}.nTE_per_IPR ;
    rm *.closest.join ;
  done ;
  count=1000 ;
  for file in *.EDTA.TEanno.gff3.mod ; do
    echo "$file";
    for i in $(seq $count); do
        echo $file $i;
        bedtools shuffle -i "$file" -g chrfile > "${file}.$i.random";
        sort -k1,1 -k4,4n -k5,5n "${file}.$i.random" > "${file}.$i.random.sorted";
        bedtools closest -k 10 -D a -a gene.bed.sorted -b "${file}.$i.random.sorted" | \
        awk '$16 >= - 500 && $16 <= 500' | sort -k4,4 | awk '$8 != - 1' > "${file}.$i.random.sorted.closest";
        rm *.random ;
        rm *.random.sorted ;
        join -1 2 -2 4 -t$'\t' -o2.9,2.12,2.15,2.16,1.2,1.3,1.4 pfam.tsv <(sort -k4 ${file}.$i.random.sorted.closest) | sort -u > ${file}.$i.random.sorted.closest.join ;
        rm *.closest ;
        sort -k6V ${file}.$i.random.sorted.closest.join | awk 'function p(){print l,c; delete a; b=c=0} NR!=1&&l!=$6{p()} ++a[$3]==1{c++} {l=$6} END{p()}' | grep -v "^ " > ${file}.$i.random.nTE_per_IPR ;
        rm *.closest.join ;
    done
  done ;
  awk -F '\t' '{print $3}' pfam.tsv | sort -u > IPR.txt ;
  awk '
  FNR==NR{
    a[$1]=$2
    next
  }
  ($1 in a){
    print $0,a[$1]
    b[$1]
    next
  }
  {
    print $1,$2 "0"
  }
  END{
    for(i in a){
      if(!(i in b)){
        print i"0"a[i]
      }
    }
  }
  ' *.mod.nTE_per_IPR IPR.txt | sort > tmp.join ;
  for file in *.random.nTE_per_IPR; do join -a1 -e '0' -o '2.2' tmp.join <(sort ${file}) > tmp.file.$file; done ;
  paste -d " " tmp.join tmp.file.* > TE_per_IPR ;
  rm tmp* ;
  rm *.nTE_per_IPR ;
  for file in *.EDTA.TEanno.gff3.mod ; do
    echo $file;
    sort -k1,1 -k4,4n -k5,5n "$file" > "$file.sorted";
    bedtools closest -k 10 -D a -a gene.bed.sorted -b "$file.sorted" | \
    awk '$16 >= - 500 && $16 <= 500' | sort -k4,4 | awk '$8 != - 1' > "${file}.sorted.closest";
    rm *.mod.sorted ;
    join -1 2 -2 4 -t$'\t' -o2.9,2.12,2.15,2.16,1.2,1.3,1.4 pfam.tsv <(sort -k4 ${file}.sorted.closest) | sort -u > ${file}.sorted.closest.join ;
    rm *.closest ;
    sort -k6V ${file}.sorted.closest.join | awk 'function p(){print l,c; delete a; b=c=0} NR!=1&&l!=$6{p()} ++a[$3]==1{c++} {l=$6} END{p()}' | grep -v "^ " > ${file}.nTE_per_IPR ;
    rm *.closest.join ;
  done ;
  cd .. ;
done
