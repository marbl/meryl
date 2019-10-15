#!/bin/bash

echo "Usage: ./spectra-hap-cn.sh <reads.meryl> <hap1.meryl> <hap2.meryl> <k> <asm1.fasta> [asm2.fasta out-prefix]"
read=$1
read_hap1=$2 # .meryl    Haplotype1 specific kmers with counts from reads
read_hap2=$3 # .meryl    Haplotype2 specific kmers with counts from reads
k=$4    # kmer
asm1_fa=$5 # .fasta    Haplotype1 assembly
asm2_fa=$6 # .fasta    Haplotype2 assembly
name=$7    # output prefix
if [ -z $name ]; then
	name="meryl"
fi

if [ -z $read_hap1 ]; then
	echo "No input provided. Exit."
	exit -1
fi

read=${read/.meryl/}		  # all read counts
read_hap1=${read_hap1/.meryl/}    # pat specific mers with read counts
read_hap2=${read_hap2/.meryl/}    # mat specific mers with read counts

hap_hist=$name.spectra-hap.hist
cn_hist=$name.spectra-hap-cn.hist

for hap in $read_hap1 $read_hap2
do
    # Get histogram from all and inherited hap-mers
    echo -e "Assembly\tkmer_multiplicity\tCount" > $hap.$hap_hist
    meryl histogram $read.meryl | awk -v asm="read-total" '{print asm"\t"$0}' >> $hap.$hap_hist
    meryl histogram $hap.meryl | awk -v asm="${hap}-inherited" '{print asm"\t"$0}' >> $hap.$hap_hist
done

for asm_fa in $asm1_fa $asm2_fa
do
    asm=`echo $asm_fa | sed 's/.fasta$//g' | sed 's/.fa$//g' | sed 's/.fasta.gz$//g' | sed 's/.fa.gz$//g'`
    if ! [[ "$(ls -A $asm.meryl 2> /dev/null )" ]]; then
        echo "Generate meryl db for $asm_fa"
        meryl count k=$k output $asm.meryl $asm_fa
    fi
    
    echo -e "Copies\tkmer_multiplicity\tCount" > $asm.$cn_hist

    # For each haplotype
    for read in $read_hap1 $read_hap2
    do
	echo -e "Start processing $read from $asm\n"

	# Convert counts to read counts
	meryl intersect output $read.$asm.meryl $read.meryl $asm.meryl
	meryl statistics $read.$asm.meryl | head -n3 | tail -n1 | awk -v asm=$asm -v read=$read '{print asm"\t"read"\t"$1"\t"$2}' >> hap.counts
	meryl statistics $read.$asm.meryl | head -n4 | tail -n1 | awk -v asm=$asm -v read=$read '{print asm"\t"read"\t"$1"\t"$2}' >> hap.counts
	meryl histogram $read.$asm.meryl  | awk -v asm=$asm '{print asm"\t"$0}' >> $read.$hap_hist
	rm -r $read.$asm.meryl

:<<'END'
	# Copy-number histogram per haplotype specific k-mers
	meryl difference output read.$read.0.meryl $read.meryl $asm.meryl
	meryl histogram read.$read.0.meryl | awk -v read=$read '{print read"-only\t"$0}' >> $asm.$cn_hist
	rm -r read.$read.0.meryl

	for i in $(seq 1 4)
	do
	    meryl intersect output read.$read.$i.meryl $read.meryl [ equal-to $i ${asm}.meryl ]
	    meryl histogram read.$read.$i.meryl | awk -v cn=$i '{print cn"\t"$0}' >> $asm.$cn_hist
	    rm -r read.$read.$i.meryl
	done

	meryl intersect output read.$read.gt$i.meryl $read.meryl [ greater-than $i ${asm}.meryl ]
	meryl histogram read.$read.gt$i.meryl | awk -v cn=">$i" '{print cn"\t"$0}' >> $asm.$cn_hist
	rm -r read.$read.gt$i.meryl
END

    done
done


