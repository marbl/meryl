#!/bin/bash

echo "Usage: ./spectra-hap-cn.sh <reads.meryl> <hap1.meryl> <hap2.meryl> <k> <asm1.fasta> [asm2.fasta] <out-prefix>"
if [[ $# -lt 6 ]]; then
	echo "Not enough arguements given."
	exit -1
fi

read=$1
read_hap1=$2 # .meryl    Haplotype1 specific kmers with counts from reads
read_hap2=$3 # .meryl    Haplotype2 specific kmers with counts from reads
k=$4    # kmer
asm1_fa=$5 # .fasta    Haplotype1 assembly
asm2_fa=$6 # .fasta    Haplotype2 assembly
name=$7    # output prefix
if [ -z $name ]; then
	name=$6
	asm2_fa=""
fi

if [ -z $read_hap1 ]; then
	echo "No input provided. Exit."
	exit -1
fi

if [ -z $scripts ]; then
    echo "Setting \$scripts path"
    scripts=`which meryl`
    scripts=`dirname $scripts`
    scripts="$scripts/../../scripts"
    echo "Done!"
fi

read=${read/.meryl/}		  # all read counts
read_hap1=${read_hap1/.meryl/}    # pat specific mers with read counts
read_hap2=${read_hap2/.meryl/}    # mat specific mers with read counts

hap_hist=$name.hapmers.hist
cn_hist=$name.spectra-hap-cn.hist

# Get histogram from all and inherited hap-mers
echo -e "Assembly\tkmer_multiplicity\tCount" > $hap_hist
meryl histogram $read.meryl | awk -v asm="read-total" '{print asm"\t"$0}' >> $hap_hist

for hap in $read_hap1 $read_hap2
do
    meryl histogram $hap.meryl | awk -v asm="${hap}" '{print asm"\t"$0}' >> $hap_hist
done

## Remove this line if R is properly installed ##
module load R                                   #
#################################################

echo "
$scripts/plot/plot_spectra_asm.R -f $hap_hist -o ${hap_hist/.hist/}"
$scripts/plot/plot_spectra_asm.R -f $hap_hist -o ${hap_hist/.hist/}

for asm_fa in $asm1_fa $asm2_fa
do
    asm=`echo $asm_fa | sed 's/.fasta$//g' | sed 's/.fa$//g' | sed 's/.fasta.gz$//g' | sed 's/.fa.gz$//g'`
    if ! [[ "$(ls -A $asm.meryl 2> /dev/null )" ]]; then
        echo "Generate meryl db for $asm_fa"
        meryl count k=$k output $asm.meryl $asm_fa
    fi
    
    # For each haplotype
    for read_hap in $read_hap1 $read_hap2
    do
	read_hap=${read_hap/.meryl/}

	echo -e "Start processing $read_hap from $asm\n"

        echo -e "Copies\tkmer_multiplicity\tCount" > $asm.$read_hap.$cn_hist

	# Convert counts to read_hap counts
	meryl intersect output $read_hap.$asm.meryl $read_hap.meryl $asm.meryl
	#meryl statistics $read_hap.$asm.meryl | head -n3 | tail -n1 | awk -v asm=$asm -v read=$read '{print asm"\t"read"\t"$1"\t"$2}' >> hap.counts
	#meryl statistics $read_hap.$asm.meryl | head -n4 | tail -n1 | awk -v asm=$asm -v read=$read '{print asm"\t"read"\t"$1"\t"$2}' >> hap.counts
	#meryl histogram $read_hap.$asm.meryl  | awk -v asm=$asm '{print asm"\t"$0}' >> $read_hap.$hap_hist
	#rm -r $read_hap.$asm.meryl

	# Copy-number histogram per haplotype specific k-mers
        meryl difference output read.$read_hap.0.meryl $read_hap.meryl $asm.meryl
        meryl histogram read.$read_hap.0.meryl | awk -v read=$read_hap '{print read"-only\t"$0}' >> $asm.$read_hap.$cn_hist
        rm -r read.$read_hap.0.meryl

	for i in $(seq 1 4)
	do
	    meryl intersect output read.$read_hap.$i.meryl $read_hap.$asm.meryl [ equal-to $i ${asm}.meryl ]
	    meryl histogram read.$read_hap.$i.meryl | awk -v cn=$i '{print cn"\t"$0}' >> $asm.$read_hap.$cn_hist
	    rm -r read.$read_hap.$i.meryl
	done

	meryl intersect output read.$read_hap.gt$i.meryl $read_hap.$asm.meryl [ greater-than $i ${asm}.meryl ]
	meryl histogram read.$read_hap.gt$i.meryl | awk -v cn=">$i" '{print cn"\t"$0}' >> $asm.$read_hap.$cn_hist
	rm -r read.$read_hap.gt$i.meryl

        echo "
        $scripts/plot/plot_spectra_cn.R -f $asm.$read_hap.$cn_hist -o $name.$asm.$read_hap"
        $scripts/plot/plot_spectra_cn.R -f $asm.$read_hap.$cn_hist -o $name.$asm.$read_hap

    done

done

