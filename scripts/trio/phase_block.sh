#!/bin/bash

echo "Usage: ./phase_block.sh <asm.fasta> <hap1.meryl> <hap2.meryl> <out>"

if [[ $# -lt 4 ]]; then
	exit -1
fi

scaff=$1
scaff=${scaff/.fasta/}

hap1=$2
hap2=$3
out=$4

if [ -z $scripts ]; then
    echo "Setting \$scripts path"
    scripts=`which meryl`
    scripts=`dirname $scripts`
    scripts="$scripts/../../scripts"
    scripts=`readlink -e $scripts`
    echo "Done!"
fi

module load bedtools
module load samtools
igvtools=$tools/IGVTools/igvtools


if [ ! -e $scaff.gaps ]; then
	echo "
	Get gaps"
	java -jar -Xmx4g $scripts/trio/fastaGetGaps.jar $scaff.fasta $scaff.gaps
fi
awk '{print $1"\t"$2"\t"$3"\tgap"}' $scaff.gaps > $scaff.gaps.bed
cat $scaff.gaps.bed > $scaff.bed

samtools faidx $scaff.fasta

if [ ! -s $out.sort.bed ]; then
	echo "
	Generate pat and mat marker sites bed"
	for hap in $hap1 $hap2
	do
		hap=${hap/.meryl/}
		echo "
		-- $hap"
		hap_short=${hap/.inherited/}
		if [ ! -s $out.$hap.bed ]; then
			meryl-lookup -dump -memory 4 -sequence $scaff.fasta -mers $hap.meryl | awk -v hap=$hap_short '$(NF-4)=="T" {print $1"\t"$(NF-5)"\t"($(NF-5)+21)"\t"hap}' > $out.$hap.bed
		fi
		cat $out.$hap.bed >> $out.bed

		if [ ! -s $out.$hap.$tdf ]; then
			$igvtools count $out.$hap.bed $out.$hap.tdf $scaff.fasta.fai
		fi
	done

	echo "
	Sort $out.bed"
	bedtools sort -i $out.bed > $out.sort.bed
fi

#$scripts/plot/plot_block.sh <in.sort.bed> <out> <num_switch> <short_range> [include_gaps] 
echo "
$scripts/plot/plot_block.sh $out.sort.bed $out 10 20000 T"
$scripts/plot/plot_block.sh $out.sort.bed $out 10 20000 T

