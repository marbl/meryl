#!/bin/bash

if [[ $# -lt 4 ]]; then
	echo "Usage: ./plot_block.sh <in.sort.bed> <out> <num_switch> <short_range> [include_gaps]"
	echo "	<in.sort.bed>:	generated with trio/phase_block.sh"
	echo "	<out>:		output prefix; automatically appends given <num_switch> and <short_range>"
	echo "	<num_switch>:	number of switches allowed in <short_range>"
	echo "	<short_range>:	interval to be determined as short-range switch (bp)"
	echo "	<include_gaps>:	T for including gaps. Only set for phased assemblies"
	echo "Arang Rhie, 2020-01-07. arrhie@gmail.com"
	exit -1
fi

bed=$1
out=$2
num_switch=$3
short_range=$4
include_gaps=$5

out=$out.${num_switch}_$short_range

if [ -z $scripts ]; then
    echo "Setting \$scripts path"
    scripts=`which meryl`
    scripts=`dirname $scripts`
    scripts="$scripts/../../scripts"
    scripts=`readlink -e $scripts`
    echo "Done!"
fi

if [ -s $out.phased_block.bed 2> /dev/null ]; then
	rm $out.phased_block.bed
	rm $out.switch.bed 2> /dev/null
fi

echo "
java -jar -Xmx1g $scripts/trio/bedMerToPhaseBlock.jar $bed $out $num_switch $short_range $include_gaps"
java -jar -Xmx1g $scripts/trio/bedMerToPhaseBlock.jar $bed $out $num_switch $short_range $include_gaps

SWITCH_ERR=`awk -v swi=0 -v tot=0 '{swi+=$(NF-1); tot+=$NF} END { print swi"\t"tot"\t"((100.0*swi)/tot)"%" }' $out.phased_block.bed`
echo "Num. switches / Total markers / Switch error rate (%): $SWITCH_ERR"

echo "
java -jar -Xmx1g $scripts/eval/bedCalcN50.jar $out.phased_block.bed | tail -n1 | awk -v out=$out -v swi=\"$SWITCH_ERR\" '{print out\"\t\"\$0\"\tswi}' - >> $out.phased_block.stats"
java -jar -Xmx1g $scripts/eval/bedCalcN50.jar $out.phased_block.bed | tail -n1 | awk -v out=$out -v swi="$SWITCH_ERR" '{print out"\t"$0"\t"swi}' - >> $out.phased_block.stats

count=$out.phased_block.counts

echo -e "Block\tRange\tMat\tPat\tSize" > $count
awk '{ swi=$(NF-1); tot=$NF; {if ($4=="mat") { mat=(tot-swi); pat=swi; } else if ($4=="pat") { mat=swi; pat=(tot-swi); }} {print $4"\t"$1"_"$2"_"$3"\t"mat"\t"pat"\t"($3-$2)}}' $out.phased_block.bed >> $count

module load R

echo "
Rscript $scripts/plot/plot_blob.R -f $count -o $out.phased_block.blob"
Rscript $scripts/plot/plot_blob.R -f $count -o $out.phased_block.blob


