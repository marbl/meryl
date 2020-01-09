#!/bin/bash

echo "Usage: ./spectra-cn.sh <read.meryl> <k> <asm1.fasta> [asm2.fasta] out-prefix"
echo -e "\t<read.meryl>\t: Generated with meryl count from i.e. illumina wgs reads"
echo -e "\t<k>\t: k-size used to build <read.meryl>"
echo -e "\t<asm1.fasta>\t: haplotype 1 assembly. gzipped or not"
echo -e "\t[asm2.fasta]\t: haplotype 2 assembly. gzipped or not"
echo -e "\t<out-prefix>: output prefix. Required."
echo
echo -e "<asm1.spectra-cn> will be generated as output."
echo -e "\tWhen asm2 and out-prefix are given,"
echo -e "\t\t<asm2.spectra-cn> and"
echo -e "\t\t<out-prefix>.spectra-hap.hist"
echo -e "\twill be generated in addition."
echo


if [[ $# -lt 4  ]]; then
    echo "No args provided. Exit."
    exit -1
fi

if [ -z $scripts ]; then
    echo "Setting \$scripts path"
    scripts=`which meryl`
    scripts=`dirname $scripts`
    scripts="$scripts/../../scripts"
    echo "Done!"
fi

read=$1
k=$2
asm1_fa=$3
asm2_fa=$4
name=$5

if [ -z $name ]; then
	name=$4
	asm2_fa=""
fi

if [ -s $name ]; then
        echo "$name already exists. Provide a different name."
        exit -1
fi


## Remove this line if R is properly installed ##
module load R					#
#################################################

asm1=`echo $asm1_fa | sed 's/.fasta$//g' | sed 's/.fa$//g' | sed 's/.fasta.gz$//g' | sed 's/.fa.gz$//g'`

echo "Generate spectra-cn.hist per assemblies"
for asm_fa in $asm1_fa $asm2_fa	# will generate only for asm1_fa if asm2_fa is empty
do
	asm=`echo $asm_fa | sed 's/.fasta$//g' | sed 's/.fa$//g' | sed 's/.fasta.gz$//g' | sed 's/.fa.gz$//g'`

	# Generate meryl db per assembly
	if [ ! -e $asm.meryl ]; then
	    echo
	    echo "Generate meryl db for $asm"
	    meryl count memory=31 k=$k output ${asm}.meryl $asm_fa
	fi

	echo
	echo "Collect read counts per asm copies"
	hist=$name.asm.$asm.spectra-cn.hist
	echo -e "Copies\tkmer_multiplicity\tCount" > $hist
	echo "
	meryl difference output read.k$k.$asm.0.meryl $read $asm.meryl"
	meryl difference output read.k$k.$asm.0.meryl $read $asm.meryl
	meryl histogram read.k$k.$asm.0.meryl | awk '{print "read-only\t"$0}' >> $hist

	for i in $(seq 1 4)
	do
	    echo "
	    meryl intersect output read.k$k.$asm.$i.meryl $read [ equal-to $i ${asm}.meryl ]"
	    meryl intersect output read.k$k.$asm.$i.meryl $read [ equal-to $i ${asm}.meryl ]
	    meryl histogram read.k$k.$asm.$i.meryl | awk -v cn=$i '{print cn"\t"$0}' >> $hist
	    rm -r read.k$k.$asm.$i.meryl
	done
	echo "
	meryl intersect output read.k$k.$asm.gt$i.meryl $read [ greater-than $i ${asm}.meryl ]"
	meryl intersect output read.k$k.$asm.gt$i.meryl $read [ greater-than $i ${asm}.meryl ]
	meryl histogram read.k$k.$asm.gt$i.meryl | awk -v cn=">$i" '{print cn"\t"$0}' >> $hist
	rm -r read.k$k.$asm.gt$i.meryl

	echo "
	meryl difference output $asm.0.meryl ${asm}.meryl $read"
	meryl difference output $asm.0.meryl ${asm}.meryl $read
	meryl statistics ${asm}.0.meryl  | head -n4 | tail -n1 | awk -v asm=$asm '{print asm"\t"$1"\t"$2}' >> $name.asmonly.txt

	echo "Plot $hist"
	echo "
	$scripts/plot/plot_spectra_cn.R -f $hist -o $name.asm.$asm"
	$scripts/plot/plot_spectra_cn.R -f $hist -o $name.asm.$asm

	ASM_ONLY=`meryl statistics ${asm}.0.meryl  | head -n4 | tail -n1 | awk '{print $2}'`
        TOTAL=`meryl statistics ${asm}.meryl  | head -n4 | tail -n1 | awk '{print $2}'`
        ERROR=`echo "$ASM_ONLY $TOTAL" | awk -v k=$k '{print (1-(1-$1/$2)^(1/k))}'`
        QV=`echo "$ASM_ONLY $TOTAL" | awk -v k=$k '{print (-10*log(1-(1-$1/$2)^(1/k))/log(10))}'`
        echo -e "$asm\t$ASM_ONLY\t$TOTAL\t$QV\t$ERROR" >> $name.qv
done

if [[ "$asm2_fa" = "" ]]; then
	echo "No asm2_fa given. Done."

	hist=$name.spectra-asm.hist

	# $asm1 only
	meryl intersect output read.k$k.$asm1.meryl $read ${asm1}.meryl

	# Write output
	echo -e "Assembly\tkmer_multiplicity\tCount" > $hist
	meryl histogram read.k$k.$asm1.0.meryl | awk '{print "read-only\t"$0}' >> $hist
	meryl histogram read.k$k.$asm1.meryl | awk -v hap="$asm1" '{print hap"\t"$0}' >> $hist

	echo "
	Plot $hist"
	echo "
	$scripts/plot/plot_spectra_asm.R -f $hist -o $name.asm"
	$scripts/plot/plot_spectra_asm.R -f $hist -o $name.asm

	echo "Clean up"
	rm -r ${asm1}.0.meryl read.k$k.$asm1.0.meryl read.k$k.$asm1.meryl


	exit 0
fi
asm2=`echo $asm2_fa | sed 's/.fasta$//g' | sed 's/.fa$//g' | sed 's/.fasta.gz$//g' | sed 's/.fa.gz$//g'`
rm -r read.k$k.$asm2.0.meryl

echo
echo "Generate spectra-cn.hist for combined $asm1_fa and $asm2_fa"
echo "\"Is my diploid assembly having k-mers in expected copy numbers?\""
hist=$name.spectra-cn.hist

echo "
Union-sum: ${asm1} + ${asm2} + shared kmer counts (asm)"
meryl union-sum output ${asm1}_${asm2}_union.meryl ${asm1}.meryl ${asm2}.meryl

echo "
0-counts in the asm; only seen in the reads"
meryl difference output read.k$k.0.meryl $read ${asm1}_${asm2}_union.meryl
echo

echo -e "Copies\tkmer_multiplicity\tCount" > $hist
meryl histogram read.k$k.0.meryl | awk '{print "read-only\t"$0}' >> $hist
for i in $(seq 1 4)
do
    echo "
    meryl intersect output read.k$k.$i.meryl $read [ equal-to $i ${asm1}_${asm2}_union.meryl ]"
    meryl intersect output read.k$k.$i.meryl $read [ equal-to $i ${asm1}_${asm2}_union.meryl ]
    meryl histogram read.k$k.$i.meryl | awk -v cn=$i '{print cn"\t"$0}' >> $hist
    rm -r read.k$k.$i.meryl
done

echo "
meryl intersect output read.k$k.gt$i.meryl $read [ greater-than $i ${asm1}_${asm2}_union.meryl ]"
meryl intersect output read.k$k.gt$i.meryl $read [ greater-than $i ${asm1}_${asm2}_union.meryl ]
meryl histogram read.k$k.gt$i.meryl | awk -v cn=">$i" '{print cn"\t"$0}' >> $hist
rm -r read.k$k.gt$i.meryl

echo "
Plot $name.spectra-cn.hist"
echo "
$scripts/plot/plot_spectra_cn.R -f $hist -o $name"
$scripts/plot/plot_spectra_cn.R -f $hist -o $name

echo "
## Count k-mers only seen in the assemblies, not in the reads ##"
echo "
0-counts in the read; only seen in the assembly"
meryl difference output ${asm1}_or_${asm2}.0.meryl ${asm1}_${asm2}_union.meryl $read
meryl intersect output ${asm1}_and_${asm2}.0.meryl $asm1.0.meryl $asm2.0.meryl

meryl statistics ${asm1}_or_${asm2}.0.meryl  | head -n4 | tail -n1 | awk '{print "both\t"$1"\t"$2}' >> $name.asmonly.txt
meryl statistics ${asm1}_and_${asm2}.0.meryl | head -n4 | tail -n1 | awk '{print "shared\t"$1"\t"$2}' >> $name.asmonly.txt

echo "
Get QV"
ASM_ONLY=`meryl statistics ${asm1}_or_${asm2}.0.meryl  | head -n4 | tail -n1 | awk '{print $2}'`
TOTAL=`meryl statistics ${asm1}_${asm2}_union.meryl | head -n4 | tail -n1 | awk '{print $2}'`
ERROR=`echo "$ASM_ONLY $TOTAL" | awk -v k=$k '{print (1-(1-$1/$2)^(1/k))}'`
QV=`echo "$ASM_ONLY $TOTAL" | awk -v k=$k '{print (-10*log(1-(1-$1/$2)^(1/k))/log(10))}'`
echo -e "Both\t$ASM_ONLY\t$TOTAL\t$QV\t$ERROR" >> $name.qv
rm -r ${asm1}_and_${asm2}.0.meryl ${asm1}_or_${asm2}.0.meryl ${asm1}_${asm2}_union.meryl ${asm1}.0.meryl ${asm2}.0.meryl

echo "
## Generate spectra-asm.hist for combined $asm1_fa and $asm2_fa ##"
echo "\"Is the assembled distinct portion bigger in one of the two assemblies?\""
hist=$name.spectra-asm.hist

# Get ${asm1} / ${asm2} / shared kmers
meryl difference output ${asm2}_only.meryl ${asm2}.meryl ${asm1}.meryl
meryl difference output ${asm1}_only.meryl ${asm1}.meryl ${asm2}.meryl
meryl intersect output ${asm1}_shrd.meryl ${asm1}.meryl ${asm2}.meryl

# $asm1 only
meryl intersect output read.k$k.$asm1.meryl $read ${asm1}_only.meryl

# $asm2 only
meryl intersect output read.k$k.$asm2.meryl $read ${asm2}_only.meryl

# shared ($asm1 and $asm2)
meryl intersect output read.k$k.shrd.meryl $read ${asm1}_shrd.meryl

# Write output
echo -e "Assembly\tkmer_multiplicity\tCount" > $hist
meryl histogram read.k$k.0.meryl | awk '{print "read-only\t"$0}' >> $hist
meryl histogram read.k$k.$asm1.meryl | awk -v hap="$asm1-only" '{print hap"\t"$0}' >> $hist
meryl histogram read.k$k.$asm2.meryl | awk -v hap="$asm2-only" '{print hap"\t"$0}' >> $hist
meryl histogram read.k$k.shrd.meryl | awk -v hap="shared" '{print hap"\t"$0}' >> $hist

echo "
Plot $hist"
echo "
$scripts/plot/plot_spectra_asm.R -f $hist -o $name.asm"
$scripts/plot/plot_spectra_asm.R -f $hist -o $name.asm

echo "
Clean up"
rm -r read.k$k.0.meryl read.k$k.$asm1.meryl read.k$k.$asm2.meryl read.k$k.shrd.meryl ${asm1}_only.meryl ${asm2}_only.meryl ${asm1}_shrd.meryl

