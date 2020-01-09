#! /bin/bash

if [[ "$#" -lt 5 ]]; then
	echo "Usage: ./cn_hap.sh <read-db.meryl> <hap1.meryl> <hap2.meryl> <k> <asm1.fasta> [asm2.fasta] <out>"
	exit -1
fi

readdb=$1
hap1=$2
hap2=$3
k=$4
asm1=$5
asm2=$6
out=$7

if [ -e $out ]; then
        echo "$out already exists. Provide a different name. (Are we missing the <out>?)"
        exit -1
fi

if [ -z $out ]; then
	out=$6
	asm2=""
fi

mkdir -p logs

echo "
Get spectra-cn plots and QV stats"
name=$out.spectra-cn
log=logs/$name.log
nohup sh "$tools/meryl/scripts/eval"/spectra-cn.sh $readdb $k $asm1 $asm2 $out > $log

echo "
Get blob plots"
name=$out.blob
log=logs/$name.log
nohup sh $tools/meryl/scripts/trio/hap_blob.sh $hap1 $hap2 $asm1 $asm2 $out > $log

echo "
Get haplotype specfic spectra-cn plots"
name=$out.spectra-hap
log=logs/$name.log
nohup sh $tools/meryl/scripts/trio/spectra-hap.sh $readdb $hap1 $hap2 $k $asm1 $asm2 $out > $log

echo "
Get phase blocks"
name=$out.phase-block
log=logs/$name.log

echo "
For $asm1" > $log
nohup sh $tools/meryl/scripts/trio/phase_block.sh $asm1 $hap1 $hap2 $out.${asm1/.fasta/}  >> $log

if [ -z $asm2 ] ; then
	exit 0
fi

echo "
For $asm2" >> $log
nohup sh $tools/meryl/scripts/trio/phase_block.sh $asm2 $hap1 $hap2 $out.${asm2/.fasta/}  >> $log

