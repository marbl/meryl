#! /bin/bash

if [[ "$#" -lt 3 ]]; then
	echo "Usage: ./_submit_cn.sh <read-db.meryl> <k> <asm1.fasta> [asm2.fasta] <out>"
	exit -1
fi

readdb=$1
k=$2
asm1=$3
asm2=$4
out=$5

args="$readdb $k $asm1 $asm2 $out"
if [ -z $out ]; then
	out="$asm2"
	args="$readdb $k $asm1 $out"
fi

cpus=32
mem=48g
name=$out.spectra-cn
script="$tools/meryl/scripts/eval/spectra-cn.sh"
partition=quick
walltime=4:00:00
path=`pwd`

mkdir -p logs
log=logs/$name.%A.log

echo "\
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args"
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args 
