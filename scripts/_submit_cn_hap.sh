#! /bin/bash

if [[ "$#" -lt 5 ]]; then
	echo "Usage: ./_submit_cn_hap.sh <read-db.meryl> <hap1.meryl> <hap2.meryl> <k> <asm1.fasta> [asm2.fasta] <out>"
	exit -1
fi

readdb=$1
hap1=$2
hap2=$3
k=$4
asm1=$5
asm2=$6
out=$7

if [ -z $out ]; then
	out=$6
	asm2=""
fi

mkdir -p logs
# All jobs are expected to finish within 4 hours
partition=quick
walltime=4:00:00
path=`pwd`
extra=""

#### Get spectra-cn plots and QV stats
cpus=32
mem=48g
name=$out.spectra-cn
script="$tools/meryl/scripts/eval/spectra-cn.sh"
args="$readdb $k $asm1 $asm2 $out"
log=logs/$name.%A.log

echo "\
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args"
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args > cn.jid
jid=`cat cn.jid`

#### Get blob plots
cpus=8
mem=10g

script="$tools/meryl/scripts/trio/hap_blob.sh"
# ./hap_blob.sh <hap1.meryl> <hap2.meryl> <asm1.fasta> [asm2.fasta] <out>
args="$hap1 $hap2 $asm1 $asm2 $out"
name=$out.blob
log=logs/$name.%A.log

echo "\
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args"
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args > blob.jid


#### Get haplotype specfic spectra-cn plots
cpus=24
mem=10g
extra="--dependency=afterok:$jid"	# Re-uses asm.meryl dbs in spectra-cn.sh.

name=$out.spectra-hap
script="$tools/meryl/scripts/trio/spectra-hap.sh"
args="$readdb $hap1 $hap2 $k $asm1 $asm2 $out"
log=logs/$name.%A.log

echo "\
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args"
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args > hap.jid

#### Get phase blocks
cpus=12
mem=24g
extra=""

script="$tools/meryl/scripts/trio/phase_block.sh"
# ./phase_block.sh <asm.fasta> <mat.meryl> <pat.meryl> <out>

if [ -z $asm2 ] ; then
	# Only one assembly given.
	args="$asm1 $hap1 $hap2 $out"
        name=$out.phase-block
else
	args="$asm1 $hap1 $hap2 $out.${asm1/.fasta/}"
	name=$out.phase-block.${asm1/.fasta/}
fi
log=logs/$name.%A.log

echo "\
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args"
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args > block1.jid

if [ -z $asm2 ] ; then
	exit 0
fi

args="$asm2 $hap1 $hap2 $out.${asm2/.fasta/}"
name=$out.phase-block.${asm2/.fasta/}
log=logs/$name.%A.log

echo "\
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args"
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args > block2.jid

