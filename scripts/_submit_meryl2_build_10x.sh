#!/bin/bash

# meryl_count=/path/to/meryl/scripts/meryl_count/
meryl_count=$tools/meryl/scripts/meryl_count

if [ -z $1 ]; then
    echo "Usage: ./_submit_meryl2_build.sh <k-size> <R1.fofn> <R2.fofn> <out_prefix> [extra]"
    echo -e "\t<k-size>: kmer size k"
    echo -e "\t<R1.fofn>: Read 1. The first 23 bases will get stripped off."
    echo -e "\t<R2.fofn>: Read 2. Will be processed as normal."
    echo -e "\t<out_prefix>: Final merged meryl db will be named as <out_prefix>.meryl"
    exit -1
fi

k=$1
R1=$2
R2=$3
out_prefix=$4
extra=$5

mkdir -p logs

# Split files >10GB
cpus=20
mem=4g
name=$out_prefix.split
partition=quick
walltime=4:00:00
path=`pwd`
log=logs/$name.%A_%a.log

wait_for=""
split=0

script=$meryl_count/split_10x.sh
LEN=`wc -l $R1 | awk '{print $1}'`

echo "R1 will be split to trim off the barcodes."
split_arrs="1-$LEN"
args="$R1"
echo "
sbatch -D $path -J $name --array=$split_arrs --partition=$partition --mem=$mem --cpus-per-task=$cpus --time=$walltime --error=$log --output=$log $script $args"
sbatch -D $path -J $name --array=$split_arrs --partition=$partition --mem=$mem --cpus-per-task=$cpus --time=$walltime --error=$log --output=$log $script $args > split_jid
split_jid=`cat split_jid`
wait_for="${wait_for}afterok:$split_jid,"


#####
echo "$R2 will be split if >12G"

script=$meryl_count/split.sh
LEN2=`wc -l $R2 | awk '{print $1}'`

split_arrs=""

for i in $(seq 1 $LEN2)
do
        fq=`sed -n ${i}p $R2`
        GB=`du -k $fq  | awk '{printf "%.0f", $1/1024/1024}'`
        if [[ $GB -lt 12 ]]; then
            echo "$fq is $GB, less than 12GB. Skip splitting."
            echo $fq >> $R2.$i
        else
            echo "$fq is $GB, over 12GB. Will split and run meryl in parallel. Split files will be in $R2.$i"
	    split_arrs="$split_arrs$i," # keep the line nums $i to split
            split=1
            echo
        fi
done

if [[ $split -eq 1 ]]; then
        split_arrs=${split_arrs%,}
	args="$R2"
        echo "
        sbatch -D $path -J $name --array=$split_arrs --partition=$partition --mem=$mem --cpus-per-task=$cpus --time=$walltime --error=$log --output=$log $script $args"
        sbatch -D $path -J $name --array=$split_arrs --partition=$partition --mem=$mem --cpus-per-task=$cpus --time=$walltime --error=$log --output=$log $script $args >> split_jid
        split_jid=`cat split_jid | tail -n1`
        wait_for="${wait_for}afterok:$split_jid,"
fi


cpus=2
mem=1g
name=$out_prefix.concat
script=$meryl_count/concat_splits.sh
args="$k $R1 $out_prefix $R2"
partition=quick
walltime=10:00
path=`pwd`
log=logs/$name.%A.log
wait_for="--dependency=${wait_for%,}"
echo "$wait_for"
echo "
sbatch -D $path -J $name --partition=$partition --mem=$mem --cpus-per-task=$cpus --time=$walltime $wait_for --error=$log --output=$log $script $args"
sbatch -D $path -J $name --partition=$partition --mem=$mem --cpus-per-task=$cpus --time=$walltime $wait_for --error=$log --output=$log $script $args

