#!/bin/bash

# export PATH=/path/to/canu1.8/bin/:$PATH
# export meryl_count=/path/to/meryl/scripts/

if [ -z $1 ]; then
    echo "Usage: ./_submit_meryl2_build.sh <k-size> <input.fofn> <out_prefix> [extra]"
    echo -e "\t<k-size>: kmer size k"
    echo -e "\t<input.fofn>: ls *.fastq.gz > input.fofn. accepts fasta, fastq, gzipped or not."
    echo -e "\t<out_prefix>: Final merged meryl db will be named as <out_prefix>.meryl"
    exit -1
fi

k=$1
input_fofn=$2
out_prefix=$3
extra=$4

LEN=`wc -l $input_fofn | awk '{print $1}'`
offset=$((LEN/1000))
leftovers=$((LEN%1000))

mkdir -p logs

cpus=32 # Max: 64 per each .meryl/ file writer
mem=110g
name=meryl_count
script=$meryl_count/meryl2_count.sh
partition=norm
walltime=1-0
path=`pwd`
log=logs/$name.%A_%a.log

if [ -e meryl_count.jid ]; then
    echo "Removing meryl_count.jid"
    cat meryl_count.jid
    rm meryl_count.jid
fi

for i in $(seq 0 $offset)
do
    args="$k $input_fofn $i"
    if [[ $i -eq $offset ]]; then
        arr_max=$leftovers
    else
        arr_max=1000
    fi
    echo "\
    sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --array=1-$arr_max --time=$walltime --error=$log --output=$log $script $args"
    sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --array=1-$arr_max --time=$walltime --error=$log --output=$log $script $args >> meryl_count.jid
done

# Wait for these jobs
WAIT="afterok:"

job_nums=`wc -l meryl_count.jid | awk '{print $1}'`

if [[ $job_nums -eq 1 ]]; then
   jid=`cat meryl_count.jid`
   WAIT=$WAIT$jid
else

for jid in $(cat meryl_count.jid)
do
    WAIT=$WAIT$jid","
done
fi

## Collect .meryl list
if [ -e  meryl_count.meryl.list ]; then
    echo "Removing meryl_count.meryl.list"
    cat meryl_count.meryl.list
    rm  meryl_count.meryl.list
fi

for line_num in $(seq 1 $LEN)
do
    input=`sed -n ${line_num}p $input_fofn`
    name=`echo $input | sed 's/.fastq.gz$//g' | sed 's/.fq.gz$//g' | sed 's/.fasta$//g' | sed 's/.fa$//g' | sed 's/.fasta.gz$//g' | sed 's/.fa.gz$//g'`
    name=`basename $name`
    echo "$name.$k.$line_num.meryl" >> meryl_count.meryl.list
done

cpus=48 # Max: 64 per each .meryl/ file writer
mem=72g
name=meryl_union_sum
script=$meryl_count/meryl2_union_sum.sh
args="$k meryl_count.meryl.list $out_prefix"
echo "\
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path --dependency=$WAIT --time=$walltime --error=$log --output=$log $script $args"
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path --dependency=$WAIT --time=$walltime --error=$log --output=$log $script $args > meryl_union_sum.jid

