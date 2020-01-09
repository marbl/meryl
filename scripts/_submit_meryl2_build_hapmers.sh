#!/bin/bash

# build=/path/to/codes/and/scripts/
build=$tools/meryl/scripts/build

if [ -z $1 ]; then
    echo "Usage: ./_submit_meryl2_build_hapmers.sh <k-size> <hap1.fofn> <hap2.fofn>"
    echo -e "\t<k-size>: kmer size k"
    echo -e "\t<hap1.fofn>: ls *.fastq.gz > hap1.fofn. accepts fasta, fastq, gzipped or not."
    echo -e "\t<hap2.fofn>: ls *.fastq.gz > hap2.fofn. accepts fasta, fastq, gzipped or not."
    echo -e "\tIntermediate kmer dbs will be named after hap1 and hap2"
    exit -1
fi

k=$1
hap1_fofn=$2
hap2_fofn=$3

hap1=${hap1_fofn/.fofn/}
hap2=${hap2_fofn/.fofn/}

if [ -e build.jid ]; then
    echo "Removing build.jid"
    cat build.jid
    rm build.jid
fi

# Build for hap1
for hap in "$hap1" "$hap2"
do
	input_fofn=$hap.fofn
	echo $input_fofn
	
	LEN=`wc -l $input_fofn | awk '{print $1}'`
	offset=$((LEN/1000))
	leftovers=$((LEN%1000))

	mkdir -p logs

	cpus=32 # Max: 64 per each .meryl/ file writer
	mem=48g
	name=build
	script=$build/meryl2_count.sh
	partition=norm
	walltime=1-0
	path=`pwd`
	log=logs/$name.%A_%a.log

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
	    sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --array=1-$arr_max --time=$walltime --error=$log --output=$log $script $args >> build.jid
	done


	# Wait for these jobs
	WAIT="afterok:"

	job_nums=`wc -l build.jid | awk '{print $1}'`

	if [[ $job_nums -eq 1 ]]; then
	   jid=`cat build.jid`
	   WAIT=$WAIT$jid
	else

	for jid in $(cat build.jid)
	do
	    WAIT=$WAIT$jid","
	done
	fi

	## Collect .meryl list
	if [ -e  build.$hap.meryl.list ]; then
	    echo "Removing build.meryl.list"
	    cat build.$hap.meryl.list
	    rm  build.$hap.meryl.list
	fi

	for line_num in $(seq 1 $LEN)
	do
	    input=`sed -n ${line_num}p $input_fofn`
	    name=`echo $input | sed 's/.fastq.gz$//g' | sed 's/.fq.gz$//g' | sed 's/.fasta$//g' | sed 's/.fa$//g' | sed 's/.fasta.gz$//g' | sed 's/.fa.gz$//g'`
	    name=`basename $name`
	    echo "$name.$k.$line_num.meryl" >> build.$hap.meryl.list
	done

	cpus=48 # Max: 64 per each .meryl/ file writer
	mem=72g
	name=meryl_union_sum
	script=$build/meryl2_union_sum.sh
	args="$k build.$hap.meryl.list $hap"
	log=logs/$name.%A_%a.log

	echo "\
	sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path --dependency=$WAIT --time=$walltime --error=$log --output=$log $script $args"
	sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path --dependency=$WAIT --time=$walltime --error=$log --output=$log $script $args > meryl_union_sum.jid
done

# Wait for these jobs
WAIT="afterok:"

job_id=meryl_union_sum.jid
job_nums=`wc -l $job_id | awk '{print $1}'`

if [[ $job_nums -eq 1 ]]; then
   jid=`cat $job_jid`
   WAIT=$WAIT$jid
else

for jid in $(cat build.jid)
do
    WAIT=$WAIT$jid","
done
fi

cpus=48 # Max: 64 per each .meryl/ file writer
mem=20g
name=meryl_diff
script=$build/meryl2_diff.sh
args="$hap1.meryl $hap2.meryl"
log=logs/$name.%A_%a.log
partition=quick

echo "\
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path --dependency=$WAIT --time=$walltime --error=$log --output=$log $script $args"
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path --dependency=$WAIT --time=$walltime --error=$log --output=$log $script $args > meryl_diff.jid


