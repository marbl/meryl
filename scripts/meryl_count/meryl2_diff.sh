#!/bin/bash

if [[ -z $1 ]] || [[ -z $2 ]]; then
	echo "Usage: meryl2_diff.sh <db1.meryl> <db2.meryl>"
	exit -1
fi

db1=$1
db2=$2

db1=${db1/.meryl/}
db2=${db2/.meryl/}

if [ ! -d ${db1}_only.meryl ]; then
	echo "\
	meryl difference output ${db1}_only.meryl $db1.meryl $db2.meryl"
	meryl difference output ${db1}_only.meryl $db1.meryl $db2.meryl
	echo
fi

if [ ! -d ${db2}_only.meryl ]; then
	echo "\
	meryl difference output ${db2}_only.meryl $db2.meryl $db1.meryl"
	meryl difference output ${db2}_only.meryl $db2.meryl $db1.meryl
	echo
fi

echo "
sh $tools/TrioBinning/meryl_count/meryl2_filt.sh ${db1}_only.meryl"
sh $tools/TrioBinning/meryl_count/meryl2_filt.sh ${db1}_only.meryl

echo "
sh $tools/TrioBinning/meryl_count/meryl2_filt.sh ${db2}_only.meryl"
sh $tools/TrioBinning/meryl_count/meryl2_filt.sh ${db2}_only.meryl

