#!/bin/bash

if [ -z $1 ]; then
	echo "Usage: ./meryl2_filt.sh <in.meryl>"
	exit -1
fi

db=$1
db=${db/.meryl/}

echo "Generate $db.hist"
meryl histogram $db.meryl > $db.hist

echo "
java -jar -Xmx1g $tools/TrioBinning/meryl_count/kmerHistToPloidyDepth.jar $db.hist
"
java -jar -Xmx1g $tools/TrioBinning/meryl_count/kmerHistToPloidyDepth.jar $db.hist > $db.hist.ploidy

cat $db.hist.ploidy

filt=`sed -n 2p $db.hist.ploidy | awk '{print $NF}'`

echo "Filter out kmers <= $filt"

echo "
meryl greater-than $filt output $db.gt$filt.meryl $db.meryl
"
meryl greater-than $filt output $db.gt$filt.meryl $db.meryl


