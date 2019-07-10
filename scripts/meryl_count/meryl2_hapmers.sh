#!/bin/sh

if [[ -z $1 ]] || [[ -z $2 ]] || [[ -z $3 ]] || [[ -z $4 ]]; then
    echo "Usage: ./meryl_count_asm.sh <hapA.mcdat> <hapB.mcdat> <asmA.fasta> <asmB.fasta>"
    exit -1
fi

hap_a=$1	# meryl hap_a.mcdat
hap_b=$2	# meryl hap_b.mcdat

a_fa=$3		# hap_a.fasta
b_fa=$4		# hap_b.fasta

hap_a=${hap_a/.mcdat/}
hap_b=${hap_b/.mcdat/}

a_name=`echo $a_fa | sed 's/.fasta$//g' | sed 's/.fa$//g'`
b_name=`echo $b_fa | sed 's/.fasta$//g' | sed 's/.fa$//g'`

echo "Generate qr dbs"
if ! [ -s $hap_a.qr ]; then
echo "\
simple-dump -m 21 -s $hap_a -e $hap_a.qr"
simple-dump -m 21 -s $hap_a -e $hap_a.qr
echo
fi

if ! [ -s $hap_b.qr ]; then
echo "\
simple-dump -m 21 -s $hap_b -e $hap_b.qr"
simple-dump -m 21 -s $hap_b -e $hap_b.qr
echo
fi

echo "Start counting"
if ! [ -s ${a_name}"_"$hap_a.counts ]; then
echo "\
simple-dump -m 21 -e $hap_a.qr -f $a_fa | awk -v name=$a_name \'{print name\\t\$0}\' > ${a_name}_$hap_a.counts"
simple-dump -m 21 -e $hap_a.qr -f $a_fa | awk -v name=$a_name '{print name"\t"$0}'   > ${a_name}"_"$hap_a.counts
echo
fi

if ! [ -s ${b_name}"_"$hap_a.counts ]; then
echo "\
simple-dump -m 21 -e $hap_a.qr -f $b_fa | awk -v name=$b_name \'{print name\\t\$0}\' > ${b_name}_$hap_a.counts"
simple-dump -m 21 -e $hap_a.qr -f $b_fa | awk -v name=$b_name '{print name"\t"$0}'   > ${b_name}"_"$hap_a.counts
echo
fi

if ! [ -s ${a_name}"_"$hap_b.counts ]; then
echo "\
simple-dump -m 21 -e $hap_b.qr -f $a_fa | awk -v name=$a_name \'{print name\\t\$0}\' > ${a_name}"_"$hap_b.counts"
simple-dump -m 21 -e $hap_b.qr -f $a_fa | awk -v name=$a_name '{print name"\t"$0}'   > ${a_name}"_"$hap_b.counts
echo
fi

if ! [ -s ${b_name}"_"$hap_b.counts ]; then
echo "\
simple-dump -m 21 -e $hap_b.qr -f $b_fa | awk -v name=$b_name \'{print name\\t\$0}\' > ${b_name}"_"$hap_b.counts"
simple-dump -m 21 -e $hap_b.qr -f $b_fa | awk -v name=$b_name '{print name"\t"$0}'   > ${b_name}"_"$hap_b.counts
echo
fi

echo "Merge"
echo "\
paste ${a_name}_$hap_a.counts ${a_name}_$hap_b.counts > $a_name.counts"
paste ${a_name}"_"$hap_a.counts ${a_name}"_"$hap_b.counts > $a_name.counts

echo "\
paste ${b_name}_$hap_a.counts ${b_name}_$hap_b.counts > $b_name.counts"
paste ${b_name}"_"$hap_a.counts ${b_name}"_"$hap_b.counts > $b_name.counts

echo "\
echo -e "Assembly\tContig\t$hap_a\t$hap_b\tTotal" > hapmers.counts"
echo -e "Assembly\tContig\t$hap_a\t$hap_b\tTotal" > hapmers.counts

echo "\
cat $a_name.counts $b_name.counts | awk '{print \$1\t\$2\t\$5\t\$NF\t\$3}' >> hapmers.counts"
cat $a_name.counts $b_name.counts | awk '{print $1"\t"$2"\t"$5"\t"$NF"\t"$3}' >> hapmers.counts
