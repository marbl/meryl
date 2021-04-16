#!/bin/sh
#
#  Test the various reports from meryl-lookup.  Length 31 has been manually verified.
#

length=393
length=31

if [ $length = 31 ] ; then
  c0=910bbe43b6e8b3ad3002efde90811ad6
  c1=910bbe43b6e8b3ad3002efde90811ad6
  c2=f9e20c9897ef5b3923a9a9461f5f46ab
  c3=c1203d8370df8206d01593364cc71629
  c4=f9e20c9897ef5b3923a9a9461f5f46ab
  c5=c1203d8370df8206d01593364cc71629
  c6=78fd9671127c21ab404ac9f786a2a70f
  c7=78fd9671127c21ab404ac9f786a2a70f
  c8=6bcb11732b1e235d593bbfe78e693f62
  c9=491f863345c5caee7ad07b6a9ba38e33
fi

if [ $length = 393 ] ; then
  c0=453d04e4b06289d53c77726289be581c
  c1=453d04e4b06289d53c77726289be581c
  c2=34a4dd6ab26d771637e4230163e0a1e0
  c3=ef17f64379a1ec6791c9df8a53f6f086
  c4=4e15fe5f80bee2d9de9bdc5e2da58766
  c5=ef17f64379a1ec6791c9df8a53f6f086
  c6=bf34073b01131a58e9cea145c6549d4f
  c7=bf34073b01131a58e9cea145c6549d4f
  c8=f23897fe79b213ba70809a13263a9f7a
  c9=5c4ade5bb958ee0daf3bcfeeae37698c
fi

#  Create base strings.
#
#  Results are a little weird if there is a spurious kmer spanning block AD or DC.

if [ ! -e A.fasta ] ; then  seqrequester generate -seed 315 -min $length -max $length -sequences 1 -ident A > seqs-A.fasta 2> /dev/null ; fi
if [ ! -e B.fasta ] ; then  seqrequester generate -seed 316 -min $length -max $length -sequences 1 -ident B > seqs-B.fasta 2> /dev/null ; fi
if [ ! -e C.fasta ] ; then  seqrequester generate -seed 211 -min $length -max $length -sequences 1 -ident C > seqs-C.fasta 2> /dev/null ; fi
if [ ! -e D.fasta ] ; then  seqrequester generate -seed 111 -min $length -max $length -sequences 1 -ident D > seqs-D.fasta 2> /dev/null ; fi

#  Prepare genome sequences.

rm -f  genome.fasta
rm -fr genome.meryl

cat       >> genome.fasta seqs-A.fasta
tail -n 1 >> genome.fasta seqs-B.fasta
tail -n 1 >> genome.fasta seqs-C.fasta

cat       >> genome.fasta seqs-A.fasta
tail -n 1 >> genome.fasta seqs-C.fasta

cat       >> genome.fasta seqs-B.fasta

cat       >> genome.fasta seqs-B.fasta

meryl k=22 memory=1g count genome.fasta output genome.meryl 2> /dev/null

#  Prepare query sequences.

rm -f seq-fwds.fasta

echo      >> seq-fwds.fasta ">ABC"
tail -n 1 >> seq-fwds.fasta seqs-A.fasta
tail -n 1 >> seq-fwds.fasta seqs-B.fasta
tail -n 1 >> seq-fwds.fasta seqs-C.fasta

echo      >> seq-fwds.fasta ">AC"
tail -n 1 >> seq-fwds.fasta seqs-A.fasta
tail -n 1 >> seq-fwds.fasta seqs-C.fasta

echo      >> seq-fwds.fasta ">ADC"
tail -n 1 >> seq-fwds.fasta seqs-A.fasta
tail -n 1 >> seq-fwds.fasta seqs-D.fasta
tail -n 1 >> seq-fwds.fasta seqs-C.fasta

echo      >> seq-fwds.fasta ">D"
tail -n 1 >> seq-fwds.fasta seqs-D.fasta

seqrequester extract -rc seq-fwds.fasta > seq-revs.fasta 2> /dev/null

#
#  Test!
#

fail=0

#  BED
meryl-lookup -bed -sequence seq-fwds.fasta -mers genome.meryl -labels genome > result-fwds.bed 2> /dev/null
meryl-lookup -bed -sequence seq-revs.fasta -mers genome.meryl -labels genome > result-revs.bed 2> /dev/null

if [ `cat result-fwds.bed | md5` != $c0 ] ; then
  fail=1
  echo "meryl-lookup -bed -sequence seq-fwds.fasta -mers genome.meryl -labels genome > result-fwds.bed"
  echo " - FAIL!"
fi

if [ `cat result-revs.bed | md5` != $c1 ] ; then
  fail=1
  echo "meryl-lookup -bed -sequence seq-revs.fasta -mers genome.meryl -labels genome > result-revs.bed"
  echo " - FAIL!"
fi

#  WIGGLE count and depth
meryl-lookup -wig-count -sequence seq-fwds.fasta -mers genome.meryl > result-fwds-count.wig 2> /dev/null
meryl-lookup -wig-depth -sequence seq-fwds.fasta -mers genome.meryl > result-fwds-depth.wig 2> /dev/null

meryl-lookup -wig-count -sequence seq-revs.fasta -mers genome.meryl > result-revs-count.wig 2> /dev/null
meryl-lookup -wig-depth -sequence seq-revs.fasta -mers genome.meryl > result-revs-depth.wig 2> /dev/null

if [ `cat result-fwds-count.wig | md5` != $c2 ] ; then
  fail=1
  echo "meryl-lookup -wig-count -sequence seq-fwds.fasta -mers genome.meryl > result-fwds-count.wig"
  echo " - FAIL!"
fi

if [ `cat result-fwds-depth.wig | md5` != $c3 ] ; then
  fail=1
  echo "meryl-lookup -wig-depth -sequence seq-fwds.fasta -mers genome.meryl > result-fwds-depth.wig"
  echo " - FAIL!"
fi

if [ `cat result-revs-count.wig | md5` != $c4 ] ; then
  fail=1
  echo "meryl-lookup -wig-count -sequence seq-revs.fasta -mers genome.meryl > result-revs-count.wig"
  echo " - FAIL!"
fi

if [ `cat result-revs-depth.wig | md5` != $c5 ] ; then
  fail=1
  echo "meryl-lookup -wig-depth -sequence seq-revs.fasta -mers genome.meryl > result-revs-depth.wig"
  echo " - FAIL!"
fi

#  Existence
meryl-lookup -existence -sequence seq-fwds.fasta -mers genome.meryl > result-fwds.existence 2> /dev/null
meryl-lookup -existence -sequence seq-revs.fasta -mers genome.meryl > result-revs.existence 2> /dev/null

if [ `cat result-fwds.existence | md5` != $c6 ] ; then
  fail=1
  echo "meryl-lookup -existence -sequence seq-fwds.fasta -mers genome.meryl > result-fwds.existence 2> /dev/null"
  echo " - FAIL!"
fi

if [ `cat result-revs.existence | md5` != $c7 ] ; then
  echo "meryl-lookup -existence -sequence seq-revs.fasta -mers genome.meryl > result-revs.existence 2> /dev/null"
  echo " - FAIL!"
fi

#  Include/Exclude
meryl-lookup -include -sequence seq-fwds.fasta -mers genome.meryl -output result-include.fasta 2> /dev/null
meryl-lookup -exclude -sequence seq-fwds.fasta -mers genome.meryl -output result-exclude.fasta 2> /dev/null

if [ `cat result-include.fasta | md5` != $c8 ] ; then
  fail=1
  echo "meryl-lookup -include -sequence seq-fwds.fasta -mers genome.meryl -output result-include.fasta 2> /dev/null"
  echo " - FAIL!"
fi

if [ `cat result-exclude.fasta | md5` != $c9 ] ; then
  fail=1
  echo "meryl-lookup -exclude -sequence seq-fwds.fasta -mers genome.meryl -output result-exclude.fasta 2> /dev/null"
  echo " - FAIL!"
fi

#  Include/Exclude of Illumina reads


#  Helpful for regenerating the checksums above.

if [ $fail -gt 0 ] ; then
  echo c0=`cat result-fwds.bed       | md5`
  echo c1=`cat result-revs.bed       | md5`
  echo c2=`cat result-fwds-count.wig | md5`
  echo c3=`cat result-fwds-depth.wig | md5`
  echo c4=`cat result-revs-count.wig | md5`
  echo c5=`cat result-revs-depth.wig | md5`
  echo c6=`cat result-fwds.existence | md5`
  echo c7=`cat result-revs.existence | md5`
  echo c8=`cat result-include.fasta  | md5`
  echo c9=`cat result-exclude.fasta  | md5`
fi

#  Bye.

exit

