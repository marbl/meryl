# Scripts

Here, bash scripts are provided to provide guidelines for the following purposes:

1.	Build k-mer dbs in batches
2.	Copy-number comparison (spectra-cn and spectra-hap)
3.	Haplotype specific mers for trio binning and making blob-plots

All scripts are tested on slurm environment for diploid genomes up to 3.3Gb.
Modify the `_submit_*` scripts to match your system.
`cpus` and `mem` are suggestive.


# 1. Build k-mer dbs

`_submit_meryl2_build.sh` is doing the following 1. Build and 2. Merge steps automatically.
`_submit_meryl2_build_10x.sh` trims the 10X barcodes from illumina reads, and runs the 1. Build 2. Merge steps.
The barcode trimming step uses `scaff_reads` from [scaff10x](https://github.com/wtsi-hpag/Scaff10X).

## 1. Build

Count k-mers in <in.sequence>, could be fastq/fasta file, gzipped or not.

```
meryl count k=$k_size output $output.meryl <in.sequence>
```

Will generate canonical k-mer counts.

`Threads`: By default, meryl uses the maximum number of cores available, up to 64.
This can be limited by defining `threads=$cores`.
If you want to see results faster, split the <in.sequence> file to multiple pieces and count each kmers and merge at the end.

`Memory`: Meryl automatically detects available memory and uses that.
To prevent overhead for allocating too much memory, use something like `memory=48` which will use up to 48G.
When the k-mer space is too big and requires more memory, meryl automatically splits countig jobs in batches and merges at the end.

Run `meryl`  with no args to see more details.

## 2. Merge

If counted from multiple fastq files, merge them at the end.
```
ulimit -Sn 32000
meryl memory=$mem union-sum output $output_meryl $meryl_dbs
```
meryl will try to make 64 x 2 files per meryl dbs to merge.

## 3. Histogram
What’s in my meryl db?
This will give you a short summary of the k-mers, unique / distinct / total.

```
meryl statistics $output_meryl | head
```

Without the head, the histogram is printed after the summary, which could go very long.

Get the histogram: This gives a traditional 2-column histogram, ready for genomescope or any other tool.

```
meryl histogram $output_meryl > $output_meryl.hist
```

Ploidy estimation including erroneous k-mer boundaries can be made by this script:
```
Java -jar -Xmx256m /path/to/this/meryl-accessories/scripts/kmerHistToPloidyDepth.jar $output_meryl.hist > $output_meryl.hist.ploidy
```

The output is a tdf file, with
```
Ploidy	Depth	Boundary
```

Depth is the depth where we see a ‘peak’ for ploidy > 0.
Boundaries are the dips in the histogram.

# 2. Copy-number comparison (spectra-cn)
This script is the meryl version spectra-cn plot, introduced and implemented by [KAT](https://github.com/TGAC/KAT).
The scripts here provide a way to draw unstacked, stacked histograms by copy numbers or specific to assemblies. No dependency to KAT.

## 1. Get kmers from the read set and assembly (assemblies)
Use the above tricks for building the meryl db for the reads.
For each assembly,
```
meryl count memory=8 output $meryl_db <asm.fasta>
```
Should work fine with most of the genomes. This is using 8G of memory (or lower), taking a few seconds to build the $meryl_db.

## 2. Get assembly-specific and shared kmers
When working with 2 genome assemblies,
such as one primary and alternate haplotype assemblies from diploid genomes,
or two haplotype assemblies generated with trio-binning,
the interest will be to see how many kmers are over / under represented compared to the expected copy numbers.

Here, we will get assembly-specific (private to one of the assemblies) kmers and the shared kmers in separate temporary meryl dbs for convenience.

hap1 and hap2 can be set as

>hap1=paternal

>hap2=maternal

Or

>hap1=pri

>hap2=alt


```
meryl difference output $name.${hap2}_only.meryl $name.${hap2}.meryl $name.${hap1}.meryl
meryl difference output $name.${hap1}_only.meryl $name.${hap1}.meryl $name.${hap2}.meryl
meryl intersect output $name.${hap1}_shrd.meryl $name.${hap1}.meryl $name.${hap2}.meryl
meryl intersect output $name.${hap2}_shrd.meryl $name.${hap2}.meryl $name.${hap1}.meryl
```

This generates 4 meryl dbs; and we need one more that sums up all the kmer counts:
```
meryl union-sum output $name.${hap1}_${hap2}_shrd.meryl $name.${hap1}_only.meryl $name.${hap1}_shrd.meryl $name.${hap2}_only.meryl $name.${hap2}_shrd.meryl
```

## 3. Intersect with the read meryl db
Now, intersect the above kmer set of interest with the kmers from the read set, and get the kmer counts.

* 0-counts in the read set; i.e. only seen in the assembly.
Most likely assembly errors or kmers difficult to get from illumina reads

```
meryl difference output $name.only.meryl $name.${hap1}_${hap2}_shrd.meryl $name.k21.meryl
```

This kmer set is a useful indicator of bp errors.
The number of k-mers are associate-able with bp, so roughly the number of distinct kmers in this set are representing number of erroneous bases in the assembly.
It is difficult to distinguish kmers not sequence-able with illumina reads though, so it might be useful for comparing overal bp error rate in different assemblies from the same individual.

```
meryl statistics $name.only.meryl | head
```

Switch the two inputs now to get:

* 0-counts in the asm; i.e. only seen in the illumina reads
```
meryl difference output $name.k21.0.meryl $name.k21.meryl $name.${hap1}_${hap2}_shrd.meryl
```

* 1~4 copy kmers
```
for i in $(seq 1 4)
do
    meryl intersect output $name.k21.$i.meryl $name.k21.meryl [ equal-to $i $name.${hap1}_${hap2}_shrd.meryl ]
done
```

* \>4 copy kmers
```
meryl intersect output $name.k21.gt$i.meryl $name.k21.meryl [ greater-than $i $name.${hap1}_${hap2}_shrd.meryl ]
```

## 4. Spectra-cn plots

The final histogram file can be obtained with:
```
echo -e "Copies\tkmer_multiplicity\tCount" > $name.spectra-cn.hist
meryl histogram $name.k21.0.meryl | awk '{print "illumina-only\t"$0}' >> $name.spectra-cn.hist
for i in $(seq 1 4)
do
    meryl histogram $name.k21.$i.meryl | awk -v cn=$i '{print cn"\t"$0}' >> $name.spectra-cn.hist
done
meryl histogram $name.k21.gt$i.meryl | awk -v cn=">$i" '{print cn"\t"$0}' >> $name.spectra-cn.hist
```

## 5. Spectra-hap plots
Since we are interested in comparing the portion of kmers each haplotype (or pseudo haplotype) assemblies have,
let’s make some plots with kmers grouped by specific / shared kmers.

```
meryl intersect output $name.k21.$hap1.meryl $name.k21.meryl $name.${hap1}_only.meryl
meryl intersect output $name.k21.$hap2.meryl $name.k21.meryl $name.${hap2}_only.meryl
meryl intersect output $name.k21.shrd.meryl $name.k21.meryl $name.${hap1}_shrd.meryl
```

Again, get the histograms by groups.
```
echo -e "kmer\tkmer_multiplicity\tCount" > $name.spectra-hap.hist
meryl histogram $name.k21.0.meryl | awk '{print "illumina-only\t"$0}' >> $name.spectra-hap.hist
meryl histogram $name.k21.$hap1.meryl | awk -v hap="$hap1-only" '{print hap"\t"$0}' >> $name.spectra-hap.hist
meryl histogram $name.k21.$hap2.meryl | awk -v hap="$hap2-only" '{print hap"\t"$0}' >> $name.spectra-hap.hist
meryl histogram $name.k21.shrd.meryl | awk -v hap="shared" '{print hap"\t"$0}' >> $name.spectra-hap.hist
```

`spectracn.R` is an example to generate plots for `$name.spectra-hap.hist`. More instructions will be posted later.


# 3. Haplotype specific mers for trio binning and making blob-plots

`_submit_meryl2_build_hapmers.sh` does the binning.
For making blob-plots, use the temporary scripts uploaded to [here]( https://github.com/arangrhie/TrioBinning/blob/master/hap_kmer_plot.R).
