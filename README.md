# Meryl

This is 'meryl', a near total rewrite of 'meryl' that appeared in both
[project kmer](http://kmer.sourceforge.net/) and
[Celera Assembler](http://wgs-assembler.sourceforge.net/).

*IMPORTANT*: Get the latest meryl code from this repo. This is not compatible with old meryl dbs built from canu 1.8 or earlier. The new meryl is significantly faster than the previous version.

meryl dbs are no longer in `.mcdat` and `.mcidx` file format. Meryl db is now designed as a DIRECTORY, containing 64 binaries + 64 indexes (128 files).

### Dependency
* gcc 4.8 or higher

### Installation
```
git clone https://github.com/marbl/meryl.git
cd meryl/src
make -j 24
export PATH=/path/to/meryl/bin:$PATH
```

# Sequence

This is 'sequence', a utility for working with sequence files.


# Scripts

Here, bash scripts are provided to give a guideline for the following use-cases:

1.	[Build k-mer dbs in batches](https://github.com/marbl/meryl/tree/master/scripts#1-build-k-mer-dbs)
2.	[Copy-number comparison (spectra-cn and spectra-hap)](https://github.com/marbl/meryl/tree/master/scripts#2-copy-number-comparison-spectra-cn)
3.	[Haplotype specific mers for trio binning and making blob-plots](https://github.com/marbl/meryl/tree/master/scripts#3-haplotype-specific-mers-for-trio-binning-and-making-blob-plots)
4.  [Include or exclude read pairs having specific kmers in one of the reads](https://github.com/marbl/meryl/blob/master/scripts/README.md#4-include-or-exclude-read-pairs-having-specific-kmers-in-one-of-the-reads)

See [scripts](https://github.com/marbl/meryl/tree/master/scripts) for more details.
