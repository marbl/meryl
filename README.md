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

Here, bash scripts are provided to provide guidelines for the following purposes:

1.	Build k-mer dbs in batches
2.	Copy-number comparison (spectra-cn and spectra-hap)
3.	Haplotype specific mers for trio binning and making blob-plots

See [scripts](https://github.com/marbl/meryl/tree/master/scripts) for more details.
