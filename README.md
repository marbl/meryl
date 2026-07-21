# Meryl

This is 'meryl', a near total rewrite of 'meryl' that appeared in both
[project kmer](http://kmer.sourceforge.net/) and
[Celera Assembler](http://wgs-assembler.sourceforge.net/).

*IMPORTANT*: Get the latest meryl code from this repo. This is not compatible with old meryl dbs built from canu 1.8 or earlier. The new meryl is significantly faster than the previous version.

### Change log (each with lots of bug fixes)

* v1.4.2: Direct kmer counting from `bam` / `cram`

* v1.4.1: One critical bug fix (bitArray) found right after v1.4 release

* v1.4: `meryl compress` kmers, experimental meryl2

* v1.3: `meryl-lookup` to lookup kmers

* v1.2: `meryl-analyze` to analyze kmers, `meryl-simple` for simpler kmer counts

* v1.1: Brand-new meryl! `.mcdat` and `.mcidx` file formats are no longer used. Meryl db is designed as a DIRECTORY, containing 64 binaries + 64 indexes (128 files)

### Dependency
* gcc 10.2+ or clang 12+

### Installation

* Release version: download a stable [release](https://github.com/marbl/meryl/releases/latest) version

  | Platform            | Filename                              |
  | ------------------- | ------------------------------------- |
  | Linux x86_64        | `meryl-1.4.2.Linux-amd64.tar.xz`      |
  | Linux arm64         | `meryl-1.4.2.Linux-arm64.tar.xz`      |
  | macOS Apple Silicon | `meryl-1.4.2.Darwin-arm64.tar.xz`     |
  | macOS Intel         | `meryl-1.4.2.Darwin-amd64.tar.xz`     |

  ```shell
  # Example for Linux-amd64
  wget https://github.com/marbl/meryl/releases/download/v1.4.2/meryl-1.4.2.Linux-amd64.tar.xz
  tar -xJf meryl-1.4.2.Linux-amd64.tar.xz
  export PATH=/path/to/meryl-1.4.2/bin:$PATH
  ```

* MacOS: after unpacking, run `xattr -dr com.apple.quarantine meryl-1.4.2.Darwin-*/`.
  On MDM-managed Macs, Gatekeeper might still block the binaries — build from source instead.

* Experimental tip (use git 2.25.1 or higher):
  ```shell
  git clone --recursive https://github.com/marbl/meryl.git
  
  # build
  cd meryl/src
  make -j 24
  export PATH=/path/to/meryl/build/bin:$PATH
  ```

## Evaluate assemblies with k-mers and more

See [Merqury](https://github.com/marbl/merqury).

## Citing Meryl

We didn't want to suffer the world with yet another k-mer counting paper, so we stuffed meryl into the merqury methods:
>Rhie, A., Walenz, B.P., Koren, S. et al. Merqury: reference-free quality, completeness, and phasing assessment for genome assemblies. Genome Biol 21, 245 (2020). https://doi.org/10.1186/s13059-020-02134-9
