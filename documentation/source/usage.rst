.. _usage:

{MOSTLY A PLACEHOLDER}

===========================
meryl2 COMMAND LINE OPTIONS
===========================

Global Options
--------------

Global options apply to every action.  Global options are expressed as
command line switches.  K-mer size, label size, memory usage, number of
threads, and simple sequence reductions are all global options.

.. code-block:: none
  :caption: meryl2 usage.
  :linenos:

  KMER and LABEL SIZE:
    -k K             Set kmer size to K (6 <= K <= 64).  Legacy format "k=K".
    -l L             Set label size to L bits (0 <= L <= 64).

  SIMPLE SEQUENCE REDUCTION:
    --compress       Homopolymer compress input sequence before forming kmers.
    --ssr ssr-spec   Apply 'ssr-spec' to input sequence before forming kmers.

  RESOURCE USAGE:
    -m M             Use up to M GB memory for counting.  Legacy format "memory=M".
    --memory M       Use up to M GB memory for counting.

    -t        T      Use T threads for counting and processing.
    --threads T      Use T threads for counting and processing.

  LOGGING:
    -V[V[V[...]]]    Increae verbosity by the length of this option.
    -Q               Be absolutely silent.
    -P               Show progress.
    -C               Show processing tree and stop.

  USAGE:
    -h               Display command line help.
    --help           Display command line help.
    help             Display command line help.

Obsolete forms of some of these are still allowed.  The kmer size can be set
with `k=<kmer-size>`, memory limits with `memory=<memory-in-gigabytes>`,
thread usage with `threads=<thread-count>`.

Obsolete option `-E` has been removed.  This used to estimate the size of an
imput that could be counted in some given memory size.

.. code-block:: none
  :caption: meryl2 debugging actions.
  :linenos:

  dumpIndex <meryl-database>
    Report the parameters used for storing data in this meryl database.

  dumpFile <meryl-database-file>
    Dump the raw data from a merylData file in a meryl database.


==================================
meryl2-lookup COMMAND LINE OPTIONS
==================================

.. code-block:: none
  :caption: meryl2-lookup usage.
  :linenos:

  usage: meryl2-lookup <report-type> \
           -sequence <input1.fasta> [<input2.fasta>] \
           -output   <output1>      [<output2>] \
           -mers     <input1.meryl> [<input2.meryl>] [...] [-estimate] \
           -labels   <input1name>   [<input2name>]   [...]

    Compare kmers in input sequences against kmers in input meryl databases.

    Input sequences (-sequence) can be FASTA or FASTQ, uncompressed, or
    compressed with gzip, xz, or bzip2.

    Report types:

    -bed:
       Generate a BED format file showing the location of kmers in
       any input database on each sequence in 'input1.fasta'.
       Each kmer is reported in a separate bed record.

    -bed-runs:
       Generate a BED format file showing the location of kmers in
       any input database on each sequence in 'input1.fasta'.
       Overlapping kmers are combined into a single bed record.

    -wig-count:
       Generate a WIGGLE format file showing the multiplicity of the
       kmer starting at each position in the sequence, if it exists in
       an input kmer database.

    -wig-depth:
       Generate a WIGGLE format file showing the number of kmers in
       any input database that cover each position in the sequence.

    -existence:
       Generate a tab-delimited line for each input sequence with the
       number of kmers in the sequence, in the database and common to both.

    -include:
    -exclude:
       Copy sequences from 'input1.fasta' (and 'input2.fasta') to the
       corresponding output file if the sequence has at least one kmer
       present (include) or no kmers present (exclude) in 'input1.meryl'.

  Run `meryl2-lookup <report-type> -help` for details on each method.

