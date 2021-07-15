.. _reference:

===================
Parameter Reference
===================



Databases
=========

A **meryl database** stores kmers, values and labels in a
compressed binary format.

The **value** of a kmer is typically the number of times -- its **count** --
it occurs in some input sequence.

The **label** of a kmer is a 64-bit binary string assigned to each kmer by
the user.  Labels are envisioned to be most useful as a collection of yes/no
or true/false flags, for example, "kmer occurred in file 1".  **Labels are
not fully implemented yet and will appear in a future version.**

The database is stored as 64 independent files, each file storing some subset
of kmers.  Each file can be processed independently, allowing meryl to use up
to 64 threads.

Details: Each kmer is stored by breaking the binary represntation into three
pieces: a file prefix, a block prefix, and a suffix.  A k-mer needs 2*k bits
to represent it.  The file prefix is 6 bits wide, representing one of the 64
possible files in a database.  Inside a file, kmers are stored in blocks,
each kmer in a block will have the same block prefix.  The suffix data for a
block is stored using Elias-Fano encoding (CITE) where each suffix is split
into two pieces.  The first piece is encoded as a unary offset from the last
piece, and the second piece is a N-log2(N)+1 binary value.  At present,
values are stored as plain 32-bit integers, stored immediately after the kmer
data.


Actions
=======

Meryl processing is built around *actions*.  An action reads kmers (or sequence) from one or many inputs
and writes kmers to a single output.  There are four basic actions:

  :count:       create a meryl database from FASTA or FASTQ inputs
  :filter:      select kmers based on their value of label from a single input
  :modify:      change kmer values or labels in a single input
  :present-in:  combine multiple inputs into one output

Numerous aliases exist for common operations.


Count
-----

The **count** action reads sequence from input files and converts to a list of
kmers and the number of times each kmer occurs in the input sequence.

Input sequences can be in either FASTA, FASTQ, raw bases, or if compiled with
Canu support, in a Canu seqStore database.  Sequence files can be gzip, bzip2
or xz compressed.

An output database must be supplied to all **count** actions.  Kmers are written
to both the output database and provided as input to other actions.

By default, **count** will count and store kmers in their canonical form.
Actions **count-forward** and **count-reverse** can be used to instead process
kmers as they occur in the sequence (**count-forward**) or as the occur in the
reverse-complement sequence (**count-reverse**).

Unless derived from another input database, **count** actions must have a kmer size
supplied on the command line with the **-k** option (or via the deprecated inline option **k=<kmer-size>**).

Counting is resource intense.  By default, meryl will use all available
resources on the host machine, or the limits as supplied by a grid manager
such as Slurm, PBS or SGE.  Lower limits can be specified with options **-m**
and **-t** (or with deprecated inline options
**memory=<max-memory-to-use-in-gigabytes>** and
**threads=<max-number-of-cpus-to-use>**).

Homopolymer-compressed kmers can be counted by supplying the **compress** flag
to the **count** action.

.. code-block:: shell

  % meryl count output kmers.meryl compress sequence.fasta

Counting Algorithms
-------------------

Two algorithms are used for counting kmers.  The algorithm that is expected
to use the least memory is used.  The choice depends mostly on kmer size, but
a little bit on the size of the input sequences.

Small Kmers
~~~~~~~~~~~
For small kmers (at most 16) meryl counts kmers directly by associating an
integer count with each possible kmer.  This has the benefit of being simple
and uses a constant amount of memory regardless of the size of the input, but
quickly exhausts memory for even moderate kmer sizes.

There are 4\ :sup:`k` kmers of size k; for k=16, there are 4,294,967,296
possible kmers.  At a minimum, counting 16-mers with this method will use 8
GB of memory, independent of input size.  As counts increase, memory is added
in chunks of 512 MB.

Details: Each integer counter is initially a 16-bit value.  Once any count
exceeds 2\ :sup:`16` = 65,535 another bit is added to all value, resulting
in 17-bit values for every kmer.  Once any count then exceeds 2\ :sup:`17` =
131,072, another bit is added, and so on.

This method uses only a single thread to read the input sequence and
increment counters in the array, but multiple threads can be used to generate
the output database.

Large Kmers
~~~~~~~~~~~

For large kmers (generally, larger than 16), meryl converts an input sequence
to a list of all kmers, duplicates included, in it.  The list of kmers is
sorted.  It is now trivial to scan the list, counting how many times each
kmer is present, and immediately writing the kmer and its value to the output
database.

If all kmers in an input sequence do not fit in memory, a partial result is
written to disk.  After all input sequences have been processed, the partial
results are combined into a single output database.  In practice, this method
requires several gigabytes of memory to minimize the overhead of writing and
merging partial results.

This algorithm can use multiple threads for every stage.

Details: Each kmer is split into a prefix and a suffix.  The prefix is used
to select a list to which the suffix is added.  A trade off is made between a
small prefix (resulting in few lists that store large suffixes) and a large
prefix (resulting in many lists where the overhead of each list could use
more space than the lists themselves).  When the (approximate) size of all
lists exceeds a user-supplied threshold, each list is sorted, the suffixes
are counted, and output to an intermediate database.  After all kmers are
processed, the intermediate databases are merged into one.

Filter
------

Modify
------

Present-In
----------
