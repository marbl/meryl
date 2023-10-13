.. note::
  The original meryl used order ACTG for a reason that turned out to be
  incorrect.  It was believed that complementing a binary sequence would be
  easier in that order, but it is just as easy in the normal order.  The
  revised order does have the appealing property that GC content can be
  computed by counting the number of low-order bits set in each base, where the
  more standard ACGT order requires additional operations.

  In the revised (ACTG) order, flipping the first bit of each two-bit encoded base will
  complement the base.  This can be done with an exclusive or against b101010.

  .. code-block:: none

    A  C  T  G
    00 01 10 11
    v| v| v| v| -- flip first bit to complement
    10 11 00 01
     T  G  A  C

    compl = kmer ^ 0b101010
    NumGC = popcount(kmer & b010101)

  In the usual (ACGT) order, complementation can be accomplished by flipping
  every bit; an exclusive-or against b111111.  GC content can be computed by
  counting bits also, but we need to count the number of two-bit encoded
  bases where the first and second bits differ.  This requires shifting the
  encoded k-mer one place, comparing the -- now overlapping -- adjacent bits
  with XOR, and finally counting the number of set bits.

  .. code-block:: none

    A  C  G  T
    00 01 10 11
    vv vv vv vv -- flip every bit to complement
    11 10 01 00
     T  G  C  A

    compl = kmer ^ 0b111111
    NumGC = popcount( (kmer>>1 ^ kmer) & 0b010101 )

  It is yet to be decided if meryl2 will also use the same order (maintaining
  compatibility with meryl1) or if it will use the more typical order
  (maintaining compatibility with the rest of the world)..
