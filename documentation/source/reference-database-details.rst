.. sidebar:: Database Implementation Details

  Each k-mer is stored by breaking the binary represntation into three
  pieces: a file prefix, a block prefix, and a suffix.  A k-mer needs 2*k bits
  to represent it.  The file prefix is 6 bits wide, representing one of the 64
  possible files in a database.  Inside a file, k-mers are stored in blocks,
  each k-mer in a block will have the same block prefix.  The suffix data for a
  block is stored using Elias-Fano encoding (CITE) where each suffix is split
  into two pieces.  The first piece is encoded as a unary offset from the last
  piece, and the second piece is a N-log2(N)+1 binary value.  At present,
  values are stored as plain 32-bit integers, stored immediately after the k-mer
  data.
