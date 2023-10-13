.. note::
  A k-mer is a short sequence of nucleotide bases.  The **forward k-mer** is
  from the supplied sequence orientation, while the **reverse-k-mer** is from
  the reverse-complemented sequence.  The **canonical k-mer** is the
  lexicographically smaller of the forward-mer and reverse-mer.

  For example, the sequence GATCTCA has five forward 3-mers: GAT, ATC, TCT,
  CTC and TCA.  The canonical k-mers are found by reverse-complementing each
  of those and picking the lexicographically smaller:

  .. code-block:: none
    :caption: Forward, reverse and canonical 3-mers.
    :linenos:

    G
    A (GAT, ATC) -> ATC
    T (ATC, GAT) -> ATC
    C (TCT, AGA) -> AGA
    T (CTC, GAG) -> CTC
    C (TCA, TGA) -> TCA
    A

  In meryl, k-mers can be up to 64 bases long and are canonical by default.
