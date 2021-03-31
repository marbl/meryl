
/******************************************************************************
 *
 *  This file is part of meryl, a genomic k-kmer counter with nice features.
 *
 *  This software is based on:
 *    'Canu' v2.0              (https://github.com/marbl/canu)
 *  which is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "runtime.H"

#include "kmers.H"
#include "system.H"
#include "sequence.H"
#include "bits.H"
#include "strings.H"
#include <fstream>
#include <iostream>


void
filter(dnaSeqFile                       *sfile1,
       dnaSeqFile                       *sfile2,
       compressedFileWriter             *ofile1,
       compressedFileWriter             *ofile2,
       std::vector<merylExactLookup *>  &klookup,
       bool                              outputIfFound) {

  //  Do nothing if there are no sequences.

  if ((sfile1 == NULL) && (sfile2 == NULL))
    return;

  //  While we load sequences from all files supplied...

  dnaSeq  seq1;
  dnaSeq  seq2;

  uint64   nReads      = 0;
  uint64   nReadsFound = 0;

  while (((sfile1 == NULL) || (sfile1->loadSequence(seq1))) &&
         ((sfile2 == NULL) || (sfile2->loadSequence(seq2)))) {
    uint32 nKmerFound = 0;

    nReads++;

    if (seq1.length() > 0) {
      kmerIterator  kiter(seq1.bases(), seq1.length());

      while (kiter.nextMer())
        if ((klookup[0]->value(kiter.fmer()) > 0) ||
            (klookup[0]->value(kiter.rmer()) > 0))
          nKmerFound++;
    }

    if (seq2.length() > 0) {
      kmerIterator  kiter(seq2.bases(), seq2.length());

      while (kiter.nextMer())
        if ((klookup[0]->value(kiter.fmer()) > 0) ||
            (klookup[0]->value(kiter.rmer()) > 0))
          nKmerFound++;
    }

    //  Report the sequence if:
    //    any kmers are found and     ifFound
    //    no  kmers are found and not ifFound

    if ((nKmerFound > 0) == outputIfFound) {
      nReadsFound++;

      if (sfile1 != NULL) {
        if (seq1.quals()[0] == 0)   fprintf(ofile1->file(), ">%s nKmers=%u\n%s\n",        seq1.ident(), nKmerFound, seq1.bases());
        else                        fprintf(ofile1->file(), "@%s nKmers=%u\n%s\n+\n%s\n", seq1.ident(), nKmerFound, seq1.bases(), seq1.quals());
      }

      if (sfile2 != NULL) {
        if (seq2.quals()[0] == 0)   fprintf(ofile2->file(), ">%s nKmers=%u\n%s\n",        seq2.ident(), nKmerFound, seq2.bases());
        else                        fprintf(ofile2->file(), "@%s nKmers=%u\n%s\n+\n%s\n", seq2.ident(), nKmerFound, seq2.bases(), seq2.quals());
      }
    }
  }

  fprintf(stderr, "\nIncluding %lu reads (or read pairs) out of %lu.\n", nReadsFound, nReads);
}

