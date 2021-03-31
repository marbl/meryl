
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

using namespace std;  //  For ifstream and string.


void
reportExistence(dnaSeqFile                       *sfile,
                compressedFileWriter             *ofile,
                std::vector<merylExactLookup *>  &klookup,
                std::vector<const char *>        &klabel) {
  dnaSeq   seq;

  const uint32 ctgn = sfile->numberOfSequences();
  const int threads = omp_get_max_threads();
  fprintf(stderr, "\nTotal of %u sequences found. Will be processed over %d threads\n", ctgn, threads);

  //  Initializing output variables
  string seqNames[ctgn];
  uint64 nKmer[ctgn];
  uint64 nKmerFound[ctgn][klookup.size()];
  for (uint32 ii=0; ii<ctgn; ii++)  {
    nKmer[ii]=0;
    for (uint32 dd=0; dd<klookup.size(); dd++) {
      nKmerFound[ii][dd] = 0;
    } 
  } 

#pragma omp parallel for private(seq)
  for (uint32 seqId = 0; seqId < ctgn; seqId++) {

#pragma omp critical
    {
      sfile->findSequence(seqId);
      sfile->loadSequence(seq);
    }

    kmerIterator  kiter(seq.bases(), seq.length());

    seqNames[seqId] = seq.ident();

    char kmerString[65];

    while (kiter.nextMer()) {
      nKmer[seqId]++;

      for (uint32 dd=0; dd<klookup.size(); dd++) {
        if ((klookup[dd]->value(kiter.fmer()) > 0) ||
            (klookup[dd]->value(kiter.rmer()) > 0))
          nKmerFound[seqId][dd]++;
      }
    }
  }

  //  Flush output
  fprintf(stderr, "Writing sequence results.\n");

  for (uint32 seqId = 0; seqId < ctgn; seqId++) {

    fprintf(ofile->file(), "%s\t%lu", seqNames[seqId].c_str(), nKmer[seqId]);
    for (uint32 dd=0; dd<klookup.size(); dd++) {
      fprintf(ofile->file(), "\t%lu\t%lu", klookup[dd]->nKmers(), nKmerFound[seqId][dd]);
    }
    fprintf(ofile->file(), "\n");
  }
}


