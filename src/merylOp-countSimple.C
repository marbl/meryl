
/******************************************************************************
 *
 *  This file is part of 'sequence' and/or 'meryl', software programs for
 *  working with DNA sequence files and k-mers contained in them.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2018-FEB-26
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.license' in the root directory of this distribution contains
 *  full conditions and disclaimers.
 */

#include "meryl.H"
#include "libbits.H"



void
merylOperation::countSimple(void) {
  uint64          bufferMax  = 1300000;
  uint64          bufferLen  = 0;
  char           *buffer     = new char     [bufferMax];
  bool            endOfSeq   = false;

  uint64          kmersLen   = 0;
  kmerTiny       *kmers      = new kmerTiny [bufferMax];

  kmerTiny        fmer;
  kmerTiny        rmer;

  uint32          kmerLoad   = 0;
  uint32          kmerValid  = fmer.merSize() - 1;
  uint32          kmerSize   = fmer.merSize();

  uint64          maxKmer    = (uint64)1 << (2 * kmerSize);

  char            str[32];
  char            strf[32];
  char            strr[32];

  fprintf(stderr, "\n");
  fprintf(stderr, "merylOp::count()-- STARTING for operation %s from inputs\n", toString(_operation));

  for (uint32 ii=0; ii<_inputs.size(); ii++)
    fprintf(stderr, "merylOp::count()--   %s\n", _inputs[ii]->_name);

  if (kmerSize == 0)
    fprintf(stderr, "ERROR: need a kmer size (-k).\n"), exit(1);

  //  Allocate memory.

  fprintf(stderr, "Allocating 8-bit storage for " F_U64 " kmers.\n",
          maxKmer);

  //  Tiny counter uses a fixed-width array (8-bit counts, for convenience)
  //  and a variable-width array to store counts.  Each kmer is explicitly stored.

  uint8        *lowBits    = new uint8    [maxKmer];
  bitArray     *highBits   = new bitArray [64];
  uint32        highBitMax = 0;

  memset(lowBits,  0, sizeof(uint8) * maxKmer);

  //  Load bases, count!

  for (uint32 ii=0; ii<_inputs.size(); ii++) {
    fprintf(stderr, "Loading kmers from '%s' into buckets.\n", _inputs[ii]->_name);

    while (_inputs[ii]->_sequence->loadBases(buffer, bufferMax, bufferLen, endOfSeq)) {
      //fprintf(stderr, "read %lu bases from '%s'\n", bufferLen, _inputs[ii]->_name);

      kmersLen = 0;

      for (uint64 bb=0; bb<bufferLen; bb++) {
        if ((buffer[bb] != 'A') && (buffer[bb] != 'a') &&   //  If not valid DNA, don't
            (buffer[bb] != 'C') && (buffer[bb] != 'c') &&   //  make a kmer, and reset
            (buffer[bb] != 'G') && (buffer[bb] != 'g') &&   //  the count until the next
            (buffer[bb] != 'T') && (buffer[bb] != 't')) {   //  valid kmer is available.
          kmerLoad = 0;
          continue;
        }

        fmer.addR(buffer[bb]);
        rmer.addL(buffer[bb]);

        if (kmerLoad < kmerValid) {   //  If not a full kmer, increase the length we've
          kmerLoad++;                 //  got loaded, and keep going.
          continue;
        }

        if      (_operation == opCount)
          kmers[kmersLen++] = (fmer < rmer) ? fmer : rmer;

        else if (_operation == opCountForward)
          kmers[kmersLen++] = fmer;

        else
          kmers[kmersLen++] = rmer;
      }

      if (endOfSeq)                   //  If the end of the sequence, clear
        kmerLoad = 0;                 //  the running kmer.

      //  Now, just pass our list of kmers to the counting engine.

      for (uint64 kk=0; kk<kmersLen; kk++) {
        uint64  kidx = (uint64)kmers[kk];
        uint32  hib  = 0;

        assert(kidx < maxKmer);

        //  If we can add one to the low bits, do it and get outta here.

        if (lowBits[kidx] < 255) {
          lowBits[kidx]++;
          continue;
        }

        //  Otherwise, we need to do some manual addition.

        lowBits[kidx] = 0;

        for (uint32 hib=0; hib < 64; hib++) {
          highBits[hib].allocate(maxKmer);

          if (highBits[hib].flipBit(kidx) == 0) {   //  If not set, set it,
            highBitMax = max(highBitMax, hib);      //  remember the possible maximum bit set,
            break;                                  //  and stop.
          }
        }
      }
    }

    //  Would like some kind of report here on the kmers loaded from this file.

    delete _inputs[ii]->_sequence;
    _inputs[ii]->_sequence = NULL;
  }

  //  Finished loading kmers.  Free up some space before dumping.

  delete [] kmers;
  delete [] buffer;

  //  Then dump.

  fprintf(stderr, "Dumping, with %d threads.\n", omp_get_max_threads());

  //  The number of blocks MUST be a power of two.

  uint32                 wPrefix    = 10;
  uint32                 wSuffix    = fmer.merSize() * 2 - wPrefix;

  uint64                 nPrefix    = ((uint64)1 << wPrefix);
  uint64                 nSuffix    = ((uint64)1 << wSuffix);

  uint64                 sMask      = ((uint64)1 << wSuffix) - 1;

  kmerCountFileWriter   *outputFile = new kmerCountFileWriter(_outputName, fmer.merSize(), wPrefix, wSuffix);

  uint32                 nThreads    = omp_get_max_threads();

#pragma omp parallel for
  for (uint64 pp=0; pp<nPrefix; pp++) {
    uint64  bStart   = pp * nSuffix;
    uint64  bEnd     = pp * nSuffix + nSuffix;

    uint64  *sBlock  = new uint64 [nSuffix];
    uint32  *cBlock  = new uint32 [nSuffix];
    uint64   nKmers  = 0;

    for (uint64 bp=bStart; bp<bEnd; bp++) {
      uint32  count = 0;

      for (uint32 aa=highBitMax+1; aa-- > 0; ) {    //  Reconstruct the count.
        count <<= 1;
        count  |= highBits[aa].getBit(bp);
      }

      count <<= 8;
      count  |= lowBits[bp];

      if (count > 0) {
        sBlock[nKmers] = bp & sMask;
        cBlock[nKmers] = count;
        nKmers++;
      }
    }

    //if (nKmers > 0)
    //  fprintf(stderr, "Dumping block pp %lu from %lu-%lu with %lu kmers.\n", pp, bStart, bEnd, nKmers);

    outputFile->addBlock(pp, nKmers, sBlock, cBlock);

    delete [] sBlock;
    delete [] cBlock;
  }

  //  Cleanup.

  delete    outputFile;

  delete [] lowBits;
  delete [] highBits;
}
