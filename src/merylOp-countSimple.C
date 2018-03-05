
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

  uint64          kmersLen   = 0;
  kmerTiny       *kmers      = new kmerTiny [bufferMax];

  kmerTiny        kmer;
  uint32          kmerLoad   = 0;
  uint32          kmerValid  = kmer.merSize() - 1;
  uint32          kmerSize   = kmer.merSize();

  uint64          maxKmer    = (uint64)1 << kmerSize;

  char            str[32];

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

  uint8        *lowBits  = new uint8    [maxKmer];
  bitArray     *highBits = new bitArray [64];

  memset(lowBits,  0, sizeof(uint8) * maxKmer);

  //  Load bases, count!

  for (uint32 ii=0; ii<_inputs.size(); ii++) {
    fprintf(stderr, "Loading kmers from '%s' into buckets.\n", _inputs[ii]->_name);

    while (_inputs[ii]->_sequence->loadBases(buffer, bufferMax, bufferLen)) {

      //fprintf(stderr, "read %lu bases from '%s'\n", bufferLen, _inputs[ii]->_name);

      //  Process the buffer of bases into a new list of kmers.
      //
      //  If not a valid base, reset the kmer size counter and skip the base.
      //
      //  Otherwise, a valid base.  Add it to the kmer, then save the kmer
      //  in the list of kmers if it is a full valid kmer.

      kmersLen = 0;

      for (uint64 bb=0; bb<bufferLen; bb++) {
        if ((buffer[bb] != 'A') && (buffer[bb] != 'a') &&
            (buffer[bb] != 'C') && (buffer[bb] != 'c') &&
            (buffer[bb] != 'G') && (buffer[bb] != 'g') &&
            (buffer[bb] != 'T') && (buffer[bb] != 't')) {
          kmerLoad = 0;
          continue;
        }

        kmer.addR(buffer[bb]);

        if (kmerLoad == kmerValid)
          kmers[kmersLen++] = kmer;
        else
          kmerLoad++;
      }

      //  If we didn't read a full buffer, the sequence ended, and we need to reset the kmer.

      if (bufferLen < bufferMax) {
        //fprintf(stderr, "END OF SEQUENCE\n");
        kmerLoad = 0;
      }

      //  Now, just pass our list of kmers to the counting engine.

      for (uint64 kk=0; kk<kmersLen; kk++) {
        uint64  kidx = (uint64)kmers[kk];
        uint32  hib  = 0;

        //  If we can add one to the low bits, do it and get outta here.

        if (lowBits[kidx] < 255) {
          lowBits[kidx]++;
          continue;
        }

        //  Otherwise, we need to do some manual addition.

        lowBits[kidx] = 0;

        for (uint32 hib=0; hib < 64; hib++) {
          highBits[hib].allocate();

          if (highBits[hib].flipBit(kidx) == 0)   //  If not set,
            break;                                //  set it and stop.
        }
      }
    }

    //  Would like some kind of report here on the kmers loaded from this file.

    delete _inputs[ii]->_sequence;
    _inputs[ii]->_sequence = NULL;
  }

  //  Finished loading kmers.  Free up some space.

  delete [] kmers;
  delete [] buffer;

  //  Now just dump!

  fprintf(stderr, "Dumping.\n");

  FILE **outputFiles = new FILE * [omp_get_max_threads()];

#if 0
  for (uint32 ii=0; ii<omp_get_max_threads(); ii++) {
    char     outputName[FILENAME_MAX+1];

    snprintf(outputName, FILENAME_MAX, "%s", _outputName);
    AS_UTL_mkdir(outputName);

    snprintf(outputName, FILENAME_MAX, "%s/%03lu.%3lu.merylData", _outputName, ii, 0);
    outputFiles[ii] = AS_UTL_openOutputFile(outputName);
  }

#pragma omp parallel for
  for (uint64 pp=0; pp<nPrefix; pp++) {
    int32  tn = omp_get_thread_num();

    data[pp].sort(pp);
    data[pp].dump(pp, tn, outputFiles[tn]);
    data[pp].clear();
  }

  for (uint32 ii=0; ii<omp_get_max_threads(); ii++) {
    AS_UTL_closeFile(outputFiles[ii]);
  }
#endif

  //  Cleanup.

  delete [] lowBits;
  delete [] highBits;
}
