
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



uint64
scaledNumber(uint64 n, uint32 div=1024) {

  if (n > 9999)   n /= div;
  if (n > 9999)   n /= div;
  if (n > 9999)   n /= div;
  if (n > 9999)   n /= div;
  if (n > 9999)   n /= div;

  return(n);
}


char
scaledUnit(uint64 n, uint32 div=1024) {
  char u = ' ';

  if (n > 9999)  {  n /= div; u = 'k';  }
  if (n > 9999)  {  n /= div; u = 'M';  }
  if (n > 9999)  {  n /= div; u = 'G';  }
  if (n > 9999)  {  n /= div; u = 'T';  }
  if (n > 9999)  {  n /= div; u = 'P';  }

  return(u);
}


uint64                                     //  Output: Estimated memory size in bytes
estimateSizes(uint64   UNUSED(maxMemory),          //  Input:  Maximum allowed memory in bytes
              uint64   nKmerEstimate,      //  Input:  Estimated number of kmers in the input
              uint32   merSize,            //  Input:  Size of kmer
              uint32  &wPrefix_,           //  Output: Number of bits in the prefix (== bucket address)
              uint64  &nPrefix_,           //  Output: Number of prefixes there are (== number of buckets)
              uint32  &wData_,             //  Output: Number of bits in kmer data
              uint64  &wDataMask_) {       //  Output: A mask to return just the data of the mer

  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "For "F_U64 " million " F_U32 "-mers:\n", nKmerEstimate / 1000000, merSize);
  fprintf(stderr, "\n");

  uint64   minMemory = UINT64_MAX;

  fprintf(stderr, "prefix     # of   struct   kmers/    segs/     data    total\n");
  fprintf(stderr, "  bits   prefix   memory   prefix   prefix   memory   memory\n");
  fprintf(stderr, "------  -------  -------  -------  -------  -------  -------\n");

  for (uint32 wp=1; wp < 2 * merSize; wp++) {
    uint64  nPrefix          = (uint64)1 << wp;                        //  Number of prefix == number of blocks of data

    uint64  kmersPerPrefix   = nKmerEstimate / nPrefix + 1;            //  Expected number of kmers we need to store per prefix

    uint64  segSize          = 8192 * 64;                              //  BITS per segment
    uint64  kmersPerSeg      = segSize / (2 * merSize - wp);           //  Kmers per segment

    uint64  segsPerPrefix    = kmersPerPrefix / kmersPerSeg + 1;       //  


    uint64  structMemory     = ((sizeof(merylCountArray) * nPrefix) +                  //  Basic structs
                                (sizeof(uint64 *)        * nPrefix * segsPerPrefix));  //  Pointers to segments

    uint64  dataMemory       = nPrefix * segsPerPrefix * segSize / 8;

    uint64  totalMemory      = structMemory + dataMemory;

    fprintf(stderr, "%6" F_U32P "  %4" F_U64P " %cP  %4" F_U64P " %cB  %4" F_U64P " %cM  %4" F_U64P " %cS  %4" F_U64P " %cB  %4" F_U64P " %cB",
            wp,
            scaledNumber(nPrefix),        scaledUnit(nPrefix),
            scaledNumber(structMemory),   scaledUnit(structMemory),
            scaledNumber(kmersPerPrefix), scaledUnit(kmersPerPrefix),
            scaledNumber(segsPerPrefix),  scaledUnit(segsPerPrefix),
            scaledNumber(dataMemory),     scaledUnit(dataMemory),
            scaledNumber(totalMemory),    scaledUnit(totalMemory));

    if (totalMemory < minMemory) {
      fprintf(stderr, "  *\n");

      minMemory  = totalMemory;
      wPrefix_   = wp;
      nPrefix_   = nPrefix;
      wData_     = 2 * merSize - wp;
      wDataMask_ = uint64MASK(wData_);

    } else {
      fprintf(stderr, "\n");
    }

    if (totalMemory > 4 * minMemory)
      break;
  }

  fprintf(stderr, "\n");
  fprintf(stderr, "minMemory   " F_U64 " %cB\n",  scaledNumber(minMemory), scaledUnit(minMemory));
  fprintf(stderr, "wPrefix     " F_U32 "\n",       wPrefix_);
  fprintf(stderr, "nPrefix     " F_U64 "\n",      nPrefix_);
  fprintf(stderr, "wData       " F_U32 "\n",       wData_);
  fprintf(stderr, "wDataMask   0x%016" F_X64P "\n", wDataMask_);
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");

  return(minMemory);
}



void
merylOperation::count(void) {
  uint64          bufferMax  = 1300000;
  uint64          bufferLen  = 0;
  char           *buffer     = new char     [bufferMax];

  uint64          kmersLen   = 0;
  kmerTiny       *kmers      = new kmerTiny [bufferMax];

  kmerTiny        kmer;
  uint32          kmerLoad   = 0;
  uint32          kmerValid  = kmer.merSize() - 1;
  uint32          kmerSize   = kmer.merSize();

  char            str[32];

  fprintf(stderr, "\n");
  fprintf(stderr, "merylOp::count()-- STARTING for operation %s from inputs\n", toString(_operation));

  for (uint32 ii=0; ii<_inputs.size(); ii++)
    fprintf(stderr, "merylOp::count()--   %s\n", _inputs[ii]->_name);

  if (kmerSize == 0)
    fprintf(stderr, "ERROR: need a kmer size (-k).\n"), exit(1);

  if (_numMers == 0)
    fprintf(stderr, "ERROR: need a number of mers (-n) estimate (for now).\n"), exit(1);

  //  Optimize memory for some expected number of kmers.

  uint32    wPrefix   = 0;
  uint64    nPrefix   = 0;
  uint32    wData     = 0;
  uint64    wDataMask = 0;

  estimateSizes(0, _numMers, kmerSize, wPrefix, nPrefix, wData, wDataMask);

  //  Allocate memory.

  merylCountArray::set(wData, 64 * 8192);

  fprintf(stderr, "Allocating " F_U64 " buckets, each with " F_U32 " bits (" F_U32 " words, " F_U32 " kmers) of storage.\n",
          nPrefix,
          merylCountArray::getSegSize_bits(),
          merylCountArray::getSegSize_bits() / 64,
          merylCountArray::getSegSize_kmers());

  merylCountArray  *data = new merylCountArray [nPrefix];

  //  Load bases, count!


  for (uint32 ii=0; ii<_inputs.size(); ii++) {
    fprintf(stderr, "Loading kmers from '%s' into buckets.\n", _inputs[ii]->_name);

    while (_inputs[ii]->_sequence->loadBases(buffer, bufferMax, bufferLen)) {

      //fprintf(stderr, "read " F_U64 " bases from '%s'\n", bufferLen, _inputs[ii]->_name);

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
        uint64  pp = (uint64)kmers[kk] >> wData;
        uint64  mm = (uint64)kmers[kk]  & wDataMask;

        assert(pp < nPrefix);

        data[pp].add(mm);
      }

      //for (uint64 kk=0; kk<kmersLen; kk++)
      //  fprintf(stderr, "%03" F_U64P " 0x%08" F_X64P " %s\n", ss, (uint64)kmers[ss], kmers[ss].toString(str));
    }

    //  Would like some kind of report here on the kmers loaded from this file.

    delete _inputs[ii]->_sequence;
    _inputs[ii]->_sequence = NULL;
  }

  //  Finished loading kmers.  Free up some space.

  delete [] kmers;
  delete [] buffer;

  //  MAke output files, one per thread.  Sort, dump and erase each block.

  fprintf(stderr, "Creating outputs.\n");

  kmerCountFileWriter *out = new kmerCountFileWriter(_outputName, 6, kmer.merSize(), wPrefix, wData);

  fprintf(stderr, "Writing outputs.\n");

#pragma omp parallel for
  for (uint64 pp=0; pp<nPrefix; pp++) {
    data[pp].sort(pp);
    data[pp].dump(pp, out);
    data[pp].clear();
  }

  delete [] out;

  //  Dump.

  delete [] data;
}
