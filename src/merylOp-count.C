
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

  uint64   minMemory = UINT64_MAX;

  fprintf(stderr, "\n");
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

    if ((wp >= 3) && (totalMemory < minMemory)) {
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
  fprintf(stderr, "Expecting to use " F_U64 " %cB memory to count " F_U64 " million " F_U32 "-mers.\n",
          scaledNumber(minMemory), scaledUnit(minMemory),
          nKmerEstimate / 1000000, merSize);
  fprintf(stderr, "\n");

  return(minMemory);
}



void
merylOperation::count(void) {
  uint64          bufferMax  = 1300000;
  uint64          bufferLen  = 0;
  char           *buffer     = new char     [bufferMax];
  bool            endOfSeq   = false;

  //uint64          kmersLen   = 0;
  //kmerTiny       *kmers      = new kmerTiny [bufferMax];

  kmerTiny        fmer;
  kmerTiny        rmer;

  uint32          kmerLoad   = 0;
  uint32          kmerValid  = fmer.merSize() - 1;
  uint32          kmerSize   = fmer.merSize();

  char            fstr[65];
  char            rstr[65];

  if (fmer.merSize() == 0)
    fprintf(stderr, "ERROR: Kmer size (-k) not supplied.\n"), exit(1);

  if (_numMers == 0)
    fprintf(stderr, "ERROR: Estimate of number of kmers (-n) not supplied.\n"), exit(1);

  fprintf(stderr, "\n");
  fprintf(stderr, "Counting %s%s%s " F_U32 "-mers from " F_SIZE_T " input file%s:\n",
          (_operation == opCount)        ? "canonical" : "",
          (_operation == opCountForward) ? "forward" : "",
          (_operation == opCountReverse) ? "reverse" : "",
          fmer.merSize(), _inputs.size(), (_inputs.size() == 1) ? "" : "s");

  for (uint32 ii=0; ii<_inputs.size(); ii++)
    fprintf(stderr, "  %s\n", _inputs[ii]->_name);

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
  fprintf(stderr, "\n");

  merylCountArray  *data = new merylCountArray [nPrefix];

  //  Load bases, count!

  for (uint32 ii=0; ii<_inputs.size(); ii++) {
    fprintf(stderr, "Loading kmers from '%s' into buckets.\n", _inputs[ii]->_name);

    while (_inputs[ii]->_sequence->loadBases(buffer, bufferMax, bufferLen, endOfSeq)) {
      //fprintf(stderr, "read " F_U64 " bases from '%s'\n", bufferLen, _inputs[ii]->_name);

      //kmersLen = 0;

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

        bool    useF = (_operation == opCountForward);
        uint64  pp   = 0;
        uint64  mm   = 0;

        if (_operation == opCount)
          useF = (fmer < rmer);


        if (useF == true) {
          pp = (uint64)fmer >> wData;
          mm = (uint64)fmer  & wDataMask;
          //fprintf(stderr, "F %s %s %u pp %lu mm %lu\n", fmer.toString(fstr), rmer.toString(rstr), fmer.merSize(), pp, mm);
        }

        else {
          pp = (uint64)rmer >> wData;
          mm = (uint64)rmer  & wDataMask;
          //fprintf(stderr, "R %s %s %u pp %lu mm %lu\n", fmer.toString(fstr), rmer.toString(rstr), rmer.merSize(), pp, mm);
        }
        
        assert(pp < nPrefix);

        data[pp].add(mm);
      }

      if (endOfSeq)                   //  If the end of the sequence, clear
        kmerLoad = 0;                 //  the running kmer.
    }

    //  Would like some kind of report here on the kmers loaded from this file.

    delete _inputs[ii]->_sequence;
    _inputs[ii]->_sequence = NULL;
  }

  //  Finished loading kmers.  Free up some space.

  //delete [] kmers;
  delete [] buffer;

  //  MAke output files, one per thread.  Sort, dump and erase each block.

  fprintf(stderr, "\n");
  fprintf(stderr, "Creating " F_U64 " output files in directory '%s', using " F_S32 " threads.\n",
          nPrefix, _outputName, omp_get_max_threads());

  kmerCountFileWriter *out = new kmerCountFileWriter(_outputName, wPrefix);

#pragma omp parallel for
  for (uint64 pp=0; pp<nPrefix; pp++) {
    data[pp].sort(pp);
    data[pp].dump(pp, out);
    data[pp].clear();
  }

  //  Cleanup.

  delete    out;
  delete [] data;

  fprintf(stderr, "\n");
  fprintf(stderr, "Finished counting.\n");

  //  Convert this to a pass-through operation.

  _inputs.clear();

  fprintf(stderr, "Convert to input for '%s'\n", _outputName);

  addInput(_outputName, new kmerCountFileReader(_outputName));

  _outputName[0]  = 0;
  _output         = NULL;

  _operation      = opUnion;

}
