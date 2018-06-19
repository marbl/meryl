
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

#include "libkmer.H"
#include "libbits.H"

#include "AS_UTL_fileIO.H"


uint32 kmerTiny::_merSize   = 0;
uint64 kmerTiny::_fullMask  = 0;
uint64 kmerTiny::_leftMask  = 0;
uint32 kmerTiny::_leftShift = 0;



kmerCountStatistics::kmerCountStatistics() {
  _numUnique     = 0;
  _numDistinct   = 0;
  _numTotal      = 0;

  _histMax       = 32 * 1024 * 1024;      //  256 MB of histogram data.
  _hist          = new uint64 [_histMax];

  for (uint64 ii=0; ii<_histMax; ii++)
    _hist[ii] = 0;

  _hbigLen       = 0;
  _hbigMax       = 0;
  _hbigCount     = NULL;
  _hbigNumber    = NULL;
}


kmerCountStatistics::~kmerCountStatistics() {
  delete [] _hist;
  delete [] _hbigCount;
  delete [] _hbigNumber;
}


void
kmerCountStatistics::dump(stuffedBits *bits) {
  bits->setBinary(64, _numUnique);
  bits->setBinary(64, _numDistinct);
  bits->setBinary(64, _numTotal);

  //  Find the last used histogram value.

  uint32  histLast = _histMax;

  while (histLast-- > 0)
    if (_hist[histLast] > 0)
      break;

  histLast++;

  bits->setBinary(32, histLast);                //  Out of order relative to struct to keep
  bits->setBinary(32, _hbigLen);                //  the arrays below word-aligned.

  bits->setBinary(64, histLast, _hist);
  bits->setBinary(64, _hbigLen, _hbigCount);
  bits->setBinary(64, _hbigLen, _hbigNumber);
}


void
kmerCountStatistics::dump(FILE        *outFile) {
  stuffedBits  *bits = new stuffedBits;

  dump(bits);

  bits->dumpToFile(outFile);

  delete bits;
}


void
kmerCountStatistics::load(stuffedBits *bits) {
  uint32  histLast;

  _numUnique   = bits->getBinary(64);
  _numDistinct = bits->getBinary(64);
  _numTotal    = bits->getBinary(64);

  histLast     = bits->getBinary(32);
  _hbigLen     = bits->getBinary(32);

  assert(_hist != NULL);

  _hist        = bits->getBinary(64, histLast, _hist);
  _hbigCount   = bits->getBinary(64, _hbigLen);
  _hbigNumber  = bits->getBinary(64, _hbigLen);
}


void
kmerCountStatistics::load(FILE        *inFile) {
  stuffedBits  *bits = new stuffedBits;

  bits->loadFromFile(inFile);

  load(bits);

  delete bits;
}






static
char *
constructBlockName(char   *prefix,
                   uint64  outIndex,
                   uint32  numFiles,
                   uint32  iteration) {
  char *name = new char [FILENAME_MAX+1];
  char  bits[67];

  bits[0] = '0';
  bits[1] = 'x';

  uint32 mask = 1;
  uint32 bp   = 2;

  for (; mask < numFiles; mask <<= 1)            //  Technically, writes the reverse of the
    bits[bp++] = (outIndex & mask) ? '1' : '0';  //  prefix, but who cares?

  bits[bp] = 0;

  snprintf(name, FILENAME_MAX, "%s/%s[%03u].merylData", prefix, bits, iteration);

  return(name);
}




kmerCountFileReader::kmerCountFileReader(const char *inputName, bool ignoreStats) {

  //  Save the input name for later use.

  strncpy(_inName, inputName, FILENAME_MAX);

  //  Initialize to nothing.

  _prefixSize      = 0;
  _suffixSize      = 0;
  _merSize         = 0;
  _numFilesBits    = 0;
  _numFiles        = 0;
  _datFiles        = NULL;

  _prefix          = 0;
  _activeMer       = 0;
  _activeFile      = 0;
  _activeIteration = 1;

  _nKmers          = 0;
  _nKmersMax       = 1024;
  _suffixes        = new uint64 [_nKmersMax];
  _counts          = new uint32 [_nKmersMax];

  //  Fail if _inName isn't a directory, or if the index file doesn't exist.

  char   N[FILENAME_MAX+1];

  snprintf(N, FILENAME_MAX, "%s/merylIndex", _inName);

  if (AS_UTL_fileExists(N) == false)
    fprintf(stderr, "ERROR: '%s' doesn't appear to be a meryl input; file '%s' doesn't exist.\n",
            _inName, N), exit(1);

  //  Open the index and load.

  stuffedBits  *indexData = new stuffedBits(N);

  uint64  m1 = indexData->getBinary(64);
  uint64  m2 = indexData->getBinary(64);

  if ((m1 != 0x646e496c7972656dllu) ||
      (m2 != 0x0000617461447865llu))
    fprintf(stderr, "ERROR: '%s' doesn't look like a meryl input; file '%s' fails magic number check.\n",
            _inName, N), exit(1);

  _prefixSize    = indexData->getBinary(32);
  _suffixSize    = indexData->getBinary(32);
  _merSize       = indexData->getBinary(32);
  _numFilesBits  = indexData->getBinary(32);
  _numFiles      = indexData->getBinary(32);
  _numIterations = indexData->getBinary(32);

  _datFiles      = new FILE                     * [_numIterations];
  _blocks        = new kmerCountFileReaderBlock * [_numIterations];

  _datFiles[0]   = NULL;    //  Zeroth iteration isn't used.
  _blocks[0]     = NULL;

  for (uint32 ii=1; ii<_numIterations; ii++) {
    _datFiles[ii] = NULL;
    _blocks[ii]   = new kmerCountFileReaderBlock();
  }

  if (ignoreStats == false)
    _stats.load(indexData);

  delete indexData;

  fprintf(stderr, "Opened '%s'.  Found prefixSize %u suffixSize %u merSize %u numFiles %u numIterations %u\n",
          _inName, _prefixSize, _suffixSize, _merSize, _numFiles, _numIterations);

  if (kmer::merSize() == 0)     //  If the global kmer size isn't set yet, set it.
    kmer::setSize(_merSize);    //  Then make sure all files are the same.

  if (kmer::merSize() != _merSize)
    fprintf(stderr, "mer size mismatch, can't process this set of files.\n"), exit(1);
}



kmerCountFileReader::~kmerCountFileReader() {

  delete [] _suffixes;
  delete [] _counts;

  for (uint32 ii=1; ii<_numIterations; ii++)
    AS_UTL_closeFile(_datFiles[ii]);

  delete [] _datFiles;

  for (uint32 ii=1; ii<_numIterations; ii++)
    delete _blocks[ii];

  delete [] _blocks;
}



void
kmerCountFileReader::openFile(uint32 idx, uint32 iteration) {

  if (_datFiles[iteration])
    return;

  char  *name = constructBlockName(_inName, idx, _numFiles, iteration);

  if (AS_UTL_fileExists(name))
    _datFiles[iteration] = AS_UTL_openInputFile(name);

  delete [] name;
}



bool
kmerCountFileReader::nextMer(void) {

  _activeMer++;

  //  If we've still got data, just update and get outta here.

  if (_activeMer < _nKmers) {
    _kmer.setPrefixSuffix(_prefix, _suffixes[_activeMer], _suffixSize);
    _count = _counts[_activeMer];
    return(true);
  }

  //  Otherwise, we need to load another block.

  //fprintf(stdout, "nextMer()-- need another block.\n");

  //  Make sure all files are opened.

 loadAgain:
  for (uint32 ii=1; ii<_numIterations; ii++)
    if (_datFiles[ii] == NULL)
      openFile(_activeFile, ii);

  //  Load blocks.

  bool loaded = false;

  for (uint32 ii=1; ii<_numIterations; ii++)
    loaded |= _blocks[ii]->loadBlock(_activeFile, ii, _datFiles[ii]);

  //  If nothing loaded. open a new file and try again.

  if (loaded == false) {
    for (uint32 ii=1; ii<_numIterations; ii++)
      AS_UTL_closeFile(_datFiles[ii]);

    _activeFile++;

    if (_numFiles <= _activeFile)
      return(false);

    goto loadAgain;
  }

  //  At least one block has loaded data.  Figure out which ones are current, decode their data,
  //  and merge those kmers into one list under out control.
  //
  //  decodeBlock() marks the block as having no data, so the next time we loadBlock() it will
  //  read more data from disk.  For blocks that don't get decoded, they retain whatever was
  //  loaded, and do not load another block in loadBlock().

  _prefix = UINT64_MAX;

#ifdef SHOW_LOAD
  for (uint32 ii=1; ii<_numIterations; ii++)
    fprintf(stdout, "block %u prefix %016lx nKmers %lu\n", ii, _blocks[ii]->prefix(), _blocks[ii]->nKmers());
#endif

  for (uint32 ii=1; ii<_numIterations; ii++)             //  Find the prefix we should
    if (_blocks[ii]->nKmers() > 0)                       //  be decoding.
      _prefix = min(_prefix, _blocks[ii]->prefix());

  uint32  totActive = 0;
  uint32  active[_numIterations] = {0};

  for (uint32 ii=1; ii<_numIterations; ii++)             //  Remember which ones are active.
    if ((_blocks[ii]->prefix() == _prefix) &&
        (_blocks[ii]->nKmers() > 0))
      active[totActive++] = ii;

  for (uint32 ii=0; ii<totActive; ii++)                  //  Decode blocks with the
    _blocks[ active[ii] ]->decodeBlock();                 //  next prefix.

  uint64  totnKmers = 0;
  uint64  maxnKmers = 0;

  for (uint32 ii=0; ii<totActive; ii++) {                //  Decide how many kmers we
    uint64  bk = _blocks[ active[ii] ]->nKmers();         //  will need to store here.

    totnKmers += bk;
    maxnKmers  = max(maxnKmers, bk);
  }

  resizeArrayPair(_suffixes, _counts, 0, _nKmersMax, totnKmers, resizeArray_doNothing);

#ifdef SHOW_LOAD
  fprintf(stdout, "nextMer()-- found _prefix 0x%016lx totActive %u totnKmers %lu maxnKmers %lu\n",
          _prefix, totActive, totnKmers, maxnKmers);
#endif

  //  If only one active, we can just copy the data from it to us.

  if (totActive == 1) {
    _nKmers = totnKmers;

    memcpy(_suffixes, _blocks[ active[0] ]->_suffixes, sizeof(uint64) * _nKmers);
    memcpy(_counts,   _blocks[ active[0] ]->_counts,   sizeof(uint32) * _nKmers);
  }

  //  Otherwise, we need to merge data.

  else {
    uint32    p[totActive] = {0};  //  Position in s[] and c[]
    uint64    l[totActive] = {0};  //  Number of entries in s[] and c[]
    uint64   *s[totActive] = {0};  //  Pointer to the suffixes for piece x
    uint32   *c[totActive] = {0};  //  Pointer to the counts   for piece x

    for (uint32 ii=0; ii<totActive; ii++) {
      p[ii] = 0;
      l[ii] = _blocks[ active[ii] ]->nKmers();
      s[ii] = _blocks[ active[ii] ]->_suffixes;
      c[ii] = _blocks[ active[ii] ]->_counts;
    }

#ifdef SHOW_LOAD
    for (uint32 ii=0; ii<totActive; ii++)
      fprintf(stdout, "merging block %u (in active[%u]) with %lu kmers\n", active[ii], ii, l[ii]);
#endif

    _nKmers = 0;

    while (1) {
      uint64  minSuffix = UINT64_MAX;
      uint32  sumCount  = 0;

      //  Find the smallest suffix over all the active inputs,
      //  and remember the sum of their counts.

      for (uint32 ii=0; ii<totActive; ii++) {
        if (s[ii] == NULL)   //  If all done, s[] is set to NULL.
          continue;

        if (minSuffix > s[ii][ p[ii] ]) {
          minSuffix = s[ii][ p[ii] ];
          sumCount  = c[ii][ p[ii] ];
        }

        else if (minSuffix == s[ii][ p[ii] ]) {
          sumCount += c[ii][ p[ii] ];
        }
      }

      //  If no counts, we're done.

      if ((minSuffix == UINT64_MAX) && (sumCount == 0))
        break;

      //  Set the suffix/count in our merged list, reallocating if needed.

      _suffixes[_nKmers] = minSuffix;
      _counts  [_nKmers] = sumCount;

      _nKmers++;

      if (_nKmers > _nKmersMax)
        fprintf(stderr, "_nKmers %lu > _nKmersMax %lu\n", _nKmers, _nKmersMax);
      assert(_nKmers <= _nKmersMax);

      //  Move to the next element of the lists we pulled data from.  If the list is
      //  exhausted, mark it as so.

      for (uint32 ii=0; ii<totActive; ii++) {
        if (s[ii] == NULL)   //  If all done, s[] is set to NULL.
          continue;

        if (minSuffix == s[ii][ p[ii] ]) {
          if (++p[ii] >= l[ii]) {
            s[ii] = NULL;
            c[ii] = NULL;
          }
        }
      }
    }

#ifdef SHOW_LOAD
    fprintf(stdout, "merged %lu kmers from %lu total and %lu max\n", _nKmers, totnKmers, maxnKmers);
#endif
  }

  //  Phew, lots of work just to do this.

  _activeMer = 0;

  _kmer.setPrefixSuffix(_prefix, _suffixes[_activeMer], _suffixSize);
  _count = _counts[_activeMer];

  return(true);
}



void
kmerCountFileWriter::initialize(uint32 prefixSize) {

  //  If our mer size is set, we're already intiialized.

  if (_initialized == true)
    return;

  //  But if the global mersize isn't set, we're hosed.

  if (kmer::merSize() == 0)
    fprintf(stderr, "kmerCountFileWriter::initialize()-- asked to initialize, but kmer::merSize() is zero!\n"), exit(1);

  //  The count operations will write data in parallel, and we defer initialization until a block of data
  //  is actually written.  Thus, we need to have some 

#pragma omp critical
  if (_initialized == false) {

    //  When counting large sets, we might run out of memory.  If so, whatever we have loaded is
    //  dumped to disk, and a new batch is counted.  This counts how many batches we've written.

    _iteration = 1;

    //  If the prefixSize is zero, set it to (arbitrary) 1/4 the kmer size.  This happens in the
    //  streaming writer (which is used when meryl does any non-count operation).  The prefixSize here
    //  just controls how often we dump blocks to the file.

    if (_prefixSize == 0)
      _prefixSize = prefixSize;

    if (_prefixSize == 0)
      _prefixSize = min((uint32)8, 2 * kmer::merSize() / 3);  //  In bits!

    //  Now that we know the kmer size, we can set up the rest of our stuff.
    //
    //  Decide how many files to write.  We can make up to 2^32 files, but will
    //  run out of file handles _well_ before that.  For now, limit to 2^6 = 64 files.

    _prefixShift   = 2 * kmer::merSize() - _prefixSize;
    _suffixMask    = uint64MASK(_prefixShift);

    _suffixSize    = 2 * kmer::merSize() - _prefixSize;
    _merSize       =     kmer::merSize();
    _numFilesBits  = (_prefixSize < 7) ? _prefixSize : 6;
    _numFiles      = (uint32)1 << _numFilesBits;
    _datFiles      = new FILE * [_numFiles];
    _locks         = new pthread_mutex_t [_numFiles];

    //  And remember that we haven't created any output files yet.

    for (uint32 ii=0; ii<_numFiles; ii++) {
      _datFiles[ii] = NULL;
      pthread_mutex_init(&_locks[ii], NULL);
    }

    //  Now we're initialized!

    fprintf(stderr, "kmerCountFileWriter()-- Creating '%s' with prefixSize %u suffixSize %u merSize %u numFiles %u\n",
            _outName, _prefixSize, _suffixSize, _merSize, _numFiles);

    _initialized = true;
  }
}



kmerCountFileWriter::kmerCountFileWriter(const char *outputName,
                                         uint32      prefixSize) {

  //  Note that we're not really initialized yet.  We could call initialize() in some cases,
  //  but the interesting one can't initialized() until the first meryl input file is opened,
  //  so we don't initialize any of them.

  _initialized   = false;

  //  Save the output directory name, and try to make it.  If we can't we'll fail quickly.

  strncpy(_outName, outputName, FILENAME_MAX);

  AS_UTL_mkdir(_outName);

  //  Set the output file name (directory) and prefix size.  If we know the kmer size already,
  //  we can finish the initialization, otherwise, it'll get done on the first write.
  //
  //  prefixSize is allowed to be zero.  If so, it'll get reset to something in initialize().

  _batchPrefix   = 0;
  _batchNumKmers = 0;
  _batchMaxKmers = 131072;
  _batchSuffixes = NULL;
  _batchCounts   = NULL;

  _prefixShift   = 0;
  _suffixMask    = 0;

  _prefixSize    = prefixSize;
  _suffixSize    = 0;
  _merSize       = 0;
  _numFilesBits  = 0;
  _numFiles      = 0;
  _datFiles      = NULL;
  _locks         = NULL;
}



kmerCountFileWriter::~kmerCountFileWriter() {

  //  Dump any last block and release the space.

  if (_batchNumKmers > 0)
    addBlock(_batchPrefix, _batchNumKmers, _batchSuffixes, _batchCounts);

  delete [] _batchSuffixes;
  delete [] _batchCounts;

  //  Close all the data files.

  for (uint32 ii=0; ii<_numFiles; ii++) {
    AS_UTL_closeFile(_datFiles[ii]);
    pthread_mutex_destroy(&_locks[ii]);
  }

  delete [] _datFiles;
  delete [] _locks;

  //  Then create and store a master index.

  stuffedBits  *indexData = new stuffedBits;

  indexData->setBinary(64, 0x646e496c7972656dllu);
  indexData->setBinary(64, 0x0000617461447865llu);
  indexData->setBinary(32, _prefixSize);
  indexData->setBinary(32, _suffixSize);
  indexData->setBinary(32, _merSize);
  indexData->setBinary(32, _numFilesBits);
  indexData->setBinary(32, _numFiles);
  indexData->setBinary(32, _iteration);

  _stats.dump(indexData);

  fprintf(stderr, "kmerCountFileWriter()-- Closing '%s' with prefixSize %u suffixSize %u merSize %u numFiles %u\n",
          _outName, _prefixSize, _suffixSize, _merSize, _numFiles);

  char     N[FILENAME_MAX+1];
  FILE    *F;

  snprintf(N, FILENAME_MAX, "%s/merylIndex", _outName);

  F = AS_UTL_openOutputFile(N);
  indexData->dumpToFile(F);
  AS_UTL_closeFile(F);

  delete indexData;
}



void
kmerCountFileWriter::addMer(kmer   k,
                            uint32 c) {

  if (_initialized == false)   //  Unfortunately, we need to know prefixShift and suffixMask, so need
    initialize();              //  to be initialized.  Sadly, we need to check for every mer.

  if (_batchSuffixes == NULL) {
    _batchSuffixes = new uint64 [_batchMaxKmers];
    _batchCounts   = new uint32 [_batchMaxKmers];
  }

  uint64  prefix = (uint64)k >> _prefixShift;
  uint64  suffix = (uint64)k  & _suffixMask;

  bool  dump1 = (_batchNumKmers >= _batchMaxKmers);
  bool  dump2 = (_batchPrefix != prefix) && (_batchNumKmers > 0);

#if 0
  fprintf(stderr, "prefixShift %u suffixMask 0x%016lx\n", _prefixShift, _suffixMask);
  fprintf(stderr, "dump1 %d _batchNumKmers %lu _batchMaxKmers %lu\n", dump1, _batchNumKmers, _batchMaxKmers);
  fprintf(stderr, "dump2 %d _batchPrefix 0x%016lx\n", dump2, _batchPrefix);
  fprintf(stderr, "              prefix 0x%016lx\n", prefix);
  fprintf(stderr, "                   k 0x%016lx\n", (uint64)k);
#endif

  if (dump1 || dump2) {
    addBlock(_batchPrefix, _batchNumKmers, _batchSuffixes, _batchCounts);

    _batchPrefix   = prefix;
    _batchNumKmers = 0;
  }

  _batchSuffixes[_batchNumKmers] = suffix;
  _batchCounts[_batchNumKmers] = c;

  _batchNumKmers++;
}



uint64
kmerCountFileWriter::firstPrefixInFile(uint32 ff) {
  uint32  nfb = logBaseTwo32(_numFiles - 1);     //  Number of prefix bits stored in each file.
  uint64  pp  = ff;

  assert(_initialized);

  pp <<= (_prefixSize - nfb);   //  The first prefix is just the file number shifted left.

  return(pp);
}



uint64
kmerCountFileWriter::lastPrefixInFile(uint32 ff) {
  uint32  nfb = logBaseTwo32(_numFiles - 1);     //  Number of prefix bits stored in each file.
  uint64  pp  = ff + 1;

  assert(_initialized);

  pp <<= (_prefixSize - nfb);   //  The first prefix of the _next_ file.
  pp  -= 1;                     //  Subtract one to get the last valid prefix for file ff.

  return(pp);
}



uint32
kmerCountFileWriter::fileNumber(uint64  prefix) {

  assert(_initialized);

  if (_initialized == false)     //  Take care of any deferred intiialization.  This one isn't
    initialize();                //  as bad as addMer(), 'cause we're working on a block of mers.

  //  Based on the prefix, decide what output file to write to.
  //  The prefix has _prefixSize bits.  We want to save the highest _numFiles bits.

  uint32  nfb = logBaseTwo32(_numFiles - 1);     //  Number of prefix bits stored in each file.
  uint32  oi  = prefix >> (_prefixSize - nfb);   //  File index we are writing to.

  if (oi >= _numFiles)
    fprintf(stderr, "addBlock()-- prefix 0x%016lx _prefixSize %u nfb %u oi %u >= numFiles %u\n",
            prefix, _prefixSize, nfb, oi, _numFiles);
  assert(oi < _numFiles);

  return(oi);
}



void
kmerCountFileWriter::addBlock(uint64  prefix,
                              uint64  nKmers,
                              uint64 *suffixes,
                              uint32 *counts) {

  //if (nKmers == 0)
  //  return;

  assert(_initialized);

  if (_initialized == false)     //  Take care of any deferred intiialization.  This one isn't
    initialize();                //  as bad as addMer(), 'cause we're working on a block of mers.

  //  Based on the prefix, decide what output file to write to.
  //  The prefix has _prefixSize bits.  We want to save the highest _numFiles bits.

  uint32  nfb = logBaseTwo32(_numFiles - 1);     //  Number of prefix bits stored in each file.
  uint32  oi  = prefix >> (_prefixSize - nfb);   //  File index we are writing to.

#if 0
  fprintf(stderr, "addBlock()-- thread %2u prefix 0x%016lx nkmers %8lu _prefixSize %u nfb %u oi %u >= numFiles %u\n",
          omp_get_thread_num(), prefix, nKmers, _prefixSize, nfb, oi, _numFiles);
#endif

  if (oi >= _numFiles)
    fprintf(stderr, "addBlock()-- prefix 0x%016lx nkmers %lu _prefixSize %u nfb %u oi %u >= numFiles %u\n",
            prefix, nKmers, _prefixSize, nfb, oi, _numFiles);
  assert(oi < _numFiles);

  //  Open a new file, if needed.

  if (_datFiles[oi] == NULL) {
    char    *name = constructBlockName(_outName, oi, _numFiles, _iteration);

    pthread_mutex_lock(&_locks[oi]);

    if (_datFiles[oi] == NULL)
      _datFiles[oi] = AS_UTL_openOutputFile(name);

    pthread_mutex_unlock(&_locks[oi]);

    delete [] name;
  }

  //  Figure out the optimal size of the Elias-Fano prefix.  It's just log2(N)-1.

  uint32  unaryBits = 0;
  uint64  unarySum  = 1;
  while (unarySum < nKmers) {
    unaryBits  += 1;
    unarySum  <<= 1;
  }

  uint32  binaryBits = _suffixSize - unaryBits;

  //fprintf(stderr, "for prefix 0x%08" F_X64P " N=" F_U64 ", unary " F_U32 "\n",
  //        prefix, nKmers, unaryBits);

  //  Dump data.

  stuffedBits   *dumpData = new stuffedBits;

  dumpData->setBinary(64, 0x7461446c7972656dllu);    //  Magic number, part 1.
  dumpData->setBinary(64, 0x0a3030656c694661llu);    //  Magic number, part 2.

  dumpData->setBinary(64, prefix);
  dumpData->setBinary(64, nKmers);
  //dumpData->setBinary(64, 0);   //  number of bits in data, not needed

  dumpData->setBinary(8,  1);                        //  Kmer coding type
  dumpData->setBinary(32, unaryBits);                //  Kmer coding parameters
  dumpData->setBinary(32, binaryBits);
  dumpData->setBinary(64, 0);

  dumpData->setBinary(8,  1);                        //  Count coding type
  dumpData->setBinary(64, 0);                        //  Count coding parameters
  dumpData->setBinary(64, 0);

  //  Split the kmer suffix into two pieces, one unary encoded offsets and one binary encoded.

  uint64  lastPrefix = 0;
  uint64  thisPrefix = 0;

  for (uint32 kk=0; kk<nKmers; kk++) {
    thisPrefix = suffixes[kk] >> binaryBits;

    dumpData->setUnary(thisPrefix - lastPrefix);
    dumpData->setBinary(binaryBits, suffixes[kk]);

    lastPrefix = thisPrefix;
  }

  //  Save the counts, too.  Eventually these will be cleverly encoded.  Really.

  uint64  lastCount = 0;
  uint64  thisCount = 0;

  for (uint32 kk=0; kk<nKmers; kk++) {
    dumpData->setBinary(32, counts[kk]);
  }

  //  Store how many bits are stored in this block of bits.

  //dumpData->setPosition(4 * 64);
  //dumpData->setBinary(64, dumpData->getPosition());

  //  Dump data to disk and cleanup.

  pthread_mutex_lock(&_locks[oi]);

  //fprintf(stderr, "           dump %lu = 8 * %lu + %lu bits data to file %u\n",
  //        dumpData->getPosition(),
  //        dumpData->getPosition() / 8, dumpData->getPosition() % 8,
  //        oi);

  dumpData->dumpToFile(_datFiles[oi]);

  pthread_mutex_unlock(&_locks[oi]);

  delete dumpData;

  //  Finally, don't forget to insert the counts into the histogram!

  for (uint32 kk=0; kk<nKmers; kk++) {
    _stats.addCount(counts[kk]);
  }
}



void
kmerCountFileWriter::incrementIteration(void) {

  for (uint32 ii=0; ii<_numFiles; ii++) {
    AS_UTL_closeFile(_datFiles[ii]);
    _datFiles[ii] = NULL;
  }

  _iteration++;
}
