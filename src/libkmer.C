
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
                   uint32  numFiles) {
  char *name = new char [FILENAME_MAX+1];
  char  bits[67];

  bits[0] = '0';
  bits[1] = 'x';

  uint32 mask = 1;
  uint32 bp   = 2;

  for (; mask < numFiles; mask <<= 1)            //  Technically, writes the reverse of the
    bits[bp++] = (outIndex & mask) ? '1' : '0';  //  prefix, but who cares?

  bits[bp] = 0;

  snprintf(name, FILENAME_MAX, "%s/%s.merylData", prefix, bits);

  return(name);
}




kmerCountFileReader::kmerCountFileReader(const char *inputName, bool ignoreStats) {

  //  Save the input name for later use.

  strncpy(_inName, inputName, FILENAME_MAX);

  //  Initialize to nothing.

  _prefixSize   = 0;
  _suffixSize   = 0;
  _merSize      = 0;
  _numFilesBits = 0;
  _numFiles     = 0;
  _datFiles     = NULL;

  _prefix      = 0;
  _activeMer   = 0;
  _activeFile  = 0;

  _nKmers      = 0;
  _nKmersMax   = 1024;
  _suffixes    = new uint64 [_nKmersMax];
  _counts      = new uint32 [_nKmersMax];

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

  _prefixSize   = indexData->getBinary(32);
  _suffixSize   = indexData->getBinary(32);
  _merSize      = indexData->getBinary(32);
  _numFilesBits = indexData->getBinary(32);
  _numFiles     = indexData->getBinary(32);
  _datFiles     = new FILE * [_numFiles];

  for (uint32 ii=0; ii<_numFiles; ii++)
    _datFiles[ii] = NULL;

  if (ignoreStats == false)
    _stats.load(indexData);

  delete indexData;

  fprintf(stderr, "Opened '%s'.  Found prefixSize %u suffixSize %u merSize %u numFiles %u\n",
          _inName, _prefixSize, _suffixSize, _merSize, _numFiles);

  //  If the global kmer size isn't set yet, set it.  Then make sure all files are the same.

  kmer  k;

  if (k.merSize() == 0)
    k.setSize(_merSize);

  if (k.merSize() != _merSize)
    fprintf(stderr, "mer size mismatch, can't process this set of files.\n"), exit(1);

  //  Remember that we haven't opened any input files yet.

  for (uint32 oi=0; oi<_numFiles; oi++)
    _datFiles[oi] = NULL;

  //  But for simplicity, open all files here.

  //for (uint32 oi=0; oi<_numFiles; oi++)
  //  openFile(oi);
}



kmerCountFileReader::~kmerCountFileReader() {

  delete [] _suffixes;
  delete [] _counts;

  for (uint32 ii=0; ii<_numFiles; ii++)
    AS_UTL_closeFile(_datFiles[ii]);

  delete [] _datFiles;
}



void
kmerCountFileReader::openFile(uint32 idx) {

  if (_datFiles[idx])
    return;

  char  *name = constructBlockName(_inName, idx, _numFiles);

  if (AS_UTL_fileExists(name))
    _datFiles[idx] = AS_UTL_openInputFile(name);

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

  stuffedBits   *dumpData = new stuffedBits();

 loadAgain:

  if (_datFiles[_activeFile] == NULL)
    openFile(_activeFile);

  bool  loaded = dumpData->loadFromFile(_datFiles[_activeFile]);    //  Try loading the next block

  if (loaded == false) {
    AS_UTL_closeFile(_datFiles[_activeFile]);                       //  If nothing loaded, close the
    _datFiles[_activeFile] = NULL;                                  //  current file.

    _activeFile++;                                                  //  Move to the next file.

    if (_numFiles <= _activeFile) {                                 //  If no more files,
      delete dumpData;                                              //  we're done.
      return(false);
    }

    goto loadAgain;                                                 //  Otherwise, try loading again.
  }

  //  We've loaded a block, so decode it.

  uint64 m1 = dumpData->getBinary(64);    //  Magic number, part 1.
  uint64 m2 = dumpData->getBinary(64);    //  Magic number, part 2.

  if ((m1 != 0x7461446c7972656dllu) ||
      (m2 != 0x0a3030656c694661llu)) {
    fprintf(stderr, "kmerCountFileReader::nextMer()-- Magic number mismatch in activeFile " F_U32 " position " F_U64 ".\n",
            _activeFile, dumpData->getPosition());
    fprintf(stderr, "kmerCountFileReader::nextMer()-- Expected 0x7461446c7972656d got 0x%016" F_X64P "\n", m1);
    fprintf(stderr, "kmerCountFileReader::nextMer()-- Expected 0x0a3030656c694661 got 0x%016" F_X64P "\n", m2);
    exit(1);
  }

  _prefix = dumpData->getBinary(64);
  _nKmers = dumpData->getBinary(64);
  //_nBits  = dumpData->getBinary(64);

  //  But if there are no kmers in this block, go back and get another block.
  //  This _should_ be caught by the writer, and the writer should not be writing empty blocks.
  //  Just in case.

  if (_nKmers == 0) {
    //fprintf(stderr, "rCountFileReader::nextMer()-- Empty block.\n");
    goto loadAgain;
  }

  uint32 kCode      = dumpData->getBinary(8);
  uint32 unaryBits  = dumpData->getBinary(32);
  uint32 binaryBits = dumpData->getBinary(32);
  uint64 k1         = dumpData->getBinary(64);

  uint32 cCode      = dumpData->getBinary(8);
  uint64 c1         = dumpData->getBinary(64);
  uint64 c2         = dumpData->getBinary(64);

  //  Resize _suffixes and _counts if needed

  resizeArrayPair(_suffixes, _counts, 0, _nKmersMax, _nKmers, resizeArray_doNothing);

  //  Decode the data!

  uint64  thisPrefix = 0;

  for (uint32 kk=0; kk<_nKmers; kk++) {
    thisPrefix += dumpData->getUnary();

    _suffixes[kk] = (thisPrefix << binaryBits) | (dumpData->getBinary(binaryBits));
  }

  for (uint32 kk=0; kk<_nKmers; kk++) {
    _counts[kk] = dumpData->getBinary(32);
  }

  //  All done.

  delete dumpData;

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



void
kmerCountFileWriter::addBlock(uint64  prefix,
                              uint64  nKmers,
                              uint64 *suffixes,
                              uint32 *counts) {

  if (nKmers == 0)
    return;

  if (_initialized == false)     //  Take care of any deferred intiialization.  This one isn't
    initialize();                //  as bad as addMer(), 'cause we're working on a block of mers.

  //  Based on the prefix, decide what output file to write to.
  //  The prefix has _prefixSize bits.  We want to save the highest _numFiles bits.

  uint32  nfb = logBaseTwo32(_numFiles - 1);     //  Number of prefix bits stored in each file.
  uint32  oi  = prefix >> (_prefixSize - nfb);   //  File index we are writing to.

  if (oi >= _numFiles)
    fprintf(stderr, "addBlock()-- prefix 0x%016lx nkmers %lu _prefixSize %u nfb %u oi %u >= numFiles %u\n",
            prefix, nKmers, _prefixSize, nfb, oi, _numFiles);
  assert(oi < _numFiles);

  //  Open a new file, if needed.

  if (_datFiles[oi] == NULL) {
    char    *name = constructBlockName(_outName, oi, _numFiles);

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


