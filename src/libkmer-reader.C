
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




kmerCountFileReader::kmerCountFileReader(const char *inputName, bool ignoreStats, bool beVerbose) {
  char   N[FILENAME_MAX+1];

  //  Save the input name for later use, but fail if
  //  the index file isn't found.

  strncpy(_inName, inputName, FILENAME_MAX);

  snprintf(N, FILENAME_MAX, "%s/merylIndex", _inName);

  if (AS_UTL_fileExists(N) == false)
    fprintf(stderr, "ERROR: '%s' doesn't appear to be a meryl input; file '%s' doesn't exist.\n",
            _inName, N), exit(1);

  //  Open the index and initialize from it.

  stuffedBits  *indexData = new stuffedBits(N);

  uint64  m1 = indexData->getBinary(64);
  uint64  m2 = indexData->getBinary(64);

  if ((m1 != 0x646e496c7972656dllu) ||
      (m2 != 0x0000617461447865llu))
    fprintf(stderr, "ERROR: '%s' doesn't look like a meryl input; file '%s' fails magic number check.\n",
            _inName, N), exit(1);

  _prefixSize    = indexData->getBinary(32);
  _suffixSize    = indexData->getBinary(32);

  _numFilesBits  = indexData->getBinary(32);
  _numBlocksBits = indexData->getBinary(32);

  _numFiles      = (uint64)1 << _numFilesBits;
  _numBlocks     = (uint64)1 << _numBlocksBits;

  _datFile       = NULL;
  _block         = new kmerCountFileReaderBlock();

  _kmer          = kmer();
  _count         = 0;

  _prefix        = 0;

  _activeMer     = 0;
  _activeFile    = 0;

  _nKmers        = 0;
  _nKmersMax     = 1024;
  _suffixes      = new uint64 [_nKmersMax];
  _counts        = new uint32 [_nKmersMax];

  if (ignoreStats == false)
    _stats.load(indexData);

  delete indexData;

  //  Check and setup the mer size if needed.

  uint32  merSize = (_prefixSize + _suffixSize) / 2;

  if (beVerbose) {
    fprintf(stderr, "Opened '%s'.\n", _inName);
    fprintf(stderr, "  prefixSize     %u\n", _prefixSize);
    fprintf(stderr, "  suffixSize     %u\n", _suffixSize);
    fprintf(stderr, "  numFilesBits   %u (%u files)\n", _numFilesBits, _numFiles);
    fprintf(stderr, "  numBlocksBits  %u (%u blocks)\n", _numBlocksBits, _numBlocks);
  }

  if (kmer::merSize() == 0)    //  If the global kmer size isn't set yet, set it.
    kmer::setSize(merSize);    //  Then make sure all files are the same.

  if (kmer::merSize() != merSize)
    fprintf(stderr, "mer size mismatch, can't process this set of files.\n"), exit(1);
}



kmerCountFileReader::~kmerCountFileReader() {

  delete [] _suffixes;
  delete [] _counts;

  AS_UTL_closeFile(_datFile);

  delete    _block;
}



bool
kmerCountFileReader::nextMer(void) {

  _activeMer++;

  //  If we've still got data, just update and get outta here.
  //  Otherwise, we need to load another block.

  if (_activeMer < _nKmers) {
    _kmer.setPrefixSuffix(_prefix, _suffixes[_activeMer], _suffixSize);
    _count = _counts[_activeMer];
    return(true);
  }

  //  Make sure all files are opened.

 loadAgain:
  if (_datFile == NULL)
    _datFile = openInputBlock(_inName, _activeFile, _numFiles);

  //  Load blocks.

  bool loaded = _block->loadBlock(_datFile, _activeFile);

  //  If nothing loaded. open a new file and try again.

  if (loaded == false) {
    AS_UTL_closeFile(_datFile);

    _activeFile++;

    if (_numFiles <= _activeFile)
      return(false);

    goto loadAgain;
  }

  //  Make sure we have space for the decoded data.

  _prefix = _block->prefix();
  _nKmers = _block->nKmers();

#ifdef SHOW_LOAD
  fprintf(stdout, "LOADED prefix %016lx nKmers %lu\n", _prefix, _nKmers);
#endif

  resizeArrayPair(_suffixes, _counts, 0, _nKmersMax, _nKmers, resizeArray_doNothing);

  //  Decode the block into _OUR_ space.
  //
  //  decodeBlock() marks the block as having no data, so the next time we loadBlock() it will
  //  read more data from disk.  For blocks that don't get decoded, they retain whatever was
  //  loaded, and do not load another block in loadBlock().

  _block->decodeBlock(_suffixes, _counts);

  //  Reset iteration, and load the first kmer.

  _activeMer = 0;

  _kmer.setPrefixSuffix(_prefix, _suffixes[_activeMer], _suffixSize);
  _count = _counts[_activeMer];

  return(true);
}
