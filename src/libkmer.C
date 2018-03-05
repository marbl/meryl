
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


uint32 kmerTiny::_merSize;
uint32 kmerTiny::_merSpan;
uint64 kmerTiny::_mask;
uint32 kmerTiny::_leftShift;





kmerCountFileReader::kmerCountFileReader(const char *inputName) {
}



kmerCountFileReader::~kmerCountFileReader() {
}
  


kmerCountFileWriter::kmerCountFileWriter(const char *outputName,
                                         uint32      splitting,
                                         uint32      merSize,
                                         uint32      prefixSize,
                                         uint32      suffixSize) {

  if (splitting > 32)
    fprintf(stderr, "kmerCountFileWriter()-- ERROR: splitting %u too high, must be at most 32.\n", splitting), exit(1);

  //  Save the output name for later use.

  strncpy(_outName, outputName, FILENAME_MAX);

  _prefixSize = prefixSize;
  _suffixSize = suffixSize;
  _merSize    = merSize;
  _splitting  = splitting;

  //  Allocate output files, one file per 'split', stored in a directory supplied by the user.
  //
  //  splitting = 0 -> 1 file.
  //  splitting = 1 -> 2 files.
  //  ...
  //  splitting >= prefixLength -> maximum files, one per prefix.

  AS_UTL_mkdir(_outName);

  _datFilesLen = (uint32)1 << splitting;
  _datFiles    = new FILE * [_datFilesLen];

  memset(_datFiles, 0, sizeof(FILE *) * _datFilesLen);
}



kmerCountFileWriter::~kmerCountFileWriter() {

  for (uint32 ii=0; ii<_datFilesLen; ii++)
    AS_UTL_closeFile(_datFiles[ii]);
}



void
kmerCountFileWriter::addMer(kmer   k,
                            uint32 c) {
}



void
kmerCountFileWriter::addBlock(uint64  prefix,
                              uint64  nKmers,
                              uint64 *suffixes,
                              uint32 *counts) {

  //  Map the prefix to an output file.
  //  If more split files than prefixes, ignore the unused ones.

  uint32  oib = _prefixSize;
  uint32  oi  =  prefix;

  if (_splitting < _prefixSize) {
    oib = _splitting;
    oi  = prefix >> (_prefixSize - _splitting);;
  }

  //  Open a new file, if needed.

  if (_datFiles[oi] == NULL) {
    kmer     k;

    //  The kmer size is fixed globally, but we want to use it for converting
    //  a prefix (oi) to an ACGT string.  We just ignore the unset stuff.

    k.setPrefixSuffix(0, oi, 0);

    char     N[FILENAME_MAX+1];
    char     M[_merSize+1];

    k.toString(M);

    snprintf(N, FILENAME_MAX, "%s/%s.merylData",
             _outName, M + _merSize - (oib + 1) / 2);

    fprintf(stderr, "mer %s -> file %s\n", M, N);

    _datFiles[oi] = AS_UTL_openOutputFile(N);
  }

  //  Dump data.

  stuffedBits   *dumpData = new stuffedBits;

  dumpData->setBinary(64, 0x7461446c7972656dllu);
  dumpData->setBinary(64, 0x0a3030656c694661llu);
  dumpData->setBinary(64, 0);  //  Number of unique mers
  dumpData->setBinary(64, 0);  //  Number of distinct mers
  dumpData->setBinary(64, 0);  //  Number of total mers
  dumpData->setBinary(64, 0);  //  Histogram position
  dumpData->setBinary(64, 0);  //  Histogram length
  dumpData->setBinary(64, 0);  //  Histogram max count

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

  //  Start of the block
  dumpData->setBinary(64, prefix);
  dumpData->setBinary(64, nKmers);

  dumpData->setBinary(8,  1);    //  Kmer coding type
  dumpData->setBinary(32, unaryBits);
  dumpData->setBinary(32, binaryBits);
  dumpData->setBinary(64, 0);

  dumpData->setBinary(8,  1);    //  Count coding type
  dumpData->setBinary(64, 0);
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

  uint64  lastCount = 0;
  uint64  thisCount = 0;

  for (uint32 kk=0; kk<nKmers; kk++) {
    dumpData->setBinary(binaryBits, counts[kk]);
  }

  dumpData->dumpToFile(_datFiles[oi]);

  delete dumpData;
}


