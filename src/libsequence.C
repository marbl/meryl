
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

#include "libsequence.H"


dnaSeqFile::dnaSeqFile(const char *filename, bool indexed) {

  _file     = new compressedFileReader(filename);
  _buffer   = new readBuffer(_file->file());

  _index    = NULL;
  _indexLen = 0;
  _indexMax = 0;

  if (indexed == false)
    return;

  if (_file->isCompressed() == true)
    fprintf(stderr, "ERROR: cannot index compressed input '%s'.\n", filename), exit(1);

  if (_file->isNormal() == false)
    fprintf(stderr, "ERROR: cannot index pipe input.\n"), exit(1);

  generateIndex();
}



dnaSeqFile::~dnaSeqFile() {
  delete    _file;
  delete    _buffer;
  delete [] _index;
}



bool
dnaSeqFile::findSequence(uint64 i) {

  if (_indexLen == 0)   return(false);
  if (_indexLen <= i)   return(false);

  _buffer->seek(_index[i]._fileOffset);

  return(true);
}



uint64
dnaSeqFile::sequenceLength(uint64 i) {

  if (_indexLen == 0)   return(UINT64_MAX);
  if (_indexLen <= i)   return(UINT64_MAX);

  return(_index[i]._sequenceLength);
}




bool
dnaSeqFile::findSequence(const char *name) {
  return(false);
}



bool
dnaSeqFile::loadIndex(void) {
  char   indexName[FILENAME_MAX+1];

  snprintf(indexName, FILENAME_MAX, "%s.index", _file->filename());

  if (AS_UTL_fileExists(indexName) == false)
    return(false);

  FILE   *indexFile = AS_UTL_openInputFile(indexName);

  AS_UTL_safeRead(indexFile, &_indexLen, "indexLen",        1, sizeof(uint64));

  _index = new dnaSeqIndexEntry [_indexLen];
  AS_UTL_safeRead(indexFile,  _index,    "index",   _indexLen, sizeof(dnaSeqIndexEntry));

  AS_UTL_closeFile(indexFile, indexName);

  return(true);
}



void
dnaSeqFile::saveIndex(void) {
  char   indexName[FILENAME_MAX+1];

  snprintf(indexName, FILENAME_MAX, "%s.index", _file->filename());

  FILE   *indexFile = AS_UTL_openOutputFile(indexName);

  AS_UTL_safeWrite(indexFile, &_indexLen, "indexLen",      1, sizeof(uint64));
  AS_UTL_safeWrite(indexFile,  _index,    "index", _indexLen, sizeof(dnaSeqIndexEntry));

  AS_UTL_closeFile(indexFile, indexName);
}



void
dnaSeqFile::generateIndex(void) {
  uint32          nameMax = 0;
  char           *name    = NULL;
  uint64          seqMax  = 0;
  char           *seq     = NULL;
  uint8          *qlt     = NULL;
  uint64          seqLen  = 0;

  if (loadIndex() == true)
    return;

  _indexLen = 0;
  _indexMax = 1048576;
  _index    = new dnaSeqIndexEntry [_indexMax];

  _index[_indexLen]._fileOffset     = _buffer->tell();
  _index[_indexLen]._sequenceLength = 0;

  //  While we read sequences:
  //    update the length of the sequence (we've already save the position)
  //    make space for more sequences
  //    save the position of the next sequence
  //
  while (loadSequence(name, nameMax, seq, qlt, seqMax, seqLen) == true) {
    _index[_indexLen]._sequenceLength = seqLen;

    increaseArray(_index, _indexLen, _indexMax, 1048576);

    _indexLen++;

    _index[_indexLen]._fileOffset     = _buffer->tell();
    _index[_indexLen]._sequenceLength = 0;
  }

  //for (uint32 ii=0; ii<_indexLen; ii++)
  //  fprintf(stderr, "%u offset %lu length %lu\n", ii, _index[ii]._fileOffset, _index[ii]._sequenceLength);

  if (_indexLen > 0)
    saveIndex();
}



uint64
dnaSeqFile::loadFASTA(char   *&name,     uint32   nameMax,
                      char   *&seq,
                      uint8  *&qlt,      uint64   seqMax) {
  uint64  nameLen = 0;
  uint64  seqLen  = 0;
  char    ch      = _buffer->read();

  assert(ch == '>');

  //  Read the header line into the name string.

  for (ch=_buffer->read(); (ch != '\n') && (ch != 0); ch=_buffer->read()) {
    if (nameLen+1 >= nameMax)
      resizeArray(name, nameLen, nameMax, 3 * nameMax / 2);
    name[nameLen++] = ch;
  }

  //  Read sequence, skipping whitespace, until we hit a new sequence (or eof).

  for (ch=_buffer->readuntil('>'); (ch != '>') && (ch != 0); ch=_buffer->readuntil('>')) {
    if (ch == '\n')
      continue;

    assert(_buffer->eof() == false);

    if (seqLen+1 >= seqMax)
      resizeArrayPair(seq, qlt, seqLen, seqMax, 3 * seqMax / 2);

    seq[seqLen] = ch;
    qlt[seqLen] = 0;

    seqLen++;
  }

  name[nameLen] = 0;
  seq[seqLen] = 0;
  qlt[seqLen] = 0;

  assert(nameLen < nameMax);
  assert(seqLen  < seqMax);

  return(seqLen);
}



uint64
dnaSeqFile::loadFASTQ(char   *&name,     uint32   nameMax,
                      char   *&seq,
                      uint8  *&qlt,      uint64   seqMax) {
  uint32  nameLen = 0;
  uint64  seqLen  = 0;
  uint64  qltLen  = 0;
  char    ch      = _buffer->read();  

  assert(ch == '@');

  //  Read the header line into the name string.

  for (ch=_buffer->read(); (ch != '\n') && (ch != 0); ch=_buffer->read()) {
    if (nameLen+1 >= nameMax)
      resizeArray(name, nameLen, nameMax, 3 * nameMax / 2);
    name[nameLen++] = ch;
  }

  //  Read sequence.

  for (ch=_buffer->read(); (ch != '\n') && (ch != 0); ch=_buffer->read()) {
    if (seqLen+1 >= seqMax)
      resizeArrayPair(seq, qlt, seqLen, seqMax, 3 * seqMax / 2);
    seq[seqLen++] = ch;
  }

  //  Skip header line

  for (ch=_buffer->read(); (ch != '\n') && (ch != 0); ch=_buffer->read()) {
    ;
  }

  //  Read qualities.

  for (ch=_buffer->read(); (ch != '\n') && (ch != 0); ch=_buffer->read()) {
    if (qltLen+1 >= seqMax)
      resizeArrayPair(seq, qlt, qltLen, seqMax, 3 * seqMax / 2);
    qlt[qltLen++] = ch;
  }

  //fprintf(stderr, "READ FASTQ name %u seq %lu qlt %lu\n", nameLen, seqLen, qltLen);

  name[nameLen] = 0;
  seq[seqLen] = 0;
  qlt[qltLen] = 0;

  assert(nameLen < nameMax);
  assert(seqLen  < seqMax);
  assert(qltLen  < seqMax);
  assert(seqLen == qltLen);

  return(seqLen);
}



bool
dnaSeqFile::loadSequence(char   *&name,     uint32   nameMax,
                         char   *&seq,
                         uint8  *&qlt,      uint64   seqMax,
                         uint64  &seqLen) {

  if (nameMax == 0)
    resizeArray(name, 0, nameMax, (uint32)1024);

  if (seqMax == 0)
    resizeArrayPair(seq, qlt, 0, seqMax, (uint64)65536);

  while (_buffer->peek() == '\n')
    _buffer->read();

  if      (_buffer->peek() == '>')
    seqLen = loadFASTA(name, nameMax,
                       seq,
                       qlt, seqMax);

  else if (_buffer->peek() == '@')
    seqLen = loadFASTQ(name, nameMax,
                       seq,
                       qlt, seqMax);

  else
    return(false);

  return(true);
}



bool
dnaSeqFile::loadBases(char    *seq,
                      uint64   maxLength,
                      uint64  &seqLength) {

  seqLength = 0;

  if (_buffer->eof() == true)
    return(false);

  //  If this is a new file, skip the first name line.

  if (_buffer->tell() == 0) {
    while (_buffer->peek() == '\n')    //  Skip whitespace before the first name line.
      _buffer->read();

#if 0
    for (char ch = _buffer->read(); (ch != '\n') && (ch != 0); ch = _buffer->read())
      ;
#else
    _buffer->skipAhead('\n');
#endif
  }

  //  Skip whitespace.

  while (_buffer->peek() == '\n')
    _buffer->read();

  //  We're now at sequence, so load until we're not in sequence or out of space.

  char    ch = _buffer->read();

  while (_buffer->eof() == false) {

    //  If we hit the next sequence, skip the header, leaving us at the start of the bases.

    if (ch == '>') {
#if 0
      for (ch = _buffer->read(); (ch != '\n') && (ch != 0); ch = _buffer->read())   //  Skip the name of the next sequence
        ;
#else
      _buffer->skipAhead('\n');
#endif
      return(true);
    }

    if (ch == '+') {
#if 0
      for (ch = _buffer->read(); (ch != '\n') && (ch != 0); ch = _buffer->read())   //  Skip the quality name line
        ;
      for (ch = _buffer->read(); (ch != '\n') && (ch != 0); ch = _buffer->read())   //  Skip the qualities
        ;
      for (ch = _buffer->read(); (ch != '\n') && (ch != 0); ch = _buffer->read())   //  Skip the name of the next sequence
        ;
#else
      _buffer->skipAhead('\n');
      _buffer->skipAhead('\n');
      _buffer->skipAhead('\n');
#endif
      return(true);
    }

    //  Otherwise, add the base and move ahead.

#if 0
    if (ch != '\n')
      seq[seqLength++] = ch;
#else
    seqLength += _buffer->copyUntil('\n', seq + seqLength, maxLength - seqLength);
#endif

    if (seqLength == maxLength)
      return(true);

    ch = _buffer->read();
  }

  return(seqLength > 0);     //  We've hit EOF.  If no bases were loaded, indicate we're all done.
}
