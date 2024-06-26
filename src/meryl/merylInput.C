
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

#include "meryl.H"



merylInput::merylInput(merylOperation *o) {
  _operation = o;
  strncpy(_name, toString(_operation->getOperation()), FILENAME_MAX);
}



merylInput::merylInput(const char *n, merylFileReader *s, uint32 threadFile) {
  _stream = s;

  if (threadFile != UINT32_MAX)
    _stream->enableThreads(threadFile);

  strncpy(_name, n, FILENAME_MAX);
}



merylInput::merylInput(const char *n, dnaSeqFile *f, bool doCompression) {
  _sequence = f;

  _homopolyCompress = doCompression;
  _lastByte         = 0;

  strncpy(_name, n, FILENAME_MAX);
}



#ifndef CANU

merylInput::merylInput(const char *n, sqStore *s, uint32 segment, uint32 segmentMax) {
}

#else

merylInput::merylInput(const char *n, sqStore *s, uint32 segment, uint32 segmentMax) {
  _store = s;

  _sqBgn = 1;                                   //  C-style, not the usual
  _sqEnd = _store->sqStore_lastReadID() + 1;   //  sqStore semantics!

  if (segmentMax > 1) {
    uint64  nBases = 0;

    for (uint32 ss=1; ss <= _store->sqStore_lastReadID(); ss++)
      nBases += _store->sqStore_getReadLength(ss);

    uint64  nBasesPerSeg = nBases / segmentMax;

    _sqBgn = 0;
    _sqEnd = 0;

    nBases = 0;

    for (uint32 ss=1; ss <= _store->sqStore_lastReadID(); ss++) {
      nBases += _store->sqStore_getReadLength(ss);

      if ((_sqBgn == 0) && ((nBases / nBasesPerSeg) == segment - 1))
        _sqBgn = ss;

      if ((_sqEnd == 0) && ((nBases / nBasesPerSeg) == segment))
        _sqEnd = ss;
    }

    if (segment == segmentMax)                      //  Annoying special case; if the last segment,
      _sqEnd = _store->sqStore_lastReadID() + 1;   //  sqEnd is set to the last read, not N+1.

    fprintf(stderr, "merylInput-- segment %u/%u picked reads %u-%u out of %u\n",
            segment, segmentMax, _sqBgn, _sqEnd, _store->sqStore_lastReadID());
  }

  _read        = new sqRead;
  _readID      = _sqBgn - 1;       //  Incremented before loading the first read
  _readPos     = UINT32_MAX;

  strncpy(_name, n, FILENAME_MAX);
}

#endif


merylInput::merylInput(const char *n) {
  _operation        = NULL;
  _stream           = NULL;
  _sequence         = NULL;
  _store            = NULL;
  _filename         = n;

  _isMultiSet       = false;

  _value            = 0;
  _valid            = false;

  _sqBgn            = 0;
  _sqEnd            = 0;

  _homopolyCompress = false;
  _lastByte         = 0;
}


merylInput::~merylInput() {
  delete _stream;
  delete _operation;
  delete _sequence;
  delete _store;
  delete _read;
}



void
merylInput::initialize(void) {
  if (_operation) {
    _operation->initialize();
    _isMultiSet = _operation->isMultiSet();
  }

  if (_stream) {
    _isMultiSet = _stream->isMultiSet();
  }
}



void
merylInput::nextMer(void) {
  char kmerString[256];

  if (_stream) {
    //fprintf(stderr, "merylIn::nextMer(%s)-- (stream)\n", _name);

    _valid = _stream->nextMer();
    _kmer  = _stream->theFMer();
    _value = _stream->theValue();
  }

  if (_operation) {
    //fprintf(stderr, "merylIn::nextMer(%s)-- (operation)\n", _name);

    _valid = _operation->nextMer();
    _kmer  = _operation->theFMer();
    _value = _operation->theValue();
  }

  //fprintf(stderr, "merylIn::nextMer(%s)-- now have valid=" F_U32 " kmer %s count " F_U64 "\n",
  //        _name, _valid, _kmer.toString(kmerString), _value);
  //fprintf(stderr, "\n");
}




#ifndef CANU

bool
merylInput::loadBasesFromCanu(char    *seq,
                              uint64   maxLength,
                              uint64  &seqLength,
                              bool    &endOfSequence) {
  return(false);
}

#else

bool
merylInput::loadBasesFromCanu(char    *seq,
                              uint64   maxLength,
                              uint64  &seqLength,
                              bool    &endOfSequence) {

  //  If no read currently loaded, load one, or return that we're done.
  //  We need to loop so we can ignore the length zero reads in seqStore
  //  that exist after correction/trimming.

  while (_readPos >= _read->sqRead_length()) {
    _readID++;

    if (_readID >= _sqEnd)  //  C-style iteration, not usual sqStore semantics.
      return(false);

    _store->sqStore_getRead(_readID, _read);
    _readPos = 0;
  }

  //  How much of the read is left to return?

  uint32  len = _read->sqRead_length() - _readPos;

  assert(len > 0);

  //  If the output space is big enough to hold the rest of the read, copy it,
  //  flagging it as the end of a sequence, and setup to load the next read.

  if (len <= maxLength) {
    memcpy(seq, _read->sqRead_sequence() + _readPos, sizeof(char) * len);

    _readPos       = _read->sqRead_length();

    seqLength      = len;
    endOfSequence  = true;
  }

  //  Otherwise, only part of the data will fit in the output space.

  else {
    memcpy(seq, _read->sqRead_sequence() + _readPos, sizeof(char) * maxLength);

    _readPos      += maxLength;

    seqLength      = maxLength;
    endOfSequence  = false;
  }

  return(true);
}

#endif



bool
merylInput::loadBases(char    *seq,
                      uint64   maxLength,
                      uint64  &seqLength,
                      bool    &endOfSequence) {
  bool  gotBases = false;

  seqLength     = 0;
  endOfSequence = true;

  if (_stream)      gotBases = false;
  if (_operation)   gotBases = false;
  if (_sequence)    gotBases = _sequence->loadBases(seq, maxLength, seqLength, endOfSequence);
  if (_store)       gotBases = loadBasesFromCanu(seq, maxLength, seqLength, endOfSequence);

  //  Homopoly compress if there are bases.
  if ((gotBases) && (_homopolyCompress))
    seqLength = homopolyCompress(seq, seqLength, seq, NULL, _lastByte);

  //  Save the last byte of the buffer.
  if ((seqLength > 0) && (endOfSequence == false))
    _lastByte = seq[seqLength - 1];
  else
    _lastByte = 0;

  return(gotBases);
}
