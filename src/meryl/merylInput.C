
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



merylInput::merylInput(merylOpTemplate *ot) {
  _template = ot;
}



merylInput::merylInput(merylOpCompute *oc) {
  _compute = oc;
}



merylInput::merylInput(merylFileReader *s/*, uint32 slice*/) {
  _stream = s;

  //if (slice != uint32max)
  //  _stream->enableThreads(slice);

  //  Grab the first kmer from the input.  Without this
  //  merylOpCompute::nextMer() doesn't load any data - the active list is
  //  empty and no inputs get refreshed.
  //
  nextMer();
}



merylInput::merylInput(dnaSeqFile *f, bool doCompression) {
  _sequence         = f;

  _homopolyCompress = doCompression;
}



#ifndef CANU

merylInput::merylInput(sqStore *s, uint32 segment, uint32 segmentMax) {
}

#else

merylInput::merylInput(sqStore *s, uint32 segment, uint32 segmentMax) {
  _store            = s;

  _sqBgn            = 1;                                  //  C-style, not the usual
  _sqEnd            = _store->sqStore_lastReadID() + 1;   //  sqStore semantics!

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
      _sqEnd = _store->sqStore_lastReadID() + 1;    //  sqEnd is set to the last read, not N+1.

    fprintf(stderr, "merylInput-- segment %u/%u picked reads %u-%u out of %u\n",
            segment, segmentMax, _sqBgn, _sqEnd, _store->sqStore_lastReadID());
  }

  _read        = new sqRead;
  _readID      = _sqBgn - 1;       //  Incremented before loading the first read
  _readPos     = uint32max;
}

#endif



merylInput::~merylInput() {
  delete _stream;
  delete _compute;
  delete _sequence;
  delete _store;
  delete _read;
}



void
merylInput::nextMer(void) {

  if (_stream) {
    _valid = _stream->nextMer();
    _kmer  = _stream->theFMer();
  }

  if (_compute) {
    _valid = _compute->nextMer();
    _kmer  = _compute->theFMer();
  }
}



char const *
merylInput::name(void) {

  if (isFromTemplate() == true) {
    return("no name template");
  }

  if (isFromOperation() == true) {
    return("no name compute");
  }

  if (isFromDatabase() == true) {
    return(_stream->filename());
  }

  if (isFromSequence() == true) {
    return(_sequence->filename());
  }

#ifdef CANU
  if (isFromStore() == true) {
    return(_store->sqStore_path());
  }
#endif

  return("no name joe");
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
  if (_compute)     gotBases = false;
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
