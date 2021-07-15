
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

#ifndef MERYLINPUT_H
#define MERYLINPUT_H

#include "meryl.H"
#include "kmers.H"
#include "sequence.H"

#ifdef CANU                   //  If we're building under
#include "sqStore.H"          //  Canu, include support
#else                         //  to read the seqStore
typedef  uint32  sqStore;     //  directly, otherwise,
typedef  uint32  sqRead;      //  add fake types to reduce
#endif                        //  the number of ifdef blocks.

class merylOpTemplate;
class merylOpCompute;


class merylInput {
public:
  merylInput(merylOpTemplate *ot);
  merylInput(merylOpCompute  *oc);
  merylInput(merylFileReader *s/*, uint32 slice=uint32max*/);
  merylInput(dnaSeqFile *f, bool doCompression);
  merylInput(sqStore *s, uint32 segment, uint32 segmentMax);
  ~merylInput();

  char const *name(void);

  char const *inputType(void) {
    if (_template)    return("meryl-template");
    if (_compute)     return("meryl-operation");
    if (_stream)      return("meryl-database");
    if (_sequence)    return("sequence-file");
    if (_store)       return("canu-seqStore");

    return("invalid-input");
  }

  //  The 'main' function - load the next kmer from the input source
  //  and save it internally so it can be returned in _kmer.
  void   nextMer(void);

private:
  bool   loadBasesFromCanu(char    *seq,
                           uint64   maxLength,
                           uint64  &seqLength,
                           bool    &endOfSequence);

public:
  bool   loadBases(char    *seq,
                   uint64   maxLength,
                   uint64  &seqLength,
                   bool    &endOfSequence);

  bool   isFromTemplate(void)     { return(_template  != nullptr); };
  bool   isFromOperation(void)    { return(_compute   != nullptr); };
  bool   isFromDatabase(void)     { return(_stream    != nullptr); };
  bool   isFromSequence(void)     { return(_sequence  != nullptr); };
  bool   isFromStore(void)        { return(_store     != nullptr); };

  bool   isCompressedFile(void)   { return((_sequence != nullptr) &&
                                           (_sequence->isCompressed())); };

  merylOpTemplate          *_template  = nullptr;
  merylOpCompute           *_compute   = nullptr;
  merylFileReader          *_stream    = nullptr;
  dnaSeqFile               *_sequence  = nullptr;
  sqStore                  *_store     = nullptr;

  //  For _compute and _stream: a copy of the 'active' kmer.

  kmer                      _kmer;
  bool                      _valid     = false;

  //  For _store: the position we're at in the store.

  uint32                    _sqBgn     = 0;
  uint32                    _sqEnd     = 0;

  sqRead                   *_read      = nullptr;
  uint32                    _readID    = 0;
  uint32                    _readPos   = uint32max;

  //  For homopoly compression: a flag, and the last base read.

  bool                      _homopolyCompress = false;
  char                      _lastByte         = 0;
};

#endif  //  MERYLINPUT_H
