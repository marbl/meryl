
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

#include "AS_global.H"
#include "AS_UTL_fileIO.H"

#include "libsequence.H"
#include "libkmer.H"
#include "libbits.H"

#include <vector>
#include <stack>
#include <algorithm>
using namespace std;


enum merylOp {
  opCount,
  opCountForward,
  opCountReverse,
  opUnion,
  opUnionMin,
  opUnionMax,
  opUnionSum,
  opIntersect,
  opIntersectMin,
  opIntersectMax,
  opIntersectSum,
  opDifference,
  opSymmetricDifference,
  opComplement,
  opNothing
};


char const *
toString(merylOp op) {
  switch (op) {
    case opCount:                return("opCount");                break;
    case opCountForward:         return("opCountForward");         break;
    case opCountReverse:         return("opCountReverse");         break;
    case opUnion:                return("opUnion");                break;
    case opUnionMin:             return("opUnionMin");             break;
    case opUnionMax:             return("opUnionMax");             break;
    case opUnionSum:             return("opUnionSum");             break;
    case opIntersect:            return("opIntersect");            break;
    case opIntersectMin:         return("opIntersectMin");         break;
    case opIntersectMax:         return("opIntersectMax");         break;
    case opIntersectSum:         return("opIntersectSum");         break;
    case opDifference:           return("opDifference");           break;
    case opSymmetricDifference:  return("opSymmetricDifference");  break;
    case opComplement:           return("opComplement");           break;
    case opNothing:              return("opNothing");              break; 
  }

  assert(0);
  return(NULL);
}




class merylOperation;


class merylInput {
public:
  merylInput(merylOperation *o);
  merylInput(char *n, kmerStreamReader *s);
  merylInput(char *n, dnaSeqFile *f);
  ~merylInput();

  void   nextMer(void);

  merylOperation     *_operation;
  kmerStreamReader  *_stream;
  dnaSeqFile         *_sequence;

  char                _name[FILENAME_MAX];

  kmer                _kmer;
  uint64              _count;
  bool                _valid;
};




class merylOperation {
public:
  merylOperation(merylOp op=opNothing, uint32 k=0, uint64 numMers=0, uint32 threads=1, uint64 memory=0);
  ~merylOperation();

  void    addInput(merylOperation *operation);
  void    addInput(char *name, kmerStreamReader *reader);
  void    addInput(char *name, dnaSeqFile *sequence);

  void    addOutput(char *name, kmerStreamWriter *writer);

  void    setOperation(merylOp op) { _operation = op;    };
  merylOp getOperation(void)       { return(_operation); };

  bool    isCounting(void) {
    return((_operation == opCount)        ||
           (_operation == opCountForward) ||
           (_operation == opCountReverse));
  };

  bool    isNormal(void) {
    return(isCounting() == false);
  };

  kmer   &theFMer(void)            { return(_kmer);   };
  uint64  theCount(void)           { return(_count);  };

  bool    nextMer(void);
  bool    validMer(void)           { return(_valid);  };

  void    count(void);

private:
  void    findMinCount(void);
  void    findMaxCount(void);
  void    findSumCount(void);

  vector<merylInput *>           _inputs;

  merylOp                        _operation;

  uint32                         _k;
  uint64                         _numMers;
  uint32                         _threads;
  uint64                         _memory;

  char                           _outputName[FILENAME_MAX];
  kmerStreamWriter             *_output;

  kmer                           _smallest;

  uint32                         _actLen;
  uint64                        *_actCount;
  uint32                        *_actIndex;

  kmer                           _kmer;
  uint64                         _count;
  bool                           _valid;
};








merylInput::merylInput(merylOperation *o) {
  _operation   = o;
  _stream      = NULL;
  _sequence    = NULL;
  _count       = 0;
  _valid       = false;

  strncpy(_name, toString(_operation->getOperation()), FILENAME_MAX);
}


merylInput::merylInput(char *n, kmerStreamReader *s) {
  _operation   = NULL;
  _stream      = s;
  _sequence    = NULL;
  _count       = 0;
  _valid       = false;

  strncpy(_name, n, FILENAME_MAX);
}


merylInput::merylInput(char *n, dnaSeqFile *f) {
  _operation   = NULL;
  _stream      = NULL;
  _sequence    = f;
  _count       = 0;
  _valid       = true;    //  Trick nextMer into doing something without a valid mer.

  strncpy(_name, n, FILENAME_MAX);
}


merylInput::~merylInput() {
  fprintf(stderr, "Destroy input %s\n", _name);
  delete _stream;
  delete _operation;
}


void
merylInput::nextMer(void) {
  char kmerString[256];

  if (_stream) {
    fprintf(stderr, "merylIn::nextMer('%s')--\n", _name);

    _valid = _stream->nextMer();
    _kmer  = _stream->theFMer();
    _count = _stream->theCount();

    fprintf(stderr, "merylIn::nextMer('%s')-- now have valid=%u kmer %s count %lu\n",
            _name, _valid, _kmer.toString(kmerString), _count);
    fprintf(stderr, "\n");
  }

  if (_operation) {
    fprintf(stderr, "merylIn::nextMer(%s)--\n", _name);

    _valid = _operation->nextMer();
    _kmer  = _operation->theFMer();
    _count = _operation->theCount();

    fprintf(stderr, "merylIn::nextMer(%s)-- now have valid=%u kmer %s count %lu\n",
            _name, _valid, _kmer.toString(kmerString), _count);
    fprintf(stderr, "\n");
  }

  if (_sequence) {
    fprintf(stderr, "merylIn::nextMer(%s)--\n", _name);
  }
}








merylOperation::merylOperation(merylOp op, uint32 k, uint64 numMers, uint32 threads, uint64 memory) {
  fprintf(stderr, "Create operation '%s'\n", toString(op));

  _operation     = op;
  _k             = k;
  _numMers       = numMers;
  _threads       = threads;
  _memory        = memory;

  _outputName[0] = 0;
  _output        = NULL;

  _actLen        = 0;
  _actCount      = new uint64 [1024];
  _actIndex      = new uint32 [1024];

  _count         = 0;
  _valid         = true;
}


merylOperation::~merylOperation() {
  fprintf(stderr, "Destroy op %s\n", toString(_operation));

  for (uint32 ii=0; ii<_inputs.size(); ii++)
    delete _inputs[ii];

  delete    _output;
  delete [] _actCount;
  delete [] _actIndex;
}




void
merylOperation::addInput(merylOperation *operation) {
  fprintf(stderr, "Adding input from operation '%s' to operation '%s'\n",
          toString(operation->_operation), toString(_operation));

  _inputs.push_back(new merylInput(operation));
  _actIndex[_actLen++] = _inputs.size() - 1;
}


void
merylOperation::addInput(char *name, kmerStreamReader *reader) {
  fprintf(stderr, "Adding input from file '%s' to operation '%s'\n",
          name, toString(_operation));

  _inputs.push_back(new merylInput(name, reader));
  _actIndex[_actLen++] = _inputs.size() - 1;
}


void
merylOperation::addInput(char *name, dnaSeqFile *sequence) {
  fprintf(stderr, "Adding input from file '%s' to operation '%s'\n",
          name, toString(_operation));

  _inputs.push_back(new merylInput(name, sequence));
  _actIndex[_actLen++] = _inputs.size() - 1;
}



void
merylOperation::addOutput(char *name, kmerStreamWriter *writer) {
  fprintf(stderr, "Adding output to file '%s' from operation '%s'\n",
          name, toString(_operation));

  strncpy(_outputName, name, FILENAME_MAX);
  _output = writer;
}




void
merylOperation::findMinCount(void) {
  _count = _actCount[0];
  for (uint32 ii=1; ii<_actLen; ii++)
    if (_actCount[ii] < _count)
      _count = _actCount[ii];
}


void
merylOperation::findMaxCount(void) {
  _count = _actCount[0];
  for (uint32 ii=1; ii<_actLen; ii++)
    if (_count < _actCount[ii])
      _count = _actCount[ii];
}


void
merylOperation::findSumCount(void) {
  _count = 0;
  for (uint32 ii=0; ii<_actLen; ii++)
    _count += _actCount[ii];
}


bool
merylOperation::nextMer(void) {

  char  kmerString[256];

  //  Find the smallest kmer in the _inputs, and save their counts in _actCount.
  //  Mark which input was used in _actIndex.

  fprintf(stderr, "\n");
  fprintf(stderr, "merylOp::nextMer()-- STARTING for operation %s\n",
          toString(_operation));

  for (uint32 ii=0; ii<_inputs.size(); ii++)
    fprintf(stderr, "merylOp::nextMer()--   CURRENT STATE: input %s kmer %s count %lu %s\n",
            _inputs[ii]->_name,
            _inputs[ii]->_kmer.toString(kmerString),
            _inputs[ii]->_count,
            _inputs[ii]->_valid ? "valid" : "INVALID");

  //  Grab the next mer for every input that was active in the last iteration.
  //  (on the first call, all inputs were 'active' last time)
  //
  for (uint32 ii=0; ii<_actLen; ii++) {
    fprintf(stderr, "merylOp::nextMer()-- CALL NEXTMER on input actIndex %u\n", _actIndex[ii]);
    _inputs[_actIndex[ii]]->nextMer();
  }

  _actLen = 0;

  //  Log.

  for (uint32 ii=0; ii<_inputs.size(); ii++)
    fprintf(stderr, "merylOp::nextMer()--   BEFORE OPERATION: input %s kmer %s count %lu %s\n",
            _inputs[ii]->_name,
            _inputs[ii]->_kmer.toString(kmerString),
            _inputs[ii]->_count,
            _inputs[ii]->_valid ? "valid" : "INVALID");


  //  Build a list of the inputs that have the smallest kmer.

  for (uint32 ii=0; ii<_inputs.size(); ii++) {
    if (_inputs[ii]->_valid == false)
      continue;

    //  If we have no active kmer, or the input kmer is smaller than the one we
    //  have, reset the list.

    if ((_actLen == 0) ||
        (_inputs[ii]->_kmer < _kmer)) {
      _actLen = 0;
      _kmer              = _inputs[ii]->_kmer;
      _actCount[_actLen] = _inputs[ii]->_count;
      _actIndex[_actLen] = ii;
      _actLen++;

      fprintf(stderr, "merylOp::nextMer()-- Active kmer %s from input %s. reset\n", _kmer.toString(kmerString), _inputs[ii]->_name);
    }

    //  Otherwise, if the input kmer is the one we have, save the count to the list.

    else if (_inputs[ii]->_kmer == _kmer) {
      //_kmer             = _inputs[ii]->_kmer;
      _actCount[_actLen] = _inputs[ii]->_count;
      _actIndex[_actLen] = ii;
      _actLen++;

      fprintf(stderr, "merylOp::nextMer()-- Active kmer %s from input %s\n", _kmer.toString(kmerString), _inputs[ii]->_name);
    }

    //  Otherwise, the input kmer comes after the one we're examining, ignore it.

    else {
    }
  }

  //  If no active kmers, we're done.

  if (_actLen == 0) {
    fprintf(stderr, "merylOp::nextMer()-- No inputs found, all done here.\n");
    fprintf(stderr, "\n");
    _valid = false;
    return(false);
  }

  //  Otherwise, active kmers!  Figure out what the count should be.

  fprintf(stderr, "merylOp::nextMer()-- op %s activeLen %u kmer %s\n", toString(_operation), _actLen, _kmer.toString(kmerString));

  //  If math-subtract gets implemented, use negative-zero to mean "don't output" and positive-zero
  //  to mean zero.  For now, count=0 means don't output.

  //  Set the count to zero, meaning "don't output the kmer".  Intersect depends on this,
  //  skipping most of it's work if all files don't have the kmer.
  _count = 0;

  switch (_operation) {
    case opCount:
    case opCountForward:
    case opCountReverse:
      count();
      fprintf(stderr, "DONE COUNTING, EXIT.\n");
      exit(1);
      break;

    case opUnion:                           //  Union, retain largest count
      _count = 1;
      break;

    case opUnionMin:                        //  Union, retain smallest count
      findMinCount();
      break;

    case opUnionMax:                        //  Union, retain largest count
      findMaxCount();
      break;

    case opUnionSum:                        //  Union, sum all counts
      findSumCount();
      break;

    case opIntersect:                       //  Intersect
      if (_actLen == _inputs.size())
        _count = 1;
      break;

    case opIntersectMin:                    //  Intersect, retain smallest count
      if (_actLen == _inputs.size())
        findMinCount();
      break;

    case opIntersectMax:                    //  Intersect, retain largest count
      if (_actLen == _inputs.size())
        findMaxCount();
      break;

    case opIntersectSum:                    //  Intersect, sum all counts
      if (_actLen == _inputs.size())
        findSumCount();
      break;

    case opDifference:
      break;

    case opSymmetricDifference:
      break;

    case opComplement:
      break;

    case opNothing:
      break;
  }  

  //  If flagged for output, output!

  fprintf(stderr, "merylOp::nextMer()-- FINISHED for operation %s with kmer %s count %lu%s\n",
          toString(_operation), _kmer.toString(kmerString), _count, ((_output != NULL) && (_count != 0)) ? " OUTPUT" : "");
  fprintf(stderr, "\n");

  if ((_output != NULL) &&
      (_count  != 0))
    _output->addMer(_kmer, _count);

  return(true);
}




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


#undef DEBUG_COUNT


class merylCountArray {
public:
  merylCountArray(void) {
    _prefix      = 0;
    _suffix      = NULL;
    _counts      = NULL;

    _segAlloc    = 16;
    _segments    = NULL;

    allocateArray(_segments, _segAlloc);

    _segments[0] = NULL;
    _segments[1] = NULL;
    _segments[2] = NULL;
    _segments[3] = NULL;

    _bitsLen     = 0;
    _bitsMax     = _segSize;
  };

  ~merylCountArray() {
    delete [] _suffix;
    delete [] _counts;

    removeSegments();
  };

  void      removeSegments(void) {
    uint32  ss = 0;
    uint64  pp = 0;

    if (_segments == NULL)
      return;

    while (pp < _bitsLen) {
      delete [] _segments[ss];
      ss += 1;
      pp += _segSize;
    }

    delete [] _segments;
    _segments = NULL;

    _bitsLen = 0;
    _bitsMax = 0;
  };

  void      addSegment(uint32 seg) {
    if (seg >= _segAlloc)
      resizeArray(_segments, _segAlloc, _segAlloc, 2 * _segAlloc, resizeArray_copyData | resizeArray_clearNew);

    assert(_segments[seg] == NULL);

    _segments[seg] = new uint64 [_segSize / 64];
  };

  //  Add a value to the table.
  //
  //  wordPos is 0 for the high bits and 63 for the bit that represents integer 1.

  void      add(uint64 value) {
    uint64  seg       = _bitsLen / _segSize;   //  Which segment are we in?
    uint64  segPos    = _bitsLen % _segSize;   //  Bit position in that segment.

    uint32  word      = segPos / 64;           //  Which word are we in=?
    uint32  wordBgn   = segPos % 64;           //  Bit position in that word.
    uint32  wordEnd   = wordBgn + _width;

#ifdef DEBUG_COUNT
    fprintf(stderr, "\n");
    fprintf(stderr, "add()-- bitsLen=%lu seg=%lu segPos=%lu word=%u wordBgn=%u wordEnd=%u\n",
            _bitsLen, seg, segPos, word, wordBgn, wordEnd);
#endif

    //  Increment the position.

    _bitsLen += _width;

    //  If the first word and the first position, we need to allocate a segment.
    //  This catches both the case when _bitsLen=0 (we've added nothing) and when
    //  _bitsLen=_segSize (we've added exactly one segment worth of kmers).

    if ((word    == 0) &&
        (wordBgn == 0))
      addSegment(seg);


    //  If there is enough space in the current word, stash the bits there.
    //  If it's a new word, we need to special case an initialization of the word.

    if      (wordBgn == 0) {
      _segments[seg][word]  = (value << (64 - wordEnd));
#ifdef DEBUG_COUNT
      fprintf(stderr, "add()-- 0 value %016lx word[%lu][%u] now 0x%016lx\n", value, seg, word, _segments[seg][word]);
#endif
    }

    else if (wordEnd <= 64) {
      _segments[seg][word] |= (value << (64 - wordEnd));
#ifdef DEBUG_COUNT
      fprintf(stderr, "add()-- 1 value %016lx word[%lu][%u] now 0x%016lx\n", value, seg, word, _segments[seg][word]);
#endif
    }

    //  Otherwise, the value spans two words.  If these can be in the same block,
    //  stash the bits there.

    else if (segPos + _width <= _segSize) {
      uint32   extraBits = wordEnd - 64;

      assert(wordEnd > 64);

      _segments[seg][word+0] |= (value >>        extraBits);
      _segments[seg][word+1]  = (value << (64 -  extraBits));
#ifdef DEBUG_COUNT
      fprintf(stderr, "add()-- 2 value %016lx word[%lu][%u] now 0x%016lx\n", value, seg, word+0, _segments[seg][word+0]);
      fprintf(stderr, "add()-- 2 value %016lx word[%lu][%u] now 0x%016lx\n", value, seg, word+1, _segments[seg][word+1]);
#endif
    }

    //  Otherwise, the word spans two segments.  At least we know it
    //  hits the last word in the first segment and the first word in
    //  the second segment.  And that the second segment hasn't been
    //  allocated yet.  And that it's just the same math as the previous
    //  case, just in different segments instead of words.

    else {
      uint32 W         = word;  //  Just to keep things pretty.  I love my optimizer!
      uint32 extraBits = wordEnd - 64;

      addSegment(seg+1);

      _segments[seg+0][W] |= (value >>        extraBits);
      _segments[seg+1][0]  = (value << (64 -  extraBits));
#ifdef DEBUG_COUNT
      fprintf(stderr, "add()-- 3 value %016lx word[%lu][%u] now 0x%016lx\n", value, seg+0, W, _segments[seg+0][W]);
      fprintf(stderr, "add()-- 3 value %016lx word[%lu][%u] now 0x%016lx\n", value, seg+1, 0, _segments[seg+1][0]);
#endif
    }
  };

  uint64    get(uint64 kk) {
    uint64  bitPos    = kk * _width;

    uint64  seg       = bitPos / _segSize;   //  Which segment are we in?
    uint64  segPos    = bitPos % _segSize;   //  Bit position in that segment.

    uint32  word      = segPos / 64;           //  Which word are we in=?
    uint32  wordBgn   = segPos % 64;           //  Bit position in that word.
    uint32  wordEnd   = wordBgn + _width;

    uint64  bits      = 0;

    //  If the bits are entirely in a single word, copy them out.

#ifdef DEBUG_COUNT
    fprintf(stderr, "\n");
    fprintf(stderr, "get()-- kk=%lu bitPos=%lu seg=%lu segPos=%lu word=%u wordBgn=%u wordEnd=%u\n",
            kk, bitPos, seg, segPos, word, wordBgn, wordEnd);
#endif

    if      (wordEnd <= 64) {
      bits = (_segments[seg][word] >> (64 - wordEnd)) & uint64MASK(_width);
#ifdef DEBUG_COUNT
      fprintf(stderr, "get()-- 1 from 0x%016lx -> 0x%016lx\n", _segments[seg][word], bits);
#endif
    }

    //  Otherwise, the value spans two words.  If these are in the same block,
    //  grab them.

    else if (segPos + _width <= _segSize) {
      //fprintf(stderr, "return 2\n");
      uint32   extraBits = wordEnd - 64;

      assert(wordEnd > 64);

      bits  = (_segments[seg][word+0] & uint64MASK(_width - extraBits)) << extraBits;
#ifdef DEBUG_COUNT
      fprintf(stderr, "get()-- 2 from 0x%016lx -> 0x%016lx\n", _segments[seg][word+0], bits);
#endif
      bits |= (_segments[seg][word+1] >> (64 - extraBits) & uint64MASK(extraBits));
#ifdef DEBUG_COUNT
      fprintf(stderr, "get()-- 2 from 0x%016lx -> 0x%016lx\n", _segments[seg][word+1], bits);
#endif
    }

    //  Otherwise, the word spans two segments.  At least we know it
    //  hits the last word in the first segment and the first word in
    //  the second segment.  And that the second segment hasn't been
    //  allocated yet.  And that it's just the same math as the previous
    //  case, just in different segments instead of words.

    else {
      //fprintf(stderr, "return 3\n");
      uint32 W         = word;  //  Just to keep things pretty.  I love my optimizer!
      uint32 extraBits = wordEnd - 64;

      bits  = (_segments[seg+0][W] & uint64MASK(_width - extraBits)) << extraBits;
#ifdef DEBUG_COUNT
      fprintf(stderr, "get()-- 3 from 0x%016lx -> 0x%016lx\n", _segments[seg+0][W], bits);
#endif
      bits |= (_segments[seg+1][0] >> (64 - extraBits) & uint64MASK(extraBits));
#ifdef DEBUG_COUNT
      fprintf(stderr, "get()-- 3 from 0x%016lx -> 0x%016lx\n", _segments[seg+1][0], bits);
#endif
    }

    return(bits);  // & uint64MASK(_width));
  };

  void      sort(uint64 prefix) {
    uint64   nValues = _bitsLen / _width;

    assert(_bitsLen % _width == 0);

    if (nValues == 0)
      return;

    uint64  *values = new uint64 [nValues];

    //  Unpack the data into _suffix.

    for (uint64 kk=0; kk<nValues; kk++)
      values[kk] = get(kk);

    //  Sort the data

    //fprintf(stderr, "sort()--  prefix=0x%016lx with %lu values\n", prefix, nValues);

#ifdef _GLIBCXX_PARALLEL
    __gnu_sequential::
#else
    std::
#endif
    sort(values, values + nValues);

    //  Count the number of distinct kmers
  
    uint64  nKmers = 1;

    for (uint64 kk=1; kk<nValues; kk++)
      if (values[kk-1] != values[kk])
        nKmers++;

    //  Allocate space for the outputs

    _suffix = new uint64 [nKmers];
    _counts = new uint32 [nKmers];

    //  And count.

    uint64  tKmers = 0;

    _counts[tKmers] = 1;
    _suffix[tKmers] = values[0];

    for (uint64 kk=1; kk<nValues; kk++) {
      if (values[kk-1] != values[kk]) {
        tKmers++;
        _counts[tKmers] = 0;
        _suffix[tKmers] = values[kk];
      }

      _counts[tKmers]++;
    }

    tKmers++;

#ifdef DEBUG_COUNT
    fprintf(stderr, "sort()-- nKmers %lu tKmers %lu\n", nKmers, tKmers);
#endif

    //  Remove all the temporary data.

    delete [] values;

    removeSegments();

    //  Report.

    _bitsLen = nKmers;

#ifdef DEBUG_COUNT
    for (uint32 kk=0; kk<nKmers; kk++)
      fprintf(stderr, "sort()-- kk %u count %u data 0x%016lx\n", kk, _counts[kk], _suffix[kk]);
#endif
  };


  void             report(uint64 prefix) {
    kmerTiny  kmer;
    char      str[64] = {0};

    for (uint32 kk=0; kk<_bitsLen; kk++) {
      kmer.setPrefixSuffix(prefix, _suffix[kk], _width);

      fprintf(stderr, "%s %u\n", kmer.toString(str), _counts[kk]);
    }
  };


  void             dump(uint64 prefix, const char *outputNamePrefix) {
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
    while (unarySum < _bitsLen) {
      unaryBits  += 1;
      unarySum  <<= 1;
    }

    uint32  binaryBits = _width - unaryBits;

    fprintf(stderr, "for prefix 0x%08lx N=%lu, unary %u\n", prefix, _bitsLen, unaryBits);

    //  Start of the block
    dumpData->setBinary(64, prefix);
    dumpData->setBinary(64, _bitsLen);

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

    for (uint32 kk=0; kk<_bitsLen; kk++) {
      thisPrefix = _suffix[kk] >> binaryBits;

      dumpData->setUnary(thisPrefix - lastPrefix);
      dumpData->setBinary(binaryBits, _suffix[kk]);

      lastPrefix = thisPrefix;
    }

    uint64  lastCount = 0;
    uint64  thisCount = 0;

    for (uint32 kk=0; kk<_bitsLen; kk++) {
      dumpData->setBinary(binaryBits, _counts[kk]);
    }

#if 1
    char     outputName[FILENAME_MAX+1];

    uint64   outDir  = prefix / 1000;
    uint64   outFile = prefix;

    snprintf(outputName, FILENAME_MAX, "%s", outputNamePrefix);
    AS_UTL_mkdir(outputName);

    snprintf(outputName, FILENAME_MAX, "%s/%03lu", outputNamePrefix, outDir);
    AS_UTL_mkdir(outputName);

    snprintf(outputName, FILENAME_MAX, "%s/%03lu/%06lu.merylData", outputNamePrefix, outDir, outFile);

    FILE *F = AS_UTL_openOutputFile(outputName);
    dumpData->dumpToFile(F);
    AS_UTL_closeFile(F);
#endif

    delete dumpData;
  };


  void             clear(void) {
    removeSegments();
    delete [] _suffix;   _suffix = NULL;
    delete [] _counts;   _counts = NULL;
  };


  static void      set(uint32 width, uint32 segsize) {
    _width   = width;
    _segSize = segsize;
  };

  static uint32    getSegSize_bits(void)    {  return(_segSize);           };
  static uint32    getSegSize_kmers(void)   {  return(_segSize / _width);  };

private:
  uint64   _prefix;
  uint64  *_suffix;
  uint32  *_counts;

  static
  uint32   _width;        //  Size of the element we're storing

  static
  uint32   _segSize;      //  Number of bits in a segment.
  uint32   _segAlloc;     //  Number of segments we _can_ allocate (size of the array below).
  uint64 **_segments;     //  Constant size blocks.

  uint64   _bitsLen;      //  Number of bits stored.      Number of kmers (after sorting).
  uint64   _bitsMax;      //  Number of bits allocated.
};

uint32  merylCountArray::_width   = 0;
uint32  merylCountArray::_segSize = 8192 * 64;






uint64                                     //  Output: Estimated memory size in bytes
estimateSizes(uint64   maxMemory,          //  Input:  Maximum allowed memory in bytes
              uint64   nKmerEstimate,      //  Input:  Estimated number of kmers in the input
              uint32   merSize,            //  Input:  Size of kmer
              uint32  &wPrefix_,           //  Output: Number of bits in the prefix (== bucket address)
              uint64  &nPrefix_,           //  Output: Number of prefixes there are (== number of buckets)
              uint32  &wData_,             //  Output: Number of bits in kmer data
              uint64  &wDataMask_) {       //  Output: A mask to return just the data of the mer

  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "For "F_U64 " million %u-mers:\n", nKmerEstimate / 1000000, merSize);
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

    fprintf(stderr, "%6u  %4lu %cP  %4lu %cB  %4lu %cM  %4lu %cS  %4lu %cB  %4lu %cB",
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
  fprintf(stderr, "minMemory   %lu %cB\n",  scaledNumber(minMemory), scaledUnit(minMemory));
  fprintf(stderr, "wPrefix     %u\n",       wPrefix_);
  fprintf(stderr, "nPrefix     %lu\n",      nPrefix_);
  fprintf(stderr, "wData       %u\n",       wData_);
  fprintf(stderr, "wDataMask   0x%016lx\n", wDataMask_);
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

  fprintf(stderr, "Allocating " F_U64 " buckets, each with %u bits (%u words, %u kmers) of storage.\n",
          nPrefix,
          merylCountArray::getSegSize_bits(),
          merylCountArray::getSegSize_bits() / 64,
          merylCountArray::getSegSize_kmers());

  merylCountArray  *data = new merylCountArray [nPrefix];

  //  Load bases, count!


  for (uint32 ii=0; ii<_inputs.size(); ii++) {
    fprintf(stderr, "Loading kmers from '%s' into buckets.\n", _inputs[ii]->_name);

    while (_inputs[ii]->_sequence->loadBases(buffer, bufferMax, bufferLen)) {

      //fprintf(stderr, "read %lu bases from '%s'\n", bufferLen, _inputs[ii]->_name);

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
      //  fprintf(stderr, "%03lu 0x%08lx %s\n", ss, (uint64)kmers[ss], kmers[ss].toString(str));
    }

    //  Would like some kind of report here on the kmers loaded from this file.

    delete _inputs[ii]->_sequence;
    _inputs[ii]->_sequence = NULL;
  }

  //  Finished loading kmers.  Free up some space.

  delete [] kmers;
  delete [] buffer;

  //  Sort each merCountArray.

  fprintf(stderr, "Sorting and dumping.\n");

#pragma omp parallel for
  for (uint64 pp=0; pp<nPrefix; pp++) {
    data[pp].sort(pp);
    //data[pp].report(pp);
    data[pp].dump(pp, _outputName);
    data[pp].clear();
  }

  //  Dump.

  delete [] data;
}



int
main(int argc, char **argv) {
  stack<merylOperation *>   opStack;

  uint32                    merSize    = 0;
  uint64                    numMers    = 0;
  uint32                    maxThreads = (uint32)1;
  uint64                    maxMem     = (uint64)2 * 1024 * 1024 * 1024;

  uint32                    outputArg  = UINT32_MAX;

  vector<char *>  err;
  for (int32 arg=1; arg < argc; arg++) {
    char    mcidx[FILENAME_MAX];
    char    mcdat[FILENAME_MAX];
    char    opt[FILENAME_MAX];
    uint32  optLen = strlen(argv[arg]);

    //uint32   creating    = 0;
    uint32   terminating = 0;

    //  Scan for options.

    if      (strcmp(argv[arg], "-h") == 0) {
      err.push_back(NULL);
      continue;
    }

    else if (strcmp(argv[arg], "-k") == 0) {
      merSize = strtouint32(argv[++arg]);
      kmerTiny::setSize(merSize);
      continue;
    }

    else if (strcmp(argv[arg], "-n") == 0) {
      numMers = strtouint64(argv[++arg]);
      continue;
    }

    else if (strcmp(argv[arg], "-m") == 0) {
      maxMem = strtouint64(argv[++arg]) * 1024 * 1024 * 1024;
      continue;
    }

    else if (strcmp(argv[arg], "-t") == 0) {
      maxThreads = strtouint32(argv[++arg]);
      continue;
    }


    //  Save a copy of the argument, just in case it's a filename, before
    //  we munge out any brackets.

    strncpy(opt, argv[arg], FILENAME_MAX);

    //  If we have [ as the first character, make a new operation (with no
    //  operation set yet) and add it to the inputs, then push on the stack.
    //
    //  We can get 0 or 1 open bracket at a time.  Seeing two in a row is an error,
    //  but it isn't caught.

    if (opStack.empty() == true) {
      //creating = true;
    }

    if (opt[0] == '[') {
      strncpy(opt, argv[arg]+1, FILENAME_MAX);
      optLen--;

      fprintf(stderr, "CREATENEW operation\n");

      //creating = true;
    }

    //  If we have a ] as the last character, strip it off and remember.
    //
    //  We can get any number of closing brackets.

    while (opt[optLen-1] == ']') {
      opt[optLen-1] = 0;
      optLen--;

      fprintf(stderr, "TERMINATE operation\n");

      terminating++;
    }

    //  Now that brackets are stripped, make meryl database names for the arg.

    snprintf(mcidx, FILENAME_MAX, "%s.mcidx", opt);
    snprintf(mcdat, FILENAME_MAX, "%s.mcdat", opt);

    merylOp             op       = opNothing;
    kmerStreamWriter  *writer   = NULL;
    kmerStreamReader  *reader   = NULL;
    dnaSeqFile         *sequence = NULL;


    if      (opt[0] == 0)
      ;  //  Got a single bracket, nothing to do here except make it not be an error.

    else if (strcmp(opt, "union") == 0)                  op = opUnion;
    else if (strcmp(opt, "union-min") == 0)              op = opUnionMin;
    else if (strcmp(opt, "union-max") == 0)              op = opUnionMax;
    else if (strcmp(opt, "union-sum") == 0)              op = opUnionSum;
    else if (strcmp(opt, "intersect") == 0)              op = opIntersect;
    else if (strcmp(opt, "intersect-min") == 0)          op = opIntersectMin;
    else if (strcmp(opt, "intersect-max") == 0)          op = opIntersectMax;
    else if (strcmp(opt, "intersect-sum") == 0)          op = opIntersectSum;
    else if (strcmp(opt, "difference") == 0)             op = opDifference;
    else if (strcmp(opt, "symmetric-difference") == 0)   op = opSymmetricDifference;
    else if (strcmp(opt, "complement") == 0)             op = opComplement;
    else if (strcmp(opt, "count") == 0)                  op = opCount;
    else if (strcmp(opt, "count-forward") == 0)          op = opCountForward;
    else if (strcmp(opt, "count-reverse") == 0)          op = opCountReverse;

    //  If we see 'output', flag the next arg as being the output name.
    //  If this arg is flagged as output, add an output using the bracket-stripped name.
    //  If this arg is a valid meryl file, make it an input.

    else if (strcmp(opt, "output") == 0)
      outputArg = arg+1;

    else if (arg == outputArg)
      writer = new kmerStreamWriter(opt, merSize, 0, merSize/2, false);

    else if ((opStack.size() > 0) &&
             (opStack.top()->isCounting() == false) &&
             (AS_UTL_fileExists(mcidx)    == true) &&
             (AS_UTL_fileExists(mcdat)    == true))
      reader = new kmerStreamReader(opt);

    else if ((opStack.size() > 0) &&
             (opStack.top()->isCounting() == true) &&
             (AS_UTL_fileExists(opt)      == true))
      sequence = new dnaSeqFile(opt);

    else {
      char *s = new char [1024];
      snprintf(s, 1024, "Don't know what to do with '%s'.\n", opt);
      err.push_back(s);
    }


    //  Create a new operation or set inputs or output.  Or do nothing.

    if (op != opNothing) {
      merylOperation *newOp = new merylOperation(op, merSize, numMers, maxThreads, maxMem);

      if (opStack.empty() == false)
        opStack.top()->addInput(newOp);

      opStack.push(newOp);
      op = opNothing;
    }

    if ((writer != NULL) &&                    //  If nothing on the stack, wait for an
        (opStack.size() > 0)) {                //  operation to show up.
      opStack.top()->addOutput(opt, writer);
      writer = NULL;
    }

    if (reader != NULL) {                      //  A reader exists only if an
      opStack.top()->addInput(opt, reader);    //  operation is on the stack.
      reader = NULL;
    }

    if (sequence != NULL) {                    //  A sequence input exists only if
      opStack.top()->addInput(opt, sequence);  //  an operation is on the stack.
      sequence = NULL;
    }

    //  Now that we're done with the operation, if it was terminating, pop it off the stack.

    for (; terminating > 0; terminating--) {
      fprintf(stderr, "TERMINATE op '%s'\n", toString(opStack.top()->getOperation()));

      opStack.pop();
    }
  }

  //  If any errors, fail.

  if ((argc == 1) ||        //  No commands
      (err.size() > 0)) {   //  Errors
    fprintf(stderr, "usage: %s ...\n", argv[0]);
    fprintf(stderr, "  -k <K>       work with mers of size K bases\n");
    fprintf(stderr, "  -n <N>       expect N mers in the input\n");
    fprintf(stderr, "  -m <M>       use no more than M gigabytes of memory\n");
    fprintf(stderr, "  -t <T>       use no more than T compute threads\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii] != NULL)
        fprintf(stderr, "%s\n", err[ii]);

    exit(1);
  }

  //  Now just walk through the kmers until nothing is left.

  merylOperation *op = opStack.top();

  fprintf(stderr, "START\n");
  fprintf(stderr, "START operation %s\n", toString(op->getOperation()));
  fprintf(stderr, "START\n");

  //  The counting operations need to be special cased.

  //if ((op->getOperation() == opCount) ||
  //    (op->getOperation() == opCountForward) ||
  //    (op->getOperation() == opCountReverse))
  //  op->count();

  //  Now just walk through the kmers until nothing is left.

  while (op->nextMer() == true)
    ;

  //  Done!

  fprintf(stderr, "DONE\n");
  fprintf(stderr, "DONE operation %s\n", toString(op->getOperation()));
  fprintf(stderr, "DONE\n");

  delete op;

  return(0);
}
