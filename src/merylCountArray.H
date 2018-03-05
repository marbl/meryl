
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

#ifndef MERYLCOUNTARRAY_H
#define MERYLCOUNTARRAY_H

#include "meryl.H"
#include "libbits.H"

class merylCountArray {
public:
  merylCountArray(void);
  ~merylCountArray();

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
    fprintf(stderr, "add()-- bitsLen=" F_U64 " seg=" F_U64 " segPos=" F_U64 " word=" F_U32 " wordBgn=" F_U32 " wordEnd=" F_U32 "\n",
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
      fprintf(stderr, "add()-- 0 value %016" F_X64P " word[" F_U64 "][" F_U32 "] now 0x%016" F_X64P "\n", value, seg, word, _segments[seg][word]);
#endif
    }

    else if (wordEnd <= 64) {
      _segments[seg][word] |= (value << (64 - wordEnd));
#ifdef DEBUG_COUNT
      fprintf(stderr, "add()-- 1 value %016" F_X64P " word[" F_U64 "][" F_U32 "] now 0x%016" F_X64P "\n", value, seg, word, _segments[seg][word]);
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
      fprintf(stderr, "add()-- 2 value %016" F_X64P " word[" F_U64 "][" F_U32 "] now 0x%016" F_X64P "\n", value, seg, word+0, _segments[seg][word+0]);
      fprintf(stderr, "add()-- 2 value %016" F_X64P " word[" F_U64 "][" F_U32 "] now 0x%016" F_X64P "\n", value, seg, word+1, _segments[seg][word+1]);
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
      fprintf(stderr, "add()-- 3 value %016" F_X64P " word[" F_U64 "][" F_U32 "] now 0x%016" F_X64P "\n", value, seg+0, W, _segments[seg+0][W]);
      fprintf(stderr, "add()-- 3 value %016" F_X64P " word[" F_U64 "][" F_U32 "] now 0x%016" F_X64P "\n", value, seg+1, 0, _segments[seg+1][0]);
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
    fprintf(stderr, "get()-- kk=" F_U64 " bitPos=" F_U64 " seg=" F_U64 " segPos=" F_U64 " word=" F_U32 " wordBgn=" F_U32 " wordEnd=" F_U32 "\n",
            kk, bitPos, seg, segPos, word, wordBgn, wordEnd);
#endif

    if      (wordEnd <= 64) {
      bits = (_segments[seg][word] >> (64 - wordEnd)) & uint64MASK(_width);
#ifdef DEBUG_COUNT
      fprintf(stderr, "get()-- 1 from 0x%016" F_X64P " -> 0x%016" F_X64P "\n", _segments[seg][word], bits);
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
      fprintf(stderr, "get()-- 2 from 0x%016" F_X64P " -> 0x%016" F_X64P "\n", _segments[seg][word+0], bits);
#endif
      bits |= (_segments[seg][word+1] >> (64 - extraBits) & uint64MASK(extraBits));
#ifdef DEBUG_COUNT
      fprintf(stderr, "get()-- 2 from 0x%016" F_X64P " -> 0x%016" F_X64P "\n", _segments[seg][word+1], bits);
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
      fprintf(stderr, "get()-- 3 from 0x%016" F_X64P " -> 0x%016" F_X64P "\n", _segments[seg+0][W], bits);
#endif
      bits |= (_segments[seg+1][0] >> (64 - extraBits) & uint64MASK(extraBits));
#ifdef DEBUG_COUNT
      fprintf(stderr, "get()-- 3 from 0x%016" F_X64P " -> 0x%016" F_X64P "\n", _segments[seg+1][0], bits);
#endif
    }

    return(bits);  // & uint64MASK(_width));
  };

  void      sort(uint64 UNUSED(prefix)) {
    uint64   nValues = _bitsLen / _width;

    assert(_bitsLen % _width == 0);

    if (nValues == 0)
      return;

    uint64  *values = new uint64 [nValues];

    //  Unpack the data into _suffix.

    for (uint64 kk=0; kk<nValues; kk++)
      values[kk] = get(kk);

    //  Sort the data

    //fprintf(stderr, "sort()--  prefix=0x%016" F_X64P " with " F_U64 " values\n", prefix, nValues);

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
    fprintf(stderr, "sort()-- nKmers " F_U64 " tKmers " F_U64 "\n", nKmers, tKmers);
#endif

    //  Remove all the temporary data.

    delete [] values;

    removeSegments();

    //  Report.

    _bitsLen = nKmers;

#ifdef DEBUG_COUNT
    for (uint32 kk=0; kk<nKmers; kk++)
      fprintf(stderr, "sort()-- kk " F_U32 " count " F_U32 " data 0x%016" F_X64P "\n", kk, _counts[kk], _suffix[kk]);
#endif
  };


  void             report(uint64 prefix) {
    kmerTiny  kmer;
    char      str[64] = {0};

    for (uint32 kk=0; kk<_bitsLen; kk++) {
      kmer.setPrefixSuffix(prefix, _suffix[kk], _width);

      fprintf(stderr, "%s " F_U32 "\n", kmer.toString(str), _counts[kk]);
    }
  };


  void             dump(uint64 prefix, kmerCountFileWriter *out) {
    out->addBlock(prefix, _bitsLen, _suffix, _counts);
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



#endif  //  MERYLCOUNTARRAY_H
