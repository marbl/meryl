
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

#include <algorithm>

#include "meryl.H"




//  A suffix-with-value.  It's whole point is to store a 128-bit kmer, a
//  32-bit count, and a 64-bit label in something that can be 32-bit aligned.
//
//  This is used to sort the bit-packed input kmers.

class swv {
public:
  void      set(kmdata suffix, kmvalu value, kmlabl label) {
    _s[0] = (uint32)(suffix >> 96);
    _s[1] = (uint32)(suffix >> 64);
    _s[2] = (uint32)(suffix >> 32);
    _s[3] = (uint32)(suffix);

    _v    = value;

    _l[0] = (uint32)(label >> 32);
    _l[1] = (uint32)(label);
  };

  kmdata    getSuffix(void) const {
    kmdata   s;

    s   = _s[0];   s <<= 32;
    s  |= _s[1];   s <<= 32;
    s  |= _s[2];   s <<= 32;
    s  |= _s[3];

    return(s);
  };

  kmvalu    getValue(void) const {
    return(_v);
  };

  kmvalu    getLabel(void) const {
    kmlabl   l;

    l  = _l[0];   l <<= 32;
    l |= _l[1];

    return(l);
  };

private:
  uint32   _s[4];   //  This bit of ugly is splitting the suffix into multiple 32-bit words, so
  kmvalu   _v;      //  an swv can be aligned to 4-byte boundaries for 32-bit kmvalu.
  uint32   _l[2];   //
};


bool
lessThan(swv const &a, swv const &b) {
  kmdata  as = a.getSuffix();   kmvalu av = a.getValue();   kmlabl al = a.getLabel();
  kmdata  bs = b.getSuffix();   kmvalu bv = b.getValue();   kmlabl bl = b.getLabel();

  if (as < bs)   return(true);
  if (as > bs)   return(false);

  if (av < bv)   return(true);
  if (av > bv)   return(false);

  if (al < bl)   return(true);
  if (al > bl)   return(false);

  return(false);
}





//  Initialize a count array by setting
//   - size of the suffix to store (value and label are set in initializeValues())
//   - the kmer prefix we're counting for
//   - how big each chunk of data should be
//   - how much storage we're currently using.
//
//  The size of the data chunk is arbitrary, is in bits, and reserves
//  8 words for malloc() information.
uint64
merylCountArray::initialize(uint64 prefix, uint32 width) {
  _sWidth       = width;

  _prefix       = prefix;

  _bitsPerPage  = getPageSize() * 8;
  _segSize      = pagesPerSegment() * _bitsPerPage - 8 * 64;

  _nBitsOldSize = usedSize();

  return(_nBitsOldSize);
}



void
merylCountArray::initializeValues(uint32 valueWidth, uint32 labelWidth) {
  _vWidth = std::min(32u, valueWidth);   _vals = new stuffedBits();
  _lWidth = std::min(64u, labelWidth);   _labs = new stuffedBits();
}



void
merylCountArray::initializeForTesting(uint32 width, uint32 nwords) {
  _sWidth       = width;

  assert(_sWidth > 0);

  _prefix       = 0;
  _suffix       = nullptr;
  _counts       = nullptr;

  _nKmers       = 0;

  _bitsPerPage  = getPageSize() * 8;
  _nReAlloc     = 0;

  _segSize      = 64 * nwords;
  _segAlloc     = 0;
  _segments     = nullptr;

  _nBits        = 0;
  _nBitsTrigger = 0;
  _nBitsOldSize = 0;

  clearStats();
}



void
merylCountArray::clearStats(void) {
#ifdef ADD_INSTRUMENT
  memset(nTests, 0, sizeof(uint64) * 6);
  memset(nStart, 0, sizeof(uint64) * 6 * 64);
#endif
}



void
merylCountArray::dumpStats(void) {
#ifdef ADD_INSTRUMENT
  fprintf(stderr, "        oneThis      twoThis     thrThis  one-plus-one one-plus-two two-plus-one\n");
  fprintf(stderr, "-- ------------ ------------ ------------ ------------ ------------ ------------\n");
  fprintf(stderr, "   %12lu %12lu %12lu %12lu %12lu %12lu\n",
          nTests[0], nTests[1], nTests[2], nTests[3], nTests[4], nTests[5]);
  fprintf(stderr, "-- ------------ ------------ ------------ ------------ ------------ ------------\n");

  for (uint32 ii=0; ii<64; ii++) {
    fprintf(stderr, "%2u %12lu %12lu %12lu %12lu %12lu %12lu\n",
            ii, nStart[0][ii], nStart[1][ii], nStart[2][ii], nStart[3][ii], nStart[4][ii], nStart[5][ii]);
  }
#endif
}



void
merylCountArray::dumpData(void) {

  for (uint32 ss=0; ss<_segAlloc; ss++) {
    if (_segments[ss] == nullptr)
      continue;
      
    fprintf(stderr, "seg[%2d]", ss);

    for (uint64 bb=0; bb<_segSize; bb += 64)
      fprintf(stderr, " 0x%016lx", _segments[ss][bb/64]);

    fprintf(stderr, "\n");
  }
}



merylCountArray::~merylCountArray() {

  removeSegments();
  removeValues();

  delete [] _suffix;
  delete [] _counts;
}



void
merylCountArray::removeSegments(void) {

  if (_segments == nullptr)               //  If no segments, then
    return;                               //  we've already removed them.

  for (uint32 ss=0; ss<_segAlloc; ss++)   //  Release the segment memory.
    delete [] _segments[ss];

  delete [] _segments;                    //  Release the list of segments...

  _nReAlloc  = 0;

  _segAlloc = 0;                          //  Don't forget to
  _segments = nullptr;                    //  forget about it.

  _nBits        = 0;                      //  Indicate that we've stored no data.
  _nBitsTrigger = 0;
  _nBitsOldSize = usedSize();
}



void
merylCountArray::removeValues(void) {
  delete _vals;    _vals  = nullptr;
  delete _labs;    _labs  = nullptr;
}



void
merylCountArray::addSegment(uint32 seg) {

  if (_segAlloc == 0) {
    resizeArray(_segments, _segAlloc, _segAlloc, 64, _raAct::copyData | _raAct::clearNew);
    _nReAlloc++;
  }
  if (seg >= _segAlloc) {
    resizeArray(_segments, _segAlloc, _segAlloc, 2 * _segAlloc, _raAct::copyData | _raAct::clearNew);
    _nReAlloc++;
  }
  assert(_segments[seg] == nullptr);

  _segments[seg] = new uint64 [_segSize / 64];

  //memset(_segments[seg], 0, sizeof(uint64) * _segSize / 64);
}




//  Unpack the suffixes and remove the data.
kmdata *
merylCountArray::unpackSuffixes(uint64 nSuffixes) {
  kmdata  *suffixes  = new kmdata [nSuffixes];

  //fprintf(stderr, "Allocate %lu suffixes, %lu bytes\n", nSuffixes, sizeof(uint64) * nSuffixes);
  //fprintf(stderr, "Sorting prefix 0x%016" F_X64P " with " F_U64 " total kmers\n", _prefix, nSuffixes);

  for (uint64 kk=0; kk<nSuffixes; kk++)
    suffixes[kk] = get(kk);

  removeSegments();

  return(suffixes);
}



swv *
merylCountArray::unpackSuffixesAndValues(uint64 nSuffixes) {
  swv  *suffixes  = new swv [nSuffixes];

  assert(_vals != nullptr);
  assert(_labs != nullptr);

  //fprintf(stderr, "Allocate %lu suffixes, %lu bytes\n", nSuffixes, sizeof(swv) * nSuffixes);
  //fprintf(stderr, "Sorting prefix 0x%016" F_X64P " with " F_U64 " total kmers\n", _prefix, nSuffixes);

  _vals->setPosition(0);
  _labs->setPosition(0);

  if      (_vWidth == 0)
    for (uint64 kk=0; kk<nSuffixes; kk++)
      suffixes[kk].set(get(kk), _vals->getEliasDelta(),    _labs->getBinary(_lWidth));
  else
    for (uint64 kk=0; kk<nSuffixes; kk++)
      suffixes[kk].set(get(kk), _vals->getBinary(_vWidth), _labs->getBinary(_lWidth));

  removeSegments();
  removeValues();

  return(suffixes);
}




//
//  Converts raw kmers listed in _segments into counted kmers listed in _suffix and _counts.
//
void
merylCountArray::countSingleKmers(void) {
  uint64   nSuffixes = _nBits / _sWidth;
  kmdata  *suffixes  = unpackSuffixes(nSuffixes);

  //  Sort the data

  std::sort(suffixes, suffixes + nSuffixes);

  //  Count the number of distinct kmers, and allocate space for them.

  uint64  nk = 1;

  for (uint64 kk=1; kk<nSuffixes; kk++)
    if (suffixes[kk-1] != suffixes[kk])
      nk++;

  _suffix = new kmdata [nk];
  _counts = new kmvalu [nk];

  //  And generate the counted kmer data.

  _nKmers = 0;

  _counts[_nKmers] = 1;
  _suffix[_nKmers] = suffixes[0];

  for (uint64 kk=1; kk<nSuffixes; kk++) {
    if (suffixes[kk-1] != suffixes[kk]) {
      _nKmers++;
      _counts[_nKmers] = 0;
      _suffix[_nKmers] = suffixes[kk];
    }

    _counts[_nKmers]++;
  }

  _nKmers++;

  //  Remove all the temporary data.

  delete [] suffixes;
};



void
merylCountArray::countSingleKmersWithValues(void) {
  uint64       nSuffixes = _nBits / _sWidth;
  swv         *suffixes  = unpackSuffixesAndValues(nSuffixes);

  //  Sort the data

  std::sort(suffixes, suffixes + nSuffixes, lessThan);

  //  Count the number of distinct kmers, and allocate space for them.

  uint64  nk = 1;

  for (uint64 kk=1; kk<nSuffixes; kk++)
    if (suffixes[kk-1].getSuffix() != suffixes[kk].getSuffix())
      nk++;

  _suffix = new kmdata [nk];
  _counts = new kmvalu [nk];
  _labels = new kmlabl [nk];

  //  And generate the counted kmer data.

  _nKmers = 0;

  _suffix[_nKmers] = suffixes[0].getSuffix();
  _counts[_nKmers] = suffixes[0].getValue();
  _labels[_nKmers] = suffixes[0].getLabel();

  for (uint64 kk=1; kk<nSuffixes; kk++) {
    if (suffixes[kk-1].getSuffix() != suffixes[kk].getSuffix()) {
      _nKmers++;
      _suffix[_nKmers] = suffixes[kk].getSuffix();
      _counts[_nKmers] = 0;
      _labels[_nKmers] = 0;
    }

    _counts[_nKmers] += suffixes[kk].getValue();
    _labels[_nKmers] += suffixes[kk].getLabel();
  }

  _nKmers++;

  //  Remove all the temporary data.

  delete [] suffixes;
};



void
merylCountArray::countMultiSetKmers(void) {
  uint64      nSuffixes = _nBits / _sWidth;
  swv        *suffixes  = unpackSuffixesAndValues(nSuffixes);

  //  Sort the data

  std::sort(suffixes, suffixes + nSuffixes, lessThan);

  //  In a multi-set, we dump each and every kmer that is loaded, no merging.

  _suffix = new kmdata [nSuffixes];
  _counts = new kmvalu [nSuffixes];

  //  And generate the counted kmer data.

  _nKmers = nSuffixes;

  for (uint64 kk=0; kk<nSuffixes; kk++) {
    _counts[kk] = suffixes[kk].getValue();
    _suffix[kk] = suffixes[kk].getSuffix();
  }

  //  Remove all the temporary data.

  delete [] suffixes;
};






void
merylCountArray::countKmers(void) {

  //fprintf(stderr, "merylCountArray::countKmers()-- _nBits %lu -- values=%c multi-set=%c\n",
  //        _nBits, (_vals == nullptr) ? 'n' : 'Y', (_multiSet == false) ? 'n' : 'Y');

  if (_nBits == 0) {    //  If no kmers stored, nothing to do, so just
    removeSegments();   //  remove the (unused) storage and return.
    return;
  }

  assert(_nBits % _sWidth == 0);

  if (_vals == nullptr)
    countSingleKmers();
  else
    if (_multiSet == false)
      countSingleKmersWithValues();
    else
      countMultiSetKmers();
}



void
merylCountArray::dumpCountedKmers(merylBlockWriter *out, kmlabl label) {
  out->addCountedBlock(_prefix, _nKmers, _suffix, _counts, _labels, label);
}



void
merylCountArray::removeCountedKmers(void) {

  delete [] _suffix;   _suffix = nullptr;
  delete [] _counts;   _counts = nullptr;

  _nKmers = 0;
}



uint64
merylCountArray::add(kmdata suffix) {

  //  Compute current position, the advance pointer to next available spot.

  uint64  nBits     = _nBits;
  uint64  seg       = nBits / _segSize;    //  Which segment are we in?
  uint64  segPos    = nBits % _segSize;    //  Bit position in that segment.

  _nBits += _sWidth;

  //  word position counts from high to low; 0 for the high bit and 63 for
  //  the bit that represents integer 1.

  uint32  word      = segPos / 64;         //  Which word are we in=?
  uint32  wordBgn   = segPos % 64;         //  Bit position in that word.
  uint32  wordEnd   = wordBgn + _sWidth;

#ifdef ADD_VERBOSE
  fprintf(stderr, "at position nBits %6u - seg %2u word %2u bgn %3u end %3u\n", nBits, seg, word, wordBgn, wordEnd);
#endif

  assert(wordEnd <= 192);

  //  If the first word and the first position, we need to allocate a segment.
  //  This catches both the case when nBits=0 (we've added nothing) and when
  //  nBits=_segSize (we've added exactly one segment worth of kmers).

  if ((word    == 0) &&
      (wordBgn == 0))
    addSegment(seg);

  //  Initialize the word if it's a new one.

  if (wordBgn == 0)
    _segments[seg][word] = 0;

  //
  //  Add the suffix bits to the list of suffixes.  It's gory.  The suffix
  //  will be at most 128 bits (actually, '128 - prefix' bits.
  //
  //  It will fit:
  //
  //    completely in a new word (so we need to initialize the whole word)
  //    completely in an existing word
  //    in 2 words in the current segment
  //    in 3 words in the current segment
  //    in 1 word in the current and 1 word in the next segment
  //    in 1 word in the current and 2 words in the next segment
  //    in 2 words in the current and 1 word in the next segment
  //  

  bool    sameSeg = (nBits + _sWidth <= seg * _segSize + _segSize);          //  suffix fits in this segment

  //  Move all kmer bits to words in this segment.

  if (sameSeg) {
    bool   oneWord =                     (wordEnd <=  64);    //  suffix fits in this word
    bool   twoWord = ((64  < wordEnd) && (wordEnd <= 128));   //  suffix fits in this and next word
    bool   thrWord =  (128 < wordEnd);                        //  suffix needs three words

#ifdef ADD_VERBOSE
    fprintf(stderr, "sameSeg one=%d two=%d thr=%d\n", oneWord, twoWord, thrWord);
#endif

#ifdef ADD_INSTRUMENT
    if (oneWord) {
      nTests[0]++;
      nStart[0][wordBgn]++;
    }

    else if (twoWord) {
      nTests[1]++;
      nStart[1][wordBgn]++;
    }

    else if (thrWord) {
      nTests[2]++;
      nStart[2][wordBgn]++;
    }

    else {
      assert(0);
    }
#endif

    if (oneWord) {
      uint64   sta = (suffix << (64 - wordEnd));

      _segments[seg][word] |= sta;
    }

    if (twoWord) {
      uint32   staBits = 64 - wordBgn;
      uint32   endBits = _sWidth - staBits;

#ifdef ADD_VERBOSE
      fprintf(stderr, "staBits %u endBits %u - sum %u\n", staBits, endBits, staBits + endBits);
#endif

      uint64   sta = (suffix >>        endBits);
      uint64   end = (suffix << (64 -  endBits));

#ifdef ADD_VERBOSE
      fprintf(stderr, "twoWord sta 0x%016lx\n", sta);
      fprintf(stderr, "twoWord end 0x%016lx\n", end);
#endif

      _segments[seg][word+0] |= sta;
      _segments[seg][word+1]  = end;
    }

    if (thrWord) {
      uint32   staBits = 64 - wordBgn;
      uint32   endBits = _sWidth - 64 - staBits;

#ifdef ADD_VERBOSE
      fprintf(stderr, "staBits %u midBits %u endBits %u - sum %u\n", staBits, 64, endBits, staBits + 64 + endBits);
#endif

      uint64   sta = (suffix >> (64 + endBits));
      uint64   mid = (suffix >>      (endBits));
      uint64   end = (suffix << (64 - endBits));

#ifdef ADD_VERBOSE
      fprintf(stderr, "thrWord sta 0x%016lx\n", sta);
      fprintf(stderr, "thrWord mid 0x%016lx\n", mid);
      fprintf(stderr, "thrWord end 0x%016lx\n", end);
#endif

      _segments[seg][word+0] |= sta;
      _segments[seg][word+1]  = mid;
      _segments[seg][word+2]  = end;
    }
  }

  //  Move all kmer bits to multiple segments.

  else {
    uint32  thisBits = seg * _segSize + _segSize - nBits;
    uint32  nextBits = _sWidth - thisBits;

#ifdef ADD_VERBOSE
    fprintf(stderr, "diffSeg - thisBits %3u nextBits %3u\n", thisBits, nextBits);
#endif

    assert(wordBgn < _segSize);

    assert(thisBits <= 128);
    assert(nextBits <= 128);

    bool   oneThis =                      (thisBits <=  64);
    bool   twoThis = ((64  < thisBits) && (thisBits <= 128));
    bool   thrThis =  (128 < thisBits);

    bool   oneNext =                      (nextBits <=  64);
    bool   twoNext = ((64  < nextBits) && (nextBits <= 128));
    bool   thrNext =  (128 < nextBits);

#ifdef ADD_VERBOSE
    fprintf(stderr, "multSeg oneThis=%d twoThis=%d thrThis=%d  oneNext=%d twoNext=%d thrNext=%d\n",
            oneThis, twoThis, thrThis,
            oneNext, twoNext, thrNext);
#endif

    assert(thrThis == false);
    assert(thrNext == false);

    //

    assert((oneThis && oneNext) ||
           (oneThis && twoNext) ||
           (twoThis && oneNext));

#ifdef ADD_INSTRUMENT
    if      (oneThis && oneNext) {
      nTests[3]++;
      nStart[3][wordBgn]++;
    }

    else if (oneThis && twoNext) {
      nTests[4]++;
      nStart[4][wordBgn]++;
    }

    else if (twoThis && oneNext) {
      nTests[5]++;
      nStart[5][wordBgn]++;
    }

    else {
      assert(0);
    }
#endif

    addSegment(seg+1);

    //  Move kmer bits to one or two words in this segment.

    if (oneThis) {
      //  Move the left-most bits of the suffix to the right-most end of the last word.
      uint64  sta = (suffix >> (_sWidth - thisBits));

      assert(word+0 == _segSize/64-1);

      _segments[seg+0][word+0] |= sta;
    }

    if (twoThis) {
      //  Move the left-most bits of the suffix to the right-most end of the second to last word,
      //  then move the middle 64-bits to the last word.
      uint64  sta = (suffix >> (nextBits + 64));
      uint64  mid = (suffix >> (nextBits));

      assert(word+1 == _segSize/64-1);

      _segments[seg+0][word+0] |= sta;
      _segments[seg+0][word+1]  = mid;
    }

    //  Move kmer bits to one or two words in the next segment.

    if (oneNext) {
      uint64  sta = (suffix << (64 - nextBits));

      _segments[seg+1][0]  = sta;
    }

    if (twoNext) {
      uint64  mid = (suffix >> (nextBits - 64));
      uint64  end = (suffix << (128 - nextBits));

      _segments[seg+1][0]  = mid;
      _segments[seg+1][1]  = end;
    }
  }

  return(usedSizeDelta());
}




uint64
merylCountArray::addValue(kmvalu value) {

  if (_vals == nullptr)
    return(0);

  if (_vWidth == 0)
    return(_vals->setEliasDelta(value));
  else
    return(_vals->setBinary(_vWidth, value));
}


uint64
merylCountArray::addLabel(kmlabl label) {

  if (_labs == nullptr)
    return(0);

  if (_lWidth == 0)
    return(0);
  else
    return(_labs->setBinary(_lWidth, label));
}






kmdata
merylCountArray::get(uint64 kk) {
  uint32  width = _sWidth;
  kmdata  bits  = 0;

  uint64  bitBgn    = kk * _sWidth;

  uint64  seg       = bitBgn / _segSize;   //  Which segment are we in?
  uint64  segPos    = bitBgn % _segSize;   //  Bit position in that segment.

  uint32  word      = segPos / 64;         //  Which word are we in=?
  uint32  wordBgn   = segPos % 64;         //  Bit position in that word.
  uint32  wordEnd   = wordBgn + _sWidth;

  //  If the bits are entirely in one word, be done.

  if (wordEnd <= 64) {
    bits = (_segments[seg][word] >> (64 - wordEnd)) & buildLowBitMask<uint64>(_sWidth);

    return(bits);
  }

  //  Otherwise, the data spans multiple words, and possibly multiple
  //  segments.
  //
  //  Iterate, copying in partial or full words until we get all the bits.

  while (width > 0) {
    seg       = bitBgn / _segSize;   //  Which segment are we in?
    segPos    = bitBgn % _segSize;   //  Bit position in that segment.

    word      = segPos / 64;         //  Which word are we in=?
    wordBgn   = segPos % 64;         //  Bit position in that word.

    //  If the bits start in the middle of a word, we know that they extend
    //  to the right side of the word.  We also know this can only happen for
    //  the first set of bits, so 'bits' should be empty already, and we don't
    //  need to make space for the new bits.
    //
    //  Tricky: if wordBgn == 0, we need to copy in a full word (as data that
    //  fit entirely inside one word is handled in the special case above) and
    //  that is done in the middle block below.

    if (wordBgn > 0) {
      uint64  w = _segments[seg][word];
      uint32  s = 64 - wordBgn;

      w <<= wordBgn;      //  Mask out the high bits we don't care about,
      w >>= wordBgn;      //  and shift back to the original location.

      assert(bits == 0);

      bits <<= s;         //  Make space in the result for the new
      bits  |= w;         //  bits and append them.

      width  -= s;
      bitBgn += s;
    }

    //  Otherwise, if we have more than a full word left to copy in,
    //  copy in the next full word of data.

    else if (width >= 64) {
      uint64  w = _segments[seg][word];
      uint32  s = 64;

      assert(wordBgn == 0);

      bits <<= s;         //  Make space in the result for the new
      bits  |= w;         //  bits and append them.

      width  -= s;
      bitBgn += s;
    }

    //  Lastly, we must have less than a full word left to copy.  We know
    //  that these bits are left aligned, so all we need to do is shift
    //  them to be right aligned.

    else {
      uint64  w = _segments[seg][word];
      uint32  s = width;

      w >>= (64 - width);

      assert(wordBgn == 0);
      assert(width < 64);

      bits <<= s;         //  Make space in the result for the new
      bits  |= w;         //  bits and append them.

      width  -= s;
      bitBgn += s;
    }
  }

  return(bits);
}



kmdata
merylCountArray::getSimple(uint64 kk) {
  kmdata  ret  = 0;
  uint64  dPos = _sWidth * kk;

  //  Slow, but easy to verify.  And, importantly, completely different from
  //  how we add words.
  //
  //  Get the dPos'th bit, shift it onto the return value.

  for (uint32 ii=0; ii < _sWidth; ii++, dPos++) {
    uint64  seg       = dPos / _segSize;     //  Which segment are we in?
    uint64  segPos    = dPos % _segSize;     //  Bit position in that segment.
    uint32  word      = segPos / 64;         //  Which word are we in=?
    uint32  wordBgn   = segPos % 64;         //  Bit position in that word.
    uint32  wordEnd   = wordBgn + _sWidth;

    assert(_segments[seg] != nullptr);

    uint64  w = _segments[seg][word];

    w <<= wordBgn;
    w >>= 63;

    ret <<= 1;
    ret  |= w;
  }

  return(ret);
}
