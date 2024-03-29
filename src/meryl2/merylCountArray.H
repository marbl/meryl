
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

#ifndef MERYLCOUNTARRAY_H
#define MERYLCOUNTARRAY_H

#include "meryl.H"

class swv;   //  Defined and only used in merylCountArray.C


//  The basic meryl counting object.  It collects a list of kmers from the
//  input, sorts them, the returns a list of distinct kmers with the number
//  of times each occurs in the list.
//
//  A kmer is split into a prefix and a suffix.  All kmers in one
//  merylCountArray object will share the same prefix (the client needs to do
//  this for us).  The merylCountArray will save the lower _sWidth bits of
//  each input kmer.
//
//  If the kmer comes with a count or color, those are saved too.  Otherwise,
//  only the kmer suffix bits are saved.
//
//  Since we don't know how many kmers will be input to us, and don't want to
//  keep reallocating data, our list is grown in segments of size _segSize
//  (segsize in initialize()) kilo bytes.

#undef ADD_VERBOSE
#define ADD_INSTRUMENT

class merylCountArray {
public:
  merylCountArray() {
  };
  ~merylCountArray();

  uint64    initialize(uint64 prefix, uint32 width);

  void      initializeValues(uint32 valueWidth, uint32 labelWidth);

  void      enableMultiSet(bool enable)    {  _multiSet = enable;  };

public:
  void      initializeForTesting(uint32 width, uint32 nwords);
  void      dumpData(void);

private:
  void      removeSegments(void);
  void      removeValues(void);
  void      addSegment(uint32 seg);

  //  Count the number of times we hit each case, and where the word starts in that case
public:
  void      clearStats(void);
  void      dumpStats(void);

#ifdef ADD_INSTRUMENT
  uint64    nTests[6];
  uint64    nStart[6][64];
#endif

  //  Add a suffix to the table.
  //
public:
  uint64    add(kmdata suffix);
  uint64    addValue(kmvalu value);
  uint64    addLabel(kmlabl label);

private:
  kmdata   *unpackSuffixes(uint64 nSuffixes);
  swv      *unpackSuffixesAndValues(uint64 nSuffixes);

  //
  //  Return the kkth kmer suffix stored in the array.  This is only used in sort(),
  //  and only to convert the bit-packed _seg data into unpacked words, so could
  //  be optimized for that case.  I don't expect much of a performance gain.
  //
public:
  kmdata    get(uint64 kk);
  kmdata    getSimple(uint64 kk);


public:
  uint64           numBits(void)        {  return(_nBits);           };
  uint64           numKmers(void)       {  return(_nBits / _sWidth); };



public:
  //  Using 1 here is probably not the most time efficient value, but the
  //  memory usage estimate seems to be the most accurate with it.
  //
  static
  uint32            pagesPerSegment(void)    { return(1); };

  //  We're doing accounting ourself instead of asking the OS for the current
  //  process size because some OSs (FreeBSD, probably MacOS) don't decrease
  //  the size and we need to when these tables are too full opposed to
  //  simply allocated.
  //
  //  If the memset in merylCountArray() is enabled, this calculation does
  //  not represent the amount of resident memory.
  //
  //  The number of pages used for data is complicated.  We're allocating in
  //  blocks of pagesPerSegment() but reserving a few words for OS overhead.
  //
  //    fSegmsUsed: the number of full segments allocated, _segSize bits in each.
  //    pPagesUsed: the leftover bits, plus one partial page used
  //
  uint64           usedSize(void) {
    uint64 memUsed    = 0;
 
    uint64 fSegmsUsed = _nBits / _segSize;
    uint64 pPagesUsed = (_nBits - fSegmsUsed * _segSize) / _bitsPerPage + 1;

    uint64 pagesUsed  = fSegmsUsed * pagesPerSegment() + pPagesUsed;

    memUsed += sizeof(merylCountArray);         //  For our metadata
    memUsed += pagesUsed * _bitsPerPage / 8;    //  For the packed kmer data
    memUsed += _segAlloc * sizeof(uint64 **);   //  For pointers to segments

    return(memUsed);
  };

  //  Returns the change in size since the last call, but sets a threshold so
  //  we don't spend a bunch of time calling usedSize().
  //
  uint64           usedSizeDelta(void) {

    if (_nBits < _nBitsTrigger)
      return(0);

    uint64  newSize   = usedSize();
    uint64  sizeDelta = newSize - _nBitsOldSize;

    _nBitsTrigger += _bitsPerPage / 16;
    _nBitsOldSize  = newSize;

    return(sizeDelta);
  };


private:
  void             countSingleKmers(void);
  void             countSingleKmersWithValues(void);
  void             countMultiSetKmers(void);
public:
  void             countKmers(void);
  void             dumpCountedKmers(merylBlockWriter *out, kmlabl label);
  void             removeCountedKmers(void);


private:
  uint32           _sWidth = 0;         //  Size of the suffix we're storing
  uint32           _vWidth = 0;         //  Size of the values we're storing
  uint32           _lWidth = 0;         //  Size of the labels we're storing

  uint64           _prefix = 0;         //  The kmer prefix we're storing data for
  kmdata          *_suffix = nullptr;   //  After sorting, the suffix of each kmer
  kmvalu          *_counts = nullptr;   //  After sorting, the number of times we've seen this kmer
  kmlabl          *_labels = nullptr;   //  After sorting, the label associated with this kmer

  uint64           _nKmers = 0;         //  Number of kmers.

  uint64           _bitsPerPage = 0;
  uint64           _nReAlloc    = 0;

private:
  uint32           _segSize  = 0;       //  Number of bits in each segment.
  uint32           _segAlloc = 0;       //  Number of segments we're allowed to allocate  (size of the array below).
  uint64         **_segments = nullptr; //  An array of blocks of data.

  stuffedBits     *_vals = nullptr;     //  An array of values.
  stuffedBits     *_labs = nullptr;     //  An array of labels.

  uint64           _nBits        = 0;   //  Number of bits stored.
  uint64           _nBitsTrigger = 0;   //  Number of bits we need to store for a size recalculation to occur.
  uint64           _nBitsOldSize = 0;   //  Last computed size.

  bool             _multiSet = false;   //  Treat the input kmers as a multiset.
};



#endif  //  MERYLCOUNTARRAY_H
