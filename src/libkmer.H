
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

#ifndef LIBKMER_H
#define LIBKMER_H

#include "AS_global.H"

//  merSize 1 NOT supported.  Fails _leftShift.


class  kmerTiny {
public:
  kmerTiny() {
    _mer = 0;
  };

  ~kmerTiny() {
  };

  static
  void        setSize(uint32 ms) {
    _merSize   = ms;
    _merSpan   = ms;
    _mask      = uint64MASK(ms * 2);
    _leftShift = ((2 * ms - 2) % 64) - 1;
    fprintf(stderr, "setSize(%u) %u %u 0x%016lx %u\n", ms, _merSize, _merSpan, _mask, _leftShift);
  };

  static
  uint32      merSize(void) { return(_merSize); };

  void        addR(char base)   { _mer = (((_mer << 2) & _mask) | ((base & 0x06) >> 1));           };
  void        addL(char base)   { _mer = (((_mer >> 2)        ) | ((base & 0x06) << _leftShift));  };

  void        mask(uint32 merSize)         { _mer &= _mask;  };

  kmerTiny   &reverseComplement(void) {
    return(*this);
  };

public:
  bool    operator!=(kmerTiny const &r) const { return(_mer != r._mer); };
  bool    operator==(kmerTiny const &r) const { return(_mer == r._mer); };
  bool    operator< (kmerTiny const &r) const { return(_mer <  r._mer); };
  bool    operator> (kmerTiny const &r) const { return(_mer >  r._mer); };
  bool    operator<=(kmerTiny const &r) const { return(_mer <= r._mer); };
  bool    operator>=(kmerTiny const &r) const { return(_mer >= r._mer); };

public:
  char    *toString(char *str) const {
    for (uint32 ii=0; ii<_merSize; ii++) {
      uint32  bb = (((_mer >> (2 * ii)) & 0x03) << 1);
      str[_merSize-ii-1] = (bb == 0x04) ? ('T') : ('A' + bb);
    }
    str[_merSize] = 0;
    return(str);
  };

  operator uint64 () const {
    return(_mer);
  };

  void     setPrefixSuffix(uint64 prefix, uint64 suffix, uint32 width) {
    _mer  = prefix << width;
    _mer |= suffix;
  };

private:
  void     operator>>=(uint32 x)  { _mer >>= x; };
  void     operator<<=(uint32 x)  { _mer <<= x; };


private:
public:
  uint64         _mer;

  static uint32  _merSize;    //  number of bases in this mer
  static uint32  _merSpan;    //  number of sequence bases this mer spans
  static uint64  _mask;       //  mask to ensure kmer has exactly _merSize bases in it
  static uint32  _leftShift;  //  how far to shift a base to append to the left of the kmer
};

uint32 kmerTiny::_merSize;
uint32 kmerTiny::_merSpan;
uint64 kmerTiny::_mask;
uint32 kmerTiny::_leftShift;


typedef kmerTiny kmer;


class kmerFileWriter {
public:
  kmerFileWriter(const char *filename, uint32 blockBits);
  ~kmerFileWriter();

  

};





class kmerStreamReader {
public:
  kmerStreamReader(const char *fn, uint32 ms=0) {
  };
  ~kmerStreamReader() {
  };

  bool    nextMer(void)   { return(false);      };
  kmer    theFMer(void)   { return(_thisMer);   };
  uint32  theCount(void)  { return(_thisCount); };

  kmer    _thisMer;
  uint32  _thisCount;
};



class kmerStreamWriter {
public:
  kmerStreamWriter(const char *filePrefix,
                   uint32 merSize,          //  In bases
                   uint32 merComp,          //  A length, bases
                   uint32 prefixSize,       //  In bits
                   bool   positionsEnabled) {
  };
  ~kmerStreamWriter() {
  };

  void    addMer(kmer k, uint32 c) {
  };

  void    addBlock(uint64  prefix,
                   uint32  length,
                   uint64 *kmers,
                   uint32 *counts);
};





#endif  //  LIBKMER
