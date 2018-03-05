
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
//#include "libbits.H"

//  merSize 1 NOT supported.  Fails _leftShift.

class  stuffedBits;


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
    fprintf(stderr, "setSize(" F_U32 ") " F_U32 " " F_U32 " 0x%016" F_X64P " " F_U32 "\n",
            ms, _merSize, _merSpan, _mask, _leftShift);
  };

  static
  uint32      merSize(void) { return(_merSize); };

  void        addR(char base)       { _mer  = (((_mer << 2) & _mask) | ((base & 0x06) >> 1));           };
  void        addL(char base)       { _mer  = (((_mer >> 2)        ) | ((base & 0x06) << _leftShift));  };

  //void        mask(uint32 merSize)  { _mer &= _mask;  };

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


typedef kmerTiny kmer;







class kmerCountStatistics {
public:
  kmerCountStatistics();
  ~kmerCountStatistics();

  void      addCount(uint64 count) {

    if (count == 0)
      return;

    if (count == 1)
      _numUnique++;

    _numDistinct += 1;
    _numTotal    += count;

    if (count < _histMax) {
      _hist[count]++;
      return;
    }
  };

  void      dump(stuffedBits *bits);
  void      dump(FILE        *outFile);

  void      load(stuffedBits *bits);
  void      load(FILE        *inFile);


private:
  uint64              _numUnique;
  uint64              _numDistinct;
  uint64              _numTotal;

  uint32              _histMax;    //  Max count that can be stored in _hist.
  uint64             *_hist;

  uint32              _hbigLen;    //  Counts bigger than _histMax are stored
  uint32              _hbigMax;    //  as unsorted arrays.
  uint64             *_hbigCount;
  uint64             *_hbigNumber;
};



class kmerCountFileReader {
public:
  kmerCountFileReader(const char *inputName);
  ~kmerCountFileReader();

  bool    nextMer(void);
  kmer    theFMer(void)   { return(_kmer);    };
  uint32  theCount(void)  { return(_count);   };

private:
  char                _inName[FILENAME_MAX+1];

  uint32              _prefixSize;
  uint32              _suffixSize;
  uint32              _merSize;
  uint32              _numFiles;

  kmerCountStatistics _stats;

  uint32              _datFilesLen;
  FILE              **_datFiles;

  kmer                _kmer;
  uint32              _count;

  uint64              _prefix;

  uint32              _activeMer;
  uint32              _activeFile;

  uint32              _nKmers;
  uint32              _nKmersMax;
  uint64             *_suffixes;
  uint32             *_counts;
};



class kmerCountFileWriter {
public:
  kmerCountFileWriter(const char *outputName,
                      uint32      splitting,
                      uint32      merSize,
                      uint32      prefixSize = 0,
                      uint32      suffixSize = 0);

  ~kmerCountFileWriter();

  void    addMer(kmer k, uint32 c);

  void    addBlock(uint64  prefix,
                   uint64  nKmers,
                   uint64 *suffixes,
                   uint32 *counts);

private:
  char                _outName[FILENAME_MAX+1];

  uint32              _prefixSize;
  uint32              _suffixSize;
  uint32              _merSize;
  uint32              _numFiles;
  FILE              **_datFiles;
  pthread_mutex_t    *_locks;

  kmerCountStatistics _stats;

};



#endif  //  LIBKMER
