
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
#include "libbits.H"

//  merSize 1 NOT supported.  Fails _leftShift.

#undef SHOW_LOAD
#undef DONT_PRINT


class  kmerTiny {
public:
  kmerTiny() {
    _mer = 0;
  };

  ~kmerTiny() {
  };

  static
  void        setSize(uint32 ms, bool beVerbose=false) {
    _merSize    = ms;

    _fullMask   = uint64MASK(ms * 2);

    _leftMask   = uint64MASK(ms * 2 - 2);
    _leftShift  = ((2 * ms - 2) % 64);

    if (beVerbose)
      fprintf(stderr, "Set global kmer size to " F_U32 " (fullMask=0x%016" F_X64P " leftMask=0x%016" F_X64P " leftShift=" F_U32 ")\n",
              _merSize, _fullMask, _leftMask, _leftShift);
  };

  static
  uint32      merSize(void) { return(_merSize); };

  //  Push an ASCII base onto the mer, shifting the mer to the right or left
  //  to make space for the new base.  Unlike the 'standard' two-bit encoding,
  //  these encode bases as A=00, C=01, G=11, T=10.
  //
  void        addR(char base)       { _mer  = (((_mer << 2) & _fullMask) | (((base >> 1) & 0x03llu)          )              );  };
  void        addL(char base)       { _mer  = (((_mer >> 2) & _leftMask) | (((base >> 1) & 0x03llu) ^ 0x02llu) << _leftShift);  };

  //  Reverse-complementation of a kmer involves complementing the bases in
  //  the mer, revesing the order of all the bases, then aligning the bases
  //  to the low-order bits of the word.
  //
  kmerTiny   &reverseComplement(void) {

    //  Complement the bases

    _mer ^= 0xaaaaaaaaaaaaaaaallu;

    //  Reverse the mer

    _mer = ((_mer >>  2) & 0x3333333333333333llu) | ((_mer <<  2) & 0xccccccccccccccccllu);
    _mer = ((_mer >>  4) & 0x0f0f0f0f0f0f0f0fllu) | ((_mer <<  4) & 0xf0f0f0f0f0f0f0f0llu);
    _mer = ((_mer >>  8) & 0x00ff00ff00ff00ffllu) | ((_mer <<  8) & 0xff00ff00ff00ff00llu);
    _mer = ((_mer >> 16) & 0x0000ffff0000ffffllu) | ((_mer << 16) & 0xffff0000ffff0000llu);
    _mer = ((_mer >> 32) & 0x00000000ffffffffllu) | ((_mer << 32) & 0xffffffff00000000llu);

    //  Shift and mask out the bases not in the mer

    _mer >>= 64 - _merSize * 2;
    _mer  &= _fullMask;

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

  static uint32  _merSize;     //  number of bases in this mer

  static uint64  _fullMask;    //  mask to ensure kmer has exactly _merSize bases in it

  static uint64  _leftMask;    //  mask out the left-most base.
  static uint32  _leftShift;   //  how far to shift a base to append to the left of the kmer
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



class kmerCountFileReaderBlock {
public:
  kmerCountFileReaderBlock() {
    _data       = NULL;

    _prefix     = 0;
    _nKmers     = 0;
    _nKmersMax  = 0;

    _kCode      = 0;
    _unaryBits  = 0;
    _binaryBits = 0;
    _k1         = 0;

    _cCode      = 0;
    _c1         = 0;
    _c2         = 0;

    _suffixes   = NULL;
    _counts     = NULL;
  };

  ~kmerCountFileReaderBlock() {
    delete    _data;
    delete [] _suffixes;
    delete [] _counts;
  };

  bool      loadBlock(uint32 activeFile, uint32 activeIteration, FILE *inFile) {

    //  If _data exists, we've already loaded the block, but haven't used it yet.

    if (_data)
      return(true);

    //  Otherwise, allocate _data, read the block from disk.  If nothing loaded,
    //  return false.

    _data = new stuffedBits();

    _prefix = UINT64_MAX;
    _nKmers = 0;

  loadAgain:
    if (_data->loadFromFile(inFile) == false) { 
#ifdef SHOW_LOAD
      fprintf(stdout, "loadBlock()-- file %u iter %u FAILED TO LOAD.\n",
              activeFile, activeIteration);
#endif
      delete _data;  _data = NULL;
      return(false);
    }

    //  Decode the header of _data, but don't process the kmers yet.

    uint64 m1 = _data->getBinary(64);    //  Magic number, part 1.
    uint64 m2 = _data->getBinary(64);    //  Magic number, part 2.

    if ((m1 != 0x7461446c7972656dllu) ||
        (m2 != 0x0a3030656c694661llu)) {
      fprintf(stderr, "kmerCountFileReader::nextMer()-- Magic number mismatch in activeFile " F_U32 " activeIteration " F_U32 " position " F_U64 ".\n",
              activeFile, activeIteration, _data->getPosition());
      fprintf(stderr, "kmerCountFileReader::nextMer()-- Expected 0x7461446c7972656d got 0x%016" F_X64P "\n", m1);
      fprintf(stderr, "kmerCountFileReader::nextMer()-- Expected 0x0a3030656c694661 got 0x%016" F_X64P "\n", m2);
      exit(1);
    }

    _prefix = _data->getBinary(64);
    _nKmers = _data->getBinary(64);

#ifdef SHOW_LOAD
    fprintf(stdout, "loadBlock()-- file %u iter %u loaded " F_U64 " kmers for prefix 0x%016lx.\n",
            activeFile, activeIteration, _nKmers, _prefix);
#endif

    //  If no kmers are in this block, go back and get another block.  This _should_ be caught
    //  by the writer, and the writer should not be writing empty blocks.  Just in case.

    if (_nKmers == 0)
      goto loadAgain;

    //  Decode the rest of the header.

    _kCode      = _data->getBinary(8);
    _unaryBits  = _data->getBinary(32);
    _binaryBits = _data->getBinary(32);
    _k1         = _data->getBinary(64);

    _cCode      = _data->getBinary(8);
    _c1         = _data->getBinary(64);
    _c2         = _data->getBinary(64);

    return(true);
  };

  //  Decode a the data into suffixes and counts.
  void      decodeBlock() {

    if (_data == NULL)
      return;

    resizeArrayPair(_suffixes, _counts, 0, _nKmersMax, _nKmers, resizeArray_doNothing);

    uint64  thisPrefix = 0;

    for (uint32 kk=0; kk<_nKmers; kk++) {
      thisPrefix += _data->getUnary();

      _suffixes[kk] = (thisPrefix << _binaryBits) | (_data->getBinary(_binaryBits));
    }

    for (uint32 kk=0; kk<_nKmers; kk++)
      _counts[kk] = _data->getBinary(32);

    delete _data;
    _data = NULL;
  };


  uint64    prefix(void)   { return(_prefix); };
  uint64    nKmers(void)   { return(_nKmers); };

  stuffedBits  *_data;

  uint64        _prefix;       //  The prefix of all kmers in this block
  uint64        _nKmers;       //  The number of kmers in this block
  uint64        _nKmersMax;    //  The number of kmers we've allocated space for in _suffixes and _counts

  uint32        _kCode;        //  Encoding type of kmer, then 128 bits of parameters
  uint32        _unaryBits;    //    bits in the unary prefix  (of the kmer suffix)
  uint32        _binaryBits;   //    bits in the binary suffix (of the kmer suffix)
  uint64        _k1;           //    unused

  uint32        _cCode;        //  Encoding type of the counts, then 128 bits of parameters
  uint64        _c1;           //    unused
  uint64        _c2;           //    unused

  uint64       *_suffixes;     //  Decoded suffixes and counts.
  uint32       *_counts;       //
};



class kmerCountFileReader {
public:
  kmerCountFileReader(const char *inputName, bool ignoreStats=true);
  ~kmerCountFileReader();

private:
  void    openFile(uint32 idx, uint32 iteration);

public:
  bool    nextMer(void);
  kmer    theFMer(void)   { return(_kmer);    };
  uint32  theCount(void)  { return(_count);   };

  char   *filename(void)  { return(_inName);  };

private:
  char                _inName[FILENAME_MAX+1];

  uint32              _prefixSize;
  uint32              _suffixSize;
  uint32              _merSize;
  uint32              _numFilesBits;
  uint32              _numFiles;
  uint32              _numIterations;

  kmerCountStatistics _stats;

  FILE                      **_datFiles;
  kmerCountFileReaderBlock  **_blocks;

  kmer                _kmer;
  uint32              _count;

  uint64              _prefix;

  uint64              _activeMer;
  uint32              _activeFile;
  uint32              _activeIteration;

  uint64              _nKmers;
  uint64              _nKmersMax;
  uint64             *_suffixes;
  uint32             *_counts;
};



class kmerCountFileWriter {
public:
  kmerCountFileWriter(const char *outputName,
                      uint32      prefixSize = 0);

  ~kmerCountFileWriter();

public:
  void    initialize(uint32 prefixSize = 0);
  char   *filename(void)  { return(_outName);  };

public:
  void    addMer(kmer k, uint32 c);

  uint32  numberOfFiles(void)           { return(_numFiles); };
  uint64  firstPrefixInFile(uint32 ff);
  uint64  lastPrefixInFile(uint32 ff);

  uint32  fileNumber(uint64 prefix);

  void    addBlock(uint64  prefix,
                   uint64  nKmers,
                   uint64 *suffixes,
                   uint32 *counts);

  void    incrementIteration(void);



private:
  bool                _initialized;

  char                _outName[FILENAME_MAX+1];

  uint32              _iteration;       //  Large datasets can be written in multiple passes.

  uint64              _batchPrefix;     //  Temporary data used when adding
  uint64              _batchNumKmers;   //  kmers one-at-a-time.
  uint64              _batchMaxKmers;
  uint64             *_batchSuffixes;
  uint32             *_batchCounts;

  uint32              _prefixShift;     //  Parameters for converting one-at-a-time
  uint64              _suffixMask;      //  kmers into prefix/suffix kmers.

  uint32              _prefixSize;      //  Parameters for writing a block of kmers.
  uint32              _suffixSize;
  uint32              _merSize;
  uint32              _numFilesBits;
  uint32              _numFiles;
  FILE              **_datFiles;
  pthread_mutex_t    *_locks;

  kmerCountStatistics _stats;
};



#endif  //  LIBKMER
