
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

#ifndef LIBBITS_H
#define LIBBITS_H

#include "AS_global.H"

#include <algorithm>

//  Define this to enable testing that the width of the data element is greater than zero.  The
//  uint64MASK() macro (bri.h) does not generate a mask for 0.  Compiler warnings are issued,
//  because you shouldn't use this in production code.
//
//  As it's expensive, emit a warning if it's enabled.
//
//#define CHECK_WIDTH
//#warning CHECK_WIDTH is EXPENSIVE


//  Writing in the middle of data is toublesome.
//  This code currently will not split an object across two allocated blocks,
//  so if you want to rewrite in the middle, you need to make sure it
//  doesn't span a previously written block.  E.g., two writes
//  of 32 bits each could be in different blocks, and trying
//  to overwrite with a single 64-bit write could result in the first
//  block getting truncated (to the current position) and the write
//  replacing the first 64-bits in the second block, not just the 32-bits expected.
//
//  And don't even think of overwriting any of the variable length data.


inline
char *
displayWord(uint64 word, char *buffer=NULL) {
  static char b[65];

  if (buffer == NULL)
    buffer = b;

  memset(buffer, 'x', 64);

  for (uint32 ii=0; ii<64; ii++)
    buffer[ii] = (word & (uint64)1 << (63 - ii)) ? '1' : '0';
  
  buffer[64] = 0;

  return(buffer);
};


//  Return bits in a word:
//    Keeping the rightmost 64-n bits (mask out the leftmost  n bits)
//    Keeping the leftmost  64-n bits (mask out the rigthmost n bits)
//
inline uint64   clearLeftBits  (uint64 v, uint32 l)  {  if (l >= 64) return(0);  return(v & (0xffffffffffffffffllu >>       l));   };
inline uint64    saveLeftBits  (uint64 v, uint32 l)  {  if (l ==  0) return(0);  return(v & (0xffffffffffffffffllu << (64 - l)));  };
inline uint64   clearRightBits (uint64 v, uint32 r)  {  if (r >= 64) return(0);  return(v & (0xffffffffffffffffllu <<       r));   };
inline uint64    saveRightBits (uint64 v, uint32 r)  {  if (r ==  0) return(0);  return(v & (0xffffffffffffffffllu >> (64 - r)));  };

inline uint64   clearMiddleBits(uint64 v, uint32 l, uint32 r) {
  return(saveRightBits(v, r) | saveLeftBits(v, l));
};

inline uint64    saveMiddleBits(uint64 v, uint32 l, uint32 r) {
  return(clearRightBits(v, r) & clearLeftBits(v, l));
};



//  Freed, Edwin E. 1983. "Binary Magic Number" Dr. Dobbs Journal Vol. 78 (April) pp. 24-37
//    Reverse the bits in a word,
//    Count the number of set bits in a words
//    and more.
//
inline
uint64
reverseBits64(uint64 x) {
  x = ((x >>  1) & uint64NUMBER(0x5555555555555555)) | ((x <<  1) & uint64NUMBER(0xaaaaaaaaaaaaaaaa));
  x = ((x >>  2) & uint64NUMBER(0x3333333333333333)) | ((x <<  2) & uint64NUMBER(0xcccccccccccccccc));
  x = ((x >>  4) & uint64NUMBER(0x0f0f0f0f0f0f0f0f)) | ((x <<  4) & uint64NUMBER(0xf0f0f0f0f0f0f0f0));
  x = ((x >>  8) & uint64NUMBER(0x00ff00ff00ff00ff)) | ((x <<  8) & uint64NUMBER(0xff00ff00ff00ff00));
  x = ((x >> 16) & uint64NUMBER(0x0000ffff0000ffff)) | ((x << 16) & uint64NUMBER(0xffff0000ffff0000));
  x = ((x >> 32) & uint64NUMBER(0x00000000ffffffff)) | ((x << 32) & uint64NUMBER(0xffffffff00000000));
  return(x);
}

inline
uint32
reverseBits32(uint32 x) {
  x = ((x >>  1) & uint32NUMBER(0x55555555)) | ((x <<  1) & uint32NUMBER(0xaaaaaaaa));
  x = ((x >>  2) & uint32NUMBER(0x33333333)) | ((x <<  2) & uint32NUMBER(0xcccccccc));
  x = ((x >>  4) & uint32NUMBER(0x0f0f0f0f)) | ((x <<  4) & uint32NUMBER(0xf0f0f0f0));
  x = ((x >>  8) & uint32NUMBER(0x00ff00ff)) | ((x <<  8) & uint32NUMBER(0xff00ff00));
  x = ((x >> 16) & uint32NUMBER(0x0000ffff)) | ((x << 16) & uint32NUMBER(0xffff0000));
  return(x);
}


inline
uint64
uint64Swap(uint64 x) {
  x = ((x >>  8) & uint64NUMBER(0x00ff00ff00ff00ff)) | ((x <<  8) & uint64NUMBER(0xff00ff00ff00ff00));
  x = ((x >> 16) & uint64NUMBER(0x0000ffff0000ffff)) | ((x << 16) & uint64NUMBER(0xffff0000ffff0000));
  x = ((x >> 32) & uint64NUMBER(0x00000000ffffffff)) | ((x << 32) & uint64NUMBER(0xffffffff00000000));
  return(x);
}

inline
uint32
uint32Swap(uint32 x) {
  x = ((x >>  8) & uint32NUMBER(0x00ff00ff)) | ((x <<  8) & uint32NUMBER(0xff00ff00));
  x = ((x >> 16) & uint32NUMBER(0x0000ffff)) | ((x << 16) & uint32NUMBER(0xffff0000));
  return(x);
}

inline
uint16
uint16Swap(uint16 x) {
  x = ((x >>  8) & 0x00ff) | ((x <<  8) & 0xff00);
  return(x);
}


inline
uint32
countNumberOfSetBits32(uint32 x) {
  x = ((x >>  1) & uint32NUMBER(0x55555555)) + (x & uint32NUMBER(0x55555555));
  x = ((x >>  2) & uint32NUMBER(0x33333333)) + (x & uint32NUMBER(0x33333333));
  x = ((x >>  4) & uint32NUMBER(0x0f0f0f0f)) + (x & uint32NUMBER(0x0f0f0f0f));
  x = ((x >>  8) & uint32NUMBER(0x00ff00ff)) + (x & uint32NUMBER(0x00ff00ff));
  x = ((x >> 16) & uint32NUMBER(0x0000ffff)) + (x & uint32NUMBER(0x0000ffff));
  return(x);
}

inline
uint64
countNumberOfSetBits64(uint64 x) {
  x = ((x >>  1) & uint64NUMBER(0x5555555555555555)) + (x & uint64NUMBER(0x5555555555555555));
  x = ((x >>  2) & uint64NUMBER(0x3333333333333333)) + (x & uint64NUMBER(0x3333333333333333));
  x = ((x >>  4) & uint64NUMBER(0x0f0f0f0f0f0f0f0f)) + (x & uint64NUMBER(0x0f0f0f0f0f0f0f0f));
  x = ((x >>  8) & uint64NUMBER(0x00ff00ff00ff00ff)) + (x & uint64NUMBER(0x00ff00ff00ff00ff));
  x = ((x >> 16) & uint64NUMBER(0x0000ffff0000ffff)) + (x & uint64NUMBER(0x0000ffff0000ffff));
  x = ((x >> 32) & uint64NUMBER(0x00000000ffffffff)) + (x & uint64NUMBER(0x00000000ffffffff));
  return(x);
}


inline
uint32
logBaseTwo32(uint32 x) {
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  return(countNumberOfSetBits32(x));
}

inline
uint64
logBaseTwo64(uint64 x) {
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  x |= x >> 32;
  return(countNumberOfSetBits64(x));
}




//#if (__GNUC__ > 3) || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4)
//#define BUILTIN_POPCOUNT
//#endif

#ifdef BUILTIN_POPCOUNT

inline
uint32
countNumberOfSetBits32(uint32 x) {
  return(__builtin_popcount(x));
}

inline
uint64
countNumberOfSetBits64(uint64 x) {
  return(__builtin_popcountll(x));
}

#endif





class bitArray {
public:
  bitArray(uint64 maxNumBits=0) {
    _maxBitSet   = 0;
    _maxBitAvail = maxNumBits;
    _bits        = NULL;

    if (maxNumBits > 0)
      allocate();
  };

  ~bitArray(void) {
    delete [] _bits;
  };

  void     allocate(void) {
    if (_bits != NULL)
      return;

    _bits = new uint64 [_maxBitAvail / 64 + 1];

    clear();
  };

  void     clear(void) {
    memset(_bits, 0, sizeof(uint64) * (_maxBitAvail / 64 + 1));
  };

  bool     getBit(uint64 position) {
    uint64   w =      (position / 64);
    uint64   b = 63 - (position % 64);

    return((_bits[w] >> b) & 0x00000001);
  };

  void     setBit(uint64 position, bool value) {
    uint64   w =      (position / 64);
    uint64   b = 63 - (position % 64);
    uint64   m = ((uint64)1) << b;

    //fprintf(stderr, "SET pos %9" F_U64P " word %2" F_U64P " bit %2" F_U64P " value %c  0x%016" F_X64P " -> ",
    //        position, w, b, (value) ? '1' : '0', _bits[w]);

    _bits[w] &= ~m;
    _bits[w] |= ((uint64)value) << b;

    //fprintf(stderr, "0x%016" F_X64P "\n", _bits[w]);
  };

  bool     flipBit(uint64 position) {
    uint64   w =      (position / 64);
    uint64   b = 63 - (position % 64);
    uint64   m = ((uint64)1) << b;
    uint64   v = _bits[w] & m;

    _bits[w] ^= v;

    return(v >> b);
  };

private:
  uint64   _maxBitSet;
  uint64   _maxBitAvail;
  uint64  *_bits;
};




class stuffedBits {
public:
  stuffedBits(uint64 nBits=16 * 1024 * 1024 * 8);
  stuffedBits(const char *inputName);
  stuffedBits(FILE *inFile);
  ~stuffedBits();

  //  Files.

  void     dumpToFile(FILE *F);
  bool     loadFromFile(FILE *F);

  //  Management of the read/write head.

  void     setPosition(uint64 position, uint64 length = 0);
  uint64   getPosition(void);

  void     byteAlign(void);

  //  UNARY CODED DATA

  uint64   getUnary(void);
  uint64  *getUnary(uint64 number, uint64 *values);

  void     setUnary(uint64 value);
  void     setUnary(uint64 number, uint64 *values);

  //  BINARY CODED DATA

  uint64   getBinary(uint32 width);
  uint64  *getBinary(uint32 width, uint64 number, uint64 *values=NULL);

  void     setBinary(uint32 width, uint64 value);
  void     setBinary(uint32 width, uint64 number, uint64 *values);

  //  GENERALIZED UNARY CODED DATA

  uint64   getGeneralizedUnary(void);
  uint64  *getGeneralizedUnary(uint64 number, uint64 *values=NULL);

  void     setGeneralizedUnary(uint64 value);
  void     setGeneralizedUnary(uint64 number, uint64 *values);



  //  Old meryl uses preDecrement() when using compressed bucket counting.
  //  Nothing else in canu uses these, and they're painful, so left unimplemented.
#if 0
  uint64   preIncrementBinary(uint64 width, uint64 position);
  uint64   postIncrementBinary(uint64 width, uint64 position);
  uint64   preDecrementBinary(uint64 width, uint64 position);
  uint64   postDecrementBinary(uint64 width, uint64 position);
#endif

  //  ELIAS GAMMA CODED DATA

  //  ELIAS DELTA CODED DATA

  //  FIBONACCI CODED DATA

private:

  //  For writing, update the length of the block to the maximum of where we're at now and the existing length.
  //
  void      updateLen(void) {
    _dataBlockLen[_dataBlk] = std::max(_dataPos, _dataBlockLen[_dataBlk]);
  };

  //  For both reading and writing, move to the next word if we're at the end of the current one.
  //
  void      updateBit(void) {
    if (_dataBit == 0) {
      _dataWrd += 1;
      _dataBit  = 64;
    }
  };

  //  For reading operations, move to the next block if we're at the end of the current one.
  //  For writing operations, this is done before the write, in ensureSpace().
  //
  //  Should be done before any reading operation.  It isn't (strictly) needed at the end
  //  of a read.  The next read will just do it at the start, and the next write
  //  only cares that there is enough space in the current block for the write.
  //
  void      updateBlk(uint64 readLength) {

    assert(_dataBit != 0);
    assert(_dataBit <= 64);

    //  If in this block, nothing to do.

    //fprintf(stderr, "updateBlk()--  at _dataPos " F_U64P "  _dataBlockLen[%u] = " F_U64P "\n", _dataPos, _dataBlk, _dataBlockLen[_dataBlk]);

    if (_dataPos + readLength <= _dataBlockLen[_dataBlk])
      return;

    //  Otherwise, we MUST be at the end of the current block.  If not,
    //  whatever we're trying to read is in the next block (or we're trying
    //  to read something longer than what is here).

    assert(_dataPos == _dataBlockLen[_dataBlk]);

    _dataBlk += 1;

    assert(_dataBlk < _dataBlocksLen);

    _dataPos  = 0;
    _data     = _dataBlocks[_dataBlk];

    _dataWrd  = 0;
    _dataBit  = 64;
  }

  //  For writing operations, make sure there is enough space for the write in this block.
  //
  void     ensureSpace(uint64 spaceNeeded) {

    assert(_dataBit != 0);
    assert(_dataBit <= 64);

    //  If enough space in the current block, just return.

    if (_dataPos + spaceNeeded < _dataBlockLenMax)
      return;

    //  Othwerwise, terminate the current block.

    _dataBlockLen[_dataBlocksLen - 1] = _dataPos;

    //  Move to the new block.

    _dataBlocksLen++;

    if (_dataBlocksLen >= _dataBlocksMax) {
      setArraySize(_dataBlocks,   _dataBlocksLen, _dataBlocksMax, _dataBlocksLen + 128);
      setArraySize(_dataBlockBgn, _dataBlocksLen, _dataBlocksMax, _dataBlocksLen + 128);
      setArraySize(_dataBlockLen, _dataBlocksLen, _dataBlocksMax, _dataBlocksLen + 128);
    }

    assert(spaceNeeded <= _dataBlockLenMax);

    _dataPos = 0;
    _data    = _dataBlocks[_dataBlocksLen - 1] = new uint64 [_dataBlockLenMax / sizeof(uint64)];

    memset(_data, 0, sizeof(uint64) * _dataBlockLenMax / 64);

    _dataBlockBgn[_dataBlocksLen - 1] = _dataBlockBgn[_dataBlocksLen - 2] + _dataBlockLen[_dataBlocksLen - 2];
    _dataBlockLen[_dataBlocksLen - 1] = 0;

    _dataBlk += 1;
    _dataWrd  = 0;
    _dataBit  = 64;
  };


  uint64   _dataBlockLenMax; //  Allocated length, of each block (in BITS).

  uint32   _dataBlocksLen;   //  Number of allocated data blocks.
  uint32   _dataBlocksMax;   //  Number of blocks we can allocate.

  uint64  *_dataBlockBgn;    //  Starting position, in the global file, of this block.
  uint64  *_dataBlockLen;    //  Length of this block.
  uint64 **_dataBlocks;      //  Just piles of bits.  Nothing interesting here.

  uint64   _dataPos;         //  Position in this block, in BITS.
  uint64  *_data;            //  Pointer to the currently active data block.

  uint64   _dataBlk;         //  Active data block.
  uint64   _dataWrd;         //  Active word in the active data block.
  uint64   _dataBit;         //  Active bit in the active word in the active data block (aka, number of bits left in this word)


  FILE    *_inputFile;
};



#endif  //  LIBBITS_H
