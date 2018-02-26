
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

#include "libbits.H"

#include "AS_UTL_fileIO.H"


stuffedBits::stuffedBits(uint64 nBits) {
  _dataBlocksLen = 1;
  _dataBlocksMax = 64;
  _dataBlocks    = new uint64 * [_dataBlocksMax];

  _dataBlockBgn  = new uint64   [_dataBlocksMax];
  _dataBlockLen  = new uint64   [_dataBlocksMax];
  _dataBlockMax  = nBits;

  _dataPos = 0;
  _data    = _dataBlocks[0] = new uint64 [_dataBlockMax / 64];

  memset(_data, 0, sizeof(uint64) * _dataBlockMax / 64);

  _dataBlockBgn[0] = 0;
  _dataBlockLen[0] = 0;

  _dataBlk = 0;
  _dataWrd = 0;
  _dataBit = 64;
};


stuffedBits::stuffedBits(char const *name) {
};


stuffedBits::~stuffedBits() {

  fprintf(stderr, "stuffedBits()--  Deleting %u blocks storing %lu bits.\n",
          _dataBlocksLen,
          _dataBlockBgn[_dataBlocksLen - 1] + _dataBlockLen[_dataBlocksLen - 1]);

  for (uint32 ii=0; ii<_dataBlocksLen; ii++)
    delete [] _dataBlocks[ii];

  delete [] _dataBlocks;
  delete [] _dataBlockBgn;
  delete [] _dataBlockLen;
};



void
stuffedBits::dumpToFile(FILE *F) {

  AS_UTL_safeWrite(F, &_dataBlocksLen, "dataBlocksLen", sizeof(uint32), 1);
  AS_UTL_safeWrite(F, &_dataBlockMax,  "dataBlocksLen", sizeof(uint64), 1);
  AS_UTL_safeWrite(F,  _dataBlockBgn,  "dataBlocksLen", sizeof(uint64), 1);
  AS_UTL_safeWrite(F,  _dataBlockLen,  "dataBlocksLen", sizeof(uint64), 1);

  for (uint32 ii=0; ii<_dataBlocksLen; ii++)
    AS_UTL_safeWrite(F, _dataBlocks[ii], "dataBlocks", sizeof(uint64), _dataBlockLen[ii] / 64 + 1);
}



void
stuffedBits::loadFromFile(FILE *F) {

  AS_UTL_safeRead(F, &_dataBlocksLen, "dataBlocksLen", sizeof(uint32), 1);
  AS_UTL_safeRead(F, &_dataBlockMax,  "dataBlocksLen", sizeof(uint64), 1);

  _dataBlocks    = new uint64 * [_dataBlocksLen];
  _dataBlockBgn  = new uint64   [_dataBlocksLen];
  _dataBlockLen  = new uint64   [_dataBlocksLen];

  AS_UTL_safeRead(F,  _dataBlockBgn,  "dataBlocksLen", sizeof(uint64), 1);
  AS_UTL_safeRead(F,  _dataBlockLen,  "dataBlocksLen", sizeof(uint64), 1);

  for (uint32 ii=0; ii<_dataBlocksLen; ii++) {
    _dataBlocks[ii] = new uint64 [_dataBlockLen[ii] / 64 + 1];

    AS_UTL_safeRead(F, _dataBlocks[ii], "dataBlocks", sizeof(uint64), _dataBlockLen[ii] / 64 + 1);
  }
}



//  Set the position of stuffedBits to 'position'.
//  Ensure that at least 'length' bits exist in the current block.
//
void
stuffedBits::setPosition(uint64 position, uint64 length) {

  _dataBlk = 0;

  while ((position < _dataBlockBgn[_dataBlk]) && (_dataBlk < _dataBlocksLen))
    _dataBlk++;

  assert(_dataBlk < _dataBlocksLen);  //  What to do if we seek to an uninitialized piece?

  _dataPos = position - _dataBlockBgn[_dataBlk];
  _data    = _dataBlocks[_dataBlk];

  _dataWrd =      _dataPos / 64;
  _dataBit = 64 - _dataPos % 64;
}

uint64
stuffedBits::getPosition(void) {
  return(_dataBlockBgn[_dataBlk] + _dataPos);
}



uint64
stuffedBits::getUnary(void) {
  uint64  value = 0;
  uint64  wrd;

  //  Ensure we're in valid data.

  updateBlk(1);

  //  Word align us first.

  wrd = _data[_dataWrd] << (64 - _dataBit);

  //  Skip entire words.  For the first word, if we're left with only zeros
  //  after the shifting, we increase the 'value' by the number of bits
  //  we could read in the word, not the full 64 bits that are zero.

  while (wrd == 0) {
    value    += _dataBit;

    _dataPos += _dataBit;
    _dataWrd += 1;

    _dataBit  = 64;
    wrd       = _data[_dataWrd];
  }

  //  Decode the last partial word.

  wrd       = 64 - logBaseTwo64(wrd);

  value    += wrd;

  _dataPos += wrd + 1;
  _dataBit -= wrd + 1;

  updateBit();

  return(value);
}

uint64 *
stuffedBits::getUnary(uint64 number, uint64 *values) {

  for (uint64 ii=0; ii<number; ii++)
    values[ii] = getUnary();

  return(values);
}




void
stuffedBits::setUnary(uint64 value) {

  ensureSpace(value+1);

  //  If we fit entirely within this word, handle it special.

  if (value + 1 < _dataBit) {
    _data[_dataWrd] = clearMiddleBits(_data[_dataWrd], 64 - _dataBit, _dataBit + value + 1);

    _dataPos += value + 1;
    _dataBit -= value + 1;

    _data[_dataWrd] |= (uint64)1 << _dataBit;

    updateLen();

    return;
  }

  //  We fit _exactly_ in this word, special again!

  if (value + 1 == _dataBit) {
    _data[_dataWrd]  = clearRightBits(_data[_dataWrd], _dataBit);
    _data[_dataWrd] |= 1;   //  ALWAYS the last bit.

    _dataPos += value + 1;  //  ALWAYS move to the next word.
    _dataWrd += 1;
    _dataBit  = 64;

    updateLen();

    return;
  }

  //  Otherwise, we span at least two words.  First, get us word aligned,
  //  by clearing the rest of this word.

  assert(value + 1 > _dataBit);

  _data[_dataWrd] = clearRightBits(_data[_dataWrd], _dataBit);

  value    -= _dataBit;

  _dataPos += _dataBit;
  _dataWrd += 1;
  _dataBit  = 64;

  assert((_dataPos % 64) == 0);

  //  Then, set entire words to zero.

  while (value >= 64) {
    _data[_dataWrd] = 0;

    value    -= 64;

    _dataPos += 64;
    _dataWrd += 1;
    _dataBit  = 64;
  }

  //  Finally, set the last partial word.  This is similar to the cases
  //  at the start, but we know the bits will always be on the left of the word.

  _data[_dataWrd] = clearLeftBits(_data[_dataWrd], value + 1);

  _dataPos += value + 1;                      //  Skip the zero bits.
  _dataBit -= value + 1;                      //  (and the sentinel)

  _data[_dataWrd] |= (uint64)1 << _dataBit;   //  Add the sentinel.

  //  Update the pointers.

  updateLen();
  updateBit();
}



void
stuffedBits::setUnary(uint64 number, uint64 *values) {

  for (uint64 ii=0; ii<number; ii++)
    setUnary(values[ii]);
}




uint64
stuffedBits::getGeneralizedUnary(void) {
  return(0);
}



uint64 *
stuffedBits::getGeneralizedUnary(uint64 number, uint64 *values) {
  return(values);
}




void
stuffedBits::setGeneralizedUnary(uint64 value) {
}



void
stuffedBits::setGeneralizedUnary(uint64 number, uint64 *values) {
}




uint64
stuffedBits::getBinary(uint32 width) {
  uint64  value = 0;

  assert(width < 65);

  //  Ensure we're in valid data.

  updateBlk(width);

  //  If we're contained in a single word, special case.

  if (width < _dataBit) {
    value = saveRightBits(_data[_dataWrd] >> (_dataBit - width), width);

    _dataPos += width;
    _dataBit -= width;
  }

  //  If we're exactly in a single word, another special case.

  else if (width == _dataBit) {
    value = saveRightBits(_data[_dataWrd], width);

    _dataPos += width;
    _dataWrd += 1;
    _dataBit  = 64;
  }

  //  Otherwise, we're spanning two words.

  else {
    uint64  w1 =         _dataBit;
    uint64  w2 = width - _dataBit;

    uint64  l  = saveRightBits(_data[_dataWrd + 0], w1) <<       w2;
    uint64  r  = saveLeftBits (_data[_dataWrd + 1], w2) >> (64 - w2);

    value = l | r;

    _dataPos += width;
    _dataWrd += 1;
    _dataBit  = 64 - w2;
  }

  return(value);
}



uint64 *
stuffedBits::getBinary(uint32 width, uint64 number, uint64 *values) {

  for (uint64 ii=0; ii<number; ii++)
    values[ii] = getBinary(width);

  return(values);
}




void
stuffedBits::setBinary(uint32 width, uint64 value) {

  assert(width < 65);

  ensureSpace(width);

  //  Mask off any pieces we're not supposed to be seeing.

  value = saveRightBits(value, width);

  //  If we fit entirely within this word, handle it special.

  if (width < _dataBit) {
    _data[_dataWrd] = clearMiddleBits(_data[_dataWrd], 64 - _dataBit, _dataBit + width) | (value << (_dataBit - width));

    _dataPos += width;
    _dataBit -= width;
  } else

  //  We fit _exactly_ in this word, special again!

  if (width == _dataBit) {
    _data[_dataWrd] = clearRightBits(_data[_dataWrd], _dataBit) | value;

    _dataPos += width;
    _dataWrd += 1;
    _dataBit  = 64;
  }

  //  Otherwise, we span two words.

  else {
    uint32  w1 =         _dataBit;
    uint32  w2 = width - _dataBit;

    _data[_dataWrd + 0] = clearRightBits(_data[_dataWrd + 0], w1) | (value >> (     w2));
    _data[_dataWrd + 1] = clearLeftBits (_data[_dataWrd + 1], w2) | (value << (64 - w2));

    _dataPos += width;
    _dataWrd += 1;
    _dataBit  = 64 - w2;
  }

  //  updateBit() isn't needed; it's handled in the special cases.

  updateLen();
}



void
stuffedBits::setBinary(uint32 width, uint64 number, uint64 *values) {

  for (uint64 ii=0; ii<number; ii++)
    setBinary(width, values[ii]);
}



