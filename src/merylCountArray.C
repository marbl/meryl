
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

#include "meryl.H"



uint32  merylCountArray::_width   = 0;
uint32  merylCountArray::_segSize = 8192 * 64;



merylCountArray::merylCountArray(void) {
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
}



merylCountArray::~merylCountArray() {
  delete [] _suffix;
  delete [] _counts;

  removeSegments();
}
