
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

#include "meryl.H"


merylFilter::merylFilter(merylFilterQuantity type, merylFilterRelation rela, char const *str)   {
  _q = type;
  _r = rela;

  memcpy(_str, str, sizeof(char) * (FILENAME_MAX+1));
}



merylFilter::~merylFilter()  {
  delete [] _presentInNum;
  delete [] _presentInIdx;
  delete [] _presentInList;
}



bool
merylFilter::isTrue(kmer k, uint32 actLen, kmer *act, uint32 *actIdx, uint32 *actRdx) const {
  bool   result = false;

  //  compare() will fail, on purpose, if _r is merylFilterRelation::isNOP.
  //  For filters that never call compare() we don't care what _r is.

  if (_q == merylFilterQuantity::isNOP) {
    assert(0);
  }


  if (_q == merylFilterQuantity::isValue) {
    kmvalu  rhs;
    kmvalu  lhs;

    //  Get the left hand value
    if      (_vIndex1 == uint32max)            lhs = _vValue1;
    else if (_vIndex1 == 0)                    lhs = k._val;
    else if (actRdx[_vIndex1-1] < uint32max)   lhs = act[actRdx[_vIndex1-1]]._val;
    else                                       return(false);

    //  Get the right hand value
    if      (_vIndex2 == uint32max)            rhs = _vValue2;
    else if (_vIndex2 == 0)                    rhs = k._val;
    else if (actRdx[_vIndex2-1] < uint32max)   rhs = act[actRdx[_vIndex2-1]]._val;
    else                                       return(false);

    result = compare(lhs, rhs);

    //fprintf(stderr, "isTrue() isValue = %s\n", result ? "true" : "false");
  }


  if (_q == merylFilterQuantity::isLabel) {
    kmlabl  rhs;
    kmlabl  lhs;

    //  Get the left hand value
    if      (_vIndex1 == uint32max)            lhs = _vLabel1;
    else if (_vIndex1 == 0)                    lhs = k._lab;
    else if (actRdx[_vIndex1-1] < uint32max)   lhs = act[actRdx[_vIndex1-1]]._lab;
    else                                       return(false);

    //  Get the right hand value
    if      (_vIndex2 == uint32max)            rhs = _vLabel2;
    else if (_vIndex2 == 0)                    rhs = k._lab;
    else if (actRdx[_vIndex2-1] < uint32max)   rhs = act[actRdx[_vIndex2-1]]._lab;
    else                                       return(false);

    result = compare(lhs, rhs);

    //fprintf(stderr, "isTrue() isLabel = %s\n", result ? "true" : "false");
  }


  if (_q == merylFilterQuantity::isBases) {             //  The two nonsense cases below
    uint32 c = (((_countA == false) ? 0 : countA(k)) +  //  should be caught by isBasesFilter().
                ((_countC == false) ? 0 : countC(k)) +  //  We'll compute the result anyway, so
                ((_countG == false) ? 0 : countG(k)) +  //  if something does change we still
                ((_countT == false) ? 0 : countT(k)));  //  provide some sensible result.

    if      ((_vIndex1 == 0) && (_vIndex2 == 0))    //  Nonsense, just comparing  ! see comment
      result = compare(_vBases1, _vBases2);         //  the two input constants!  ! above

    else if ((_vIndex1 == 0) && (_vIndex2 != 0))    //  Useful: kmer OP constant
      result = compare(c, _vBases2);                //

    else if ((_vIndex1 != 0) && (_vIndex2 == 0))    //  Useful: constant OP kmer
      result = compare(_vBases1, c);                //

    else                                            //  Nonsense, comparing self  ! see comment
      result = compare(c, c);                       //  to self!                  ! above

    //fprintf(stderr, "isTrue() isBases = %s\n", result ? "true" : "false");
  }


  //  The index filter checks two things:
  //    does the kmer appear in the correct number of inputs
  //    does the kmer appear in ALL specified inputs
  //
  //  The first is a simple array lookup.
  //
  //  The second needs to iterate over _presentInIdx and check that the kmer
  //  is present in each of those inputs.  actRdx[] maps an input index to
  //  an act[] index, or uint32max if that input didn't supply a kmer.

  if (_q == merylFilterQuantity::isIndex) {
    result = _presentInNum[actLen];

    for (uint32 pp=0; pp<_presentInLen; pp++) {
      uint32  ai = actRdx[ _presentInList[pp] ];

      if (ai == uint32max)
        result = false;
    }
  }

  return(result);
}



void
merylFilter::finalizeFilterInputs(uint32 nInputs, std::vector<char const *> &err) {

  //  Allocate space for the lookup tables.

  _presentInNum = new bool [nInputs + 1];   //  +1 because we actually access [nInputs].
  _presentInIdx = new bool [nInputs];

  _presentInList = new uint32 [nInputs];

  //  Initialize defaults.  If nothing was specified, default to allowing
  //  'any' number of inputs, then set the state of the lookup table to
  //  'true' if 'any' is specified or 'false' otherwise.  Finally, if 'all'
  //  was specified, now that we know the number of inputs, we can set that
  //  to 'true'.

  if ((_input_num.size() == 0) && (_input_num_all == false))
    _input_num_any = true;

  for (uint32 ii=0; ii<=nInputs; ii++)
    _presentInNum[ii] = (_input_num_any == true) ? true : false;

  if (_input_num_all == true)
    _presentInNum[nInputs] = true;

  //  Check each of the entries in _input_num and set the flag for each to true.

  for (uint32 ii=0; ii<_input_num.size(); ii++) {
    uint32  a = _input_num[ii];

    if (a == 0)
      addError(err, "filter '%s' invalid; there is no 0th input database.\n", _str);
    else if (a > nInputs)
      addError(err, "filter '%s' invalid: cannot occur in %u inputs; there are only %u inputs.\n", _str, a, nInputs);
    else
      _presentInNum[a] = true;
  }

  //  Build a lookup table for the databases that the kmer must be present in.
  //   - initialize everything to false.
  //   - if '@a-all' was supplied, set all those to true.
  //   - set any explicitly specified index to true.

  for (uint32 ii=0; ii<nInputs; ii++)
    _presentInIdx[ii] = false;

  for (uint32 ii=0; ii<_input_idx.size(); ii++) {
    uint32 a = _input_idx[ii];

    if (a == 0)
      addError(err, "filter '%s' invalid; there is no 0th input database.\n", _str);
    else if (a > nInputs)
      addError(err, "filter '%s' invalid: input %u does not exist; there are only %u inputs.\n", _str, a, nInputs);
    else {
      _presentInIdx[a-1] = true;
      _presentInList[_presentInLen++] = a-1;
    }
  }

  for (uint32 ii=_input_num_at_least; ii<=nInputs; ii++) {
    _presentInIdx[ii-1] = true;
    _presentInList[_presentInLen++] = ii-1;
  }

  std::sort(_presentInList, _presentInList + _presentInLen);
}



char *
merylFilter::describe(char *str) {
  char const   *rType = nullptr;
  char          lhs[256] = {0};
  char          rhs[256] = {0};

  switch (_r) {
    case merylFilterRelation::isNOP:    rType = "unspecified";        break;
    case merylFilterRelation::isEq:     rType = "is-equal-to";        break;
    case merylFilterRelation::isNeq:    rType = "is-not-equal-to";    break;
    case merylFilterRelation::isLeq:    rType = "is-less-or-equal";   break;
    case merylFilterRelation::isGeq:    rType = "is-more-or-equal";   break;
    case merylFilterRelation::isLt:     rType = "is-less-than";       break;
    case merylFilterRelation::isGt:     rType = "is-more-than";       break;
    default:                            rType = "unspecified";        break;
  }

  //  Lots of duplication with isTrue().

  if      (_q == merylFilterQuantity::isValue) {
    if      (_vIndex1 == uint32max)            sprintf(lhs, "constant value %s", toDec(_vValue1));
    else if (_vIndex1 == 0)                    sprintf(lhs, "output kmer value");
    else                                       sprintf(lhs, "kmer value from input %u", _vIndex1);

    if      (_vIndex2 == uint32max)            sprintf(rhs, "constant value %s", toDec(_vValue2));
    else if (_vIndex2 == 0)                    sprintf(rhs, "output kmer value");
    else                                       sprintf(rhs, "kmer value from input %u", _vIndex2);

    sprintf(str, "EMIT if %s %s %s\n", lhs, rType, rhs);
  }

  else if (_q == merylFilterQuantity::isLabel) {
    if      (_vIndex1 == uint32max)            sprintf(lhs, "constant label %s", toDec(_vLabel1));
    else if (_vIndex1 == 0)                    sprintf(lhs, "output kmer label");
    else                                       sprintf(lhs, "kmer label from input %u", _vIndex1);

    if      (_vIndex2 == uint32max)            sprintf(rhs, "constant label %s", toDec(_vLabel2));
    else if (_vIndex2 == 0)                    sprintf(rhs, "output kmer label");
    else                                       sprintf(rhs, "kmer label from input %u", _vIndex2);

    sprintf(str, "EMIT if %s %s %s\n", lhs, rType, rhs);
  }

  else if (_q == merylFilterQuantity::isBases) {
    char  acgt[5] = {0};
    int32 acgtlen = 0;

    if (_countA == true)   acgt[acgtlen++] = 'A';
    if (_countC == true)   acgt[acgtlen++] = 'C';
    if (_countT == true)   acgt[acgtlen++] = 'T';
    if (_countG == true)   acgt[acgtlen++] = 'G';

    //  vIndex is 0 if nothing was supplied to the filter, which tells us
    //  to get the value from the output kmer.

    if      (_vIndex1 == 0)
      sprintf(lhs, "number of %s in the kmer", acgt);
    else
      sprintf(lhs, "%s", toDec(_vBases1));

    if      (_vIndex2 == 0)
      sprintf(rhs, "number of %s in the kmer", acgt);
    else
      sprintf(rhs, "%s", toDec(_vBases2));

    sprintf(str, "EMIT if %s %s %s\n", lhs, rType, rhs);
  }

  else if (_q == merylFilterQuantity::isIndex) {
    sprintf(str, "EMIT if <index filter not described>\n");
  }

  else {
    sprintf(str, "EMIT if <empty filter>\n");
  }


  return(str);
}
