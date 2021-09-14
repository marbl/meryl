
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


//  Should we fail if an alias is detected for an operation being set up?


bool
merylCommandBuilder::isAlias(void) {
  merylOpTemplate  *op = getCurrent();

  //  Check for _needsValue and _needsConstant.  If set, assume this word is
  //  an integer constant, and copy the constant to the filter (_needsValue)
  //  or operation (_needsConstant).  For both, the word is always consumed,
  //  but we emit an error if it's not understood.
  //
  //  _needsValue is used by filter aliases:
  //    less-than,
  //    greater-than
  //    at-least
  //    at-most
  //    equal-to
  //    not-equal-to
  //
  //  _needsConstant is used by selection aliases:
  //    increase
  //    decrease
  //    multiply
  //    divide
  //    divide-round
  //    modulo

  if ((_needsValue    == true) ||
      (_needsConstant == true)) {
    char   *optstr   = _optString;
    uint32  constant = 0;

    if      ((_needsValue == true) && (strncmp(_optString, "distinct=", 9) == 0))
      op->getLastFilter()._vValue2Distinct = strtodouble(optstr + 9);

    else if ((_needsValue == true) && (strncmp(_optString, "word-freq=", 10) == 0))
      op->getLastFilter()._vValue2WordFreq = strtodouble(optstr + 10);

    else if ((_needsValue == true) && (strncmp(_optString, "word-frequency=", 15) == 0))
      op->getLastFilter()._vValue2WordFreq = strtodouble(optstr + 15);

    else if ((_needsValue == true) && (strncmp(_optString, "threshold=", 10) == 0))
      op->getLastFilter()._vValue2 = decodeInteger(optstr, 10, 0, constant, _errors);

    else if (_needsValue == true)
      op->getLastFilter()._vValue2 = decodeInteger(optstr, 0, 0, constant, _errors);

    else if (_needsConstant == true)
      op->_valueConstant = decodeInteger(optstr, 0, 0, constant, _errors);

    else
      assert(0);

    _needsValue    = false;   //  We've processed the requested argument
    _needsConstant = false;   //  for either of these just by getting here.

    return(true);
  }

  //  Decode the alias and set up the operation.

  if      (strcmp(_optString, "union") == 0) {
    op->_valueSelect    = merylModifyValue::valueCount;
    op->_labelSelect    = merylModifyLabel::labelOr;

    merylFilter  f(merylFilterQuantity::isIndex,
                   merylFilterRelation::isNOP, false, "union");

    f._input_num_any    = true;

    op->addFilterToProduct(f);
  }

  else if (strcmp(_optString, "union-min") == 0) {
    op->_valueSelect    = merylModifyValue::valueMin;
    op->_valueConstant  = kmvalumax;
    op->_labelSelect    = merylModifyLabel::labelSelected;

    merylFilter  f(merylFilterQuantity::isIndex,
                   merylFilterRelation::isNOP, false, "union");

    f._input_num_any    = true;

    op->addFilterToProduct(f);
  }

  else if (strcmp(_optString, "union-max") == 0) {
    op->_valueSelect    = merylModifyValue::valueMax;
    op->_valueConstant  = kmvalumin;
    op->_labelSelect    = merylModifyLabel::labelSelected;

    merylFilter  f(merylFilterQuantity::isIndex,
                   merylFilterRelation::isNOP, false, "union");

    f._input_num_any    = true;

    op->addFilterToProduct(f);
  }

  else if (strcmp(_optString, "union-sum") == 0) {
    op->_valueSelect    = merylModifyValue::valueAdd;
    op->_labelSelect    = merylModifyLabel::labelOr;

    merylFilter  f(merylFilterQuantity::isIndex,
                   merylFilterRelation::isNOP, false, "union");

    f._input_num_any    = true;

    op->addFilterToProduct(f);
  }


  else if (strcmp(_optString, "intersect") == 0) {
    op->_valueSelect    = merylModifyValue::valueFirst;
    op->_labelSelect    = merylModifyLabel::labelAnd;

    merylFilter  f(merylFilterQuantity::isIndex,
                   merylFilterRelation::isNOP, false, "intersect");

    f._input_num_all    = true;

    op->addFilterToProduct(f);
  }

  else if (strcmp(_optString, "intersect-min") == 0) {
    op->_valueSelect    = merylModifyValue::valueMin;
    op->_valueConstant  = kmvalumax;
    op->_labelSelect    = merylModifyLabel::labelSelected;

    merylFilter  f(merylFilterQuantity::isIndex,
                   merylFilterRelation::isNOP, false, "intersect");

    f._input_num_all    = true;

    op->addFilterToProduct(f);
  }

  else if (strcmp(_optString, "intersect-max") == 0) {
    op->_valueSelect    = merylModifyValue::valueMax;
    op->_valueConstant  = kmvalumin;
    op->_labelSelect    = merylModifyLabel::labelSelected;

    merylFilter  f(merylFilterQuantity::isIndex,
                   merylFilterRelation::isNOP, false, "intersect");

    f._input_num_all    = true;

    op->addFilterToProduct(f);
  }

  else if (strcmp(_optString, "intersect-sum") == 0) {
    op->_valueSelect    = merylModifyValue::valueAdd;
    op->_labelSelect    = merylModifyLabel::labelAnd;

    merylFilter  f(merylFilterQuantity::isIndex,
                   merylFilterRelation::isNOP, false, "intersect");

    f._input_num_all    = true;

    op->addFilterToProduct(f);
  }


  //  kmer occurs in the first input, subtract value of all other inputs
  else if (strcmp(_optString, "subtract") == 0) {
    op->_valueSelect    = merylModifyValue::valueSub;
    op->_labelSelect    = merylModifyLabel::labelDifference;

    merylFilter  f(merylFilterQuantity::isIndex,
                   merylFilterRelation::isNOP, false, "subtract");

    f._input_num_any    = true;
    f._input_idx.push_back(1);

    op->addFilterToProduct(f);
  }

  //  kmer must occur only in the first input
  else if (strcmp(_optString, "difference") == 0) {
    op->_valueSelect    = merylModifyValue::valueFirst;
    op->_labelSelect    = merylModifyLabel::labelFirst;

    merylFilter  f(merylFilterQuantity::isIndex,
                   merylFilterRelation::isNOP, false, "difference");

    f._input_num.push_back(1);
    f._input_idx.push_back(1);

    op->addFilterToProduct(f);
  }

  else if (strcmp(_optString, "distinct") == 0) {
    assert(0);
  }

  else if (strcmp(_optString, "word-frequency") == 0) {
    assert(0);
  }

  else if (strcmp(_optString, "threshold") == 0) {
    assert(0);
  }

  else if (strcmp(_optString, "count-suffix") == 0) {
  }


  else if ((strcmp(_optString, "less-than")    == 0) ||
           (strcmp(_optString, "greater-than") == 0) ||
           (strcmp(_optString, "at-least")     == 0) ||
           (strcmp(_optString, "at-most")      == 0) ||
           (strcmp(_optString, "equal-to")     == 0) ||
           (strcmp(_optString, "not-equal-to") == 0)) {
    op->_inputsMin      = 1;
    op->_inputsMax      = 1;

    op->_valueSelect    = merylModifyValue::valueFirst;
    op->_labelSelect    = merylModifyLabel::labelFirst;

    merylFilterRelation   rel = merylFilterRelation::isNOP;

    if (strcmp(_optString, "less-than")    == 0)   rel = merylFilterRelation::isLt;
    if (strcmp(_optString, "greater-than") == 0)   rel = merylFilterRelation::isGt;
    if (strcmp(_optString, "at-least")     == 0)   rel = merylFilterRelation::isGeq;
    if (strcmp(_optString, "at-most")      == 0)   rel = merylFilterRelation::isLeq;
    if (strcmp(_optString, "equal-to")     == 0)   rel = merylFilterRelation::isEq;
    if (strcmp(_optString, "not-equal-to") == 0)   rel = merylFilterRelation::isNeq;

    merylFilter  f(merylFilterQuantity::isValue, rel, false, _optString);

    f._vIndex1 = 0;
    f._vValue2 = 0;

    op->addFilterToProduct(f);

    _needsValue = true;
  }


  else if ((strcmp(_optString, "increase")     == 0) ||
           (strcmp(_optString, "decrease")     == 0) ||
           (strcmp(_optString, "multiply")     == 0) ||
           (strcmp(_optString, "divide")       == 0) ||
           (strcmp(_optString, "divide-round") == 0) ||
           (strcmp(_optString, "modulo")       == 0)) {
    op->_inputsMin      = 1;
    op->_inputsMax      = 1;

    op->_valueSelect    = merylModifyValue::valueNOP;
    op->_valueConstant  = 0;
    op->_labelSelect    = merylModifyLabel::labelFirst;

    if (strcmp(_optString, "increase")     == 0)   op->_valueSelect = merylModifyValue::valueAdd;
    if (strcmp(_optString, "decrease")     == 0)   op->_valueSelect = merylModifyValue::valueSub;
    if (strcmp(_optString, "multiply")     == 0)   op->_valueSelect = merylModifyValue::valueMul;
    if (strcmp(_optString, "divide")       == 0)   op->_valueSelect = merylModifyValue::valueDiv;
    if (strcmp(_optString, "divide-round") == 0)   op->_valueSelect = merylModifyValue::valueDivZ;
    if (strcmp(_optString, "modulo")       == 0)   op->_valueSelect = merylModifyValue::valueMod;

    _needsConstant = true;
  }

  else {             //  Not an alias word.
    return(false);   //
  }

  return(true);
}


