
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





//  If the previous word was an alias and it needs a value,
//  Check for _needsValue and _needsConstant.  If set, assume this word is
//  an integer constant, and copy the constant to the selector (_needsValue)
//  or operation (_needsConstant).  For both, the word is always consumed,
//  but we emit an error if it's not understood.
//
//  _needsValue is used by selector aliases:
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
bool
merylCommandBuilder::isAliasConstant(void) {
  merylOpTemplate  *op = getCurrent();
  uint32            cn = 0;

  if      (_needsConstant) {
    op->_valueConstant = decodeInteger(_optString, 0, 0, cn, _errors);
  }

  else if (_needsValue) {
    if     (strncmp(_optString, "distinct=", 9) == 0)
      op->getLastSelector()._vValue2Distinct = strtodouble(_optString + 9);

    else if (strncmp(_optString, "word-freq=", 10) == 0)
      op->getLastSelector()._vValue2WordFreq = strtodouble(_optString + 10);

    else if (strncmp(_optString, "word-frequency=", 15) == 0)
      op->getLastSelector()._vValue2WordFreq = strtodouble(_optString + 15);

    else if (strncmp(_optString, "threshold=", 10) == 0)
      op->getLastSelector()._vValue2 = decodeInteger(_optString, 10, 0, cn, _errors);

    else
      op->getLastSelector()._vValue2 = decodeInteger(_optString, 0, 0, cn, _errors);
  }

  else {            //  Doesn't need a constant
    return false;   //  or a value.
  }

  _needsValue    = false;   //  We've processed the requested argument (or failed)
  _needsConstant = false;   //  either way, the constant has been handled.

  return true;
}


bool
merylCommandBuilder::isAliasWord(void) {
  merylOpTemplate  *op = getCurrent();


  //  Decode the alias and set up the operation.

  if      (strcmp(_optString, "union") == 0) {
    op->_valueAssign    = merylAssignValue::valueCount;
    op->_labelAssign    = merylAssignLabel::labelOr;

    merylSelector  f(merylSelectorQuantity::isIndex,
                     merylSelectorRelation::isNOP, false, "union");

    f._input_num_any    = true;

    op->addSelectorToProduct(f);
  }

  else if (strcmp(_optString, "union-min") == 0) {
    op->_valueAssign    = merylAssignValue::valueMin;
    op->_valueConstant  = kmvalumax;
    op->_labelAssign    = merylAssignLabel::labelSelected;

    merylSelector  f(merylSelectorQuantity::isIndex,
                     merylSelectorRelation::isNOP, false, "union");

    f._input_num_any    = true;

    op->addSelectorToProduct(f);
  }

  else if (strcmp(_optString, "union-max") == 0) {
    op->_valueAssign    = merylAssignValue::valueMax;
    op->_valueConstant  = kmvalumin;
    op->_labelAssign    = merylAssignLabel::labelSelected;

    merylSelector  f(merylSelectorQuantity::isIndex,
                     merylSelectorRelation::isNOP, false, "union");

    f._input_num_any    = true;

    op->addSelectorToProduct(f);
  }

  else if (strcmp(_optString, "union-sum") == 0) {
    op->_valueAssign    = merylAssignValue::valueAdd;
    op->_labelAssign    = merylAssignLabel::labelOr;

    merylSelector  f(merylSelectorQuantity::isIndex,
                     merylSelectorRelation::isNOP, false, "union");

    f._input_num_any    = true;

    op->addSelectorToProduct(f);
  }


  else if (strcmp(_optString, "intersect") == 0) {
    op->_valueAssign    = merylAssignValue::valueFirst;
    op->_labelAssign    = merylAssignLabel::labelAnd;

    merylSelector  f(merylSelectorQuantity::isIndex,
                     merylSelectorRelation::isNOP, false, "intersect");

    f._input_num_all    = true;

    op->addSelectorToProduct(f);
  }

  else if (strcmp(_optString, "intersect-min") == 0) {
    op->_valueAssign    = merylAssignValue::valueMin;
    op->_valueConstant  = kmvalumax;
    op->_labelAssign    = merylAssignLabel::labelSelected;

    merylSelector  f(merylSelectorQuantity::isIndex,
                     merylSelectorRelation::isNOP, false, "intersect");

    f._input_num_all    = true;

    op->addSelectorToProduct(f);
  }

  else if (strcmp(_optString, "intersect-max") == 0) {
    op->_valueAssign    = merylAssignValue::valueMax;
    op->_valueConstant  = kmvalumin;
    op->_labelAssign    = merylAssignLabel::labelSelected;

    merylSelector  f(merylSelectorQuantity::isIndex,
                     merylSelectorRelation::isNOP, false, "intersect");

    f._input_num_all    = true;

    op->addSelectorToProduct(f);
  }

  else if (strcmp(_optString, "intersect-sum") == 0) {
    op->_valueAssign    = merylAssignValue::valueAdd;
    op->_labelAssign    = merylAssignLabel::labelAnd;

    merylSelector  f(merylSelectorQuantity::isIndex,
                     merylSelectorRelation::isNOP, false, "intersect");

    f._input_num_all    = true;

    op->addSelectorToProduct(f);
  }


  //  kmer occurs in the first input, subtract value of all other inputs
  else if (strcmp(_optString, "subtract") == 0) {
    op->_valueAssign    = merylAssignValue::valueSub;
    op->_labelAssign    = merylAssignLabel::labelDifference;

    merylSelector  f(merylSelectorQuantity::isIndex,
                     merylSelectorRelation::isNOP, false, "subtract");

    f._input_num_any    = true;
    f._input_idx.push_back(1);

    op->addSelectorToProduct(f);
  }

  //  kmer must occur only in the first input
  else if (strcmp(_optString, "difference") == 0) {
    op->_valueAssign    = merylAssignValue::valueFirst;
    op->_labelAssign    = merylAssignLabel::labelFirst;

    merylSelector  f(merylSelectorQuantity::isIndex,
                     merylSelectorRelation::isNOP, false, "difference");

    f._input_num.push_back(1);
    f._input_idx.push_back(1);

    op->addSelectorToProduct(f);
  }

  else if ((strcmp(_optString, "less-than")    == 0) ||
           (strcmp(_optString, "greater-than") == 0) ||
           (strcmp(_optString, "at-least")     == 0) ||
           (strcmp(_optString, "at-most")      == 0) ||
           (strcmp(_optString, "equal-to")     == 0) ||
           (strcmp(_optString, "not-equal-to") == 0)) {
    op->_inputsMin      = 1;
    op->_inputsMax      = 1;

    op->_valueAssign    = merylAssignValue::valueFirst;
    op->_labelAssign    = merylAssignLabel::labelFirst;

    merylSelectorRelation   rel = merylSelectorRelation::isNOP;

    if (strcmp(_optString, "less-than")    == 0)   rel = merylSelectorRelation::isLt;
    if (strcmp(_optString, "greater-than") == 0)   rel = merylSelectorRelation::isGt;
    if (strcmp(_optString, "at-least")     == 0)   rel = merylSelectorRelation::isGeq;
    if (strcmp(_optString, "at-most")      == 0)   rel = merylSelectorRelation::isLeq;
    if (strcmp(_optString, "equal-to")     == 0)   rel = merylSelectorRelation::isEq;
    if (strcmp(_optString, "not-equal-to") == 0)   rel = merylSelectorRelation::isNeq;

    merylSelector  f(merylSelectorQuantity::isValue, rel, false, _optString);

    f._vIndex1 = 0;
    f._vValue2 = 0;

    op->addSelectorToProduct(f);

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

    op->_valueAssign    = merylAssignValue::valueNOP;
    op->_valueConstant  = 0;
    op->_labelAssign    = merylAssignLabel::labelFirst;

    if (strcmp(_optString, "increase")     == 0)   op->_valueAssign = merylAssignValue::valueAdd;
    if (strcmp(_optString, "decrease")     == 0)   op->_valueAssign = merylAssignValue::valueSub;
    if (strcmp(_optString, "multiply")     == 0)   op->_valueAssign = merylAssignValue::valueMul;
    if (strcmp(_optString, "divide")       == 0)   op->_valueAssign = merylAssignValue::valueDiv;
    if (strcmp(_optString, "divide-round") == 0)   op->_valueAssign = merylAssignValue::valueDivZ;
    if (strcmp(_optString, "modulo")       == 0)   op->_valueAssign = merylAssignValue::valueMod;

    _needsConstant = true;
  }

  else {             //  Not an alias word.
    return(false);   //
  }

  return(true);
}


