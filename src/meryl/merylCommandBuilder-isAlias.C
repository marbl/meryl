
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

  //  Decide if this is an alias or not.  If not, just get out of here.
  //  We'll have to do the test again later to decide how to set up
  //  the operation.

  if ((strcmp(_optString, "union")         != 0) &&
      (strcmp(_optString, "union-min")     != 0) &&
      (strcmp(_optString, "union-max")     != 0) &&
      (strcmp(_optString, "union-sum")     != 0) &&
      (strcmp(_optString, "intersect")     != 0) &&
      (strcmp(_optString, "intersect-min") != 0) &&
      (strcmp(_optString, "intersect-max") != 0) &&
      (strcmp(_optString, "intersect-sum") != 0) &&
      (strcmp(_optString, "subtract")      != 0) &&
      (strcmp(_optString, "difference")    != 0))
    return(false);

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

  else {         //  We detected an alias, but didn't decode it.
    assert(0);   //  Somebody screwed up the code.
  }

  return(true);
}


