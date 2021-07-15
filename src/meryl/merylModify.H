
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

#ifndef MERYLMODIFY_H
#define MERYLMODIFY_H

#include "meryl.H"


//  Each filter action needs to know how to set the value/label given the
//  various input values and labels.  These describe how to do so.


enum class merylModifyValue {
  valueNOP,         //  Do absoluately nothing.

  valueSet,         //  Set to _valueConstant

  valueSelected,    //  Set to the value of the kmer selected by the label filter
  valueFirst,       //  Set to the value of the kmer in the first input file

  valueMin,         //  Set to the minimum of all values and _valueConstant
  valueMax,         //  Set to the maximum of all values and _valueConstant

  valueAdd,         //  Set to the selected value + all others and _valueConstant
  valueSub,         //  Set to the selected value - all others and _valueConstant
  valueMul,         //  Set to the selected value * all others and _valueConstant
  valueDiv,         //  Set to the selected value / all others and _valueConstant
  valueDivZ,        //  Set to the selected value / all others and _valueConstant, round 0 to 1
  valueMod,         //  Set to the remainder of div

  valueCount,       //  Set to the number of files with this kmer.
};


enum class merylModifyLabel {
  labelNOP,         //  Do absoluately nothing.

  labelSet,         //  Set to _labelConstant.

  labelSelected,    //  Set to the label of the selected input file.
  labelFirst,       //  Set to the label of the kmer in the first input file.

  labelMin,         //  Set to the label of the kmer with the minimum value (e.g., valueMin).
  labelMax,         //  Set to the label of the kmer with the maximum value (e.g., valueMax).

  labelAnd,         //  AND all labels together along with the labelConstant
  labelOr,          //  OR  all labels together along with the labelConstant
  labelXor,         //  XOR all labels together along with the labelConstant
  labelDifference,  //  Unset any bits set in kmers 2..N (like opDifference)

  labelLightest,    //  Set label to lowest  bit weight label in active set
  labelHeaviest,    //  Set label to highest bit weight label in active set

  labelInvert,      //  Invert the existing label bits

  labelShiftLeft,   //  Shift the label left  _labelConstant bits
  labelShiftRight,  //  Shift the label right _labelConstant bits

  labelRotateLeft,  //  Shift the label left  _labelConstant bits
  labelRotateRight, //  Shift the label right _labelConstant bits
};


#endif  //  MERYLMODIFY_H
