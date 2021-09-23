
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

#ifndef MERYLINCLUDE
#error "Do not use merylModify.H, use meryl.H instead."
#endif


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

inline
char const *
toString(merylModifyValue m) {
  switch (m) {
    case merylModifyValue::valueNOP:       return("valueNOP");       break;
    case merylModifyValue::valueSet:       return("valueSet");       break;
    case merylModifyValue::valueSelected:  return("valueSelected");  break;
    case merylModifyValue::valueFirst:     return("valueFirst");     break;
    case merylModifyValue::valueMin:       return("valueMin");       break;
    case merylModifyValue::valueMax:       return("valueMax");       break;
    case merylModifyValue::valueAdd:       return("valueAdd");       break;
    case merylModifyValue::valueSub:       return("valueSub");       break;
    case merylModifyValue::valueMul:       return("valueMul");       break;
    case merylModifyValue::valueDiv:       return("valueDiv");       break;
    case merylModifyValue::valueDivZ:      return("valueDivZ");      break;
    case merylModifyValue::valueMod:       return("valueMod");       break;
    case merylModifyValue::valueCount:     return("valueCount");     break;
    default:
      break;
  }

  assert(0);
  return(nullptr);
}


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

inline
char const *
toString(merylModifyLabel m) {
  switch (m) {
    case merylModifyLabel::labelNOP:          return("labelNOP");          break;
    case merylModifyLabel::labelSet:          return("labelSet");          break;
    case merylModifyLabel::labelSelected:     return("labelSelected");     break;
    case merylModifyLabel::labelFirst:        return("labelFirst");        break;
    case merylModifyLabel::labelMin:          return("labelMin");          break;
    case merylModifyLabel::labelMax:          return("labelMax");          break;
    case merylModifyLabel::labelAnd:          return("labelAnd");          break;
    case merylModifyLabel::labelOr:           return("labelOr");           break;
    case merylModifyLabel::labelXor:          return("labelXor");          break;
    case merylModifyLabel::labelDifference:   return("labelDifference");   break;
    case merylModifyLabel::labelLightest:     return("labelLightest");     break;
    case merylModifyLabel::labelHeaviest:     return("labelHeaviest");     break;
    case merylModifyLabel::labelInvert:       return("labelInvert");       break;
    case merylModifyLabel::labelShiftLeft:    return("labelShiftLeft");    break;
    case merylModifyLabel::labelShiftRight:   return("labelShiftRight");   break;
    case merylModifyLabel::labelRotateLeft:   return("labelRotateLeft");   break;
    case merylModifyLabel::labelRotateRight:  return("labelRotateRight");  break;
    default:
      break;
  }

  assert(0);
  return(nullptr);
}


#endif  //  MERYLMODIFY_H
