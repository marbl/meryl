
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

#ifndef MERYLASSIGN_H
#define MERYLASSIGN_H

#ifndef MERYLINCLUDE
#error "Do not use merylAssign.H, use meryl.H instead."
#endif


//  Each action needs to know how to set the value/label given the
//  various input values and labels.  These describe how to do so.


enum class merylAssignValue {
  valueNOP,         //  Do absoluately nothing.

  valueSet,         //  Set to _valueConstant

  valueSelected,    //  Set to the value of the kmer selected by the label selector
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
toString(merylAssignValue m) {
  switch (m) {
    case merylAssignValue::valueNOP:       return("valueNOP");       break;
    case merylAssignValue::valueSet:       return("valueSet");       break;
    case merylAssignValue::valueSelected:  return("valueSelected");  break;
    case merylAssignValue::valueFirst:     return("valueFirst");     break;
    case merylAssignValue::valueMin:       return("valueMin");       break;
    case merylAssignValue::valueMax:       return("valueMax");       break;
    case merylAssignValue::valueAdd:       return("valueAdd");       break;
    case merylAssignValue::valueSub:       return("valueSub");       break;
    case merylAssignValue::valueMul:       return("valueMul");       break;
    case merylAssignValue::valueDiv:       return("valueDiv");       break;
    case merylAssignValue::valueDivZ:      return("valueDivZ");      break;
    case merylAssignValue::valueMod:       return("valueMod");       break;
    case merylAssignValue::valueCount:     return("valueCount");     break;
    default:
      break;
  }

  assert(0);
  return(nullptr);
}


enum class merylAssignLabel {
  labelNOP,         //  Do absoluately nothing.

  labelSet,         //  Set to _labelConstant.

  labelSelected,    //  Set to the label of the kmer selected by the value selector.
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
toString(merylAssignLabel m) {
  switch (m) {
    case merylAssignLabel::labelNOP:          return("labelNOP");          break;
    case merylAssignLabel::labelSet:          return("labelSet");          break;
    case merylAssignLabel::labelSelected:     return("labelSelected");     break;
    case merylAssignLabel::labelFirst:        return("labelFirst");        break;
    case merylAssignLabel::labelMin:          return("labelMin");          break;
    case merylAssignLabel::labelMax:          return("labelMax");          break;
    case merylAssignLabel::labelAnd:          return("labelAnd");          break;
    case merylAssignLabel::labelOr:           return("labelOr");           break;
    case merylAssignLabel::labelXor:          return("labelXor");          break;
    case merylAssignLabel::labelDifference:   return("labelDifference");   break;
    case merylAssignLabel::labelLightest:     return("labelLightest");     break;
    case merylAssignLabel::labelHeaviest:     return("labelHeaviest");     break;
    case merylAssignLabel::labelInvert:       return("labelInvert");       break;
    case merylAssignLabel::labelShiftLeft:    return("labelShiftLeft");    break;
    case merylAssignLabel::labelShiftRight:   return("labelShiftRight");   break;
    case merylAssignLabel::labelRotateLeft:   return("labelRotateLeft");   break;
    case merylAssignLabel::labelRotateRight:  return("labelRotateRight");  break;
    default:
      break;
  }

  assert(0);
  return(nullptr);
}


#endif  //  MERYLASSIGN_H
