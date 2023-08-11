
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

//
//  Handle value and label assignment methods.
//
//  Unlike 'output', the _curParam MUST be valid in all cases.  The word must
//  be of form 'value=something' or 'label=something'.  It would be trivial
//  to allow 'value= something'....it just looks weird.
//


bool
merylCommandBuilder::isAssignValue(char const *str) {
  merylModifyValue     select   = merylModifyValue::valueNOP;
  kmvalu               constant = kmvalumax;

  if (_curPname != opPname::pnValue)
    return false;

  //  Check for modifiers with constants at the end.

  if      (strncmp(str, "#", 1)           == 0)   { select = merylModifyValue::valueSet;        constant = decodeInteger(_optString, 6+1, 0, constant, _errors); }

  else if (strncmp(str, "min#", 4)        == 0)   { select = merylModifyValue::valueMin;        constant = decodeInteger(_optString, 6+4, 0, constant, _errors); }
  else if (strncmp(str, "max#", 4)        == 0)   { select = merylModifyValue::valueMax;        constant = decodeInteger(_optString, 6+4, 0, constant, _errors); }

  else if (strncmp(str, "add#", 4)        == 0)   { select = merylModifyValue::valueAdd;        constant = decodeInteger(_optString, 6+4, 0, constant, _errors); }
  else if (strncmp(str, "sum#", 4)        == 0)   { select = merylModifyValue::valueAdd;        constant = decodeInteger(_optString, 6+4, 0, constant, _errors); }

  else if (strncmp(str, "sub#", 4)        == 0)   { select = merylModifyValue::valueSub;        constant = decodeInteger(_optString, 6+4, 0, constant, _errors); }
  else if (strncmp(str, "dif#", 4)        == 0)   { select = merylModifyValue::valueSub;        constant = decodeInteger(_optString, 6+3, 0, constant, _errors); }

  else if (strncmp(str, "mul#", 4)        == 0)   { select = merylModifyValue::valueMul;        constant = decodeInteger(_optString, 6+4, 0, constant, _errors); }

  else if (strncmp(str, "div#", 4)        == 0)   { select = merylModifyValue::valueDiv;        constant = decodeInteger(_optString, 6+4, 0, constant, _errors); }
  else if (strncmp(str, "divzero#", 8)    == 0)   { select = merylModifyValue::valueDivZ;       constant = decodeInteger(_optString, 6+8, 0, constant, _errors); }

  else if (strncmp(str, "mod#", 4)        == 0)   { select = merylModifyValue::valueMod;        constant = decodeInteger(_optString, 6+4, 0, constant, _errors); }
  else if (strncmp(str, "rem#", 4)        == 0)   { select = merylModifyValue::valueMod;        constant = decodeInteger(_optString, 6+4, 0, constant, _errors); }

  //  Check for modifiers without constants.  Set the constant to whatever
  //  the identity is for the given modifier.

  else if (strncmp(str, "first",  6)      == 0)   { select = merylModifyValue::valueSelected;   constant = 0;         }
  else if (strncmp(str, "selected",  9)   == 0)   { select = merylModifyValue::valueSelected;   constant = kmvalumax; }

  else if (strncmp(str, "min",  4)        == 0)   { select = merylModifyValue::valueMin;        constant = kmvalumax; }
  else if (strncmp(str, "max",  4)        == 0)   { select = merylModifyValue::valueMax;        constant = 0;         }

  else if (strncmp(str, "add",  4)        == 0)   { select = merylModifyValue::valueAdd;        constant = 0;         }
  else if (strncmp(str, "sum",  4)        == 0)   { select = merylModifyValue::valueAdd;        constant = 0;         }

  else if (strncmp(str, "sub",  4)        == 0)   { select = merylModifyValue::valueSub;        constant = 0;         }
  else if (strncmp(str, "dif",  4)        == 0)   { select = merylModifyValue::valueSub;        constant = 0;         }

  else if (strncmp(str, "mul",  4)        == 0)   { select = merylModifyValue::valueMul;        constant = 1;         }

  else if (strncmp(str, "div",  4)        == 0)   { select = merylModifyValue::valueDiv;        constant = 1;         }
  else if (strncmp(str, "divzero",  8)    == 0)   { select = merylModifyValue::valueDivZ;       constant = 1;         }

  else if (strncmp(str, "mod",  4)        == 0)   { select = merylModifyValue::valueMod;        constant = 0;         }
  else if (strncmp(str, "rem",  4)        == 0)   { select = merylModifyValue::valueMod;        constant = 0;         }

  else if (strncmp(str, "count",  6)      == 0)   { select = merylModifyValue::valueCount;      constant = 0;         }

  //  Nope, don't know what this is.

  else {
    sprintf(_errors, "Unknown value selection modifier '%s'.", _optString);
  }

  //  Copy the decoded modifier to the select.

  getCurrent()->_valueSelect   = select;
  getCurrent()->_valueConstant = constant;

  return(true);
}


bool
merylCommandBuilder::isAssignLabel(char const *str) {
  merylModifyLabel     select   = merylModifyLabel::labelNOP;
  kmlabl               constant = kmlablmax;

  if (_curPname != opPname::pnLabel)
    return false;

  //  Check for modifiers with constants at the end.

  if      (strncmp(str, "#", 1)              == 0)   { select = merylModifyLabel::labelSet;          constant = decodeInteger(_optString, 6+1,  0, constant, _errors); }

  else if (strncmp(str, "min#", 4)           == 0)   { select = merylModifyLabel::labelMin;          constant = decodeInteger(_optString, 6+4,  0, constant, _errors); }
  else if (strncmp(str, "max#", 4)           == 0)   { select = merylModifyLabel::labelMax;          constant = decodeInteger(_optString, 6+3,  0, constant, _errors); }

  else if (strncmp(str, "and#", 4)           == 0)   { select = merylModifyLabel::labelAnd;          constant = decodeInteger(_optString, 6+4,  0, constant, _errors); }
  else if (strncmp(str, "or#", 3)            == 0)   { select = merylModifyLabel::labelOr;           constant = decodeInteger(_optString, 6+3,  0, constant, _errors); }
  else if (strncmp(str, "xor#", 4)           == 0)   { select = merylModifyLabel::labelXor;          constant = decodeInteger(_optString, 6+4,  0, constant, _errors); }

  else if (strncmp(str, "difference#", 11)   == 0)   { select = merylModifyLabel::labelDifference;   constant = decodeInteger(_optString, 6+11, 0, constant, _errors); }

  else if (strncmp(str, "lightest#", 9)      == 0)   { select = merylModifyLabel::labelLightest;     constant = decodeInteger(_optString, 6+9,  0, constant, _errors); }
  else if (strncmp(str, "heaviest#", 9)      == 0)   { select = merylModifyLabel::labelHeaviest;     constant = decodeInteger(_optString, 6+9,  0, constant, _errors); }

  else if (strncmp(str, "invert#", 7)        == 0)   { select = merylModifyLabel::labelHeaviest;     constant = decodeInteger(_optString, 6+7,  0, constant, _errors); }

  else if (strncmp(str, "shift-left#", 11)   == 0)   { select = merylModifyLabel::labelShiftLeft;    constant = decodeInteger(_optString, 6+11, 0, constant, _errors); }
  else if (strncmp(str, "shift-right#", 12)  == 0)   { select = merylModifyLabel::labelShiftRight;   constant = decodeInteger(_optString, 6+12, 0, constant, _errors); }

  else if (strncmp(str, "rotate-left#", 12)  == 0)   { select = merylModifyLabel::labelRotateLeft;   constant = decodeInteger(_optString, 6+12, 0, constant, _errors); }
  else if (strncmp(str, "rotate-right#", 13) == 0)   { select = merylModifyLabel::labelRotateRight;  constant = decodeInteger(_optString, 6+13, 0, constant, _errors); }

  //  Check for modifiers without constants.  Set the constant to whatever
  //  the identity is for the given modifier.

  else if (strncmp(str, "first", 6)          == 0)   { select = merylModifyLabel::labelFirst;        constant = 0;          }
  else if (strncmp(str, "selected", 9)       == 0)   { select = merylModifyLabel::labelFirst;        constant = 0;          }

  else if (strncmp(str, "min", 4)            == 0)   { select = merylModifyLabel::labelMin;          constant = 0;          }
  else if (strncmp(str, "max", 4)            == 0)   { select = merylModifyLabel::labelMax;          constant = 0;          }

  else if (strncmp(str, "and", 4)            == 0)   { select = merylModifyLabel::labelAnd;          constant = kmlablmax;  }
  else if (strncmp(str, "or", 3)             == 0)   { select = merylModifyLabel::labelOr;           constant = 0;          }
  else if (strncmp(str, "xor", 4)            == 0)   { select = merylModifyLabel::labelXor;          constant = kmlablmax;  }

  else if (strncmp(str, "difference", 11)    == 0)   { select = merylModifyLabel::labelDifference;   constant = 0;          }

  else if (strncmp(str, "lightest", 9)       == 0)   { select = merylModifyLabel::labelLightest;     constant = kmlablmax;  }
  else if (strncmp(str, "heaviest", 9)       == 0)   { select = merylModifyLabel::labelHeaviest;     constant = 0;          }

  else if (strncmp(str, "invert", 7)         == 0)   { select = merylModifyLabel::labelHeaviest;     constant = kmlablmax;  }

  else if (strncmp(str, "shift-left", 11)    == 0)   { select = merylModifyLabel::labelShiftLeft;    constant = 1;          }
  else if (strncmp(str, "shift-right", 12)   == 0)   { select = merylModifyLabel::labelShiftRight;   constant = 1;          }

  else if (strncmp(str, "rotate-left", 12)   == 0)   { select = merylModifyLabel::labelRotateLeft;   constant = 1;          }
  else if (strncmp(str, "rotate-right", 13)  == 0)   { select = merylModifyLabel::labelRotateRight;  constant = 1;          }

  //  Nope, don't know what this is.

  else {
    sprintf(_errors, "Unknown label selection modifier '%s'.", _optString);
  }

  //  Copy the decoded modifier to the select.

  getCurrent()->_labelSelect   = select;
  getCurrent()->_labelConstant = constant;

  return(true);
}



bool
merylCommandBuilder::isAssign(void) {
  char const  *str = _optString + 6;

  if (_curClass != opClass::clAssign)      //  Silently ignore non-assign classes.
    return false;

  if ((_curPname != opPname::pnValue) &&   //  Noisly complain about unknown Pnames.
      (_curPname != opPname::pnLabel)) {
#warning error
    fprintf(stderr, "expecting value or label parameter name\n");
    return false;
  }

  if ((_curParam == nullptr) ||            //  Noisly complain about missing parameters.
      (_curParam[0] == 0)) {
#warning error
    fprintf(stderr, "expecting parameter\n");
    return false;
  }

  return ((isAssignValue(str) == true) ||     //  Process the parameter.
          (isAssignLabel(str) == true));
}
