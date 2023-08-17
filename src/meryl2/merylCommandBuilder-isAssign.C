
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
merylCommandBuilder::isAssignValue(void) {
  merylAssignValue     assign   = merylAssignValue::valueNOP;
  kmvalu               constant = kmvalumax;

  if (_curPname != opPname::pnValue)
    return false;

  fprintf(stderr, "isAssignValue()-- '%s' with param '%s'\n", _optString, _curParam);

  //  Check for modifiers with constants at the end.

  if      (strncmp(_curParam, "#", 1)           == 0)   { assign = merylAssignValue::valueSet;        constant = decodeInteger(_curParam, 1, 0, constant, _errors); }

  else if (strncmp(_curParam, "min#", 4)        == 0)   { assign = merylAssignValue::valueMin;        constant = decodeInteger(_curParam, 4, 0, constant, _errors); }
  else if (strncmp(_curParam, "max#", 4)        == 0)   { assign = merylAssignValue::valueMax;        constant = decodeInteger(_curParam, 4, 0, constant, _errors); }

  else if (strncmp(_curParam, "add#", 4)        == 0)   { assign = merylAssignValue::valueAdd;        constant = decodeInteger(_curParam, 4, 0, constant, _errors); }
  else if (strncmp(_curParam, "sum#", 4)        == 0)   { assign = merylAssignValue::valueAdd;        constant = decodeInteger(_curParam, 4, 0, constant, _errors); }

  else if (strncmp(_curParam, "sub#", 4)        == 0)   { assign = merylAssignValue::valueSub;        constant = decodeInteger(_curParam, 4, 0, constant, _errors); }
  else if (strncmp(_curParam, "dif#", 4)        == 0)   { assign = merylAssignValue::valueSub;        constant = decodeInteger(_curParam, 3, 0, constant, _errors); }

  else if (strncmp(_curParam, "mul#", 4)        == 0)   { assign = merylAssignValue::valueMul;        constant = decodeInteger(_curParam, 4, 0, constant, _errors); }

  else if (strncmp(_curParam, "div#", 4)        == 0)   { assign = merylAssignValue::valueDiv;        constant = decodeInteger(_curParam, 4, 0, constant, _errors); }
  else if (strncmp(_curParam, "divzero#", 8)    == 0)   { assign = merylAssignValue::valueDivZ;       constant = decodeInteger(_curParam, 8, 0, constant, _errors); }

  else if (strncmp(_curParam, "mod#", 4)        == 0)   { assign = merylAssignValue::valueMod;        constant = decodeInteger(_curParam, 4, 0, constant, _errors); }
  else if (strncmp(_curParam, "rem#", 4)        == 0)   { assign = merylAssignValue::valueMod;        constant = decodeInteger(_curParam, 4, 0, constant, _errors); }

  //  Check for modifiers without constants.  Set the constant to whatever
  //  the identity is for the given modifier.

  else if (strncmp(_curParam, "first",  6)      == 0)   { assign = merylAssignValue::valueFirst;      constant = 0;         }
  else if (strncmp(_curParam, "selected",  9)   == 0)   { assign = merylAssignValue::valueSelected;   constant = kmvalumax; }

  else if (strncmp(_curParam, "min",  4)        == 0)   { assign = merylAssignValue::valueMin;        constant = kmvalumax; }
  else if (strncmp(_curParam, "max",  4)        == 0)   { assign = merylAssignValue::valueMax;        constant = 0;         }

  else if (strncmp(_curParam, "add",  4)        == 0)   { assign = merylAssignValue::valueAdd;        constant = 0;         }
  else if (strncmp(_curParam, "sum",  4)        == 0)   { assign = merylAssignValue::valueAdd;        constant = 0;         }

  else if (strncmp(_curParam, "sub",  4)        == 0)   { assign = merylAssignValue::valueSub;        constant = 0;         }
  else if (strncmp(_curParam, "dif",  4)        == 0)   { assign = merylAssignValue::valueSub;        constant = 0;         }

  else if (strncmp(_curParam, "mul",  4)        == 0)   { assign = merylAssignValue::valueMul;        constant = 1;         }

  else if (strncmp(_curParam, "div",  4)        == 0)   { assign = merylAssignValue::valueDiv;        constant = 1;         }
  else if (strncmp(_curParam, "divzero",  8)    == 0)   { assign = merylAssignValue::valueDivZ;       constant = 1;         }

  else if (strncmp(_curParam, "mod",  4)        == 0)   { assign = merylAssignValue::valueMod;        constant = 0;         }
  else if (strncmp(_curParam, "rem",  4)        == 0)   { assign = merylAssignValue::valueMod;        constant = 0;         }

  else if (strncmp(_curParam, "count",  6)      == 0)   { assign = merylAssignValue::valueCount;      constant = 0;         }

  //  Nope, don't know what this is.

  else {
    sprintf(_errors, "Unknown value assign:value=<parameter> in '%s'.", _optString);
    return false;
  }

  getCurrent()->_valueAssign   = assign;
  getCurrent()->_valueConstant = constant;

  resetClass();

  return true;
}


bool
merylCommandBuilder::isAssignLabel(void) {
  merylAssignLabel     assign   = merylAssignLabel::labelNOP;
  kmlabl               constant = kmlablmax;

  if (_curPname != opPname::pnLabel)
    return false;

  fprintf(stderr, "isAssignLabel()-- '%s' with param '%s'\n", _optString, _curParam);

  //  Check for modifiers with constants at the end.

  if      (strncmp(_curParam, "#", 1)              == 0)   { assign = merylAssignLabel::labelSet;          constant = decodeInteger(_curParam, 1,  0, constant, _errors); }

  else if (strncmp(_curParam, "min#", 4)           == 0)   { assign = merylAssignLabel::labelMin;          constant = decodeInteger(_curParam, 4,  0, constant, _errors); }
  else if (strncmp(_curParam, "max#", 4)           == 0)   { assign = merylAssignLabel::labelMax;          constant = decodeInteger(_curParam, 3,  0, constant, _errors); }

  else if (strncmp(_curParam, "and#", 4)           == 0)   { assign = merylAssignLabel::labelAnd;          constant = decodeInteger(_curParam, 4,  0, constant, _errors); }
  else if (strncmp(_curParam, "or#", 3)            == 0)   { assign = merylAssignLabel::labelOr;           constant = decodeInteger(_curParam, 3,  0, constant, _errors); }
  else if (strncmp(_curParam, "xor#", 4)           == 0)   { assign = merylAssignLabel::labelXor;          constant = decodeInteger(_curParam, 4,  0, constant, _errors); }

  else if (strncmp(_curParam, "difference#", 11)   == 0)   { assign = merylAssignLabel::labelDifference;   constant = decodeInteger(_curParam, 11, 0, constant, _errors); }

  else if (strncmp(_curParam, "lightest#", 9)      == 0)   { assign = merylAssignLabel::labelLightest;     constant = decodeInteger(_curParam, 9,  0, constant, _errors); }
  else if (strncmp(_curParam, "heaviest#", 9)      == 0)   { assign = merylAssignLabel::labelHeaviest;     constant = decodeInteger(_curParam, 9,  0, constant, _errors); }

  else if (strncmp(_curParam, "invert#", 7)        == 0)   { assign = merylAssignLabel::labelHeaviest;     constant = decodeInteger(_curParam, 7,  0, constant, _errors); }

  else if (strncmp(_curParam, "shift-left#", 11)   == 0)   { assign = merylAssignLabel::labelShiftLeft;    constant = decodeInteger(_curParam, 11, 0, constant, _errors); }
  else if (strncmp(_curParam, "shift-right#", 12)  == 0)   { assign = merylAssignLabel::labelShiftRight;   constant = decodeInteger(_curParam, 12, 0, constant, _errors); }

  else if (strncmp(_curParam, "rotate-left#", 12)  == 0)   { assign = merylAssignLabel::labelRotateLeft;   constant = decodeInteger(_curParam, 12, 0, constant, _errors); }
  else if (strncmp(_curParam, "rotate-right#", 13) == 0)   { assign = merylAssignLabel::labelRotateRight;  constant = decodeInteger(_curParam, 13, 0, constant, _errors); }

  //  Check for modifiers without constants.  Set the constant to whatever
  //  the identity is for the given modifier.

  else if (strncmp(_curParam, "first", 6)          == 0)   { assign = merylAssignLabel::labelFirst;        constant = 0;          }
  else if (strncmp(_curParam, "selected", 9)       == 0)   { assign = merylAssignLabel::labelSelected;     constant = 0;          }

  else if (strncmp(_curParam, "min", 4)            == 0)   { assign = merylAssignLabel::labelMin;          constant = 0;          }
  else if (strncmp(_curParam, "max", 4)            == 0)   { assign = merylAssignLabel::labelMax;          constant = 0;          }

  else if (strncmp(_curParam, "and", 4)            == 0)   { assign = merylAssignLabel::labelAnd;          constant = kmlablmax;  }
  else if (strncmp(_curParam, "or", 3)             == 0)   { assign = merylAssignLabel::labelOr;           constant = 0;          }
  else if (strncmp(_curParam, "xor", 4)            == 0)   { assign = merylAssignLabel::labelXor;          constant = kmlablmax;  }

  else if (strncmp(_curParam, "difference", 11)    == 0)   { assign = merylAssignLabel::labelDifference;   constant = 0;          }

  else if (strncmp(_curParam, "lightest", 9)       == 0)   { assign = merylAssignLabel::labelLightest;     constant = kmlablmax;  }
  else if (strncmp(_curParam, "heaviest", 9)       == 0)   { assign = merylAssignLabel::labelHeaviest;     constant = 0;          }

  else if (strncmp(_curParam, "invert", 7)         == 0)   { assign = merylAssignLabel::labelHeaviest;     constant = kmlablmax;  }

  else if (strncmp(_curParam, "shift-left", 11)    == 0)   { assign = merylAssignLabel::labelShiftLeft;    constant = 1;          }
  else if (strncmp(_curParam, "shift-right", 12)   == 0)   { assign = merylAssignLabel::labelShiftRight;   constant = 1;          }

  else if (strncmp(_curParam, "rotate-left", 12)   == 0)   { assign = merylAssignLabel::labelRotateLeft;   constant = 1;          }
  else if (strncmp(_curParam, "rotate-right", 13)  == 0)   { assign = merylAssignLabel::labelRotateRight;  constant = 1;          }

  //  Nope, don't know what this is.

  else {
    sprintf(_errors, "Unknown value assign:label=<parameter> in '%s'.", _optString);
    return false;
  }

  getCurrent()->_labelAssign   = assign;
  getCurrent()->_labelConstant = constant;

  resetClass();

  return true;
}



bool
merylCommandBuilder::isAssign(void) {

  if (_curClass != opClass::clAssign)      //  Silently ignore non-assign classes.
    return false;

  if ((_curPname != opPname::pnValue) &&   //  Noisly complain about unknown Pnames.
      (_curPname != opPname::pnLabel)) {
    sprintf(_errors, "expecting value or label parameter name in '%s'\n", _optString);
    return false;
  }

  if ((_curParam == nullptr) ||            //  Noisly complain about missing parameters.
      (_curParam[0] == 0)) {
    sprintf(_errors, "expecting parameters in '%s'\n", _optString);
    return false;
  }

  return ((isAssignValue() == true) ||     //  Process the parameter.
          (isAssignLabel() == true));
}
