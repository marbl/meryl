
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


//  Handle value and label selection methods.
//    value=min#C
//    label=#C

bool
merylCommandBuilder::isValueSelect(void) {
  merylModifyValue     select   = merylModifyValue::valueNOP;
  kmvalu               constant = kmvalumax;
  char const          *str      = _optString + 6;

  //  If not a 'value=<something>' string, do nothing.  The optString isn't for us.

  if (strncmp(_optString, "value=", 6) != 0)
    return(false);

  //  Check for modifiers with constants at the end.

  if      (strncmp(str, "#", 1)           == 0)   { select = merylModifyValue::valueSet;        constant = decodeInteger(_optString, 6+1, 0, constant, _errors); }

  else if (strncmp(str, "first#", 6)      == 0)   { select = merylModifyValue::valueSelected;   constant = decodeInteger(_optString, 6+6, 0, constant, _errors); }
  else if (strncmp(str, "selected#", 9)   == 0)   { select = merylModifyValue::valueSelected;   constant = decodeInteger(_optString, 6+9, 0, constant, _errors); }

  else if (strncmp(str, "min#", 4)        == 0)   { select = merylModifyValue::valueMin;        constant = decodeInteger(_optString, 6+4, 0, constant, _errors); }
  else if (strncmp(str, "max#", 4)        == 0)   { select = merylModifyValue::valueMax;        constant = decodeInteger(_optString, 6+4, 0, constant, _errors); }

  else if (strncmp(str, "add#", 4)        == 0)   { select = merylModifyValue::valueAdd;        constant = decodeInteger(_optString, 6+4, 0, constant, _errors); }
  else if (strncmp(str, "sum#", 4)        == 0)   { select = merylModifyValue::valueAdd;        constant = decodeInteger(_optString, 6+4, 0, constant, _errors); }

  else if (strncmp(str, "sub#", 4)        == 0)   { select = merylModifyValue::valueSub;        constant = decodeInteger(_optString, 6+4, 0, constant, _errors); }
  else if (strncmp(str, "dif#", 3)        == 0)   { select = merylModifyValue::valueSub;        constant = decodeInteger(_optString, 6+3, 0, constant, _errors); }

  else if (strncmp(str, "mul#", 4)        == 0)   { select = merylModifyValue::valueMul;        constant = decodeInteger(_optString, 6+4, 0, constant, _errors); }

  else if (strncmp(str, "div#", 4)        == 0)   { select = merylModifyValue::valueDiv;        constant = decodeInteger(_optString, 6+4, 0, constant, _errors); }
  else if (strncmp(str, "divzero#", 8)    == 0)   { select = merylModifyValue::valueDivZ;       constant = decodeInteger(_optString, 6+8, 0, constant, _errors); }

  else if (strncmp(str, "mod#", 4)        == 0)   { select = merylModifyValue::valueMod;        constant = decodeInteger(_optString, 6+4, 0, constant, _errors); }
  else if (strncmp(str, "rem#", 4)        == 0)   { select = merylModifyValue::valueMod;        constant = decodeInteger(_optString, 6+4, 0, constant, _errors); }

  //  Check for modifiers without constants.  Set the constant to whatever
  //  the identity is for the given modifier.

  else if (strncmp(str, "first",  5)      == 0)   { select = merylModifyValue::valueSelected;   constant = 0;         }
  else if (strncmp(str, "selected",  8)   == 0)   { select = merylModifyValue::valueSelected;   constant = kmvalumax; }

  else if (strncmp(str, "min",  3)        == 0)   { select = merylModifyValue::valueMin;        constant = kmvalumax; }
  else if (strncmp(str, "max",  3)        == 0)   { select = merylModifyValue::valueMax;        constant = 0;         }

  else if (strncmp(str, "add",  3)        == 0)   { select = merylModifyValue::valueAdd;        constant = 0;         }
  else if (strncmp(str, "sum",  3)        == 0)   { select = merylModifyValue::valueAdd;        constant = 0;         }

  else if (strncmp(str, "sub",  3)        == 0)   { select = merylModifyValue::valueSub;        constant = 0;         }
  else if (strncmp(str, "dif",  2)        == 0)   { select = merylModifyValue::valueSub;        constant = 0;         }

  else if (strncmp(str, "mul",  3)        == 0)   { select = merylModifyValue::valueMul;        constant = 1;         }

  else if (strncmp(str, "div",  3)        == 0)   { select = merylModifyValue::valueDiv;        constant = 0;         }
  else if (strncmp(str, "divzero",  7)    == 0)   { select = merylModifyValue::valueDivZ;       constant = 0;         }

  else if (strncmp(str, "mod",  3)        == 0)   { select = merylModifyValue::valueMod;        constant = 0;         }
  else if (strncmp(str, "rem",  3)        == 0)   { select = merylModifyValue::valueMod;        constant = 0;         }

  else if (strncmp(str, "count",  5)      == 0)   { select = merylModifyValue::valueCount;      constant = 0;         }

  //  Nope, don't know what this is.

  else {
    addError("Unknown value selection modifier '%s'.\n", _optString);
  }

  //  Copy the decoded modifier to the select.

  getCurrent()->_valueSelect   = select;
  getCurrent()->_valueConstant = constant;

  return(true);
}
