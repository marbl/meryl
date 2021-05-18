
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


//  In merylOp-count.C
uint64
findMaxInputSizeForMemorySize(uint32 kMerSize, uint64 memorySize);



//  Everything is initialized in the declaration.  Nothing really to do here.
merylCommandBuilder::merylCommandBuilder() {
  _allowedThreads = getMaxThreadsAllowed();   //  Absolute maximum limits on
  _allowedMemory  = getMaxMemoryAllowed();    //  memory= and threads= values
}



merylCommandBuilder::~merylCommandBuilder() {

  for (uint32 ii=0; ii<_opRoot.size(); ii++) {
    uint32  rr = _opRoot[ii];

    for (uint32 tt=0; tt<64; tt++)   //  Destroy threads first.
      delete _thList[tt][rr];

    delete _opList[rr];              //  Then destroy the master.
  }

  for (uint32 tt=0; tt<64; tt++)
    delete [] _thList[tt];
}



void
merylCommandBuilder::terminateOperation(void) {

  for (; _terminating > 0; _terminating--) {
    if (merylOperation::_verbosity >= sayConstruction)
      fprintf(stderr, "terminate()-  Pop operation '%s' from top of stack.\n", toString(_opStack.top()->getOperation()));
    _opStack.pop();
  }
}



void
merylCommandBuilder::setAllowedMemory(char const *memstr) {
  _allowedMemory  = (uint64)(strtodouble(memstr) * 1024.0 * 1024.0 * 1024.0);
}



void
merylCommandBuilder::setAllowedThreads(char const *thrstr) {
  _allowedThreads = strtouint64(thrstr);

  omp_set_num_threads(_allowedThreads);
}



//  Initialize, either for a new operation or just for a new option string.
//   - Strip off any leading '['s, and count the number of closing ']'s.
//   - Save some copies of the stripped string, converted to file paths.
//   - Make sure there is an operation on the stack.
//
void
merylCommandBuilder::initialize(char *opt) {

  if (merylOperation::_verbosity >= sayConstruction)
    fprintf(stderr, "initialize()-\n");

  //  Process any existing left-over termination requests.

  terminateOperation();

  //  Save a copy of the string.

  _optStringLen = 0;

  while ((_optStringLen < FILENAME_MAX) && (opt[_optStringLen] != 0)) {
    _optString[_optStringLen] = opt[_optStringLen];
    _optStringLen++;
  }

  _optString[_optStringLen] = 0;

  //  Ignore '[' at the start of the string.  Their purpose is to visually
  //  match the ']' which tells us to stop adding inputs to the current
  //  command.  There should only be one opening bracket.

  if (_optString[0] == '[') {
    for (uint32 ii=0; ii<_optStringLen; ii++)
      _optString[ii] = _optString[ii+1];

    _optStringLen--;
  }

  //  If we have a ']' as the last character, strip it off and remember that
  //  we need to close the command on the stack after we process this arg.
  //  We can get any number of closing brackets.

  while ((_optStringLen > 0) &&
         (_optString[_optStringLen-1] == ']')) {
    if (merylOperation::_verbosity >= sayConstruction)
      fprintf(stderr, "initialize()- Found terminator.\n");

    _optString[_optStringLen-1] = 0;
    _optStringLen--;

    _terminating++;
  }

  //  Save a few copies of the command line word.

  strncpy(_inoutName, _optString, FILENAME_MAX + 1);

  snprintf(_indexName, FILENAME_MAX, "%s/merylIndex", _optString);
  snprintf(_sqInfName, FILENAME_MAX, "%s/info",       _optString);
  snprintf(_sqRdsName, FILENAME_MAX, "%s/reads",      _optString);

  //  If the stack is empty, push on a new operation.

  if (_opStack.size() == 0) {
    merylOperation  *op = new merylOperation;

    if (merylOperation::_verbosity >= sayConstruction)
      fprintf(stderr, "initialize()- Add new empty root operation at position %zu.\n", _opList.size());

    _opStack.push(op);
    _opList.push_back(op);
    _opRoot.push_back(_opList.size() - 1);
  }

  if (merylOperation::_verbosity >= sayConstruction)
    fprintf(stderr, "initialize()- Got option '%s' for '%s' at stack level %zu.\n",
            _optString,
            toString(_opStack.top()->getOperation()),
            _opStack.size());
}




bool
merylCommandBuilder::processLabel(char const *val, bool failIfBad) {
  merylLabelAction  action = merylLabelAction::labelNothing;
  merylLabelFilter  filter = merylLabelFilter::labelNothing;
  kmlabl            label  = 0;

  //  Decide what label action this is.  

  if      (strncmp(val, "#", 1)           == 0)   { action = merylLabelAction::labelConstant;     label = 0;          }
  else if (strncmp(val, "@", 1)           == 0)   { action = merylLabelAction::labelSelected;     label = kmvalumax;  }

  else if (strncmp(val, "first", 5)       == 0)   { action = merylLabelAction::labelFirst;        label = 0;          }

  else if (strncmp(val, "and", 3)         == 0)   { action = merylLabelAction::labelAnd;          label = kmvalumax;  }
  else if (strncmp(val, "or", 2)          == 0)   { action = merylLabelAction::labelOr;           label = 0;          }
  else if (strncmp(val, "xor", 3)         == 0)   { action = merylLabelAction::labelXor;          label = 0;          }
  else if (strncmp(val, "difference", 4)  == 0)   { action = merylLabelAction::labelDifference;   label = 0;          }
  else if (strncmp(val, "lightest", 4)    == 0)   { action = merylLabelAction::labelLightest;     label = kmvalumax;  }
  else if (strncmp(val, "heaviest", 4)    == 0)   { action = merylLabelAction::labelHeaviest;     label = 0;          }

  else if (strncmp(val, "invert", 6)      == 0)   { action = merylLabelAction::labelInvert;       label = kmvalumax;  }
  else if (strncmp(val, "shift-left", 7)  == 0)   { action = merylLabelAction::labelShiftLeft;    label = 0;          }
  else if (strncmp(val, "left-shift", 7)  == 0)   { action = merylLabelAction::labelShiftLeft;    label = 0;          }
  else if (strncmp(val, "shift-right", 7) == 0)   { action = merylLabelAction::labelShiftRight;   label = 0;          }
  else if (strncmp(val, "right-shift", 7) == 0)   { action = merylLabelAction::labelShiftRight;   label = 0;          }

  else {
    if (failIfBad == true) {
      fprintf(stderr, "ERROR: unknown label action '%s' from option '%s'.\n", val, _optString);
      exit(1);
    } else {
      return(false);
    }
  }

  //  Scan the val for any constant and the type of constant.
  //  If there is a constant, decode it.
  //
  //  The constant can be either from a labelSet operation or tacked onto
  //  the others as a parameter, e.g., and#010010b

  uint32 valLen        = strlen(val);
  uint32 constPos      = 0;
  uint32 shiftVal      = 0;
  bool   invalidNumber = false;

  if (val[valLen-1] == 'b')   shiftVal = 1;
  if (val[valLen-1] == 'o')   shiftVal = 3;
  if (val[valLen-1] == 'd')   shiftVal = 10;
  if (val[valLen-1] == 'h')   shiftVal = 4;

  if (('0' <= val[valLen-1]) &&
      (val[valLen-1] <= '9'))
    shiftVal = 10;

  for (uint32 vv=0; vv<valLen; vv++)
    if (val[vv] == '#')
      constPos = vv+1;

  //  If we're expecting a value, decode it (or complain).

  if ((constPos > 0) && (shiftVal == 0)) {
    fprintf(stderr, "ERROR: expected number in option '%s' but didn't find type of number ('b', 'o' or 'h').\n", _optString);
    exit(1);
  }

  if ((constPos == 0) && (shiftVal > 0)) {
    fprintf(stderr, "ERROR: expected number in option '%s' but didn't find '#' marking start of number.\n", _optString);
    exit(1);
  }

  if ((constPos > 0) && (shiftVal == 10)) {
    label = strtouint64(val + constPos);
    shiftVal = 0;
  }

  if ((constPos > 0) && (shiftVal > 0)) {
    uint8  decode[256];

    for (uint32 qq=0; qq<256; qq++)
      decode[qq] = 0xff;

    if (shiftVal >= 1) {
      decode['0'] = 0x00;
      decode['1'] = 0x01;
    }

    if (shiftVal >= 3) {
      decode['2'] = 0x02;
      decode['3'] = 0x03;
      decode['4'] = 0x04;
      decode['5'] = 0x05;
      decode['6'] = 0x06;
      decode['7'] = 0x07;
    }

    if (shiftVal >= 4) {
      decode['8'] = 0x08;
      decode['9'] = 0x09;
      decode['a'] = 0x0a;   decode['A'] = 0x0a;
      decode['b'] = 0x0b;   decode['B'] = 0x0b;
      decode['c'] = 0x0c;   decode['C'] = 0x0c;
      decode['d'] = 0x0d;   decode['D'] = 0x0d;
      decode['e'] = 0x0e;   decode['E'] = 0x0e;
      decode['f'] = 0x0f;   decode['F'] = 0x0f;
    }

    for (uint32 vv=constPos; vv<valLen-1; vv++) {
      uint8 v = decode[val[vv]];

      if (v == 0xff)
        invalidNumber = true;

      label <<= shiftVal;
      label  |= v;
    }
  }

  //  Set the label action, or complain it's bogus.

  if (invalidNumber == true) {
    fprintf(stderr, "ERROR: invalid number digits in option '%s'.\n", _optString);
    exit(1);
  }

  _opStack.top()->setLabelAction(action, filter, label);

  return(true);
}


bool
merylCommandBuilder::processOptions(void) {

  //  If the string is entirely a number, treat it as either a threshold or a
  //  constant, depending on the operation.  This is used for things like
  //  "greater-than 45" and "divide 2".
  //
  //  If there is no operation, or it doesn't want a number, we fall through
  //  and return 'false' when key/val is checked below.

  bool  isNum = isDecNumber(_optString, 0);

  if ((_opStack.top()->needsThreshold() == true) && (isNum == true)) {
    _opStack.top()->setThreshold(strtouint64(_optString));
    return(true);
  }

  if ((_opStack.top()->needsConstant() == true) && (isNum == true)) {
    _opStack.top()->setConstant(strtouint64(_optString));
    return(true);
  }

  //  Some options have no value, unfortunately.

  if (strcmp(_optString, "compress") == 0) {
    _doCompression = true;
    return(true);
  }

  //  If we're a label operation, check for a label, but don't fail if we
  //  fail to find one.

  //if ((_opStack.top()->getOperation() == merylOp::opLabel) &&
  //    (_opStack.top()->_labelAction == labelNothing)) {
  //  if (processLabel(_optString, true) == true)
  //    return(true);
  //}

  //  If we get something that looks like a number, fail.

  if (_optString[0] == '#') {
    fprintf(stderr, "ERROR: suspicious option '%s' encountered; is this part of a 'label' action?\n", _optString);
    exit(1);
  }

  //  The rest should be in a key=value structure.  If we don't find that,
  //  just return.

  KeyAndValue   kv(_optString);
  char const   *key = kv.key();
  char const   *val = kv.value();

  if ((key == NULL) ||
      (val == NULL))
    return(false);

  uint32 valLen = strlen(val);
  uint32 val32  = strtouint32(val);
  uint64 val64  = strtouint64(val);
  double valDB  = strtodouble(val);


#if 0
  //  Kmer size.

  ////  This is from meryl master.  git merge wanted to put it here, but doesn't
  ////  seem to belong.

  if (strcmp(key, "k") == 0) {
    if ((kmerTiny::merSize() != 0) &&
        (kmerTiny::merSize() != val32)) {
      fprintf(stderr, "ERROR: kmer size mismatch: %u != %u\n", kmerTiny::merSize(), val32);
      exit(1);
    }
    kmerTiny::setSize(val32);
#endif


  //  Label processing
  //
  //  The label constant can be supplied in multiple forms, or not even
  //  present:
  //
  //    label=#00101010b   (or [0-7]o or [0-9A-G]h)
  //    label=and#00111b
  //    label=and

  if (strcmp(key, "label") == 0) {
    processLabel(val, true);
    return(true);
  }

  //  Number of kmers expected for counting.
  if (strcmp(key, "n") == 0) {
    _opStack.top()->setExpectedNumberOfKmers(val64);
    return(true);
  }

  //  A suffix to filter kmers by when counting.
  if (strcmp(key, "count-suffix") == 0) {
    _opStack.top()->setCountSuffix(val);
    return(true);
  }

  //  Threshold values for less-than, etc, specifed as a fraction of the
  //  total distinct kmers, or as a word-frequency, or as an absolute count.

  if ((strcmp(key, "d") == 0) ||
      (strcmp(key, "distinct") == 0)) {
    _opStack.top()->setFractionDistinct(valDB);
    return(true);
  }

  if ((strcmp(key, "f") == 0) ||
      (strcmp(key, "word-frequency") == 0)) {
    _opStack.top()->setWordFrequency(valDB);
    return(true);
  }

  if ((strcmp(key, "t") == 0) ||            //  See above for special case of this!
      (strcmp(key, "threshold") == 0)) {
    _opStack.top()->setThreshold(val64);
    return(true);
  }

  //  Segment of input, for counting from seqStore.  Useless otherwise.
  if ((strcmp(key, "segment") == 0) &&
      (isDecNumber(val, '/'))) {
    decodeRange(val, _segment, _segmentMax);
#ifndef CANU
    fprintf(stderr, "WARNING: option '%s' ignored, available only with Canu support.\n", _optString);
#endif
    return(true);
  }

  //  If nothing triggered, we don't recognize it as an option.  Maybe the
  //  filename had an '=' in it?

  return(false);
}



//  The present-in flag is a bit complicated.  It can accept any number three
//  different formats, separated by colons.
//
//    present-in:@<digits>[:...]
//              :#<digits>[:...]
//              :#<digits>-#<digits>[:...]
//
//  The parsing uses a state machine to keep track of what it's looking for.
//  The only non-obvious bit is using the next state variable to decide
//  if we've got one or two # values - basically checking why we're
//  here - and either setting a single value or a range of values.
//
//  This function parses the option string and sets merylOperation members
//    _selectedFile
//    _presentInLen
//    _presentInMax
//    _presentIn
//

void
setPresentIn(uint32 lo, uint32 hi, uint32 *&pIn, uint32 &pInLen, uint32 &pInMax) {
  if (pInLen <= hi)
    resizeArray(pIn, pInLen, pInMax, hi+1, _raAct::copyData | _raAct::clearNew);

  for (uint32 ii=lo; ii<=hi; ii++)
    pIn[ii] = 1;
}

void
merylCommandBuilder::parsePresentIn(void) {
  uint32 constexpr  EXPECT_NEW = 0x001;
  uint32      EXPECT_HASH = 0x002;
  uint32      EXPECT_AT   = 0x004;
  uint32      EXPECT_DASH = 0x008;

  char   const     *str        = _optString + 10;
  uint32            expectNext = EXPECT_NEW;

  uint32            presentBgn = 0;

  delete [] _presentIn;

#if 0
  _presentInLen = 0;
  _presentInMax = 0;
  _presentIn    = nullptr;

  _selectedFile = UINT32_MAX;
#endif

  if (0 != strncmp(_optString, "present-in:", 11))
    assert(0);

  while (*str) {
    if ((expectNext & EXPECT_NEW)  && (*str == ':')) {
      str++;
      expectNext = EXPECT_HASH | EXPECT_AT;
    }

    else if ((expectNext & EXPECT_AT) && (*str == '@')) {
      str = strtonumber(str+1, _selectedFile);
      expectNext = EXPECT_NEW;
    }

    else if ((expectNext & EXPECT_HASH) && (*str == '#')) {
      uint32 number = 0;

      str = strtonumber(str+1, number);

      if (expectNext == EXPECT_HASH) {
        setPresentIn(presentBgn, number, _presentIn, _presentInLen, _presentInMax);
      } else {
        presentBgn = number;
        setPresentIn(number, number, _presentIn, _presentInLen, _presentInMax);
      }

      expectNext = EXPECT_DASH | EXPECT_NEW;
    }

    else if ((expectNext & EXPECT_DASH) && (*str == '-')) {
      str++;
      expectNext = EXPECT_HASH;
    }

    else {
      char const  *errStr = "(the unexpected)";

      if (expectNext == (EXPECT_HASH | EXPECT_AT))   errStr = "'#' or '@'";
      if (expectNext == (EXPECT_NEW))                errStr = "':'";
      if (expectNext == (EXPECT_DASH | EXPECT_NEW))  errStr = "'-' or ':'";
      if (expectNext == (EXPECT_HASH))               errStr = "'#'";

      fprintf(stderr, "\n");
      fprintf(stderr, "Unexpected character '%c' found in '%s'; expecting %s.\n", *str, _optString, errStr);
      fprintf(stderr, "                                   %*s^\n", (int)(str - _optString), "");
      exit(1);
    }
  }
}




bool
merylCommandBuilder::processOperation(void) {
  merylOp           nOp  = merylOp::opNothing;

  merylValueAction  vAct = merylValueAction::valueNothing;
  merylValueFilter  vFil = merylValueFilter::valueNothing;
  kmvalu            valu = kmvaluzero;

  merylLabelAction  lAct = merylLabelAction::labelNothing;
  merylLabelFilter  lFil = merylLabelFilter::labelNothing;
  kmlabl            labl = kmlablzero;

  uint32            presentInLen   = 0;
  uint32           *presentIn      = nullptr;
  uint32            foundIn        = UINT32_MAX;
  bool              presentInAny   = false;         //  Kmer should be present in any input.
  bool              presentInAll   = false;         //  Kmer should be present in all inputs.
  bool              presentInFirst = false;         //  Kmer should be present in the first and any other inputs.
  bool              presentInOnly  = false;         //  Kmer should be present only in the first input.

  assert(_opStack.size() > 0);

  //  If the string is of length zero, explicitly do nothing and return
  //  success.  This was (probably) just an isolated terminating bracket.

  if      (0 == _optStringLen)
    return(true);

  //  Check for an operation string.

  //  Input from fasta, kmer is computed, count is computed, label can be set to a constant.
  if      (0 == strcmp(_optString, "count"))                  { nOp = merylOp::opCount;        }
  else if (0 == strcmp(_optString, "count-forward"))          { nOp = merylOp::opCountForward; }
  else if (0 == strcmp(_optString, "count-reverse"))          { nOp = merylOp::opCountReverse; }

  //  Input from single source, modify the value or count, or filter based on value or count.
  else if (0 == strcmp(_optString, "modify"))                 { nOp = merylOp::opPresentIn; }
  else if (0 == strcmp(_optString, "filter"))                 { nOp = merylOp::opPresentIn; }

  //  The generic operation.
  //    present-in:#1-#4:#6:@1:first
  //
  else if (0 == strncmp(_optString, "present-in:", 11)) {
    parsePresentIn();
  }

  //  Input from multiple sources.  count and label pick from one of the sources.
  else if (0 == strcmp(_optString, "union"))                  { nOp = merylOp::opPresentIn;  presentInAny = true;  vAct = merylValueAction::valueCount;  lAct = merylLabelAction::labelFirst; }
  else if (0 == strcmp(_optString, "union-min"))              { nOp = merylOp::opPresentIn;  presentInAny = true;  vAct = merylValueAction::valueMin;    lAct = merylLabelAction::labelMin;   }
  else if (0 == strcmp(_optString, "union-max"))              { nOp = merylOp::opPresentIn;  presentInAny = true;  vAct = merylValueAction::valueMax;    lAct = merylLabelAction::labelMax;   }
  else if (0 == strcmp(_optString, "union-sum"))              { nOp = merylOp::opPresentIn;  presentInAny = true;  vAct = merylValueAction::valueSum;    lAct = merylLabelAction::labelOr;    }

  else if (0 == strcmp(_optString, "intersect"))              { nOp = merylOp::opPresentIn;  presentInAll = true;  vAct = merylValueAction::valueFirst;  lAct = merylLabelAction::labelFirst; }
  else if (0 == strcmp(_optString, "intersect-min"))          { nOp = merylOp::opPresentIn;  presentInAll = true;  vAct = merylValueAction::valueMin;    lAct = merylLabelAction::labelMin;   }
  else if (0 == strcmp(_optString, "intersect-max"))          { nOp = merylOp::opPresentIn;  presentInAll = true;  vAct = merylValueAction::valueMax;    lAct = merylLabelAction::labelMax;   }
  else if (0 == strcmp(_optString, "intersect-sum"))          { nOp = merylOp::opPresentIn;  presentInAll = true;  vAct = merylValueAction::valueSum;    lAct = merylLabelAction::labelAnd;   }





  //  Input from a single source, modify the value.
  else if (0 == strcmp(_optString, "increase"))               { nOp = merylOp::opPresentIn;  presentInOnly = true;  vAct = merylValueAction::valueIncrease;     }
  else if (0 == strcmp(_optString, "decrease"))               { nOp = merylOp::opPresentIn;  presentInOnly = true;  vAct = merylValueAction::valueDecrease;     }
  else if (0 == strcmp(_optString, "multiply"))               { nOp = merylOp::opPresentIn;  presentInOnly = true;  vAct = merylValueAction::valueMultiply;     }
  else if (0 == strcmp(_optString, "divide"))                 { nOp = merylOp::opPresentIn;  presentInOnly = true;  vAct = merylValueAction::valueDivideTrunc;  }
  else if (0 == strcmp(_optString, "divide-round"))           { nOp = merylOp::opPresentIn;  presentInOnly = true;  vAct = merylValueAction::valueDivideRound;  }
  else if (0 == strcmp(_optString, "remainder"))              { nOp = merylOp::opPresentIn;  presentInOnly = true;  vAct = merylValueAction::valueRemainder;    }

  //  Input from a single source, filter based on value.
  else if (0 == strcmp(_optString, "less-than"))              { nOp = merylOp::opPresentIn;  presentInOnly = true;  vFil = merylValueFilter::valueLessThan;     }
  else if (0 == strcmp(_optString, "greater-than"))           { nOp = merylOp::opPresentIn;  presentInOnly = true;  vFil = merylValueFilter::valueGreaterThan;  }
  else if (0 == strcmp(_optString, "at-least"))               { nOp = merylOp::opPresentIn;  presentInOnly = true;  vFil = merylValueFilter::valueAtLeast;      }
  else if (0 == strcmp(_optString, "at-most"))                { nOp = merylOp::opPresentIn;  presentInOnly = true;  vFil = merylValueFilter::valueAtMost;       }
  else if (0 == strcmp(_optString, "equal-to"))               { nOp = merylOp::opPresentIn;  presentInOnly = true;  vFil = merylValueFilter::valueEqualTo;      }
  else if (0 == strcmp(_optString, "not-equal-to"))           { nOp = merylOp::opPresentIn;  presentInOnly = true;  vFil = merylValueFilter::valueNotEqualTo;   }



  //  subtract   - subtracts [2..n] from the first file, threshold to zero.
  //  difference - 
  else if (0 == strcmp(_optString, "subtract"))               { }

  else if (0 == strcmp(_optString, "difference"))             { }
  else if (0 == strcmp(_optString, "symmetric-difference"))   { }

  //  Input from a single source, output is descriptive statistics not a meryl db.
  else if (0 == strcmp(_optString, "histogram"))              { nOp = merylOp::opHistogram;  }
  else if (0 == strcmp(_optString, "statistics"))             { nOp = merylOp::opStatistics; }

  else if (0 == strcmp(_optString, "compare"))                { nOp = merylOp::opCompare;    }

  else return(false);   //  optString is not an operation.

  //  If the top-of-stack command is counting, pop it off and possibly make a
  //  new command.  Counting operations cannot take input from a command.
  //  We're also guaranteed to have either an empty stack, or a valid
  //  non-counting operation on the stack after this.

  if (_opStack.top()->isCounting() == true) {
    _opStack.pop();

    if (_opStack.size() == 0) {
      merylOperation   *op = new merylOperation();

      if (merylOperation::_verbosity >= sayConstruction)
        fprintf(stderr, "processOp()-  Add new empty root operation at position %zu.\n", _opList.size());

      _opStack.push(op);
      _opList.push_back(op);
      _opRoot.push_back(_opList.size() - 1);
    }

    assert(_opStack.top()->isCounting() == false);
  }

  //  If there is a valid command on the stack, push a new command onto the
  //  stack, and add it to the inputs list.

  if (_opStack.top()->getOperation() != merylOp::opNothing) {
    merylOperation   *op = new merylOperation();

    _opStack.top()->addInputFromOp(op);

    _opStack.push(op);
    _opList.push_back(op);

    if (merylOperation::_verbosity >= sayConstruction)
      fprintf(stderr, "processOp()-  Operation '%s' added to stack at level %zu\n",
              toString(nOp), _opStack.size());
  } else {
    if (merylOperation::_verbosity >= sayConstruction)
      fprintf(stderr, "processOp()-  Operation '%s' replaces '%s' at level %zu\n",
              toString(nOp), toString(_opStack.top()->getOperation()), _opStack.size());
  }

  //  Setup this operation.

  _opStack.top()->setOperation(nOp);
  _opStack.top()->setValueAction(vAct, vFil, valu);
  _opStack.top()->setLabelAction(lAct, lFil, labl);

  return(true);  //  optString was an operation.
}




bool
merylCommandBuilder::isOutput(void) {

  //  If we see 'output', flag the next arg as the output name.
  if (strcmp(_optString, "output") == 0) {
    _isOutput = true;
    return(true);
  }

  //  If the flag isn't set, this wasn't 'output' nor is it an output name.
  if (_isOutput == false)
    return(false);

  //  Must be an output name.  Reset the flag, and add the output
  //  to whatever operation is current.

  _isOutput   = false;

  _opStack.top()->addOutput(_inoutName);

  return(true);
}



bool
merylCommandBuilder::isPrinter(void) {

  //  If we see 'print' or 'printACGT', flag the next arg as the output name.
  if (strcmp(_optString, "print") == 0) {
    _printACGTorder = false;
    _isPrint        = true;
    return(true);
  }

  if (strcmp(_optString, "printACGT") == 0) {
    _printACGTorder = true;
    _isPrint        = true;
    return(true);
  }

  //  If the flag isn't set, this isn't a printer output name.
  if (_isPrint == false)
    return(false);

  //  Must be a printer name.  Reset the flag.  Add the name, unless it looks
  //  like a meryl database; in this case, the user requested 'print
  //  db.meryl' and output should go to stdout.  We need to add the printer
  //  to this operation, and return false ("this is not a printer option") so
  //  we can properly handle the database.

  _isPrint = false;

  if (fileExists(_indexName) == true) {
    _opStack.top()->addPrinter(nullptr, _printACGTorder);
    return(false);
  }

  _opStack.top()->addPrinter(_inoutName, _printACGTorder);
  return(true);
}




bool
merylCommandBuilder::isMerylInput(void) {

  if (fileExists(_indexName) == false)
    return(false);

  _opStack.top()->addInputFromDB(_inoutName);

  return(true);
}

bool
merylCommandBuilder::isCanuInput(std::vector<char *> &err) {

  if ((fileExists(_sqInfName) == false) ||
      (fileExists(_sqRdsName) == false))
    return(false);

#ifndef CANU
  char *s = new char [FILENAME_MAX + 129];
  snprintf(s, FILENAME_MAX + 129, "Detected seqStore input '%s', but no Canu support available.", _inoutName);
  err.push_back(s);
#endif

  _opStack.top()->addInputFromCanu(_inoutName, _segment, _segmentMax);

  _segment    = 1;
  _segmentMax = 1;

  return(true);
}

bool
merylCommandBuilder::isSequenceInput(void) {

  if (fileExists(_inoutName) == false)
    return(false);

  _opStack.top()->addInputFromSeq(_inoutName, _doCompression);

  return(true);
}





void
merylCommandBuilder::finalize(void) {

  //  Finish processing the last option string.
  terminateOperation();

  //  If no operation supplied, assume it's a print, and make it print all kmers.
#if 0
  if ((_opStack.size() > 0) &&
      (_opStack.top()->getOperation() == opNothing)) {
    if (merylOperation::_verbosity >= sayConstruction)
      fprintf(stderr, "finalize()- Change opNothing to opLessThan at stack level %zu.\n", _opStack.size());
    _opStack.top()->setOperation(opLessThan);
  }
#endif

  //  Clear the stack.
  while (_opStack.size() > 0)
    _opStack.pop();

  //  Update memory and threads for everything.
  for (uint32 oo=0; oo<_opList.size(); oo++) {
    _opList[oo]->_maxMemory  = _allowedMemory;
    _opList[oo]->_maxThreads = _allowedThreads;
  }
}



void
merylCommandBuilder::printTree(merylOperation *op, uint32 indent) {

  fprintf(stderr, "%*s%-s\n", indent, "", toString(op->getOperation()));

  if (op->_valueAction != merylValueAction::valueNothing)
    fprintf(stderr, "%*slabel=%s constant=%x\n", indent+2, "", toString(op->_valueAction), op->_valueConstant);

  if (op->_valueFilter != merylValueFilter::valueNothing)
    fprintf(stderr, "%*slabel=%s constant=%x\n", indent+2, "", toString(op->_valueFilter), op->_valueConstant);

  if (op->_labelAction != merylLabelAction::labelNothing)
    fprintf(stderr, "%*slabel=%s constant=%lx\n", indent+2, "", toString(op->_labelAction), op->_labelConstant);

  if (op->_labelFilter != merylLabelFilter::labelNothing)
    fprintf(stderr, "%*slabel=%s constant=%lx\n", indent+2, "", toString(op->_labelFilter), op->_labelConstant);


#if 0
  if (op->_mathConstant > 0)
    fprintf(stderr, "%*sconstant=%lu\n", indent+2, "", op->_mathConstant);

  if (op->_threshold != UINT64_MAX)
    fprintf(stderr, "%*sthreshold=%lu\n", indent+2, "", op->_threshold);

  if (op->_fracDist != DBL_MAX)
    fprintf(stderr, "%*sfraction-distinct=%f\n", indent+2, "", op->_fracDist);

  if (op->_wordFreq != DBL_MAX)
    fprintf(stderr, "%*sword-frequenct=%f\n", indent+2, "", op->_wordFreq);
#endif

  for (uint32 ii=0; ii<op->_inputs.size(); ii++) {
    merylInput  *in = op->_inputs[ii];

    if (in->isFromOperation() == true) {
      printTree(in->_operation, indent+2);
    }

    if (in->isFromDatabase() == true) {
      fprintf(stderr, "%*s%s\n", indent+2, "", in->_name);
    }

    if (in->isFromSequence() == true) {
      fprintf(stderr, "%*s%s%s\n", indent+2, "", in->_name, in->_homopolyCompress ? " (homopoly compressed)" : "");
    }

    if (in->isFromStore() == true) {
      fprintf(stderr, "%*s%s (reads %u through %u)\n", indent+2, "", in->_name, in->_sqBgn, in->_sqEnd);
    }
  }

  if (op->_outputO) {
    fprintf(stderr, "%*soutput to %s\n", indent+2, "", op->_outputO->filename());
  }

  if (op->_printerName) {
    fprintf(stderr, "%*sprint to %s\n", indent+2, "", op->_printerName);
  }
}



//  Clone the command tree(s) into thread-specific copies, one tree per thread.
//
//
void
merylCommandBuilder::spawnThreads(void) {
  uint32  indent = 0;

  omp_set_num_threads(_allowedThreads);

  for (uint32 tt=0; tt<64; tt++) {

    //  Construct a list of operations for each thread.
    _thList[tt] = new merylOperation * [_opList.size()];

    //  Copy operations from the main list to our thread list.
    for (uint32 oo=0; oo<_opList.size(); oo++)
      _thList[tt][oo] = new merylOperation(_opList[oo],
                                           tt,
                                           _opList[oo]->_inputs.size(),
                                           _allowedThreads, _allowedMemory);

    //  Update all the input/output files to be per-thread.
    for (uint32 oo=0; oo<_opList.size(); oo++) {
      merylOperation  *op = _thList[tt][oo];   //  The per-thread operation we're fixing up.
      merylOperation  *OP = _opList[oo];       //  The master operation we're using as a template.

      for (uint32 ii=0; ii<OP->_inputs.size(); ii++) {
        merylInput  *IN = OP->_inputs[ii];     //  The template input for this operation.

        //  If the template input is from an operation, we need to search for
        //  that operation in the master list of operations, then set the
        //  per-thread input to be the corresponding operation.
        if (IN->isFromOperation() == true) {
          uint32  inop = UINT32_MAX;

          for (uint32 xx=0; xx<_opList.size(); xx++)   //  Search all template operations for
            if (IN->_operation == _opList[xx])         //  the one that is our input.
              inop = xx;

          if (inop == UINT32_MAX)
            fprintf(stderr, "Failed to find corresponding operation.\n"), exit(1);

          op->addInputFromOp(_thList[tt][inop]);       //  Add input from the same op in our list.
        }

        //  If the template input is from a database, make a new input for
        //  just the piece we're processing in this thread (done implicitly
        //  in addInputFromDB()).
        if (IN->isFromDatabase() == true) {
          op->addInputFromDB(IN->_name);
        }

        //  We should never get inputs from a sequence file.
        if (IN->isFromSequence() == true) {
          assert(0);
          continue;
        }

        //  We should never get inputs from a Canu seqStore file.
        if (IN->isFromStore() == true) {
          assert(0);
          continue;
        }
      }
    }
  }
}
