
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

char const *begin   = "({c}";     //  Words start with a capture group.
char const *end     = ")[]\v]";   //  And end with a word-separator or an action close operator.

char const *col     = "\\s*:\\s*";
char const *equ     = "\\s*=\\s*";

char const *base2   = "({c}[+-]?[0-1]+[bB])";
char const *base8   = "({c}[+-]?[0-7]+[oO])";
char const *baseDd  = "({c}[+-]?[0-9]+[dD]?)";
char const *baseDs  = "({c}[+-]?[0-9]+[kKmMgGtTpPeE][iI]?)";
char const *baseH   = "({c}[+-]?[0-9a-fA-F]+[hH])";
char const *baseF   = "({c}[+-]?[0-9]*[.]?[0-9]+([eE][+-]?[0-9]+)?)";

//  baseF doesn't allow '0.'

char const *integer = "({c}([+-]?[01]+[bB])|"
                          "([+-]?[01234567]+[oO])|"
                          "([+-]?[0-9]+[dD]?)|"
                          "([+-]?[0-9]+[kKmMgGtTpPeE][iI]?)|"
                          "([+-]?[0-9a-fA-F]+[hH])"
                      ")";
char const *number  = "({c}([+-]?[01]+[bB])|"
                          "([+-]?[01234567]+[oO])|"
                          "([+-]?[0-9]+[dD]?)|"
                          "([+-]?[0-9]+[kKmMgGtTpPeE][iI]?)|"
                          "([+-]?[0-9a-fA-F]+[hH])|"
                          "([+-]?[0-9]*[.]?[0-9]+([eE][+-]?[0-9]+)?)"
                      ")";

char const *path    = "({c}[ \t[:print:]]+)";   //  Pathnames <SPACE>, <TAB> and all the printable characters.


void
merylCommandBuilder::compileRegExes(void) {

  _regex[0x00].compile("({c}\v+)");   //  Match word-seps, MUST be an extra capture group to fit with the rest.

  _regex[0x01].compile(begin, "\\[", end, nullptr);  //  Match a single '[' to start a new operation.
  _regex[0x02].compile(begin, "\\]", end, nullptr);  //  Match a single ']' to complete an operation.

  //egex[0x03];
  //egex[0x04];

  //  COUNTING and COUNT OPTIONS

  _regex[0x05].compile(begin, "count(-({c}({p}forward)|({p}reverse)))?(\v+n", equ, integer, ")?", end, nullptr);
  _regex[0x06].compile(begin, "count-suffix", equ, "([ACGTacgt]+)", end, nullptr);

  //egex[0x07];
  //egex[0x08];

  //  COMBINING ALIASES
  //   - subtract,   union,     union-min,     union-max,     union-sum
  //   - difference, intersect, intersect-min, intersect-max, intersect-sum
  //
  _regex[0x09].compile(begin, "union"                 "(-({c}min|max|sum))?", end, nullptr);
  _regex[0x0a].compile(begin, "(isect|in({p}tersect))""(-({c}min|max|sum))?", end, nullptr);

  _regex[0x0b].compile(begin, "sub({p}tract)?",          end, nullptr);
  _regex[0x0c].compile(begin, "diff({p}erence)?",        end, nullptr);

  //egex[0x0d];
  //egex[0x0e];
  //egex[0x0f];

  //  FILTERING ALIASES
  //   - less-than     <tdw>
  //   - greater-than  <tdw>
  //   - at-least      <tdw>
  //   - at-most       <tdw>
  //   - equal-to      <tdw>
  //   - not-equal-to  <tdw>
  //
  //  tdw matches:
  //   - <integer>
  //   - threshold=<integer>
  //   - distinct=<float>
  //   - word-frequency=<float>  (and wf=<float>)
  //
  char tdw[1024] = {0};
  sprintf(tdw, "(%s|({c,p}threshold)%s%s|({c,p}distinct)%s%s|({c,p}w[of]rd-frequency)%s%s)",
          integer, equ, integer, equ, baseF, equ, baseF);

  _regex[0x10].compile(begin, "less-than\\s+",    tdw, end, nullptr);
  _regex[0x11].compile(begin, "greater-than\\s+", tdw, end, nullptr);
  _regex[0x12].compile(begin, "at-least\\s+",     tdw, end, nullptr);
  _regex[0x13].compile(begin, "at-most\\s+",      tdw, end, nullptr);
  _regex[0x14].compile(begin, "equal-to\\s+",     tdw, end, nullptr);
  _regex[0x15].compile(begin, "not-equal-to\\s+", tdw, end, nullptr);

  //egex[0x16];
  //egex[0x17];

  //  MODIFYING ALIASES
  //   - increase      <integer>
  //   - decrease      <integer>
  //   - multiply      <float>
  //   - divide        <float>
  //   - divide-round  <float>
  //   - modulo        <integer>
  //
  _regex[0x18].compile(begin, "increase\\s+",     integer, end, nullptr);
  _regex[0x19].compile(begin, "decrease\\s+",     integer, end, nullptr);
  _regex[0x1a].compile(begin, "multiply\\s+",     number,  end, nullptr);
  _regex[0x1b].compile(begin, "divide\\s+",       number,  end, nullptr);
  _regex[0x1c].compile(begin, "divide-round\\s+", number,  end, nullptr);
  _regex[0x1d].compile(begin, "modulo\\s+",       integer, end, nullptr);

  //egex[0x1e];
  //egex[0x1f];

  //  Output parameters.
  //
  //  Compatibility mode 'output <path>' must be last so it doesn't match to 'output : db = file'.

  _regex[0x20].compile(begin, "({p}output)", col, "({p}d[ab]tabase)",        equ, path,       end, nullptr);
  _regex[0x21].compile(begin, "({p}output)", col, "({p}list)",               equ, path,       end, nullptr);
  _regex[0x22].compile(begin, "({p}output)", col, "({p}show)",                                end, nullptr);
  _regex[0x23].compile(begin, "({p}output)", col, "({p}pipe)",               equ, path,       end, nullptr);
  _regex[0x24].compile(begin, "({p}output)", col, "({p}histogram)",     "(", equ, path, ")?", end, nullptr);
  _regex[0x25].compile(begin, "({p}output)", col, "({p}stat[is]stics)", "(", equ, path, ")?", end, nullptr);

  _regex[0x26].compile(begin,    "(output)\\s",                                 path, end, nullptr);
  //egex[0x27];

  //  Input parameters.

  _regex[0x28].compile(begin, "({p}input)", col, "({p}database)", equ, path, end, nullptr);
  _regex[0x29].compile(begin, "({p}input)", col, "({p}list)",     equ, path, end, nullptr);
  _regex[0x2a].compile(begin, "({p}input)", col, "({p}pipe)",     equ, path, end, nullptr);
  _regex[0x2b].compile(begin, "({p}input)", col, "({p}action)",   equ, path, end, nullptr);

  _regex[0x2c].compile(begin, "({p}segment)\v+", baseDd, "/", baseDd, end, nullptr);  //  for canu
  //egex[0x2f];

  //  Assign.

  _regex[0x30].compile(begin, "(({p}assign)|set)", col, "({p}value)" "", equ, "([^\v]*)", end, nullptr);
  _regex[0x40].compile(begin, "(({p}assign)|set)", col, "({p}label)" "", equ, "([^\v]*)", end, nullptr);

  //  Select.
  //   - 'and' 'or' 'not'
  //
  _regex[0x50].compile(begin, "(({p}select)|get)", col, "({p}value)" "", col, "([^\v]*)", end, nullptr);
  _regex[0x60].compile(begin, "(({p}select)|get)", col, "({p}label)" "", col, "([^\v]*)", end, nullptr);
  _regex[0x70].compile(begin, "(({p}select)|get)", col, "({p}bases)" "", col, "([^\v]*)", end, nullptr);
  _regex[0x80].compile(begin, "(({p}select)|get)", col, "({p}input)" "", col, "([^\v]*)", end, nullptr);


  //  segment n/m for canu
#ifndef CANU
  //sprintf(_errors, "option '%s' available only with Canu support.", _optString);
#endif

  _regex[0xff].compile(begin, path, end, nullptr);

  _regexLen = 0x100;
}


//  Search all regexes for a match, return the first one.
//  All state is saved in _regex; all we need to know is
//  which one fired.
uint32
merylCommandBuilder::matchProgramText(uint64 pp) {
  uint32 rr=0;

  if (globals.showDetails() == true) {
    fprintf(stderr, "Search for command word in text starting at position pp=%lu:\n", pp);
    fprintf(stderr, "  %s\n", displayString(_pTxt + pp));
  }

  for (rr=0; rr<_regexLen; rr++)
    if (_regex[rr].match(_pTxt + pp) == true)
      break;

#if 0
  if (rr == 0)             //  If matching a word-separator, silently return success.
    return 0;

  if ((rr == 0) ||         //  If matching a word-separator (rr=0), a (potential)
      (rr == 0xff) ||      //
      (rr == _regexLen))   //  input file or nothing, return without dumping the match.
    return 0;

  //if (rr == 0xff) {        //  If a file match, require the file acually exist.
  //}
#endif

  if (globals.showDetails() == true) {
    for (uint32 ii=0; ii<_regex[rr].numCaptures(); ii++)
      fprintf(stderr, "  0x%02x/%02u: %s %03lu-%03lu '%s'\n", rr, ii,
              _regex[rr].isCaptureValid(ii) ? "valid" : "inval",
              _regex[rr].getCaptureBgn(ii),
              _regex[rr].getCaptureEnd(ii),
              displayString(_regex[rr].getCapture(ii)));
    fprintf(stderr, "\n");
  }

  return rr;
}


struct regexmatch {
  char const *_m = nullptr;
  uint32      _b = 0;
  uint32      _e = 0;
  uint32      _l = 0;
  bool        _v = false;
};


void
merylCommandBuilder::processProgramText(void) {
  uint64 dLen = 0;         //  For display of program text in errors.
  uint64 dMax = 0;         //
  char  *d    = nullptr;   //

  uint64      v64 = 0;
  double      vD  = 0;

  if (globals.showConstruction() == true)
    fprintf(stderr, "processProgramText()-\n");

  if (_opStack.size() == 0)    //  Add an empty operation onto the stack
    addNewOperation();         //  so we have space to work.


  for (uint64 pp=0; pp<_pTxtLen; pp++) {
    uint32 rr = matchProgramText(pp);

    //  For (potential) input matches, try to add the word as an input file.
    //  If it succeeds, we're done (we just need to do the rest of the loop
    //  do advance pp to the end of the word).  If it fais, reset rr
    //  to indicate a no-match word and fall through to reporting an error.
    //
    if ((rr == 0xff) && (isInput(_regex[rr].getCapture(1)) == false))
      rr = _regexLen;

    //  Report non-matching words as errors.
    //
    if (rr >= _regexLen) {
      char *str = new char [_pTxtLen];

      for (uint32 oo=0; (pp < _pTxtLen) && (_pTxt[pp] != '\v'); pp++) {
        char const *sym = displayLetter(_pTxt[pp]);

        while (*sym)
          str[oo++] = *sym++;

        str[oo] = 0;
      }

      sprintf(_errors, "ERROR: word '%s' not recognized.", str);
      continue;
    }

    //  Otherwise, process the match.
    //    [-] - the whole match, including any word separators (not copied)
    //    [0] - the whole match, without word separators
    //    [1] - 1st option
    //    [2] - 2nd option
    //    [3] - 3rd ...

    merylOpTemplate  *op = getCurrent();

    char const *Cm[8] = { nullptr };

    for (uint32 ii=1; ii<_regex[rr].numCaptures(); ii++)
      Cm[ii-1] = _regex[rr].getCapture(ii);

    //  Make sense of the match.

    switch (rr) {
      case 0x00:    //  Do nothing for word-separators at the start of the string.
        break;

      case 0x01: {  //  [
        terminateOperations(1);   //  Making a new operation.  Close the existing one,
        addNewOperation();        //  and make a new one.
      } break;
      case 0x02: {  //  ]
        if (terminateOperations(1, true) == false)
          sprintf(_errors, "processWord()- Extra ']' encountered in command line.");
      } break;

      case 0x03: {
      } break;
      case 0x04: {
      } break;

      case 0x05: {   //  count, count-forward, count-reverse n=<integer>
        if      (op->_isCounting)
          sprintf(_errors, "ERROR: operation is already a counting operation.\n");

        else if (op->_isSelector)
          sprintf(_errors, "ERROR: operation is a select operation, cannot also be a counting operation.\n");

        else {
          op->_isCounting = true;
          op->_counting   = new merylOpCounting(Cm[1][0]);
        }
      } break;
      case 0x06: {  //  count-suffix=[ACGTacgt]
        if (op->_isCounting == true)
          op->_counting->setCountSuffix(Cm[1]);
        else
          sprintf(_errors, "option '%s' encountered for non-counting operation.", Cm[0]);
      } break;

      case 0x07: {
      } break;
      case 0x08: {
      } break;


      case 0x09: {  //  union
        op->_valueAssign    = merylAssignValue::valueCount;
        op->_labelAssign    = merylAssignLabel::labelOr;

        merylSelector  f(merylSelectorQuantity::isIndex,
                         merylSelectorRelation::isNOP, false, "union");

        f._input_num_any    = true;

        op->addSelectorToProduct(f);

        if      ((Cm[1][0] == 'm') && (Cm[1][1] == 'i') && (Cm[1][2] == 'n')) {
          op->_valueAssign    = merylAssignValue::valueMin;
          op->_valueConstant  = kmvalumax;
          op->_labelAssign    = merylAssignLabel::labelSelected;
        }
        else if ((Cm[1][0] == 'm') && (Cm[1][1] == 'a') && (Cm[1][2] == 'x')) {
          op->_valueAssign    = merylAssignValue::valueMax;
          op->_valueConstant  = kmvalumin;
          op->_labelAssign    = merylAssignLabel::labelSelected;
        }
        else if ((Cm[1][0] == 's') && (Cm[1][1] == 'u') && (Cm[1][2] == 'm')) {
          op->_valueAssign    = merylAssignValue::valueAdd;
          op->_labelAssign    = merylAssignLabel::labelOr;
        }
        else if  (Cm[1][0] ==  0) {
          ;  //  Do nothing.
        }
        else {
          sprintf(_errors, "unknown union '%s' encountered.", Cm[0]);   //  Can't actually happen unles regex is busted.
        }
      } break;


      case 0x0a: {  //  intersect
        op->_valueAssign    = merylAssignValue::valueFirst;
        op->_labelAssign    = merylAssignLabel::labelAnd;

        merylSelector  f(merylSelectorQuantity::isIndex,
                         merylSelectorRelation::isNOP, false, "intersect");

        f._input_num_all    = true;

        op->addSelectorToProduct(f);

        if      ((Cm[1][0] == 'm') && (Cm[1][1] == 'i') && (Cm[1][2] == 'n')) {
          op->_valueAssign    = merylAssignValue::valueMin;
          op->_valueConstant  = kmvalumax;
          op->_labelAssign    = merylAssignLabel::labelSelected;
        }
        else if ((Cm[1][0] == 'm') && (Cm[1][1] == 'a') && (Cm[1][2] == 'x')) {
          op->_valueAssign    = merylAssignValue::valueMax;
          op->_valueConstant  = kmvalumin;
          op->_labelAssign    = merylAssignLabel::labelSelected;
        }
        else if ((Cm[1][0] == 's') && (Cm[1][1] == 'u') && (Cm[1][2] == 'm')) {
          op->_valueAssign    = merylAssignValue::valueAdd;
          op->_labelAssign    = merylAssignLabel::labelAnd;
        }
        else if  (Cm[1][0] ==  0) {
          ;  //  Do nothing.
        }
        else {
          sprintf(_errors, "unknown union '%s' encountered.", Cm[0]);   //  Can't actually happen unles regex is busted.
        }
      } break;


      case 0x0b: {   //  subtract
        op->_valueAssign    = merylAssignValue::valueSub;
        op->_labelAssign    = merylAssignLabel::labelDifference;

        merylSelector  f(merylSelectorQuantity::isIndex,
                         merylSelectorRelation::isNOP, false, "subtract");

        f._input_num_any    = true;
        f._input_idx.push_back(1);

        op->addSelectorToProduct(f);
      } break;

      case 0x0c: {  //  difference
        op->_valueAssign    = merylAssignValue::valueFirst;
        op->_labelAssign    = merylAssignLabel::labelFirst;

        merylSelector  f(merylSelectorQuantity::isIndex,
                         merylSelectorRelation::isNOP, false, "difference");

        f._input_num.push_back(1);
        f._input_idx.push_back(1);

        op->addSelectorToProduct(f);
      } break;

      case 0x0d: {
      } break;
      case 0x0e: {
      } break;
      case 0x0f: {
      } break;


      case 0x10:    //  less-than
      case 0x11:    //  greater-than
      case 0x12:    //  at-least
      case 0x13:    //  at-most
      case 0x14:    //  equal-to
      case 0x15: {  //  not-equal-to
        op->_inputsMin      = 1;
        op->_inputsMax      = 1;

        op->_valueAssign    = merylAssignValue::valueFirst;
        op->_labelAssign    = merylAssignLabel::labelFirst;

        merylSelectorRelation   rel = merylSelectorRelation::isNOP;

        if (rr == 0x10)   rel = merylSelectorRelation::isLt;
        if (rr == 0x11)   rel = merylSelectorRelation::isGt;
        if (rr == 0x12)   rel = merylSelectorRelation::isGeq;
        if (rr == 0x13)   rel = merylSelectorRelation::isLeq;
        if (rr == 0x14)   rel = merylSelectorRelation::isEq;
        if (rr == 0x15)   rel = merylSelectorRelation::isNeq;

        merylSelector  f(merylSelectorQuantity::isValue, rel, false, _optString);

        f._vIndex1 = 0;   //  Use value from output kmer.
        f._vValue2 = 0;   //  Reset below, or to distinct or word-freq.

        if      (Cm[1][0] == 'd')   f._vValue2Distinct = strtodouble(Cm[2]);
        else if (Cm[1][0] == 'w')   f._vValue2WordFreq = strtodouble(Cm[2]);
        else if (Cm[1][0] == 't')   f._vValue2         = decodeInteger(Cm[2], 0, 0, v64, _errors);
        else                        f._vValue2         = decodeInteger(Cm[1], 0, 0, v64, _errors);

        fprintf(stderr, "values %lu %f %f\n", f._vValue2, f._vValue2Distinct, f._vValue2WordFreq);

        op->addSelectorToProduct(f);
      } break;

      case 0x16: {
      } break;
      case 0x17: {
      } break;

      case 0x18:    //  increase
      case 0x19:    //  decrease
      case 0x1a:    //  multiply
      case 0x1b:    //  divide
      case 0x1c:    //  divide-round
      case 0x1d: {  //  modulo
        op->_inputsMin      = 1;
        op->_inputsMax      = 1;

        op->_valueAssign    = merylAssignValue::valueNOP;   //  Reset below.
        op->_valueConstant  = 0;                            //  This too.
        op->_labelAssign    = merylAssignLabel::labelFirst;

        if (rr == 0x18)    op->_valueAssign = merylAssignValue::valueAdd;
        if (rr == 0x19)    op->_valueAssign = merylAssignValue::valueSub;
        if (rr == 0x1a)    op->_valueAssign = merylAssignValue::valueMul;
        if (rr == 0x1b)    op->_valueAssign = merylAssignValue::valueDiv;
        if (rr == 0x1c)    op->_valueAssign = merylAssignValue::valueDivZ;
        if (rr == 0x1d)    op->_valueAssign = merylAssignValue::valueMod;

        //  Mul, Div and DivZ allow floats!
        op->_valueConstant = decodeInteger(Cm[1], 0, 0, v64, _errors);
      } break;

      case 0x1e: {
      } break;
      case 0x1f: {
      } break;


      case 0x20:  op->addOutputToDB  (Cm[1],        _errors);   break;  //  output:data=path
      case 0x21:  op->addOutputToList(Cm[1], false, _errors);   break;  //  output:list=path
      case 0x22:  op->addOutputToList(Cm[1], false, _errors);   break;  //  output:show
      case 0x23:  op->addOutputToPipe(Cm[1],        _errors);   break;  //  output:pipe=path
      case 0x24:  op->addHistoOutput (Cm[1],        _errors);   break;  //  output:histo(=path)
      case 0x25:  op->addStatsOutput (Cm[1],        _errors);   break;  //  output:stats(=path)
      case 0x26:  op->addOutputToDB  (Cm[1],        _errors);   break;  //  legacy output

      case 0x27: {
      } break;

      case 0x28:
      case 0x29:
      case 0x2a: {
        merylInput *in = new merylInput;
        bool        sx = false;

        if (rr == 0x27)   sx = in->registerMerylDB   (Cm[1], _errors);
        if (rr == 0x28)   sx = in->registerMerylList (Cm[1], _errors);
        if (rr == 0x29)   sx = in->registerActionPipe(Cm[1], _errors);
        // (rr == 0x2a)   sx = in->registerTemplate  (nullptr,   _errors);

        if (sx)
          op->addInput(in);
        else
          delete in;
      } break;

      case 0xff: {
      } break;

    }

    if (globals.showConstruction() == true)
      fprintf(stderr, "----------\n");

    pp += _regex[rr].getCaptureEnd(1) - 1;   //  -1 because the for loops adds one by default.
  }

  delete [] d;
}
