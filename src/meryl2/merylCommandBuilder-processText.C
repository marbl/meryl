
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





char const *begin   = "\n*({c}";
char const *end     = ")\n*\\]?";

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

char const *path    = "[\t[:print:]]+";    //  Pathnames are <TAB> plus all the printable characters.





void
merylCommandBuilder::compileRegExes(void) {

  //for (uint32 rr=0; rr<0xff; rr++)
  //  _regex[rr].enableVerbose();

  //egex[0x00] is reserved for 'no match'.
  _regex[0x01].compile(begin, "\\[", end, nullptr);
  _regex[0x02].compile(begin, "\\]", end, nullptr);

  //  COUNTING

  _regex[0x10].compile(begin, "count(-({c}({p}forward)|({p}reverse)))?", end, nullptr);

  _regex[0x11].compile(begin, "n=(", integer, ")", end, nullptr);
  _regex[0x12].compile(begin, "count-suffix", equ, "([ACGTacgt]+)", end, nullptr);

  //  ALIASES
  //
  //  union,     union-min,     union-max,     union-sum
  //  intersect, intersect-min, intersect-max, intersect-sum
  //
  _regex[0x20].compile(begin, "union"                 "(-({c}min|max|sum))?", end, nullptr);
  _regex[0x21].compile(begin, "(isect|in({p}tersect))""(-({c}min|max|sum))?", end, nullptr);

  _regex[0x22].compile(begin, "sub({p}tract)?",          end, nullptr);
  _regex[0x23].compile(begin, "diff({p}erence)?",        end, nullptr);

  //  subtract
  //  difference
  //
  //  less-than        Z | threshold=Z | distinct=F | word-frequency=F | wf=F
  //  greater-than     Z | threshold=Z | distinct=F | word-frequency=F | wf=F
  //  at-least         Z | threshold=Z | distinct=F | word-frequency=F | wf=F
  //  at-most          Z | threshold=Z | distinct=F | word-frequency=F | wf=F
  //  equal-to         Z | threshold=Z | distinct=F | word-frequency=F | wf=F
  //  not-equal-to     Z | threshold=Z | distinct=F | word-frequency=F | wf=F
  //
  //  increase         Z | threshold=Z | distinct=F | word-frequency=F | wf=F
  //  decrease         Z | threshold=Z | distinct=F | word-frequency=F | wf=F
  //  multiply     F | Z | threshold=Z | distinct=F | word-frequency=F | wf=F
  //  divide       F | Z | threshold=Z | distinct=F | word-frequency=F | wf=F
  //  divide-round F | Z | threshold=Z | distinct=F | word-frequency=F | wf=F
  //  modulo           Z | threshold=Z | distinct=F | word-frequency=F | wf=F
  //


  char tdw[1024] = {0};
  sprintf(tdw, "(%s|({c,p}threshold)%s%s|({c,p}distinct)%s%s|({c,p}w[of]rd-frequency)%s%s)",
          integer, equ, integer, equ, baseF, equ, baseF);

  //_regex[0x24].enableVerbose();
  _regex[0x24].compile(begin, "less-than\\s+",    tdw, end, nullptr);
  _regex[0x25].compile(begin, "greater-than\\s+", tdw, end, nullptr);
  _regex[0x26].compile(begin, "at-least\\s+",     tdw, end, nullptr);
  _regex[0x27].compile(begin, "at-most\\s+",      tdw, end, nullptr);
  _regex[0x28].compile(begin, "equal-to\\s+",     tdw, end, nullptr);
  _regex[0x29].compile(begin, "not-equal-to\\s+", tdw, end, nullptr);

  _regex[0x2a].compile(begin, "increase\\s+",     integer, end, nullptr);
  _regex[0x2b].compile(begin, "decrease\\s+",     integer, end, nullptr);

  _regex[0x2c].compile(begin, "multiply\\s+",     number, end, nullptr);
  _regex[0x2d].compile(begin, "divide\\s+",       number, end, nullptr);
  _regex[0x2e].compile(begin, "divide-round\\s+", number, end, nullptr);

  _regex[0x2f].compile(begin, "modulo\\s+",       integer, end, nullptr);

  //  Output parameters.
  //
  //  Compatibility mode 'output <path>' must be last so it doesn't match to 'output :db=file'.

  _regex[0x31].compile(begin, "({p}output)", col, "({p}d[ab]tabase)"   "", equ, "({c}", path, ")", end, nullptr);
  _regex[0x32].compile(begin, "({p}output)", col, "({p}list)"          "", equ, "({c}", path, ")", end, nullptr);
  _regex[0x33].compile(begin, "({p}output)", col, "({p}show)",                                     end, nullptr);
  _regex[0x34].compile(begin, "({p}output)", col, "(pipe)"             "", equ, "({c}", path, ")", end, nullptr);
  _regex[0x35].compile(begin, "({p}output)", col, "(histogram)"        "", equ, "({c}", path, ")", end, nullptr);
  _regex[0x36].compile(begin, "({p}output)", col, "(histogram)",                                   end, nullptr);
  _regex[0x37].compile(begin, "({p}output)", col, "({p}stat[is]stics)" "", equ, "({c}", path, ")", end, nullptr);
  _regex[0x38].compile(begin, "({p}output)", col, "({p}stat[is]stics)",                            end, nullptr);
  _regex[0x30].compile(begin,    "(output)\\s+({c}", path, ")",                                    end, nullptr);

  //  Assign.

  _regex[0x40].compile(begin, "(({p}assign)|set)", col, "({p}value)" "", equ, "([^\n]*)", end, nullptr);
  _regex[0x50].compile(begin, "(({p}assign)|set)", col, "({p}label)" "", equ, "([^\n]*)", end, nullptr);

  //  Select.

  _regex[0x60].compile(begin, "(({p}select)|get)", col, "({p}value)" "", col, "([^\n]*)", end, nullptr);
  _regex[0x61].compile(begin, "(({p}select)|get)", col, "({p}label)" "", col, "([^\n]*)", end, nullptr);
  _regex[0x62].compile(begin, "(({p}select)|get)", col, "({p}bases)" "", col, "([^\n]*)", end, nullptr);
  _regex[0x63].compile(begin, "(({p}select)|get)", col, "({p}input)" "", col, "([^\n]*)", end, nullptr);

  //  Input parameters.

  _regex[0x70].compile(begin, "({p}input)", col, "({p}database)"       "", equ, "({c}", path, ")", end, nullptr);
  _regex[0x71].compile(begin, "({p}input)", col, "({p}list)"           "", equ, "({c}", path, ")", end, nullptr);
  _regex[0x72].compile(begin, "({p}input)", col, "({p}pipe)"           "", equ, "({c}", path, ")", end, nullptr);
  _regex[0x73].compile(begin, "({p}input)", col, "({p}action)"         "", equ, "({c}", path, ")", end, nullptr);

  _regexLen = 0x74;
}


uint32
merylCommandBuilder::matchProgramText(uint64 pp) {
  uint64 dLen = 0;         //  For display of program text in errors.
  uint64 dMax = 0;         //
  char  *d    = nullptr;   //

  fprintf(stderr, "pTxt/%02lu: '%s'\n", pp, displayString(_pTxt + pp, d, dLen, dMax));

  //  Search all regexes for a match, return the first one.
  //  All state is saved in _regex; all we need to know is
  //  which one fired.

  for (uint32 rr=0; rr<_regexLen; rr++)
    if (_regex[rr].match(_pTxt + pp) == true) {
      for (uint32 ii=0; ii<_regex[rr].numCaptures(); ii++)
        if (_regex[rr].getCaptureValid(ii))
          fprintf(stderr, "0x%02x/%02u: %03lu-%03lu '%s'\n", rr, ii,
                  _regex[rr].getCaptureBgn(ii),
                  _regex[rr].getCaptureEnd(ii),
                  displayString(_regex[rr].getCapture(ii), d, dLen, dMax));
      fprintf(stderr, "\n");
      return rr;
    }

  return 0;
}


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

    //  If no match, emit an error and move to the next letter.
    if (rr == 0) {
      continue;
    }

    //  Otherwise, process the match.

    merylOpTemplate  *op = getCurrent();

    char const *C[8];
    for (uint32 ii=0; ii<_regex[rr].numCaptures(); ii++)
      C[ii] = _regex[rr].getCapture(ii);


    switch (rr) {
      case 0x01:   //  [
        terminateOperations(1);   //  Making a new operation.  Close the existing one,
        addNewOperation();        //  and make a new one.
        break;

      case 0x02:   //  ]
        if (terminateOperations(1, true) == false)
          sprintf(_errors, "processWord()- Extra ']' encountered in command line.");
        break;

      case 0x10:   //  count, count-forward, count-reverse
        if      (C[2][0] == 'f') {
          fprintf(stderr, "COUNT FORWARD\n");
        }
        else if (C[2][0] == 'r') {
          fprintf(stderr, "COUNT REVERSE\n");
        }
        else {
          fprintf(stderr, "COUNT\n");
        }
        break;

      case 0x11:   //  n=<integer>
        if (op->_isCounting == true)
          op->_counting->setExpectedNumberOfKmers(decodeInteger(C[2], 0, 0, v64, _errors));
        else
          sprintf(_errors, "option '%s' encountered for non-counting operation.", C[1]);
       
        break;

      case 0x12:   //  count-suffix=[ACGTacgt]
        if (op->_isCounting == true)
          op->_counting->setCountSuffix(C[2]);
        else
          sprintf(_errors, "option '%s' encountered for non-counting operation.", C[1]);
        break;

      case 0x20:   //  union
        break;
      case 0x21:   //  intersect
        break;
      case 0x22:   //  subtract
        break;
      case 0x23:   //  difference
        break;
      case 0x24:   //  less-than
        break;
      case 0x25:   //  greater-than
        break;
      case 0x26:   //  at-least
        break;
      case 0x27:   //  at-most
        break;
      case 0x28:   //  equal-to
        break;
      case 0x29:   //  not-equal-to
        break;
      case 0x2a:   //  increase
        break;
      case 0x2b:   //  decrease
        break;
      case 0x2c:   //  multiply
        break;
      case 0x2d:   //  divide
        break;
      case 0x2e:   //  divide-round
        break;
      case 0x2f:   //  modulo
        break;

      case 0x30:
        break;
      case 0x31:
        break;
      case 0x32:
        break;
      case 0x33:
        break;
      case 0x34:
        break;
      case 0x35:
        break;
      case 0x36:
        break;
      case 0x37:
        break;
      case 0x38:
        break;
      case 0x39:
        break;
      case 0x3a:
        break;
      case 0x3b:
        break;
      case 0x3c:
        break;
      case 0x3d:
        break;
      case 0x3e:
        break;
      case 0x3f:
        break;


      case 0x50:
        break;
      case 0x51:
        break;
      case 0x52:
        break;
      case 0x53:
        break;
      case 0x54:
        break;
      case 0x55:
        break;
      case 0x56:
        break;
      case 0x57:
        break;
      case 0x58:
        break;
      case 0x59:
        break;
      case 0x5a:
        break;
      case 0x5b:
        break;
      case 0x5c:
        break;
      case 0x5d:
        break;
      case 0x5e:
        break;
      case 0x5f:
        break;


      case 0x60:
        break;
      case 0x61:
        break;
      case 0x62:
        break;
      case 0x63:
        break;
      case 0x64:
        break;
      case 0x65:
        break;
      case 0x66:
        break;
      case 0x67:
        break;
      case 0x68:
        break;
      case 0x69:
        break;
      case 0x6a:
        break;
      case 0x6b:
        break;
      case 0x6c:
        break;
      case 0x6d:
        break;
      case 0x6e:
        break;
      case 0x6f:
        break;


    }

    if (globals.showConstruction() == true)
      fprintf(stderr, "----------\n");

    pp += _regex[rr].getCaptureEnd(1) - 1;   //  -1 because the for loops adds one by default.
  }

  delete [] d;
}
