
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

#include "system.H"

using namespace merylutil;

uint64  _pTxtLen = 0;
uint64  _pTxtMax = 0;
uint64  _pTxtPos = 0;
char   *_pTxt    = nullptr;



//  Distinguishing filename "file[1]" from command "[alias file]" is tricky.
void
appendProgramWord(char const *w) {

  if (w == nullptr)
    return;

  while (*w) {                                           //  Copy word to program text.
    increaseArray(_pTxt, _pTxtLen+2, _pTxtMax, 16384);   //  Get space for this letter and
    _pTxt[_pTxtLen++] = *w++;                            //  the '\n\0' terminating letters,
  }                                                      //  grabbing 16KB more as needed.

  _pTxt[_pTxtLen++] = '\n';
  _pTxt[_pTxtLen]   = '\0';
}



char const *begin   = "\n*({c}";
char const *end     = ")\n*\\]?";

char const *base2   = "({c}-?[01]+[bB])";
char const *base8   = "({c}-?[01234567]+[oO])";
char const *baseDd  = "({c}-?[0123456789]+[dD]?)";
char const *baseDs  = "({c}-?[0123456789]+[kKmMgGtTpPeE][iI]?)";
char const *baseH   = "({c}-?[0123456789abcdefABCDEF]+[hH])";
char const *baseF   = "({c}-?[0123456789]*[.]?[0123456789]+))";

char const *integer = "({c}(-?[01]+[bB])|"
                          "(-?[01234567]+[oO])|"
                          "(-?[0123456789]+[dD]?)|"
                          "(-?[0123456789]+[kKmMgGtTpPeE][iI]?)|"
                          "(-?[0123456789abcdefABCDEF]+[hH])"
                      ")";
char const *number  = "({c}(-?[01]+[bB])|"
                          "(-?[01234567]+[oO])|"
                          "(-?[0123456789]+[dD]?)|"
                          "(-?[0123456789]+[kKmMgGtTpPeE][iI]?)|"
                          "(-?[0123456789abcdefABCDEF]+[hH])|"
                          "(-?[0123456789]*[.]?[0123456789]+)"
                      ")";

char const *path    = "[\t[:print:]]+";    //  Pathnames are <TAB> plus all the printable characters.


int
main(int argc, char **argv) {
  v2::regEx  re[256];

  //for (uint32 rr=0; rr<0xff; rr++)
  //  re[rr].enableVerbose();

  re[0x01].compile(begin, "\\[", end, nullptr);
  re[0x02].compile(begin, "\\]", end, nullptr);

  //  COUNTING

  re[0x10].compile(begin, "(count)(-(({c,p}forward)|({c,p}reverse)))?", end, nullptr);

  re[0x11].compile(begin, "(n)\\s*=\\s*(", integer, ")", end, nullptr);
  re[0x12].compile(begin, "(count-suffix)\\s*=\\s*([ACGTacgt]+)", end, nullptr);

  //  ALIASES
  //
  //  union,     union-min,     union-max,     union-sum
  //  intersect, intersect-min, intersect-max, intersect-sum
  //
  //  subtract
  //  difference
  //
  //  less-than        Z | threshold=Z | distinct=F | word-frequency=F
  //  greater-than     Z | threshold=Z | distinct=F | word-frequency=F
  //  at-least         Z | threshold=Z | distinct=F | word-frequency=F
  //  at-most          Z | threshold=Z | distinct=F | word-frequency=F
  //  equal-to         Z | threshold=Z | distinct=F | word-frequency=F
  //  not-equal-to     Z | threshold=Z | distinct=F | word-frequency=F
  //
  //  increase         Z | threshold=Z | distinct=F | word-frequency=F
  //  decrease         Z | threshold=Z | distinct=F | word-frequency=F
  //  multiply     F | Z | threshold=Z | distinct=F | word-frequency=F
  //  divide       F | Z | threshold=Z | distinct=F | word-frequency=F
  //  divide-round F | Z | threshold=Z | distinct=F | word-frequency=F
  //  modulo           Z | threshold=Z | distinct=F | word-frequency=F
  //

  re[0x20].compile(begin, "union(-({c}min|max|sum))?", end, nullptr);

  re[0x21].compile(begin, "(isect|in({p}tersect))(-({c}(min|max|sum)))?", end, nullptr);

  re[0x22].compile(begin, "sub({p}tract)?",          end, nullptr);
  re[0x23].compile(begin, "diff({p}erence)?",        end, nullptr);

  re[0x24].compile(begin, "less-than\\s+",    integer, end, nullptr);
  re[0x25].compile(begin, "greater-than\\s+", integer, end, nullptr);
  re[0x26].compile(begin, "at-least\\s+",     integer, end, nullptr);
  re[0x27].compile(begin, "at-most\\s+",      integer, end, nullptr);
  re[0x28].compile(begin, "equal-to\\s+",     integer, end, nullptr);
  re[0x29].compile(begin, "not-equal-to\\s+", integer, end, nullptr);

  re[0x2a].compile(begin, "increase\\s+",     integer, end, nullptr);
  re[0x2b].compile(begin, "decrease\\s+",     integer, end, nullptr);

  re[0x2c].compile(begin, "multiply\\s+",     number, end, nullptr);
  re[0x2d].compile(begin, "divide\\s+",       number, end, nullptr);
  re[0x2e].compile(begin, "divide-round\\s+", number, end, nullptr);

  re[0x2f].compile(begin, "modulo\\s+",       integer, end, nullptr);

  //

  re[0x30].compile(begin, "({p}output)\\s*:\\s*({p}d[ab]tabase)"   "\\s*=\\s*([^\n]*)", end, nullptr);
  re[0x31].compile(begin, "({p}output)\\s*:\\s*({p}list)"          "\\s*=\\s*([^\n]*)", end, nullptr);
  re[0x32].compile(begin, "({p}output)\\s*:\\s*({p}show)",                              end, nullptr);
  re[0x33].compile(begin, "({p}output)\\s*:\\s*(pipe)"             "\\s*=\\s*([^\n]*)", end, nullptr);
  re[0x34].compile(begin, "({p}output)\\s*:\\s*(histogram)"        "\\s*=\\s*([^\n]*)", end, nullptr);
  re[0x35].compile(begin, "({p}output)\\s*:\\s*(histogram)",                            end, nullptr);
  re[0x36].compile(begin, "({p}output)\\s*:\\s*({p}stat[is]stics)" "\\s*=\\s*([^\n]*)", end, nullptr);
  re[0x37].compile(begin, "({p}output)\\s*:\\s*({p}stat[is]stics)",                     end, nullptr);

  re[0x40].compile(begin, "(({p}assign)|set)\\s*:\\s*({p}value)" "\\s*=\\s*([^\n]*)", end, nullptr);
  re[0x41].compile(begin, "(({p}assign)|set)\\s*:\\s*({p}label)" "\\s*=\\s*([^\n]*)", end, nullptr);

  re[0x50].compile(begin, "(({p}select)|get)\\s*:\\s*({p}value)" "\\s*:\\s*([^\n]*)", end, nullptr);
  re[0x51].compile(begin, "(({p}select)|get)\\s*:\\s*({p}label)" "\\s*:\\s*([^\n]*)", end, nullptr);
  re[0x52].compile(begin, "(({p}select)|get)\\s*:\\s*({p}bases)" "\\s*:\\s*([^\n]*)", end, nullptr);
  re[0x53].compile(begin, "(({p}select)|get)\\s*:\\s*({p}input)" "\\s*:\\s*([^\n]*)", end, nullptr);

  re[0x60].compile(begin, "({p}input)\\s*:\\s*({p}database)"       "\\s*=\\s*([^\n]*)", end, nullptr);
  re[0x61].compile(begin, "({p}input)\\s*:\\s*({p}list)"           "\\s*=\\s*([^\n]*)", end, nullptr);
  re[0x62].compile(begin, "({p}input)\\s*:\\s*({p}pipe)"           "\\s*=\\s*([^\n]*)", end, nullptr);
  re[0x63].compile(begin, "({p}input)\\s*:\\s*({p}action)"         "\\s*=\\s*([^\n]*)", end, nullptr);



  for (int arg=1; arg<argc; arg++)
    appendProgramWord(argv[arg]);

  uint64 dLen = 0;
  uint64 dMax = 0;
  char  *d    = nullptr;

  while (_pTxtPos < _pTxtLen) {
    bool   matched = false;

    fprintf(stderr, "pTxt:    '%s'\n", displayString(_pTxt + _pTxtPos, d, dLen, dMax));


    for (uint32 rr=0; rr<256; rr++) {
      if (re[rr].match(_pTxt + _pTxtPos) == false)
        continue;

      matched = true;

      for (uint32 ii=0; ii<re[rr].numCaptures(); ii++)
        fprintf(stderr, "0x%02x/%02u: '%s'\n", rr, ii, displayString(re[rr].getCapture(ii), d, dLen, dMax));
      fprintf(stderr, "\n");

      switch (rr) {
        case 0x00:
          break;
        case 0x01:
          break;

        case 0x10:
          break;
        case 0x11:
          break;
        case 0x12:
          break;
      }

      _pTxtPos += re[rr].getCaptureEnd(1);
      break;   //  We could just continue, no need to go back to rr=0 and start again...
    }

    if (matched == false)
      _pTxtPos++;
  }

  return 0;
}

