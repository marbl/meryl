
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

#ifndef MERYL_H
#define MERYL_H

#include "runtime.H"
#include "system.H"
#include "strings.H"
#include "bits.H"

#include "kmers.H"
#include "sequence.H"

#include <stack>
#include <vector>


#include <stdarg.h>
inline
void
addError(std::vector<char const *> &ev, char const *fmt, ...) {
  char    *str = new char [1024];
  va_list  ap;

  va_start(ap, fmt);
  vsnprintf(str, 1024, fmt, ap);
  va_end(ap);

  ev.push_back(str);
}



//  A global (variable 'verbosity') to tell what messages to output.
//
class merylVerbosity {
public:
  merylVerbosity()  {};
  ~merylVerbosity() {};

  void   beQuiet(void)               { _v = 0; }
  void   increaseVerbosity(void)     { _v++;   }

  void   enableProgressReport(void)  { _p = true; }

  bool   showNothing     (void)   { return(_v == 0); }
  bool   showStandard    (void)   { return(_v >= 1); }
  bool   showProgress    (void)   { return(_v >= 2); }
  bool   showConstruction(void)   { return(_v >= 3); }
  bool   showDetails     (void)   { return(_v >= 4); }
  bool   showEverything  (void)   { return(_v >= 5); }

private:
  uint32 _v = 1;       //  Verbosity level.
  bool   _p = false;   //  Show progress, does nothing yet.
};

extern merylVerbosity verbosity;


//  Each node in the tree is assigned a type.  Most of the nodes will be
//  opFilter.  The special cases are to make counts, output a histogram, and
//  so on.
//
//  opPrint is used for commands like 'meryl print input.meryl', where the
//  action never gets set up because we never see any action specification.
//
enum class merylOpType {
  opNothing,
  opCounting,
  opStatistics,
  opHistogram,
  opPrint,
  opFilter
};

inline
char const *
toString(merylOpType t) {
  switch (t) {
    case merylOpType::opNothing:     return("opNothing");     break;
    case merylOpType::opCounting:    return("opCounting");    break;
    case merylOpType::opStatistics:  return("opStatistics");  break;
    case merylOpType::opHistogram:   return("opHistogram");   break;
    case merylOpType::opPrint:       return("opPrint");       break;
    case merylOpType::opFilter:      return("opFilter");      break;
  }
  assert(0);
  return(nullptr);
}



#define MERYLINCLUDE

#include "merylInput.H"

#include "merylModify.H"
#include "merylFilter.H"

#include "merylOpCounting.H"
#include "merylOpTemplate.H"
#include "merylOpCompute.H"

#include "merylCountArray.H"
#include "merylCommandBuilder.H"

#undef  MERYLINCLUDE


#endif  //  MERYL_H
