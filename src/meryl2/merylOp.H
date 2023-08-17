
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
