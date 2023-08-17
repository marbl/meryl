
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

#ifndef MERYL_COMMAND_BUILDER_H
#define MERYL_COMMAND_BUILDER_H

#include "meryl.H"


enum class opClass { clNone, clOutput, clAssign, clSelect, clInput };
enum class opPname { pnNone, pnDB, pnList, pnShow, pnPipe, pnAction, pnHisto, pnStats, pnValue, pnLabel, pnBases, pnInput };

class merylOpTemplate;
class merylOpCompute;

inline
char const *toString(opClass c) {
  char const *r = "invalid";

  switch (c) {
    case opClass::clNone:    r = "clNone";    break;
    case opClass::clOutput:  r = "clOutput";  break;
    case opClass::clAssign:  r = "clAssign";  break;
    case opClass::clSelect:  r = "clSelect";  break;
    case opClass::clInput:   r = "clInput";   break;
  }

  return r;
}

inline
char const *toString(opPname n) {
  char const *r = "invalid";

  switch (n) {
    case opPname::pnNone:    r = "pnNone";    break;
    case opPname::pnDB:      r = "pnDB";      break;
    case opPname::pnList:    r = "pnList";    break;
    case opPname::pnShow:    r = "pnShow";    break;
    case opPname::pnPipe:    r = "pnPipe";    break;
    case opPname::pnAction:  r = "pnAction";  break;
    case opPname::pnHisto:   r = "pnHisto";   break;
    case opPname::pnStats:   r = "pnStats";   break;
    case opPname::pnValue:   r = "pnValue";   break;
    case opPname::pnLabel:   r = "pnLabel";   break;
    case opPname::pnBases:   r = "pnBases";   break;
    case opPname::pnInput:   r = "pnInput";   break;
  }

  return r;
}


class merylCommandBuilder {
public:
  merylCommandBuilder()  {
  }
  ~merylCommandBuilder() {
    delete [] _optString;
    for (uint32 ss=0; ss<64; ss++)   //  Delete the arrays of pointers to compute
      delete [] _thList[ss];         //  objects; the objects are deleted elsewhere
  }

private:
  void    addNewOperation(void);
  bool    terminateOperations(uint32 nter=1, bool pop=false);

private:
  bool    processValueModifier(char const *val, bool failIfBad);
  bool    processLabelModifier(char const *val, bool failIfBad);
  bool    processBasesModifier(char const *val, bool failIfBad);

  uint32                 isRelation(uint32 bgn);
  merylSelectorRelation  decodeRelation(uint32 bgn);

  void                   decodeSelector(uint32 bgn,
                                        merylSelector &f);

public:
  void    printTree(uint32 tID, uint32 indent=0);
  void    printTree(merylOpTemplate *op, uint32 tID, uint32 iID, uint32 indent);
  bool    displayTreesAndErrors(void);

private:
  opClass     _curClass = opClass::clNone;
  opPname     _curPname = opPname::pnNone;
  char const *_curParam = nullptr;

  bool    decodeWord(char const *opt);
  void    resetClass(void) {
    _curClass = opClass::clNone;
    _curPname = opPname::pnNone;
    _curParam = nullptr;
  }

private:
  uint32  findCreateTerminate(void);
public:
  void    processWord(char const *opt);         //  Add a command line word to the builder
  void    buildTrees(void);                     //  Notify us that all words have been added
  void    performCounting(uint64 allowedMemory,
                          uint32 allowedThreads);
  void    spawnThreads(uint32 allowedThreads);  //  Build 64 trees for processing kmers.
  void    runThreads(uint32 allowedThreads);

private:
  bool    isEmpty(void)      { return(_optString[0] == 0); }
  bool    isOption(void);
  bool    isAlias(void);

  bool    isInOut(void);
  bool    isInput(void);

  bool    isAssignValue(void);
  bool    isAssignLabel(void);
  bool    isAssign(void);
  bool    isSelect(void);

  bool    isValueSelector(void);
  bool    isLabelSelector(void);
  bool    isBasesSelector(void);
  bool    isInputSelector(void);

  bool    isSelectorConnective(void);

  bool    isCount(void);

  bool    isHistogram(void);
  bool    isStatistics(void);



public:
  uint32           numOperations(void)           { return(_opList.size()); };
  merylOpTemplate *getOperation(uint32 i)        { return(_opList[i]);     };

  uint32           numTrees(void)                { return(_opTree.size()); };
  merylOpTemplate *getTree(uint32 r)             { return(_opList   [ _opTree[r] ]);  };
  merylOpCompute  *getTree(uint32 r, uint32 t)   { return(_thList[t][ _opTree[r] ]);  };

  merylOpTemplate *getCurrent(void)              { assert(_opStack.size() > 0); return(_opStack.top()); };

  //  Collecting any errors encountered when building the operation tree.
public:
  uint32                      numErrors(void)   {  return(_errors.size());  }
  std::vector<char const *>  &getErrors(void)   {  return(_errors);         }

  std::vector<char const *>  _errors;

  //  Storage and flags we use for processing a command line word.
private:
  uint32    _optStringLen = 0;              //  A copy of the command line word,
  uint32    _optStringMax = 0;              //  with any '[' and ']' removed.
  char     *_optString = nullptr;           //

  uint32    _segment       = 1;             //  for Canu seqStore inputs
  uint32    _segmentMax    = 1;             //

  bool      _doCompression = false;         //  for FASTA or FASTQ file inputs

  bool      _isOutput      = false;         //  enable output to merylDB

  bool      _isPrint        = false;        //  enable output to ASCII dump
  bool      _printACGTorder = false;        //  fix canonical ordering

  bool      _isHistogram    = false;        //  enable histogram output
  bool      _isStatistics   = false;        //  enable statistics output

  //  If _needsValue or _needsConstant is set, the current selector or current
  //  action needs an integer parameter in the next word.

  bool      _needsValue     = false;
  bool      _needsConstant  = false;

  //  If _invertNextSelector is set, the sense of the next selector encountered
  //  is inverted.  Example: 's1 and not s2'.

  bool      _invertNextSelector = false;

  //  _opStack is a stack of operations, used when constructing the tree of operations.
  //  _opList is a list of operations.

  std::stack <merylOpTemplate *>    _opStack;
  std::vector<merylOpTemplate *>    _opList;
  merylOpCompute                  **_thList[64] = { nullptr };   //  Mirrors opList

  std::vector<uint32>               _opTree;
};

#endif
