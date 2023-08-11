
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

  uint32               isRelation(uint32 bgn);
  merylFilterRelation  decodeRelation(uint32 bgn);

  void                 decodeFilter(uint32 bgn,
                                    merylFilter &f);

public:
  void    printTree(merylOpTemplate *op, uint32 inputNum, uint32 indent);
  bool    displayTreesAndErrors(void);

private:
  enum class opClass { clNone, clOutput, clAssign, clSelect, clInput };
  enum class opPname { pnNone, pnDB, pnList, pnShow, pnPipe, pnAction, pnHisto, pnStats, pnValue, pnLabel, pnBases, pnInput };

  opClass     _curClass = opClass::clNone;
  opPname     _curPname = opPname::pnNone;
  char const *_curParam = nullptr;

  bool    decodeWord(char const *opt);

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

  bool    isAssignValue(char const *str);
  bool    isAssignLabel(char const *str);
  bool    isAssign(void);
  bool    isSelect(void);

  bool    isValueFilter(void);
  bool    isLabelFilter(void);
  bool    isBasesFilter(void);
  bool    isInputFilter(void);

  bool    isFilterConnective(void);

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

  //  If _needsValue or _needsConstant is set, the current filter or current
  //  action needs an integer parameter in the next word.

  bool      _needsValue     = false;
  bool      _needsConstant  = false;

  //  If _invertNextFilter is set, the sense of the next filter encountered
  //  is inverted.  Example: 'filter1 and not filter2'.

  bool      _invertNextFilter = false;

  //  _opStack is a stack of operations, used when constructing the tree of operations.
  //  _opList is a list of operations.

  std::stack <merylOpTemplate *>    _opStack;
  std::vector<merylOpTemplate *>    _opList;
  merylOpCompute                  **_thList[64] = { nullptr };   //  Mirrors opList

  std::vector<uint32>               _opTree;
};

#endif
