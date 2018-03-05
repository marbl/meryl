
/******************************************************************************
 *
 *  This file is part of 'sequence' and/or 'meryl', software programs for
 *  working with DNA sequence files and k-mers contained in them.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2018-FEB-26
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.license' in the root directory of this distribution contains
 *  full conditions and disclaimers.
 */

#ifndef MERYLOP_H
#define MERYLOP_H

#include "meryl.H"


enum merylOp {
  opCount,
  opCountForward,
  opCountReverse,
  opUnion,
  opUnionMin,
  opUnionMax,
  opUnionSum,
  opIntersect,
  opIntersectMin,
  opIntersectMax,
  opIntersectSum,
  opDifference,
  opSymmetricDifference,
  opComplement,
  opPrint,
  opNothing
};



class merylOperation {
public:
  merylOperation(merylOp op=opNothing, uint32 k=0, uint64 numMers=0, uint32 threads=1, uint64 memory=0);
  ~merylOperation();

  void    addInput(merylOperation *operation);
  void    addInput(char *name, kmerCountFileReader *reader);
  void    addInput(char *name, dnaSeqFile *sequence);

  void    addOutput(char *name, kmerCountFileWriter *writer);

  void    setOperation(merylOp op) { _operation = op;    };
  merylOp getOperation(void)       { return(_operation); };

  bool    isCounting(void) {
    return((_operation == opCount)        ||
           (_operation == opCountForward) ||
           (_operation == opCountReverse));
  };

  bool    isNormal(void) {
    return(isCounting() == false);
  };

  kmer   &theFMer(void)            { return(_kmer);   };
  uint64  theCount(void)           { return(_count);  };

  bool    nextMer(void);
  bool    validMer(void)           { return(_valid);  };

  void    count(void);
  void    countSimple(void);

private:
  void    findMinCount(void);
  void    findMaxCount(void);
  void    findSumCount(void);

  vector<merylInput *>           _inputs;

  merylOp                        _operation;

  uint32                         _k;
  uint64                         _numMers;
  uint32                         _threads;
  uint64                         _memory;

  char                           _outputName[FILENAME_MAX];
  kmerCountFileWriter           *_output;

  kmer                           _smallest;

  uint32                         _actLen;
  uint64                        *_actCount;
  uint32                        *_actIndex;

  kmer                           _kmer;
  uint64                         _count;
  bool                           _valid;
};


char const *toString(merylOp op);


#endif  //  MERYLOP_H
