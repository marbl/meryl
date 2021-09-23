
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

#ifndef MERYLOPCOMPUTE_H
#define MERYLOPCOMPUTE_H

#ifndef MERYLINCLUDE
#error "Do not use merylOpCompute.H, use meryl.H instead."
#endif





//  This class is used for performing the compute.
//
//  It is called in merylCommandBuilder::spawnThreads().  The
//  compute is done using these objects.
//
class merylOpCompute {
public:
  merylOpCompute(merylOpTemplate *op, uint32 slice, uint32 nInputs);
  ~merylOpCompute();

  //  Methods used to run this operation.
public:
  bool    nextMer(void);

private:
  void    findOutputKmer(void);
  void    findOutputValue(void);
  void    findOutputLabel(void);

public:
  kmer    theFMer(void)            { return(_kmer);   };

public:
  void    addInputFromOp  (merylOpCompute *operation);
  void    addInputFromDB  (char const *dbName, uint32 slice);

  void    addOutput    (merylOpTemplate *ot, uint32 slice);
  void    addPrinter   (merylOpTemplate *ot, uint32 slice);
  void    addStatistics(merylOpTemplate *ot, uint32 slice);

  //
  //  Pointer back to the template, so we can get other useful bits instead
  //  of copying everything here.
  //

  merylOpTemplate               *_ot   = nullptr;

  //
  //  Inputs.  A list of the inputs to this operation, and the specific slice
  //  we're working on.  The input objects are all new, different from those
  //  in the template, so we can limit each one to a specific database slice.
  //

  std::vector<merylInput *>      _inputs;
  uint32                         _inputSlice = 0;

  //
  //  Outputs.
  //

  merylHistogram                *_stats          = nullptr;    //  Per-thread stats
  merylStreamWriter             *_writerSlice    = nullptr;    //  Per-thread output
  compressedFileWriter          *_printerSlice   = nullptr;    //  Per-thread output
  compressedFileWriter          *_printer        = nullptr;    //  Either _printerSlice or ot->_printer

  //
  //  Data used to compute the action.
  //

  kmer                           _kmer;               //  The kmer we're computing.

  uint32                         _actLen = 0;         //  Number of inputs with that kmer sequence.
  kmer                          *_act    = nullptr;   //  Kmer sequence/value/label from inputs.
  uint32                        *_actIdx = nullptr;   //  Index into _inputs each kmer in _act.
  uint32                        *_actRdx = nullptr;   //  Reverse index mapping database index to _act index (or uint32max).

  //
  //  Filters.
  //

  bool   isKmerFilteredOut(void);

  //  std::vector< std::vector<merylFilter> >  _filter;

  friend class merylCommandBuilder;
};


#endif  //  MERYLOPCOMPUTE_H
