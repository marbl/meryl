
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

#ifndef MERYLOP_H
#define MERYLOP_H

#include "meryl.H"

#include <vector>

//  These would (maybe) be handy, but they change the order of kmers, so we need
//  to sort again.
//
//  opReverse,             //  Reverse the kmer.     Accepts a single input only.
//  opComplement,          //  Complement the kmer.  Accepts a single input only.
//  opReverseComplement,   //  Rev-Compl the kmer.   Accepts a single input only.

enum class merylOp {
  opNothing,

  opCount,
  opCountForward,
  opCountReverse,

  opPresentIn,

  opValue,
  opLabel,

  opHistogram,
  opStatistics,

  opCompare,
};

enum class merylValueAction {
  valueNothing,     //  Do absoluately nothing.

  valueConstant,    //  Set to _valueConstant.
  valueSelected,    //  Set to the value of the selected input file.

  valueCount,       //  Set to the number of files with this kmer.
  valueFirst,       //  Set to the value of the kmer in the first input file.

  valueMin,         //  Set to the minimum of all values.
  valueMax,         //  Set to the maximum of all values.
  valueSum,         //  Set to the sum of all values.
  valueSub,         //  Set to the first value minus all the rest.

  valueIncrease,    //  Set to _value + _valueConstant (or valumax
  valueDecrease,    //  Set to _value - _valueConstant (or zero)
  valueMultiply,    //  Set to _value * _valueConstant (or valumax)
  valueDivideTrunc, //  Set to _value / _valueConstant, integer division; truncate
  valueDivideRound, //  Set to _value / _valueConstant, real    division; round
  valueRemainder,   //  Set to _value % _valueConstant
};

enum class merylValueFilter {
  valueNothing,     //  Do absoluately nothing.

  valueLessThan,    //  Output if the value is less than _valueConstant
  valueGreaterThan, //  (and so on)
  valueAtLeast,     //  
  valueAtMost,      //  
  valueEqualTo,     //  
  valueNotEqualTo,  //  
};

enum class merylLabelAction {
  labelNothing,     //  Do absoluately nothing.

  labelConstant,    //  Set to _labelConstant.
  labelSelected,    //  Set to the label of the selected input file.

  labelFirst,       //  Set to the label of the kmer in the first input file.

  labelMin,         //  Set to the label of the kmer with the minimum value (e.g., valueMin).
  labelMax,         //  Set to the label of the kmer with the maximum value (e.g., valueMax).

  labelAnd,         //  AND all labels together along with the labelConstant
  labelOr,          //  OR  all labels together along with the labelConstant
  labelXor,         //  XOR all labels together along with the labelConstant
  labelDifference,  //  Unset any bits set in kmers 2..N (like opDifference)

  labelLightest,    //  Set label to lowest  bit weight label in active set
  labelHeaviest,    //  Set label to highest bit weight label in active set

  labelInvert,      //  Invert the existing label bits
  labelShiftLeft,   //  Shift the label left  _labelConstant bits
  labelShiftRight,  //  Shift the label right _labelConstant bits
};

enum class merylLabelFilter {
  labelNothing,     //  Do absoluately nothing.
};




enum merylVerbosity {
  sayNothing,
  sayStandard,
  sayConstruction,
  sayDetails,
  sayEverything
};


//  For the simple counting algorithm, how bit of a word to use for initial counts.
typedef uint16 lowBits_t;


class merylOperation {
public:
  merylOperation();
  merylOperation(merylOperation *op,
                 uint32 fileNum,
                 uint32 nInputs,
                 uint32 threads, uint64 memory);
  ~merylOperation();

private:
  void    clearInputs(void);
  void    checkInputs(char const *name);

  uint64  guesstimateNumberOfkmersInInput_dnaSeqFile(dnaSeqFile *sequence);
  uint64  guesstimateNumberOfkmersInInput_sqStore(sqStore *store, uint32 bgnID, uint32 endID);
  uint64  guesstimateNumberOfkmersInInput(void);

  //  Returns true
  void    configureCounting(uint64   memoryAllowed,      //  Input:  Maximum allowed memory in bytes
                            bool    &useSimple_,         //  Output: algorithm to use
                            uint32  &wPrefix_,           //  Output: Number of bits in the prefix (== bucket address)
                            uint64  &nPrefix_,           //  Output: Number of prefixes there are (== number of buckets)
                            uint32  &wData_,             //  Output: Number of bits in kmer data
                            kmdata  &wDataMask_);        //  Output: A mask to return just the data of the mer);

private:
  void    nextMer_findOutputValue(void);
  void    nextMer_filterOutputValue(void);

  void    nextMer_findOutputLabel(void);
  void    nextMer_filterOutputLabel(void);

  void    nextMer_outputKmer(void);
  void    nextMer_printKmer(void);


public:
  void    addInputFromOp  (merylOperation *operation);
  void    addInputFromDB  (char const *dbName);
  void    addInputFromCanu(char const *stName, uint32 segment, uint32 segmentMax);
  void    addInputFromSeq (char const *sqName, bool doCompression);

  void    addOutput(char const *wrName);
  void    addPrinter(char const *prName, bool ACGTorder);



  void    finalize(void);

  void    setOperation(merylOp op) { _operation = op;    };
  merylOp getOperation(void)       { return(_operation); };

  void    setValueAction(merylValueAction act, merylValueFilter fil, kmvalu cst) {
    _valueAction   = act;
    _valueFilter   = fil;
    _valueConstant = cst;
  };

  void    setLabelAction(merylLabelAction act, merylLabelFilter fil, kmlabl cst) {
    _labelAction   = act;
    _labelFilter   = fil;
    _labelConstant = cst;
  };

  void    setConstant(uint64 p)                { _mathConstant = p; };
  void    setThreshold(uint64 p)               { _threshold    = p; };
  void    setFractionDistinct(double p)        { _fracDist     = p; };
  void    setWordFrequency(double p)           { _wordFreq     = p; };

  void    setExpectedNumberOfKmers(uint64 n)   { _expNumKmers  = n; };

  void    setCountSuffix(char const *s)              {
    _countSuffixLength = strlen(s);

    for (uint32 ii=0; ii<_countSuffixLength; ii++) {
      _countSuffixString[ii] = s[ii];

      _countSuffix.addR(s[ii]);
    }
  };

  void    setMemoryLimit(uint64 m)             { _maxMemory    = m; };
  void    setThreadLimit(uint32 t)             { _maxThreads   = t; };

  uint64  getMemoryLimit(void)                 { return(_maxMemory);  };
  uint32  getThreadLimit(void)                 { return(_maxThreads); };

  bool    isCounting(void) {
    return((_operation == merylOp::opCount)        ||
           (_operation == merylOp::opCountForward) ||
           (_operation == merylOp::opCountReverse));
  };

  bool    isNormal(void) {
    return(isCounting() == false);
  };

  bool    needsThreshold(void) {
    return((_valueFilter == merylValueFilter::valueLessThan)     ||
           (_valueFilter == merylValueFilter::valueGreaterThan)  ||
           (_valueFilter == merylValueFilter::valueAtLeast)      ||
           (_valueFilter == merylValueFilter::valueAtMost)       ||
           (_valueFilter == merylValueFilter::valueEqualTo)      ||
           (_valueFilter == merylValueFilter::valueNotEqualTo));
  };

  bool    needsConstant(void) {
    return((_valueAction == merylValueAction::valueIncrease)     ||
           (_valueAction == merylValueAction::valueDecrease)     ||
           (_valueAction == merylValueAction::valueMultiply)     ||
           (_valueAction == merylValueAction::valueDivideTrunc)  ||
           (_valueAction == merylValueAction::valueDivideRound)  ||
           (_valueAction == merylValueAction::valueRemainder));
  };

  //bool    needsRealParameter(void) {
  //  return((_operation == opMostFrequent));
  //};

  bool    isMultiSet(void)         { return(_isMultiSet); };

  kmer   &theFMer(void)            { return(_kmer);       };
  kmvalu  theValue(void)           { return(_value);      };
  kmlabl  theLabel(void)           { return(_label);      };

public:
  void    doCounting(void);

public:
  void    initializeThreshold(void);
  bool    initialize(bool isMasterTree=false);

private:
  void    nextMer_findSmallestNormal(void);
  void    nextMer_findSmallestMultiSet(void);
  bool    nextMer_finish(void);

public:
  bool    nextMer(void);

  bool    validMer(void)           { return(_valid);  };

  void    countSimple(void);
  void    countThreads(uint32  wPrefix,
                       uint64  nPrefix,
                       uint32  wData,
                       kmdata  wDataMask);
  void    count(uint32  wPrefix,
                uint64  nPrefix,
                uint32  wData,
                kmdata  wDataMask);

  void    reportHistogram(void);
  void    reportStatistics(void);

  static
  void    onlyConfigure(void)               { _onlyConfig   = true;       };
  static
  void    showProgress(void)                { _showProgress = true;       };
  static
  void    increaseVerbosity(void) {
    if      (_verbosity == sayNothing)        _verbosity    = sayStandard;
    else if (_verbosity == sayStandard)       _verbosity    = sayConstruction;
    else if (_verbosity == sayConstruction)   _verbosity    = sayDetails;
    else if (_verbosity == sayDetails)        _verbosity    = sayEverything;
    else                                      _verbosity    = sayEverything;
  };
  static
  void    beQuiet(void)                     { _verbosity    = sayNothing; };


private:
  void    findMinCount(void);
  void    findMaxCount(void);
  void    findSumCount(void);
  void    subtractCount(void);

  //  Describing the inputs to this operation.

  std::vector<merylInput *>      _inputs;
  bool                           _isMultiSet = false;

  //  Describing the operation.
  //
  //    presentIn   - the number of files the kmer needs to be present in.
  //                  presentIn[nFound] == false --> do not output the kmer
  //                  presentIn[nFound] == true  --> do     output the kmer
  //
  //    valueAction - how to combine the values into a single output value.
  //
  //    labelAction - how to combine the labels into a single output label.

  uint32                         _presentInLen  = 0;
  uint32                         _presentInMax  = 0;
  uint8                         *_presentIn     = nullptr;


  merylValueAction               _valueAction   = merylValueAction::valueNothing;
  merylValueFilter               _valueFilter   = merylValueFilter::valueNothing;

  kmvalu                         _valueConstant = 0;
  kmvalu                         _valueMin      = 0;            //  Acceptable ranges of values
  kmvalu                         _valueMax      = kmvalumax;    //  to output/print.


  merylLabelAction               _labelAction   = merylLabelAction::labelNothing;
  merylLabelFilter               _labelFilter   = merylLabelFilter::labelNothing;
  kmlabl                         _labelConstant = 0;

  //  Where is the kmer present?

  uint32                         _actLen   = 0;
  kmvalu                        *_actValue = nullptr;   //  Value      of the ith active kmer
  kmlabl                        *_actLabel = nullptr;   //  Label      of the ith active kmer
  uint32                        *_actIndex = nullptr;   //  File index of the ith active kmer

  //  If _writeToOutput is not set after various operations, we've decided
  //  that this kmer should not be output, and we'll skip it.

  bool                           _writeToOutput;

  //  If the kmer is found in input file _selectedFile, these hold the value and
  //  label of that specific kmer.  If no file is selected by the user,
  //  _selectedFile will by uint32max, _selected will be false, and _selValue and
  //  _selLabel will never be used.

  uint32                         _selectedFile = uint32max;
  bool                           _selected     = false;
  kmvalu                         _selValue     = 0;
  kmlabl                         _selLabel     = 0;

  //
  //
  
  merylOp                        _operation   = merylOp::opNothing;

  uint64                         _mathConstant = 0;
  uint64                         _threshold    = uint64max;
  double                         _fracDist     = DBL_MAX;
  double                         _wordFreq     = DBL_MAX;

  char                           _countSuffixString[65] = {0};
  uint32                         _countSuffixLength     =  0;
  kmer                           _countSuffix;

  uint64                         _expNumKmers = 0;

  uint32                         _maxThreads = 0;
  uint64                         _maxMemory  = 0;

  merylHistogram                *_stats = nullptr;

  merylFileWriter               *_outputO = nullptr;    //  Main output object (this pointer owns the object)
  merylFileWriter               *_outputP = nullptr;    //  Main output object
  merylStreamWriter             *_writer  = nullptr;    //  Per-thread output

  FILE                          *_printer        = nullptr;
  char                          *_printerName    = nullptr;
  bool                           _printACGTorder = false;

  char                           _kmerString[65] = {0};

  uint32                         _fileNumber = uint32max;


  kmer                           _kmer;
  kmvalu                         _value = 0;
  kmlabl                         _label = 0;

  bool                           _valid = true;

  static bool                    _onlyConfig;
  static bool                    _showProgress;
  static merylVerbosity          _verbosity;

  friend class merylCommandBuilder;
};


char const *toString(merylOp op);
char const *toString(merylValueAction op);
char const *toString(merylValueFilter op);
char const *toString(merylLabelAction op);
char const *toString(merylLabelFilter op);


#endif  //  MERYLOP_H
