
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
#include <cmath>


void
merylOperation::findMinCount(void) {
  _value = _actCount[0];
  for (uint32 ii=1; ii<_actLen; ii++)
    if (_actCount[ii] < _value)
      _value = _actCount[ii];
}



void
merylOperation::findMaxCount(void) {
  _value = _actCount[0];
  for (uint32 ii=1; ii<_actLen; ii++)
    if (_value < _actCount[ii])
      _value = _actCount[ii];
}



void
merylOperation::findSumCount(void) {
  _value = 0;
  for (uint32 ii=0; ii<_actLen; ii++)
    _value += _actCount[ii];
}


void
merylOperation::subtractCount(void) {
  _value = _actCount[0];
  for (uint32 ii=1; ii<_actLen; ii++) {
    if ( _value > _actCount[ii] )
      _value -= _actCount[ii];
    else {
      _value = 0;
      return;
    }
  }
}


void
merylOperation::initializeThreshold(void) {

  //  If no thresholds to set, nothing to do.

  if ((_fracDist == DBL_MAX) &&
      (_wordFreq == DBL_MAX))
    return;

  //  The problem with using more than one database is that the number of
  //  distinct kmers is not known.

  if (_inputs.size() != 1) {
    fprintf(stderr, "ERROR: operation most-frequent can work with only one meryl database.\n");
    exit(1);
  }

  //  Later, we could allow streaming operations, and construction of statistics on
  //  the fly.

  bool    allDatabase = true;

  for (uint32 ii=0; ii<_inputs.size(); ii++) {
    if (_inputs[ii]->isFromDatabase() == false) {
      fprintf(stderr, "ERROR: input '%s' to operation most-frequent is not a meryl database.\n",
              _inputs[ii]->_name);
      allDatabase = false;
    }
  }

  if (allDatabase == false)
    exit(1);

  for (uint32 ii=0; ii<_inputs.size(); ii++)
    _inputs[ii]->_stream->loadStatistics();

  merylHistogram  *stats = _inputs[0]->_stream->stats();

  if (_fracDist < DBL_MAX) {
    uint64  nKmers       = 0;
    uint64  nKmersTarget = _fracDist * stats->numDistinct();

    for (uint64 ii=0; ii<stats->histogramLength(); ii++) {
      nKmers += stats->histogramOccurrences(ii);

      if (nKmers >= nKmersTarget) {
        _threshold = stats->histogramValue(ii);
        break;
      }
    }
  }

  if (_wordFreq < DBL_MAX) {
    _threshold = _wordFreq * stats->numTotal();
  }

  //  Cleanup.

  for (uint32 ii=0; ii<_inputs.size(); ii++)
    _inputs[ii]->_stream->dropStatistics();
}



bool
merylOperation::initialize(bool isMasterTree) {
  bool  proceed = true;

  //  Initialize all the inputs this operation might have.

  for (uint32 ii=0; ii<_inputs.size(); ii++)
    _inputs[ii]->initialize();

  //  Decide if we're processing a multi-set.

  for (uint32 ii=0; ii<_inputs.size(); ii++)
    _isMultiSet |= _inputs[ii]->isMultiSet();

  //  If this is a counting operation. STOP NOW!  The output cannot be
  //  initialized until after we figure out the correct prefix size (not 0 as
  //  below).

  if (isCounting() == true)
    return(true);

  //  Initialize outputs.
  //    If the master tree, initialize the merylFileWriter.
  //    If a thread tree, create a merylStreamWriter.

  if (_outputO)
    _outputO->initialize(0, isMultiSet());

  if (_outputP)
    _writer = _outputP->getStreamWriter(_fileNumber);

  //  The threshold operations need to decide on a threshold based on the histogram.

  if (isMasterTree)
    initializeThreshold();

  //  If configuring, or if the operation is pass-through with no output,
  //  don't stream the mers.  This only matters for the root node; the return
  //  value for all other nodes is ignored (those are called above, when
  //  initializing the inputs to this node).

  if (_onlyConfig == true)
    proceed = false;

  if ((_operation == opPassThrough) &&   //  'meryl print DATABASE' uses opPassThrough.
      (_printer   == NULL))              //  but has _printer set.
    proceed = false;                     //  Counting operations do not set _printer.

  return(proceed);
}



//  Perform the counting operation, then close the output.
//
void
merylOperation::doCounting(void) {
  char    name[FILENAME_MAX + 1] = { 0 };
  bool    doSimple   = false;
  bool    doThreaded = true;
  uint32  wPrefix    = 0;
  uint64  nPrefix    = 0;
  uint32  wData      = 0;
  kmdata  wDataMask  = 0;

  configureCounting(_maxMemory,
                    doSimple,
                    wPrefix,
                    nPrefix,
                    wData,
                    wDataMask);

  setNumThreads(_maxThreads);

  if (doSimple) {
    fprintf(stderr, "Start counting with SIMPLE method.\n");
    countSimple();
  }

  else if (doThreaded) {
    fprintf(stderr, "Start counting with THREADED method.\n");
    countThreads(wPrefix, nPrefix, wData, wDataMask);
  }

  else {
    fprintf(stderr, "Start counting with SEQUENTIAL method.\n");
    count(wPrefix, nPrefix, wData, wDataMask);
  }

  //  Fiddle with the operation.
  //   - remove the output; it's already been written.
  //   - remove all the inputs
  //   - convert the operation to a simple 'pass through'
  //   - add the counted output as an input

  if (_outputO)
    strncpy(name, _outputO->filename(), FILENAME_MAX + 1);   //  know which input to open later.

  delete _outputO;
  _outputO = NULL;

  clearInputs();

  _operation = opPassThrough;

  if ((name[0]     == 0) ||         //  If there is no name, we've
      (_onlyConfig == true))        //  been asked to only configure.
    return;

  addInputFromDB(name);
}



//  Build a list of the inputs that have the smallest kmer, saving their
//  counts in _actCount, and the input that it is from in _actIndex.
void
merylOperation::nextMer_findSmallestNormal(void) {

  _actLen = 0;                                       //  Reset to nothing on the list.

  for (uint32 ii=0; ii<_inputs.size(); ii++) {
    if (_inputs[ii]->_valid == false)
      continue;

    if ((_actLen == 0) ||                            //  If we have no active kmer, or the input kmer is
        (_inputs[ii]->_kmer < _kmer)) {              //  smaller than the one we have, reset the list.
      _kmer        = _inputs[ii]->_kmer;
      _actCount[0] = _inputs[ii]->_value;
      _actIndex[0] = ii;
      _actLen      = 1;

      if (_verbosity >= sayDetails) {
        char  kmerString[256];
        fprintf(stderr, "merylOp::nextMer()-- Active kmer %s from input %s. reset\n", _kmer.toString(kmerString), _inputs[ii]->_name);
      }
    }

    else if (_inputs[ii]->_kmer == _kmer) {          //  Otherwise, if the input kmer is the one we
      _actCount[_actLen] = _inputs[ii]->_value;      //  have, save the count and input to the lists.
      _actIndex[_actLen] = ii;
      _actLen++;

      if (_verbosity >= sayDetails) {
        char  kmerString[256];
        fprintf(stderr, "merylOp::nextMer()-- Active kmer %s from input %s\n", _kmer.toString(kmerString), _inputs[ii]->_name);
      }
    }

    else {                                           //  Otherwise, the input kmer comes after the
    }                                                //  one we're examining, ignore it.
  }
}



//  For multi-set operation, the list we build must have exactly one item in it.

//  THIS IS WRONG, it needs to build a list with all the stuff with the same kmer AND value.
//  It's up to the operation to decide what to do.
//  So we probably need a new operation 'merge' or something.

//  The action depends on the operations.
//    Intersect -- treat non-multiset as wildcard; add to list with same value
//    Union     -- treat non-multiset as multiset

void
merylOperation::nextMer_findSmallestMultiSet(void) {

  _actLen = 0;

  //
  //  A first pass to handle the inputs that are multi-sets.
  //    If a union operation, save only the first lowest kmer/value pair.
  //    Otherwise, add it to the list if it is the same or lower.
  //

  for (uint32 ii=0; ii<_inputs.size(); ii++) {
    if ((_inputs[ii]->_valid == false) ||
        (_inputs[ii]->isMultiSet() == false))
      continue;

    if ((_operation == opUnion) ||
        (_operation == opUnionMin) ||
        (_operation == opUnionMax) ||
        (_operation == opUnionSum)) {
      if (((_actLen == 0)) ||
          ((_inputs[ii]->_kmer  < _kmer)) ||
          ((_inputs[ii]->_kmer == _kmer) && (_inputs[ii]->_value < _actCount[0]))) {
        _kmer        = _inputs[ii]->_kmer;
        _actCount[0] = _inputs[ii]->_value;
        _actIndex[0] = ii;
        _actLen      = 1;

        if (_verbosity >= sayDetails) {
          char  kmerString[256];
          fprintf(stderr, "merylOp::nextMer()-- Active kmer %s from input %s. reset\n", _kmer.toString(kmerString), _inputs[ii]->_name);
        }
      }
    }

    else if (((_actLen == 0)) ||
             ((_inputs[ii]->_kmer  < _kmer)) ||
             ((_inputs[ii]->_kmer == _kmer) && (_inputs[ii]->_value < _actCount[0]))) {
      _kmer              = _inputs[ii]->_kmer;
      _actCount[_actLen] = _inputs[ii]->_value;
      _actIndex[_actLen] = ii;
      _actLen++;
    }
  }

  //
  //  A second pass to handle the inputs that are not multi-sets.
  //  - If a union operation, reset the list if anything is strictly smaller.
  //  - Otherwise, treat as a wildcard and add it to the list if the kmer is equal, resetting the value.
  //
  //    In both cases, add if the list is empty.
  //    - For union:      we want to add these kmers/values to the output.
  //    - For intersect:  the intersect logic will reject it, and we'll then load a new kmer for this input.
  //    - For difference: if this is not the first input, it works as expected, no output
  //                      if this is     the first input, it works as expected, the kmer is output (with original value).
  //

  for (uint32 ii=0; ii<_inputs.size(); ii++) {
    if ((_inputs[ii]->_valid == false) ||
        (_inputs[ii]->isMultiSet() == false))
      continue;

    if ((_operation == opUnion) ||
        (_operation == opUnionMin) ||
        (_operation == opUnionMax) ||
        (_operation == opUnionSum)) {
      if (((_actLen == 0)) ||
          ((_inputs[ii]->_kmer  < _kmer)) ||
          ((_inputs[ii]->_kmer == _kmer) && (_inputs[ii]->_value < _actCount[0]))) {
        _kmer        = _inputs[ii]->_kmer;
        _actCount[0] = _inputs[ii]->_value;
        _actIndex[0] = ii;
        _actLen      = 1;

        if (_verbosity >= sayDetails) {
          char  kmerString[256];
          fprintf(stderr, "merylOp::nextMer()-- Active kmer %s from input %s. reset\n", _kmer.toString(kmerString), _inputs[ii]->_name);
        }
      }
    }

    else if (((_actLen == 0)) ||
             ((_inputs[ii]->_kmer  < _kmer)) ||
             ((_inputs[ii]->_kmer == _kmer) && (_inputs[ii]->_value < _actCount[0]))) {
      _kmer              = _inputs[ii]->_kmer;
      _actCount[_actLen] = _inputs[ii]->_value;
      _actIndex[_actLen] = ii;
      _actLen++;
    }
  }
}



//  If no active kmers, we're done.  Several bits of housekeeping need to be done:
//   - Histogram operations need to finish up and report the histogram now.
//     Alternatively, it could be done in the destructor.
//   - Any outputs need to call finishIteration() to rename and/or merge
//     their intermediate outputs.
bool
merylOperation::nextMer_finish(void) {

  if (_verbosity >= sayDetails) {
    fprintf(stderr, "merylOp::nextMer()-- No inputs found, all done here.\n");
    fprintf(stderr, "\n");
  }

  _valid = false;

  if (_operation == opHistogram)
    reportHistogram();

  if (_operation == opStatistics)
    reportStatistics();

  delete _writer;
  _writer = NULL;

  return(false);
}



bool
merylOperation::nextMer(void) {

 nextMerAgain:

  //  Get some logging out of the way.

  if (_verbosity >= sayDetails) {
    fprintf(stderr, "\n");
    fprintf(stderr, "merylOp::nextMer()-- STARTING for operation %s\n",
            toString(_operation));

    if (_verbosity >= sayEverything) {
      char  kmerString[256];

      for (uint32 ii=0; ii<_inputs.size(); ii++)
        fprintf(stderr, "merylOp::nextMer()--   CURRENT STATE: input %s kmer %s count %u %s\n",
                _inputs[ii]->_name,
                _inputs[ii]->_kmer.toString(kmerString),
                _inputs[ii]->_value,
                _inputs[ii]->_valid ? "valid" : "INVALID");
    }
  }

  //  Grab the next mer for every input that was active in the last iteration.
  //  (on the first call, all inputs were 'active' last time)

  for (uint32 ii=0; ii<_actLen; ii++) {
    if (_verbosity >= sayDetails)
      fprintf(stderr, "merylOp::nextMer()-- CALL NEXTMER on input actIndex " F_U32 "\n", _actIndex[ii]);
    _inputs[_actIndex[ii]]->nextMer();
  }

  //  Build a list of the inputs that have the smallest kmer, saving their
  //  counts in _actCount, and the input that it is from in _actIndex.

  if (isMultiSet() == false)
    nextMer_findSmallestNormal();
  else
    nextMer_findSmallestMultiSet();

  //  If no active kmers, we're done.

  if (_actLen == 0)
    return(nextMer_finish());

  //  Otherwise, active kmers!  Figure out what the count should be.

  if (_verbosity >= sayDetails) {
    char  kmerString[256];
    fprintf(stderr, "merylOp::nextMer()-- op %s activeLen " F_U32 " kmer %s\n", toString(_operation), _actLen, _kmer.toString(kmerString));
  }

  //  If math-subtract gets implemented, use negative-zero to mean "don't output" and positive-zero
  //  to mean zero.  For now, count=0 means don't output.

  //  Set the count to zero, meaning "don't output the kmer".  Intersect depends on this,
  //  skipping most of it's work if all files don't have the kmer.
  _value = 0;

  switch (_operation) {
    case opCount:
    case opCountForward:
    case opCountReverse:
      fprintf(stderr, "ERROR: got %s, but shouldn't have.\n", toString(_operation));
      assert(0);
      break;

    case opPassThrough:                     //  Result of counting kmers.  Guaranteed to have
      _value = _actCount[0];                //  exactly one input file.  Also the operation that
      break;                                //  'print' of a database has.

    case opLessThan:
      _value = (_actCount[0]  < _threshold) ? _actCount[0] : 0;
      break;

    case opGreaterThan:
      _value = (_actCount[0]  > _threshold) ? _actCount[0] : 0;
      break;

    case opAtLeast:
      _value = (_actCount[0] >= _threshold) ? _actCount[0] : 0;
      break;

    case opAtMost:
      _value = (_actCount[0] <= _threshold) ? _actCount[0] : 0;
      break;

    case opEqualTo:
      _value = (_actCount[0] == _threshold) ? _actCount[0] : 0;
      break;

    case opNotEqualTo:
      _value = (_actCount[0] != _threshold) ? _actCount[0] : 0;
      break;

    case opIncrease:
      if (UINT64_MAX - _actCount[0] < _mathConstant)
        _value = ~((kmvalu)0);  //  OVERFLOW!
      else
        _value = _actCount[0] + _mathConstant;
      break;

    case opDecrease:
      if (_actCount[0] < _mathConstant)
        _value = 0;             //  UNDERFLOW!
      else
        _value = _actCount[0] - _mathConstant;
      break;

    case opMultiply:
      if (UINT64_MAX / _actCount[0] < _mathConstant)
        _value = ~((kmvalu)0);   //  OVERFLOW!
      else
        _value = _actCount[0] * _mathConstant;
      break;

    case opDivide:
      if (_mathConstant == 0)
        _value = 0;             //  DIVIDE BY ZERO!
      else
        _value = _actCount[0] / _mathConstant;
      break;
    case opDivideRound:
      if (_mathConstant == 0)
        _value = 0;             //  DIVIDE BY ZERO!
      else {
        if (_actCount[0] < _mathConstant)
          _value = 1;
        else
          _value = round (_actCount[0] / (double) _mathConstant);
      }
      break;

    case opModulo:
      if (_mathConstant == 0)
        _value = 0;             //  DIVIDE BY ZERO!
      else
        _value = _actCount[0] % _mathConstant;
      break;

    case opUnion:                           //  Union
      _value = _actLen;
      break;

    case opUnionMin:                        //  Union, retain smallest count
      findMinCount();
      break;

    case opUnionMax:                        //  Union, retain largest count
      findMaxCount();
      break;

    case opUnionSum:                        //  Union, sum all counts
      findSumCount();
      break;

    case opIntersect:                       //  Intersect
      if (_actLen == _inputs.size())
        _value = _actCount[0];
      break;

    case opIntersectMin:                    //  Intersect, retain smallest count
      if (_actLen == _inputs.size())
        findMinCount();
      break;

    case opIntersectMax:                    //  Intersect, retain largest count
      if (_actLen == _inputs.size())
        findMaxCount();
      break;

    case opIntersectSum:                    //  Intersect, sum all counts
      if (_actLen == _inputs.size())
        findSumCount();
      break;

    case opSubtract:
      if (_actIndex[0] == 0) {
        if (_actLen == 1)
          _value = _actCount[0];
        else if (_actLen > 1)
          subtractCount();
      }
      break;

    case opDifference:
      if ((_actLen == 1) && (_actIndex[0] == 0))
        _value = _actCount[0];
      break;

    case opSymmetricDifference:
      if (_actLen == 1)
        _value = _actCount[0];
      break;

    case opPloidy:
      break;

    case opCompare:
      if       (_actLen == 1) {
        char  str[65];

        fprintf(stdout, "kmer %s only in input %u\n",
                _kmer.toString(str), _actIndex[0]);
      }
      else if ((_actLen == 2) && (_actCount[0] != _actCount[1])) {
        char  str[65];

        fprintf(stdout, "kmer %s has value %u in input 1 != value %u in input 2\n",
                _kmer.toString(str), _actCount[0], _actCount[1]);
      }
      else {
      }

    case opHistogram:
      break;

    case opStatistics:
      break;

    case opNothing:
      break;
  }

  //  If the count is zero, skip this kmer and get another one.

  if (_value == 0)
    goto nextMerAgain;

  //  And if not zero, output it, print it, and return it.

  if (_verbosity >= sayDetails) {
    char  kmerString[256];
    fprintf(stderr, "merylOp::nextMer()-- FINISHED for operation %s with kmer %s count %u%s\n",
            toString(_operation), _kmer.toString(kmerString), _value, ((_outputP != NULL) && (_value != 0)) ? " OUTPUT" : "");
    fprintf(stderr, "\n");
  }

  //  If flagged for output, output!

  if (_writer != nullptr) {
    _writer->addMer(_kmer, _value);
  }

  //  If flagged for printing, print!

  if (_printer != nullptr) {
    kmer  pk = _kmer;

    if (_printACGTorder == true)
      pk.recanonicalizeACGTorder();

#pragma omp critical (printLock)
    {
      fputs(pk.toString(_kmerString), _printer);
      fputc('\t', _printer);
      fputs(toDec(_value), _printer);
      fputc('\n', _printer);
    }
  }

  //  Now just return and let the client query us to get the kmer and value.

  return(true);
}
