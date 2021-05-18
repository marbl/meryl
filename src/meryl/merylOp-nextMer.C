
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

  if ((_operation == merylOp::opPresentIn) &&   //  'meryl print DATABASE' uses opPresentIn.
      (_printer   == NULL))                     //  but has _printer set.
    proceed = false;                            //  Counting operations do not set _printer.

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

  omp_set_num_threads(_maxThreads);

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

  _operation = merylOp::opPresentIn;

  if ((name        == NULL) ||      //  If there is no name, we've
      (name[0]     == 0) ||         //  been asked to only configure.
      (_onlyConfig == true))
    return;

  addInputFromDB(name);
}



//  Build a list of the inputs that have the smallest kmer, saving their
//  counts in _actValue, and the input that it is from in _actIndex.
void
merylOperation::nextMer_findSmallestNormal(void) {

  //  Scan all the inputs.  Remember the smallest kmer we've seen, and
  //  create a list of all the inputs that kmer is in.
  //
  //   - If we have no active kmer, or the input kmer is smaller than the one
  //     we have, reset the list.
  //   - Otherwise, if the input kmer is the one we have, save the count and
  //     input to the lists.
  //
  //  Check if the kmer is present in the requested (@1) input file.

  _actLen        = 0;
  _writeToOutput = false;
  _selected      = false;
  _selValue      = 0;
  _selLabel      = 0;

  for (uint32 ii=0; ii<_inputs.size(); ii++) {

    //  Skip the file if there are no more kmers in it.

    if (_inputs[ii]->_valid == false) {
      continue;
    }

    //  If the kmer in the file is higher than the kmer we've picked,
    //  skip it.

    if ((_actLen > 0) && (_inputs[ii]->_kmer > _kmer)) {
      continue;
    }

    //  If the kmer in this file is lower than the kmer we've picked,
    //  forget all the kmers we've saved.

    if ((_actLen > 0) && (_inputs[ii]->_kmer < _kmer)) {
      _actLen        = 0;

      _writeToOutput = false;
      _selected      = false;
      _selValue      = 0;
      _selLabel      = 0;
    }

    //  If we have no kmer picked, or the kmer in this file is the same as
    //  the kmer we've picked, add the kmer.

    if (_actLen == 0)
      _kmer = _inputs[ii]->_kmer;

    if (_inputs[ii]->_kmer == _kmer) {
      _actValue[_actLen] = _inputs[ii]->_value;
      _actLabel[_actLen] = _inputs[ii]->_label;
      _actIndex[_actLen] = ii;
      _actLen++;
    }

    //  If this is the file the user selected, or if there is no selected
    //  file, declare we can output the kmer and save the value and label.

    if (_selectedFile == uint32max)
      _writeToOutput = true;

    if (_selectedFile == ii) {
      _writeToOutput = true;
      _selected      = true;
      _selValue      = _inputs[ii]->_value;
      _selLabel      = _inputs[ii]->_label;
    }
  }

  //  Now that we've found the smallest kmer and counted the number of inputs
  //  it is found in, remove the _writeToOutput flag if the kmer isn't present
  //  in the desired number of files.

  assert(_actLen < _presentInLen);

  if (_presentIn[_actLen] == false)
    _writeToOutput = false;
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

  if (_operation == merylOp::opHistogram)
    reportHistogram();

  if (_operation == merylOp::opStatistics)
    reportStatistics();

  delete _writer;
  _writer = NULL;

  return(false);
}



//  Write the kmer, value and label to a meryl database output.
//  Nothing magic here.
void
merylOperation::nextMer_outputKmer(void) {

  if (_writer == nullptr)
    return;

  _writer->addMer(_kmer, _value, _label);
}



//  Write the kmer, value and label to an ASCII output file.
//   - If requested, change the kmer from ACTG order to ACGT
//     order, this is needed so that canonical kmers are
//     correct.
//   - The mer is easy, just call its toString() method, then
//     advance to the end of line.
//   - The count is easy too.  This one will return the end
//     of line.
//   - The label is not so easy.  We should print as either
//     hex or binary, and need print only _labelWidth bits.
//   - Append a newline, then write to the output file.
//   - fputs() is NOT thread-safe, and MUST be guarded in
//     a critical section.  The other stuff done in this
//     function is all thread-save (or local).
void
merylOperation::nextMer_printKmer(void) {
  char   outstr[256] = {0};
  char  *outptr = outstr;
  kmer  pk      = _kmer;

  if (_printer == nullptr)
    return;

  if (_printACGTorder == true)
    pk.recanonicalizeACGTorder();

  pk.toString(outptr);
  while (*outptr)
    outptr++;

  *outptr++ = '\t';
  outptr = toDec(_value, outptr);

  if (kmer::labelSize() > 0) {
    *outptr++ = '\t';
    outptr = toBin(_label, outptr, kmer::labelSize());
  }

  *outptr++ = '\n';
  *outptr++ = 0;

#pragma omp critical (printLock)
  fputs(outstr, _printer);
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

  //  Grab the next mer for every input that was active in the last
  //  iteration.  (on the first call, all inputs were 'active' last time)

  for (uint32 ii=0; ii<_actLen; ii++) {
    if (_verbosity >= sayDetails)
      fprintf(stderr, "merylOp::nextMer()-- CALL NEXTMER on input actIndex " F_U32 "\n", _actIndex[ii]);
    _inputs[_actIndex[ii]]->nextMer();
  }

  //  Find the smallest kmer in any input, and remember the values and labels
  //  of the kmer in each input file.

  if (isMultiSet() == false)
    nextMer_findSmallestNormal();
  else {
    assert(0);
    nextMer_findSmallestMultiSet();
  }

  //  If no active kmers, we're done.  There are no kmers left in any input.

  if (_actLen == 0)
    return(nextMer_finish());

  //  The findSmallest() functions will decide, based on present-in, if we
  //  should skip this kmer.  If we're told to skip it, do it right now;
  //  there is no need to do any processing on values or labels.

  if (_verbosity >= sayDetails) {
    char  kmerString[256];
    fprintf(stderr, "merylOp::nextMer()-- op %s activeLen " F_U32 " kmer %s\n", toString(_operation), _actLen, _kmer.toString(kmerString));
  }

  if (_writeToOutput == false)
    goto nextMerAgain;

  //  Figure out the value of the output kmer.

  nextMer_findOutputValue();
  nextMer_filterOutputValue();

  if (_writeToOutput == false)
    goto nextMerAgain;

  //  Figure out the label of the output kmer.

  nextMer_findOutputLabel();
  nextMer_filterOutputLabel();

  if (_writeToOutput == false)
    goto nextMerAgain;

  //  Output the fully processed kmer to whatever outputs exist.

  if (_verbosity >= sayDetails) {
    char  kmerString[256];
    fprintf(stderr, "merylOp::nextMer()-- FINISHED for operation %s with kmer %s count %u label %016lx%s\n",
            toString(_operation), _kmer.toString(kmerString), _value, _label, ((_outputP != NULL) && (_value != 0)) ? " OUTPUT" : "");
    fprintf(stderr, "\n");
  }

  nextMer_outputKmer();
  nextMer_printKmer();

  //  We have now loaded the nextMer; return and let our clients query us.

  return(true);
}



void
merylOperation::nextMer_findOutputValue(void) {

  switch (_valueAction) {
    case merylValueAction::valueNothing:
      break;

    case merylValueAction::valueConstant:
      _value = _valueConstant;
      break;

    case merylValueAction::valueSelected:
      _value = _selValue;
      break;

    case merylValueAction::valueMin:
      _value = _valueConstant;
      for (uint32 ii=0; ii<_actLen; ii++)
        _value = std::min(_value, _actValue[ii]);
      break;

    case merylValueAction::valueMax:
      _value = _valueConstant;
      for (uint32 ii=0; ii<_actLen; ii++)
        _value = std::max(_value, _actValue[ii]);
      break;

    case merylValueAction::valueSum:
      _value = _valueConstant;
      for (uint32 ii=0; ii<_actLen; ii++)
        _value += _actValue[ii];
      break;

    case merylValueAction::valueSub:
      _value = _selValue + _selValue;
      for (uint32 ii=0; ii<_actLen; ii++) {
        if (_value > _actValue[ii])
          _value -= _actValue[ii];
        else
          _value = 0;
      }
      break;

    case merylValueAction::valueIncrease:
      if (kmvalumax - _value < _valueConstant)
        _value = kmvalumax;
      else
        _value = _value + _valueConstant;
      break;

    case merylValueAction::valueDecrease:
      if (_value < _valueConstant)
        _value = 0;
      else
        _value = _value - _valueConstant;
      break;

    case merylValueAction::valueMultiply:
      if (kmvalumax / _value < _valueConstant)
        _value = kmvalumax;
      else
        _value = _value * _valueConstant;
      break;

    case merylValueAction::valueDivideTrunc:
      if (_valueConstant == 0)
        _value = 0;
      else
        _value = _value / _valueConstant;
      break;

    case merylValueAction::valueDivideRound:
      if (_valueConstant == 0)
        _value = 0;
      else if (_value < _valueConstant)
        _value = 1;
      else
        _value = (kmvalu)round(_value / (double)_valueConstant);
      break;

    case merylValueAction::valueRemainder:
      if (_valueConstant == 0)
        _value = 0;
      else
        _value = _value % _valueConstant;
      break;
  }
}



//  Decide, based on the value, if we should output a kmer.
//  Kmers with value zero are never output.
//
void
merylOperation::nextMer_filterOutputValue(void) {

  assert(_writeToOutput == true);

  if (_value == 0) {
    _writeToOutput = false;
    return;
  }

  switch (_valueFilter) {
    case merylValueFilter::valueNothing:
      break;

    case merylValueFilter::valueLessThan:
      if (not(_value < _valueConstant))
        _writeToOutput = true;
      break;

    case merylValueFilter::valueGreaterThan:
      if (not(_value > _valueConstant))
        _writeToOutput = true;
      break;

    case merylValueFilter::valueAtLeast:
      if (not(_value >= _valueConstant))
        _writeToOutput = true;
      break;

    case merylValueFilter::valueAtMost:
      if (not(_value <= _valueConstant))
        _writeToOutput = true;
      break;

    case merylValueFilter::valueEqualTo:
      if (not(_value == _valueConstant))
        _writeToOutput = true;
      break;

    case merylValueFilter::valueNotEqualTo:
      if (not(_value != _valueConstant))
        _writeToOutput = true;
      break;
  }
}



void
merylOperation::nextMer_findOutputLabel(void) {

  switch (_labelAction) {
    case merylLabelAction::labelNothing:
      //assert(_actLen == 1);
      _label = _actLabel[0];
      break;

    case merylLabelAction::labelConstant:
      _label = _labelConstant;
      break;

    case merylLabelAction::labelSelected:
      _label = _selLabel;
      break;

    case merylLabelAction::labelAnd:
      _label = _labelConstant;
      for (uint32 ll=0; ll<_actLen; ll++)
        _label &= _actLabel[ll];
      break;

    case merylLabelAction::labelOr:
      _label = _labelConstant;
      for (uint32 ll=0; ll<_actLen; ll++)
        _label |= _actLabel[ll];
      break;

    case merylLabelAction::labelXor:
      _label = _labelConstant;
      for (uint32 ll=0; ll<_actLen; ll++)
        _label ^= _actLabel[ll];
      break;

    case merylLabelAction::labelDifference:
      _label = _actLabel[0] & ~_labelConstant;
      for (uint32 ll=0; ll<_actLen; ll++)
        _label &= ~_actLabel[ll];
      break;

    case merylLabelAction::labelFirst:
      assert(_actLen >= 1);
      _label = _actLabel[0];
      break;

#if 0
    case merylLabelAction::labelLast:
      assert(_actLen >= 1);
      _label = _actLabel[_actLen-1];
      break;
#endif

    case merylLabelAction::labelLightest:
      _label = _labelConstant;
      for (uint32 ll=0; ll<_actLen; ll++)
        if (countNumberOfSetBits64(_actLabel[ll]) < countNumberOfSetBits64(_label))
          _label = _actLabel[ll];
      break;

    case merylLabelAction::labelHeaviest:
      _label = _labelConstant;
      for (uint32 ll=0; ll<_actLen; ll++)
        if (countNumberOfSetBits64(_actLabel[ll]) > countNumberOfSetBits64(_label))
          _label = _actLabel[ll];
      break;

    case merylLabelAction::labelInvert:
      assert(_actLen == 1);
      _label = ~_actLabel[0];
      break;

    case merylLabelAction::labelShiftLeft:
      assert(_actLen == 1);
      _label = _actLabel[0] >> _labelConstant;
      break;

    case merylLabelAction::labelShiftRight:
      assert(_actLen == 1);
      _label = _actLabel[0] << _labelConstant;
      break;
  }
}



void
merylOperation::nextMer_filterOutputLabel(void) {
}

