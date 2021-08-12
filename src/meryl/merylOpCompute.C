
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



merylOpCompute::merylOpCompute(merylOpTemplate *ot, uint32 dbSlice, uint32 nInputs) {
  _ot            = ot;

  _inputSlice    = dbSlice;

  _act           = new kmer   [nInputs];   //  We know how many inputs there are,
  _actIdx        = new uint32 [nInputs];   //  and can allocate these now.  The actual
  _actRdx        = new uint32 [nInputs];   //  inputs are added later.
}



merylOpCompute::~merylOpCompute() {

  for (uint32 ii=0; ii<_inputs.size(); ii++)   //  Destroy the inputs, which will recursively destroy
    delete _inputs[ii];                        //  any merylOpCompute objects we have as inputs.

  _inputs.clear();

  delete    _stats;                            //  Close the outputs.
  delete    _writerSlice;
  delete    _printerSlice;

  delete [] _act;                              //  Release any scratch memory we used.
  delete [] _actIdx;
  delete [] _actRdx;
}



////////////////////////////////////////
//
//  INPUTS.
//
void
merylOpCompute::addInputFromOp(merylOpCompute *ocin) {
  _inputs.push_back(new merylInput(ocin));
}



void
merylOpCompute::addInputFromDB(char const *dbName, uint32 slice) {
  _inputs.push_back(new merylInput(new merylFileReader(dbName, slice)));
}



////////////////////////////////////////
//
//  OUTPUT to database.
//
void
merylOpCompute::addOutput(merylOpTemplate *ot, uint32 slice) {
  if (ot->_writer != nullptr)
    _writerSlice = ot->_writer->getStreamWriter(slice);
}

void
merylOpCompute::outputKmer(void) {
  if (_writerSlice != nullptr)
    _writerSlice->addMer(_kmer);
}



////////////////////////////////////////
//
//  PRINTING to ASCII files.
//
void
merylOpCompute::addPrinter(merylOpTemplate *ot, uint32 slice) {
  char  T[FILENAME_MAX+1] = { 0 };
  char  N[FILENAME_MAX+1] = { 0 };

  if (ot->_printerName == nullptr)    //  If no name defined, no
    return;                           //  print output requested.

  if (ot->_printer != nullptr)        //  If defined, we're to a single file
    return;                           //  and use _ot->_printer.

  else {                              //  Otherwise, open a per-thread output.
    strncpy(T, ot->_printerName, FILENAME_MAX);

    char   *pre = T;
    char   *suf = strchr(T, '#');
    uint32  len = 0;

    while ((suf) && (*suf == '#')) {
      *suf = 0;
      len++;
      suf++;
    }

    snprintf(N, FILENAME_MAX, "%s%0*d%s", pre, len, slice, suf);

    _printerSlice = new compressedFileWriter(N);
  }
}



void
merylOpCompute::printKmer(void) {
  char   outstr[256] = {0};
  char  *outptr = outstr;
  kmer   pk     = _kmer;
  FILE   *f     = nullptr;

#warning this should be cached
  if      (_ot->_printer != nullptr)       //  Only one of these is ever set, but
    f = _ot->_printer->file();             //  we still need to figure out which
  else if (_printerSlice != nullptr)       //  one to use.
    f = _printerSlice->file();

  if (f == nullptr)                        //  If neither is set, no print is
    return;                                //  requested for this action.

  if (_ot->_printACGTorder == true)        //  If requested, recompute the canonical
    pk.recanonicalizeACGTorder();          //  mer in ACGT order.  Yuck.

  pk.toString(outptr);                     //  Convert the kmer to ASCII, then
  while (*outptr)                          //  advance to the end of the string.
    outptr++;

  *outptr++ = '\t';                        //  Add the value.  There is always
  outptr = toDec(_kmer._val, outptr);      //  a value to add.

#warning need to print label as binary or hex, user supplied
  if (kmer::labelSize() > 0) {             //  If a label exists, add it too.
    *outptr++ = '\t';
    outptr = toBin(_kmer._lab, outptr, kmer::labelSize());
  }

  *outptr++ = '\n';                        //  Terminate the string and
  *outptr++ = 0;                           //  emit it.

#pragma omp critical (printLock)           //  fputs() is not thread safe and will
  fputs(outstr, f);                        //  happily intermix on e.g. Linux.
}



////////////////////////////////////////
//
//  STATISTICS of kmer values.
//
//  Each thread will allocate a histogram, and this can get rather large.
//  Storing the first million values in an array will take
//    64 threads * 8 bytes/val * 1048576 values = 512 MB
//
void
merylOpCompute::addStatistics(merylOpTemplate *ot, uint32 slice) {

  if ((ot->_statsFile == nullptr) &&
      (ot->_histoFile == nullptr))
    return;

  _stats = new merylHistogram(1048576);
}



////////////////////////////////////////
//
//  COMPUTING the kmer/value/label to output.
//
//   - findOutputKmer() scans all the inputs to find the smallest kmer, then
//     creates a list of the inputs with that kmer.
//
//     _kmer._mer is not valid if _actLen is zero after this function.
//
//   - findOutputValue() and findOutputLabel() are reasonably straight
//     forward, but long, and compute an output value/label based on the
//     action specified.
//
void
merylOpCompute::findOutputKmer(void) {

  _actLen = 0;

  for (uint32 ii=0; ii<_inputs.size(); ii++) {
    kmdata kmer = _inputs[ii]->_kmer._mer;

    if (_inputs[ii]->_valid == false)      //  No more kmers in the file,
      continue;                            //  skip it.

    if ((_actLen > 0) &&                   //  If we've picked a kmer already,
        (kmer > _kmer._mer))               //  and this one is bigger,
      continue;                            //  skip this one.

    if ((_actLen > 0) &&                   //  If we've picked a kmer already,
        (kmer < _kmer._mer))               //  but this one is smaller,
      _actLen = 0;                         //  forget everything we've done.

    if (_actLen == 0)                      //  Pick this one if nothing picked yet.
      _kmer._mer = kmer;

    _act[_actLen]    = _inputs[ii]->_kmer; //  Copy the kmer/value/label to the list.
    _actIdx[_actLen] = ii;

    _actLen++;
  }
}



void
merylOpCompute::findOutputValue(void) {

  switch (_ot->_valueSelect) {
    case merylModifyValue::valueNOP:
      break;

    case merylModifyValue::valueSet:
      _kmer._val = _ot->_valueConstant;
      break;

    case merylModifyValue::valueSelected:
#warning wrong
      _kmer._val = _act[0]._val;  //  Need to figure out which input to select
      break;

    case merylModifyValue::valueFirst:
#warning wrong
      _kmer._val = _act[0]._val;  //  Or do we need to verify that actIdx[0] is 0?
      break;

    case merylModifyValue::valueMin:
      _kmer._val = _ot->_valueConstant;
      for (uint32 ii=0; ii<_actLen; ii++)
        _kmer._val = std::min(_kmer._val, _act[ii]._val);
      break;

    case merylModifyValue::valueMax:
      _kmer._val = _ot->_valueConstant;
      for (uint32 ii=0; ii<_actLen; ii++)
        _kmer._val = std::max(_kmer._val, _act[ii]._val);
      break;

    case merylModifyValue::valueAdd:
      _kmer._val = _ot->_valueConstant;
      for (uint32 ii=0; ii<_actLen; ii++)
        if (kmvalumax - _kmer._val < _act[ii]._val)
          _kmer._val = kmvalumax;
        else
          _kmer._val = _kmer._val + _act[ii]._val;
      break;

    case merylModifyValue::valueSub:
#warning wrong
      _kmer._val = _ot->_valueConstant + _ot->_valueConstant;
      for (uint32 ii=0; ii<_actLen; ii++)
        if (_kmer._val > _act[ii]._val)
          _kmer._val -= _act[ii]._val;
        else
          _kmer._val = 0;
      break;

    case merylModifyValue::valueMul:
      _kmer._val = _ot->_valueConstant;
      for (uint32 ii=0; ii<_actLen; ii++)
        if (kmvalumax / _kmer._val < _act[ii]._val)
          _kmer._val = kmvalumax;
        else
          _kmer._val = _kmer._val * _act[ii]._val;
      break;


#warning wrong
    case merylModifyValue::valueDiv:
      if (_ot->_valueConstant == 0)
        _kmer._val = 0;
      else
        _kmer._val = _kmer._val / _ot->_valueConstant;
      break;

    case merylModifyValue::valueDivZ:
      if (_ot->_valueConstant == 0)
        _kmer._val = 0;
      else if (_kmer._val < _ot->_valueConstant)
        _kmer._val = 1;
      else
        _kmer._val = (kmvalu)round(_kmer._val / (double)_ot->_valueConstant);
      break;

    case merylModifyValue::valueMod:
      if (_ot->_valueConstant == 0)
        _kmer._val = 0;
      else
        _kmer._val = _kmer._val % _ot->_valueConstant;
      break;

    case merylModifyValue::valueCount:
      _kmer._val = _actLen;
      break;
  }
}



void
merylOpCompute::findOutputLabel(void) {

  switch (_ot->_labelSelect) {
    case merylModifyLabel::labelNOP:
      //assert(_actLen == 1);
      _kmer._lab = _act[0]._lab;
      break;

    case merylModifyLabel::labelSet:
      _kmer._lab = _ot->_labelConstant;
      break;

#warning wrong
    case merylModifyLabel::labelSelected:
      _kmer._lab = _act[0]._lab;
      break;

#warning wrong
    case merylModifyLabel::labelFirst:
      _kmer._lab = _act[0]._lab;
      break;

    case merylModifyLabel::labelMin:
      _kmer._lab = _ot->_labelConstant;
      for (uint32 ll=0; ll<_actLen; ll++)
        _kmer._lab = std::min(_kmer._lab, _act[ll]._lab);
      break;
    case merylModifyLabel::labelMax:
      _kmer._lab = _ot->_labelConstant;
      for (uint32 ll=0; ll<_actLen; ll++)
        _kmer._lab = std::max(_kmer._lab, _act[ll]._lab);
      break;


    case merylModifyLabel::labelAnd:
      _kmer._lab = _ot->_labelConstant;
      for (uint32 ll=0; ll<_actLen; ll++)
        _kmer._lab &= _act[ll]._lab;
      break;

    case merylModifyLabel::labelOr:
      _kmer._lab = _ot->_labelConstant;
      for (uint32 ll=0; ll<_actLen; ll++)
        _kmer._lab |= _act[ll]._lab;
      break;

    case merylModifyLabel::labelXor:
      _kmer._lab = _ot->_labelConstant;
      for (uint32 ll=0; ll<_actLen; ll++)
        _kmer._lab ^= _act[ll]._lab;
      break;

    case merylModifyLabel::labelDifference:
#warning check this
      _kmer._lab = _act[0]._lab & ~_ot->_labelConstant;
      for (uint32 ll=1; ll<_actLen; ll++)
        _kmer._lab &= ~_act[ll]._lab;
      break;

    case merylModifyLabel::labelLightest:
      _kmer._lab = _ot->_labelConstant;
      for (uint32 ll=0; ll<_actLen; ll++)
        if (countNumberOfSetBits64(_act[ll]._lab) < countNumberOfSetBits64(_kmer._lab))
          _kmer._lab = _act[ll]._lab;
      break;

    case merylModifyLabel::labelHeaviest:
      _kmer._lab = _ot->_labelConstant;
      for (uint32 ll=0; ll<_actLen; ll++)
        if (countNumberOfSetBits64(_act[ll]._lab) > countNumberOfSetBits64(_kmer._lab))
          _kmer._lab = _act[ll]._lab;
      break;

    case merylModifyLabel::labelInvert:
      assert(_actLen == 1);
      _kmer._lab = ~_act[0]._lab;
      break;

    case merylModifyLabel::labelShiftLeft:
      assert(_actLen == 1);
      _kmer._lab = _act[0]._lab >> _ot->_labelConstant;
      break;

    case merylModifyLabel::labelShiftRight:
      assert(_actLen == 1);
      _kmer._lab = _act[0]._lab << _ot->_labelConstant;
      break;

#warning wrong
    case merylModifyLabel::labelRotateLeft:
      assert(_actLen == 1);
      _kmer._lab = _act[0]._lab >> _ot->_labelConstant;
      break;

#warning wrong
    case merylModifyLabel::labelRotateRight:
      assert(_actLen == 1);
      _kmer._lab = _act[0]._lab << _ot->_labelConstant;
      break;
  }
}
