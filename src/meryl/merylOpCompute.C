
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

  //  Initialize _act so that each input is listed.  This forces
  //  merylOpCompute::nextMer() to call nextMer() on the first iteration.

  for (uint32 ii=0; ii<nInputs; ii++) {
    _actIdx[ii] = ii;
    _actRdx[ii] = ii;
  }

  _actLen = nInputs;
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

  if (ot->_printer != nullptr) {      //  If defined, we're to a single file
    _printer = ot->_printer;          //  and use _ot->_printer.
    return;
  }

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
    _printer      = _printerSlice;
  }
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

    _actRdx[ii] = uint32max;               //  Wipe the reverse map.

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

    _act   [_actLen] = _inputs[ii]->_kmer; //  Copy the kmer/value/label to the list.
    _actIdx[_actLen] = ii;
    _actRdx[ii]      = _actLen;

    _actLen++;
  }
}



void
merylOpCompute::findOutputValue(void) {
  kmvalu  q = 0;
  kmvalu  r = 0;

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
      _kmer._val = _act[0]._val;

      for (uint32 ii=1; ii<_actLen; ii++)
        if (_kmer._val > _act[ii]._val)
          _kmer._val -= _act[ii]._val;
        else
          _kmer._val = 0;

      if (_kmer._val > _ot->_valueConstant)
        _kmer._val -= _ot->_valueConstant;
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


    case merylModifyValue::valueDiv:
      _kmer._val = _act[0]._val;

      for (uint32 ii=1; ii<_actLen; ii++)
        if (_act[ii]._val > 0)
          _kmer._val /= _act[ii]._val;
        else
          _kmer._val = 0;

      if (_ot->_valueConstant > 0)
        _kmer._val /= _ot->_valueConstant;
      else
        _kmer._val = 0;

      break;


      //  Division, but with rounding instead of truncation.
      //  Additionally, values between 0 and 0.5 are rounded up to 1.
      //
      //  However, division by zero results in a 0 output.
    case merylModifyValue::valueDivZ:
      _kmer._val = _act[0]._val;

      for (uint32 ii=1; ii<_actLen; ii++)
        if (_act[ii]._val == 0)
          _kmer._val = 0;
        else if (_kmer._val < _act[ii]._val)
          _kmer._val = 1;
        else
          _kmer._val = round(_kmer._val / (double)_act[ii]._val);

      if (_ot->_valueConstant == 0)
        _kmer._val = 0;
      else if (_kmer._val < _ot->_valueConstant)
        _kmer._val = 1;
      else
        _kmer._val = round(_kmer._val / (double)_ot->_valueConstant);

      break;


    case merylModifyValue::valueMod:
      q = _act[0]._val;
      r = 0;

      for (uint32 ii=1; ii<_actLen; ii++)
        if (_act[ii]._val > 0) {
          kmvalu  qt = q / _act[ii]._val;

          r += q - qt * _act[ii]._val;
          q  = qt;
        } else {
          r += q;
          q  = 0;
        }

      if (_ot->_valueConstant > 0) {
        kmvalu  qt = q / _ot->_valueConstant;

        r += q - qt * _ot->_valueConstant;
        q  = qt;
      } else {
        r += q;
        q  = 0;
      }

      _kmer._val = r;

      break;


    case merylModifyValue::valueCount:
      _kmer._val = _actLen;
      break;
  }
}



void
merylOpCompute::findOutputLabel(void) {
  kmlabl  l;
  kmvalu  v;

  switch (_ot->_labelSelect) {
    case merylModifyLabel::labelNOP:
      break;

    case merylModifyLabel::labelSet:
      _kmer._lab = _ot->_labelConstant;
      break;

#warning wrong
    case merylModifyLabel::labelSelected:
      _kmer._lab = _act[0]._lab;
      break;

    case merylModifyLabel::labelFirst:
      _kmer._lab = _act[0]._lab;
      break;

    case merylModifyLabel::labelMin:
      l = _ot->_labelConstant;
      v = kmvalumax;

      for (uint32 ll=0; ll<_actLen; ll++)
        if (_act[ll]._val < v) {
          l = _act[ll]._lab;
          v = _act[ll]._val;
        }

      _kmer._lab = l;
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
