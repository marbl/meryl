
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

//
//  Handles input AND output options.
//
//  Expects _curClass == opClass::clInput  OR  Expects _curClass == opClass::clOutput.
//  Expects _curPname == (as below)            Expects _curPname == (as below)
//    input:database   opPname::pnDB             output:database   opPname::pnDB
//    input:list       opPname::pnList           output:list       opPname::pnList
//    input:pipe       opPname::pnPipe           output:show       opPname::pnShow
//    input:action     opPname::pnAction         output:pipe       opPname::pnPipe
//
//  1) If _curParam is defined and not empty, the option word was of the form
//     'output:database=some.merylDB', and we can do all work in this single
//     call.  At the end of this function, all three _cur members are
//     undefined because we're completely finished.
//
//  2) If _curParam is defined but empty (*curParam == 0) then we have an
//     option word that only describes the parameter class/name, e.g.,
//     'output:database='.  This is an error.
//
//  3) If _curParam is nullptr, then the option word was 'bare'; 'output:database'.
//     This is an error, except for pnShow, pnHisto and pnStats.
//

#if 0
bool
merylCommandBuilder::isInOut(void) {
  merylOpTemplate  *op = getCurrent();

  if ((_curClass != opClass::clOutput) &&
      (_curClass != opClass::clInput))
    return false;

  //
  //  Check for error conditions.
  //

  if      ((_curPname == opPname::pnShow) &&    //  pnShow will never have a filename.
           (_curParam != nullptr)) {
    sprintf(_errors, "Operation #%u specified a file path for output:show in '%s'.\n", op->_ident, _optString);
    goto finish;
  }
  else if ( (_curClass == opClass::clOutput) &&  //  pnHisto and pnStats are allowed to have or
           ((_curPname == opPname::pnHisto) ||   //  not have a filename, ignore them for now.
            (_curPname == opPname::pnStats))) {
  }
  else if ((_curParam == nullptr) ||            //  Everything else must have a filename.
           (_curParam[0] == 0)) {
    sprintf(_errors, "Operation #%u %s:%s has no output path name supplied.\n",
            op->_ident, toString(_curClass), toString(_curPname));
    goto finish;
  }

  //
  //  With errors out of the way, process as normal.  This does, however, still create
  //  the outputs with error conditions - pnShow could have a filename - which will
  //  be a problem if we ever stop stopping on error conditions.
  //

  if (_curClass == opClass::clOutput) {
    switch (_curPname) {
      case opPname::pnDB:      op->addOutputToDB  (_curParam,        _errors);  break;
      case opPname::pnList:    op->addOutputToList(_curParam, false, _errors);  break;
      case opPname::pnShow:    op->addOutputToList(_curParam, false, _errors);  break;
      case opPname::pnPipe:    op->addOutputToPipe(_curParam,        _errors);  break;
      case opPname::pnHisto:   op->addHistoOutput (_curParam,        _errors);  break;
      case opPname::pnStats:   op->addStatsOutput (_curParam,        _errors);  break;
      default:
        sprintf(_errors, "Got Pname '%s' for clOutput.\n", toString(_curPname));
        break;
    }
  }

  if (_curClass == opClass::clInput) {
    merylInput *in = new merylInput;
    bool        sx = false;

    switch (_curPname) {
      case opPname::pnDB:      sx = in->registerMerylDB   (_curParam, _errors);   break;
      case opPname::pnList:    sx = in->registerMerylList (_curParam, _errors);   break;
      case opPname::pnPipe:    sx = in->registerActionPipe(_curParam, _errors);   break;
      //se opPname::pnAction:  sx = in->registerTemplate  (nullptr,   _errors);   break;
      default:
        sprintf(_errors, "Got Pname '%s' for clInput.\n", toString(_curPname));
        break;
      }

    if (sx)
      op->addInput(in);
    else
      delete in;
  }

 finish:
  _curClass = opClass::clNone;
  _curPname = opPname::pnNone;
  _curParam = nullptr;

  return true;
}

#endif




//
//  Handle input files that do not have a class/pname supplied; that is,
//  command words that are a meryl database, a canu seqStore or a sequence file.
//

bool
merylCommandBuilder::isInput(char const *word) {

  merylOpTemplate  *op = getCurrent();
  merylInput       *in = new merylInput;

  bool as = op->acceptsSequenceInputs();

  bool im =                  (as == false) && (in->registerMerylDB (word,                        _errors, false));
  bool ic = (im == false) && (as == true)  && (in->registerSeqStore(word, _segment, _segmentMax, _errors, false));
  bool is = (ic == false) && (as == true)  && (in->registerSeqFile (word, _doCompression,        _errors, false));

  if ((im == false) &&
      (ic == false) &&
      (is == false)) {
    delete in;
    return false;
  }

  op->addInput(in);

  //_segment    = 1;   //  These only apply to Canu inputs, but no harm
  //_segmentMax = 1;   //  in resetting them for every input.

  return true;
}
