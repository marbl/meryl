
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
//     'output:database'.  All we can do here is remember remember that we
//     still need to get the parameter value from the next word.  This is
//     done by setting _curParam to nullptr and leaving _curClass and
//     _curPname set.  The next time through processWord() we'll handle the
//     left-over _curClass/_curPname BEFORE we try to decode it from the
//     current word, and end up back here for the last case...
//
//  3) ...where we handle the continuation of the above.  _curParam is
//     nullptr, and the filename is in _optString.
//

bool
merylCommandBuilder::isInOut(void) {
  merylOpTemplate  *op = getCurrent();

  if ((_curClass != opClass::clOutput) &&
      (_curClass != opClass::clInput))
    return false;

  char const *fpath = (_curParam == nullptr) ? _optString :            //  Case 3
                      (_curParam[0] != 0)    ? _curParam  : nullptr;   //  Case 1 : Case 2

  fprintf(stderr, "isInOut() enter with '%s' '%s' '%s' fpath '%s'\n", toString(_curClass), toString(_curPname), _curParam, fpath);

  if ((fpath != nullptr) && (fpath[0] == 0))  //  Blow up if the name is empty.
    sprintf(_errors, "Operation #%u has no output path name supplied.\n", op->_ident);

  //  For pnShow, there will never be a filename.  Unlike the other Pnames,
  //  we need to accept and complete this one immediately.

  if ((_curClass == opClass::clOutput) &&
      (_curPname == opPname::pnShow)) {
    op->addOutputToList(nullptr, false, _errors);

    if (fpath != nullptr)
      sprintf(_errors, "Operation #%u specified a file path for output:display in '%s'.\n", op->_ident, _optString);

    _curClass = opClass::clNone;
    _curPname = opPname::pnNone;
    _curParam = nullptr;

    return true;
  }

  //  For cases 1 and 3, we have a filepath, and can finish setting up the
  //  output for this operation.  At the end, forget all about the output
  //  parameter, because we're completely done with it.

  if (fpath) {
    if (_curClass == opClass::clOutput) {
      switch (_curPname) {
        case opPname::pnDB:      op->addOutputToDB  (fpath,          _errors);  break;
        case opPname::pnList:    op->addOutputToList(fpath,   false, _errors);  break;
        case opPname::pnShow:    op->addOutputToList(nullptr, false, _errors);  break;  //  Should never happen!
        case opPname::pnPipe:    op->addOutputToPipe(fpath,          _errors);  break;
        case opPname::pnHisto:   op->addHistoOutput (fpath,          _errors);  break;
        case opPname::pnStats:   op->addStatsOutput (fpath,          _errors);  break;
        default:
          fprintf(stderr, "Got Pname '%s' for clOutput.\n", toString(_curPname));
          assert(0);
          break;
      }
    }

    if (_curClass == opClass::clInput) {
      merylInput *in = new merylInput;
      bool        sx = false;

      switch (_curPname) {
        case opPname::pnDB:      sx = in->registerMerylDB   (fpath,   _errors);   break;
        case opPname::pnList:    sx = in->registerMerylList (fpath,   _errors);   break;
        case opPname::pnPipe:    sx = in->registerActionPipe(fpath,   _errors);   break;
        //se opPname::pnAction:  sx = in->registerTemplate  (nullptr, _errors);   break;
        default:
          fprintf(stderr, "Got Pname '%s' for clInput.\n", toString(_curPname));
          assert(0);
          break;
      }

      if (sx)
        op->addInput(in);
      else
        delete in;
    }

    _curClass = opClass::clNone;
    _curPname = opPname::pnNone;
    _curParam = nullptr;

    fprintf(stderr, "isInOut() consumed1 '%s'\n", _optString);

    return true;
  }

  //  For case 2, leave _curClass and _curPname alone, but forget all about
  //  the (empty) _curParam.

  _curParam = nullptr;

  fprintf(stderr, "isInOut() consumed2 '%s'\n", _optString);

  return true;
}



//
//  Handle input files that do not have a class/pname supplied; that is,
//  command words that are a meryl database, a canu seqStore or a sequence file.
//

bool
merylCommandBuilder::isInput(void) {

  if (_curClass != opClass::clNone)    //  This function only handles stray input
    return false;                      //  files; all others are handled above.

  merylOpTemplate  *op = getCurrent();
  merylInput       *in = new merylInput;

  bool as = op->acceptsSequenceInputs();

  bool im =                  (as == false) && (in->registerMerylDB (_optString,                        _errors, false));
  bool ic = (im == false) && (as == true)  && (in->registerSeqStore(_optString, _segment, _segmentMax, _errors, false));
  bool is = (ic == false) && (as == true)  && (in->registerSeqFile (_optString, _doCompression,        _errors, false));

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
