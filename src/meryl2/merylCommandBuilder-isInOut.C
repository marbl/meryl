
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

  if ((fpath != nullptr) && (fpath[0] == 0))  //  Blow up if the name is empty.
    sprintf(_errors, "Operation #%u has no output path name supplied.\n", op->_ident);

  //  For cases 1 and 3, we have a filepath, and can finish setting up the
  //  output for this operation.  At the end, forget all about the output
  //  parameter, because we're completely done with it.

  if (fpath) {
    if (_curClass == opClass::clOutput)
      switch (_curPname) {
        case opPname::pnDB:      op->addOutputToDB  (fpath,          _errors);  break;
        case opPname::pnList:    op->addOutputToList(fpath,   false, _errors);  break;
        case opPname::pnShow:    op->addOutputToList(nullptr, false, _errors);  break;
        case opPname::pnPipe:    op->addOutputToPipe(fpath,          _errors);  break;
        case opPname::pnHisto:   op->addHistoOutput (fpath,          _errors);  break;
        case opPname::pnStats:   op->addStatsOutput (fpath,          _errors);  break;
        default:
          fprintf(stderr, "Got Pname '%s' for clOutput.\n", toString(_curPname));
          assert(0);
          break;
      }

    if (_curClass == opClass::clInput)
      switch (_curPname) {
        case opPname::pnDB:      op->addInputFromDB  (fpath,   _errors);   break;
        case opPname::pnList:    op->addInputFromList(fpath,   _errors);   break;
        case opPname::pnPipe:    op->addInputFromPipe(fpath,   _errors);   break;
        case opPname::pnAction:  op->addInputFromOp  (nullptr, _errors);   break;
        default:
          fprintf(stderr, "Got Pname '%s' for clInput.\n", toString(_curPname));
          assert(0);
          break;
      }

    _curClass = opClass::clNone;
    _curPname = opPname::pnNone;
    _curParam = nullptr;

    return true;
  }

  //  For case 2, leave _curClass and _curPname alone, but forget all about
  //  the (empty) _curParam.

  _curParam = nullptr;

  return true;
}



//
//  Handle input files that do not have a class/pname supplied; that is,
//  command words that are a meryl database, a canu seqStore or a sequence file.
//

bool
merylCommandBuilder::isInput(void) {
  merylOpTemplate  *op = getCurrent();

  bool   merylInput = (fileExists(_optString, '/', "merylIndex"));
  bool   canuInput  = (fileExists(_optString, '/', "info") &&
                       fileExists(_optString, '/', "reads"));
  bool   seqInput   = (fileExists(_optString));

  if      (merylInput)
    op->addInputFromDB(_optString, _errors);

  else if (canuInput)
    op->addInputFromCanu(_optString, _segment, _segmentMax, _errors);

  else if (seqInput)
    op->addInputFromSeq(_optString, _doCompression, _errors);

  else
    return false;

  _segment    = 1;   //  These only apply to Canu inputs, but no harm
  _segmentMax = 1;   //  in resetting them for every input.

  return true;
}
