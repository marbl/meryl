
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
#include "matchToken.H"


//
//  Parameter Class:     Parameter Name:
//    output               database, db
//                         list
//                         show, stdout
//                         pipe
//
//    assign, set          value
//                         label
//
//    select, get          value
//                         label
//
//    input                database, db
//                         list
//                         pipe
//                         action
//
//

//
//  Detect and strip action create '[' and terminate ']' symbols from the word.
//
uint32
merylCommandBuilder::findCreateTerminate(void) {
  uint32  t = 0;

  if (_optString[0] == '[') {                              //  If told to make a new action, close the
    terminateOperations(1);                                //  existing action, make a new one, and
    addNewOperation();                                     //  remove the '['.

    for (uint32 ii=0; ii<_optStringLen; ii++)              //  There should only ever by a single '[' and
      _optString[ii] = _optString[ii+1];                   //  never a space after it, unless the user is
    _optStringLen--;                                       //  being a jerk, in which case, we'll fail later
  }                                                        //  when we try to parse '[action' or ' action'.

  for (t=0; ((_optStringLen > 0) &&                        //  Remove any ']' at the end and keep track
             (_optString[_optStringLen-1] == ']')); t++)   //  of how many so we can terminate actions
    _optString[--_optStringLen] = 0;                       //  after this word is processed.

  return t;
}


bool
merylCommandBuilder::decodeWord(char const *opt) {
  char const *pn = nullptr;   //  Eventual location of parameter name in opt.

  if (_curClass != opClass::clNone)   //  If we already have a class set,
    return false;                     //  do not attempt to decode another one!

  _curClass = opClass::clNone;
  _curPname = opPname::pnNone;
  _curParam = nullptr;

  if      (matchToken(opt, pn, "output:")     == true)  _curClass = opClass::clOutput;
  else if (matchToken(opt, pn, "assign:")     == true)  _curClass = opClass::clAssign;
  else if (matchToken(opt, pn, "set:", true)  == true)  _curClass = opClass::clAssign;
  else if (matchToken(opt, pn, "select:")     == true)  _curClass = opClass::clSelect;
  else if (matchToken(opt, pn, "get:", true)  == true)  _curClass = opClass::clSelect;
  else if (matchToken(opt, pn, "input:")      == true)  _curClass = opClass::clInput;
  else
    return false;


  if (_curClass == opClass::clOutput) {
    if      (matchToken(pn, _curParam, "database")     == true)  _curPname = opPname::pnDB;
    else if (matchToken(pn, _curParam, "db", true)     == true)  _curPname = opPname::pnDB;
    else if (matchToken(pn, _curParam, "list")         == true)  _curPname = opPname::pnList;
    else if (matchToken(pn, _curParam, "show")         == true)  _curPname = opPname::pnShow;
    else if (matchToken(pn, _curParam, "stdout", true) == true)  _curPname = opPname::pnShow;
    else if (matchToken(pn, _curParam, "pipe")         == true)  _curPname = opPname::pnPipe;
    else if (matchToken(pn, _curParam, "histogram")    == true)  _curPname = opPname::pnHisto;
    else if (matchToken(pn, _curParam, "statistics")   == true)  _curPname = opPname::pnStats;
    else
      return false;
  }

  if (_curClass == opClass::clAssign) {
    if      (matchToken(pn, _curParam, "value=")    == true)  _curPname = opPname::pnValue;
    else if (matchToken(pn, _curParam, "label=")    == true)  _curPname = opPname::pnLabel;
    else
      return false;
  }

  if (_curClass == opClass::clSelect) {
    if      (matchToken(pn, _curParam, "value:")    == true)  _curPname = opPname::pnValue;
    else if (matchToken(pn, _curParam, "label:")    == true)  _curPname = opPname::pnLabel;
    else if (matchToken(pn, _curParam, "bases:")    == true)  _curPname = opPname::pnBases;
    else if (matchToken(pn, _curParam, "input:")    == true)  _curPname = opPname::pnInput;
    else
      return false;
  }

  if (_curClass == opClass::clInput) {
    if      (matchToken(pn, _curParam, "database")  == true)  _curPname = opPname::pnDB;
    else if (matchToken(pn, _curParam, "list")      == true)  _curPname = opPname::pnList;
    else if (matchToken(pn, _curParam, "pipe")      == true)  _curPname = opPname::pnPipe;
    else if (matchToken(pn, _curParam, "action")    == true)  _curPname = opPname::pnAction;
    else
      return false;
  }

  return true;
}




void
merylCommandBuilder::processWord(char const *rawWord) {

  if (globals.showConstruction() == true)
    fprintf(stderr, "processWord()- arg '%s'\n", rawWord);

  duplicateString(rawWord,    //  Make a copy of the word and count the length.
                  _optString, _optStringLen, _optStringMax);

  uint32 terminate = findCreateTerminate();

  if ((globals.showConstruction() == true) && (terminate > 0))
    fprintf(stderr, "processWord()- Found %s terminator%s.\n",
            (terminate == 1) ? "a" : toDec(terminate),
            (terminate == 1) ? ""  : "s");

  if (_opStack.size() == 0)             //  If stack is still empty, tsk, tsk, user didn't
    addNewOperation();                  //  explicitly make an action, so make one for them.

  if ((isEmpty()  == true) ||           //  If the word is now empty (it was a single '[' or
      (isOption() == true))             //  multiple ']'), or if the word is a recognized option,
    goto finishWord;                    //  we're done.

  if (_curClass != opClass::clNone) {   //  If we already have a class/pname set, we're waiting
    if ((isInOut()  == false) &&        //  for another argument for a previous word - usually
        (isAssign() == false) &&        //  an input or output path - so take care of it now.
        (isSelect() == false) &&        //  If nobody returns true, then we have failed to
        (isInput()  == false))          //  parse an expected parameter.
      fprintf(stderr, "Failed to decode second word '%s'.\n", rawWord);
    goto finishWord;
  }

#warning o:stats and o:histo need special case here when no file is supplied

  if (isInput() == true)                //  Handle merylDBs, canu seqStores and sequence
    goto finishWord;                    //  files.  This MUST come after _curClass != clNone!

  if ((isAlias()      == true) ||       //  Handle aliases and count operations.  Ideally,
      (isCount()      == true))         //  these are only encountered with a completely empty
    goto finishWord;                    //  operation on the stack, but that isn't enforced.

  if (decodeWord(_optString) == true) { //  If we decode a word, try to parse it.
    if ((isInOut()  == false) &&        //  If all fail, report an error but still
        (isAssign() == false) &&        //  consume the word (since we decoded
        (isSelect() == false))          //  successfully).
      fprintf(stderr, "Failed to parse decoded word '%s'.\n", rawWord);
    goto finishWord;
  }

  sprintf(_errors, "Can't interpret '%s': not a meryl command, option, or recognized input file.", rawWord);

 finishWord:
  if (terminateOperations(terminate, true) == false)
    sprintf(_errors, "processWord()- Extra ']' encountered in command line.");

  if (globals.showConstruction() == true)
    fprintf(stderr, "----------\n");
}



