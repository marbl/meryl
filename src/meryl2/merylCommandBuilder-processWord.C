
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
//  Detect and strip action create '[' and terminate ']' symbols from the
//  word, returning the number of terminate symbols found.
//
#if 0
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

  return t;                                                //  Return the number of terminate encountered.
}



//
//  Decode a word of the form 'class:name[:=]value'.
//
bool
merylCommandBuilder::decodeWord(char const *opt) {
  char const *pn = nullptr;             //  Eventual location of parameter name in opt.
  char const *pv = nullptr;             //  Eventual location of parameter value in opt.

  //fprintf(stderr, "decodeWord() enter with '%s' '%s' '%s'\n", toString(_curClass), toString(_curPname), _curParam);

  if (_curClass != opClass::clNone)     //  If we already have a class set, do not attempt
    return false;                       //  to decode the word; we're waiting for the value.

  _curClass = opClass::clNone;
  _curPname = opPname::pnNone;
  _curParam = nullptr;

  //  Handle legacy 'output' usage.

#if 0
  if ((strcmp(opt, "output") == 0) ||   //  If we see just 'output' or 'out', instead of
      (strcmp(opt, "out")    == 0)) {   //  output:database, assume this is a database
    _curClass = opClass::clOutput;      //  output and a well-trained meryl-v1 user.
    _curPname = opPname::pnDB;
    _curParam = (opt[6] == 0) ? (opt + 6)   //  _curParam MUST be NUL so that isInOut()
                              : (opt + 3);  //  knows that the parameter doens't exist yet.
    return true;
  }

  if (matchToken(opt, pn, "output=") == true) {
    _curClass = opClass::clOutput;
    _curPname = opPname::pnDB;
    _curParam = opt + 7;

    return true;
  }
#endif

  //  Decode the parameter class.

  if      (matchToken(opt, pn, "output:")     == true)  _curClass = opClass::clOutput;
  else if (matchToken(opt, pn, "assign:")     == true)  _curClass = opClass::clAssign;
  else if (matchToken(opt, pn, "set:", true)  == true)  _curClass = opClass::clAssign;
  else if (matchToken(opt, pn, "select:")     == true)  _curClass = opClass::clSelect;
  else if (matchToken(opt, pn, "get:", true)  == true)  _curClass = opClass::clSelect;
  else if (matchToken(opt, pn, "input:")      == true)  _curClass = opClass::clInput;
  else
    goto returnfailure;

#warning improve error handling here

  if ((pn == nullptr) ||                      //  If we matched something without a separator,
      (*pn == 0)) {                           //  or with a missing Pname, warn and fail.
    sprintf(_errors, "Expecting class:name in parameter '%s'.", opt);
    goto returnfailure;
  }

  //  Decode the OUTPUT parameter name.

  if (_curClass == opClass::clOutput) {
    if      (matchToken(pn, pv, "database=")    == true)  _curPname = opPname::pnDB;
    else if (matchToken(pn, pv, "db=", true)    == true)  _curPname = opPname::pnDB;
    else if (matchToken(pn, pv, "list=")        == true)  _curPname = opPname::pnList;
    else if (matchToken(pn, pv, "show")         == true)  _curPname = opPname::pnShow;
    else if (matchToken(pn, pv, "stdout", true) == true)  _curPname = opPname::pnShow;
    else if (matchToken(pn, pv, "pipe=")        == true)  _curPname = opPname::pnPipe;
    else if (matchToken(pn, pv, "histogram=")   == true)  _curPname = opPname::pnHisto;
    else if (matchToken(pn, pv, "statistics=")  == true)  _curPname = opPname::pnStats;
    else if (matchToken(pn, pv, "stats=", true) == true)  _curPname = opPname::pnStats;
    else {
      sprintf(_errors, "Expecting output:database=, output:list=, output:show, output:pipe=,");
      sprintf(_errors, "  output:histogram[=] or output:statistics[=] in parameter '%s'.", opt);
      goto returnfailure;
    }

    _curParam = pv;
  }

  //  Decode the ASSIGN parameter name.

  if (_curClass == opClass::clAssign) {
    if ((matchToken(pn, pv, "value:") == true) ||
        (matchToken(pn, pv, "label:") == true) ||
        (matchToken(pn, pv, "bases:") == true) ||
        (matchToken(pn, pv, "input:") == true)) {
      sprintf(_errors, "Expecting assign:value= or assign:label= in paramter '%s'.", opt);
      sprintf(_errors, "  Note '=' instead of ':'.", opt);
      goto returnfailure;
    }

    if      (matchToken(pn, pv, "value=")       == true)  _curPname = opPname::pnValue;
    else if (matchToken(pn, pv, "label=")       == true)  _curPname = opPname::pnLabel;
    else {
      sprintf(_errors, "Expecting assign:value= or assign:label= in parameter '%s'.", opt);
      goto returnfailure;
    }

    _curParam = pv;
  }

  //  Decode the SELECT parameter name.

  if (_curClass == opClass::clSelect) {
    if ((matchToken(pn, pv, "value=") == true) ||
        (matchToken(pn, pv, "label=") == true) ||
        (matchToken(pn, pv, "bases=") == true) ||
        (matchToken(pn, pv, "input=") == true)) {
      sprintf(_errors, "Expecting select:value:, select:label:, select:bases: or select:input:");
      sprintf(_errors, "  in parameter '%s'.  Note ':' instead of '='.", opt);
      goto returnfailure;
    }

    if      (matchToken(pn, pv, "value:")       == true)  _curPname = opPname::pnValue;
    else if (matchToken(pn, pv, "label:")       == true)  _curPname = opPname::pnLabel;
    else if (matchToken(pn, pv, "bases:")       == true)  _curPname = opPname::pnBases;
    else if (matchToken(pn, pv, "input:")       == true)  _curPname = opPname::pnInput;
    else {
      sprintf(_errors, "Expecting select:value:, select:label:, select:bases: or select:input:");
      sprintf(_errors, "  in parameter '%s'.\n", opt);
      goto returnfailure;
    }

    _curParam = pv;
  }

  //  Decode the INPUT parameter name.

  if (_curClass == opClass::clInput) {
    if      (matchToken(pn, pv, "database=")    == true)  _curPname = opPname::pnDB;
    else if (matchToken(pn, pv, "db=", true)    == true)  _curPname = opPname::pnDB;
    else if (matchToken(pn, pv, "list=")        == true)  _curPname = opPname::pnList;
    else if (matchToken(pn, pv, "stdin", true)  == true)  _curPname = opPname::pnList;
    else if (matchToken(pn, pv, "pipe=")        == true)  _curPname = opPname::pnPipe;
    else if (matchToken(pn, pv, "action=")      == true)  _curPname = opPname::pnAction;
    else {
      sprintf(_errors, "Expecting input:database=, input:list=, input:stdin, input:pipe=");
      sprintf(_errors, "  or input:action= in paramter '%s'\n", opt);
      goto returnfailure;
    }

    _curParam = pv;
  }

  if (globals.showConstruction() == true)
    fprintf(stderr, "decodeWord()- '%s' -> %s::%s with param '%s'\n", opt, toString(_curClass), toString(_curPname), _curParam);

  return true;

 returnfailure:
  _curClass = opClass::clNone;
  _curPname = opPname::pnNone;
  _curParam = nullptr;

  return false;
}
#endif


//
//  
#if 0
void
merylCommandBuilder::processWord(char const *rawWord) {

  if (globals.showConstruction() == true)
    fprintf(stderr, "processWord()- '%s'\n", rawWord);

#if 1
  duplicateString(rawWord,              //  Make a copy of the word and count the length.
                  _optString, _optStringLen, _optStringMax);

  uint32 terminate = findCreateTerminate();

  if ((globals.showConstruction() == true) && (terminate > 0))
    fprintf(stderr, "processWord()- Found %s terminator%s.\n",
            (terminate == 1) ? "a" : toDec(terminate),
            (terminate == 1) ? ""  : "s");

  if (_opStack.size() == 0)             //  If stack is still empty, tsk, tsk, user didn't
    addNewOperation();                  //  explicitly make an action, so make one for them.

  //if ((isEmptyWord()  == true) ||       //  If the word is now empty (it was a single '[' or
  //    (isOptionWord() == true))         //  multiple ']'), or if the word is a recognized option,
  //  goto finishWord;                    //  we're done.

  //if ((isAliasConstant() == true))      //  If an alias constant -- or _expecting_ a constant
  //  goto finishWord;                    //  but one not found, we're done.

  if (_curClass != opClass::clNone) {   //  If we already have a class/pname set, we're waiting
    if ((isInOut()  == false) &&        //  for another argument for a previous word - usually
        (isAssign() == false) &&        //  an input or output path - so take care of it now.
        (isSelect() == false))          //  Bare inputs are handled in isInput() well below.
      fprintf(stderr, "Failed to decode second word '%s'.\n", rawWord);
    goto finishWord;
  }

  //if ((isAliasWord()    == true) ||     //  Handle aliases and count operations.  Ideally,
  //    (isCountingWord() == true))       //  these are only encountered with a completely empty
  //  goto finishWord;                    //  operation on the stack, but that isn't enforced.

  if (decodeWord(_optString) == true) { //  If we decode a word, try to parse it.
    if ((isInOut()   == false) &&       //  If all fail, report an error but still
        (isAssign()  == false) &&       //  consume the word (since we decoded
        (isSelect()  == false))         //  successfully).
      fprintf(stderr, "Failed to parse decoded word '%s'.\n", rawWord);
    goto finishWord;
  }

  if ((isSelectConnective() == true))   //  Handle 'and', 'or' and 'not'.
    goto finishWord;

  if ((isInput() == true))              //  This is last so that files such as 'o:d' are
    goto finishWord;                    //  treated as a parameter instead of an input.

  sprintf(_errors, "Can't interpret '%s': not a meryl command, option, or recognized input file.", rawWord);

 finishWord:
  if (terminateOperations(terminate, true) == false)
    sprintf(_errors, "processWord()- Extra ']' encountered in command line.");

  if (globals.showConstruction() == true)
    fprintf(stderr, "----------\n");
#endif
}


#endif
