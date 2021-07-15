
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



void
merylCommandBuilder::processWord(char const *opt) {
  uint32  terminate = 0;

  if (verbosity.showConstruction() == true)
    fprintf(stderr, "processWord()- arg '%s'\n", opt);

  //  Save a copy of the string.

  _optStringLen = 0;

  while ((_optStringLen < FILENAME_MAX) && (opt[_optStringLen] != 0)) {
    _optString[_optStringLen] = opt[_optStringLen];
    _optStringLen++;
  }

  _optString[_optStringLen] = 0;

  //  If there is a '[' at the start of the string, push on a new operation and
  //  remove the bracket from the string.

  if (_optString[0] == '[') {
    addNewOperation();

    for (uint32 ii=0; ii<_optStringLen; ii++)
      _optString[ii] = _optString[ii+1];

    _optStringLen--;
  }

  //  If the stack is still empty, tsk tsk, the dear user didn't supply
  //  an initial command creating '[' -- which usually happens with the
  //  legacy usage, e.g., 'meryl output x.meryl union ....'.

  if (_opStack.size() == 0)
    addNewOperation();

  //  If there are ']' at the end of the string, strip them off and remember that
  //  we need to close the command on the stack after we process this arg.
  //  We can get any number of closing brackets.

  while ((_optStringLen > 0) &&
         (_optString[_optStringLen-1] == ']')) {
    _optString[_optStringLen-1] = 0;
    _optStringLen--;

    terminate++;
  }

  if (verbosity.showConstruction() == true) {
    if (terminate == 1)
      fprintf(stderr, "processWord()- Found a terminator.\n");
    else if (terminate > 1)
      fprintf(stderr, "processWord()- Found %u terminators.\n", terminate);
  }

  //  Save a few copies of the command line word expanded to useful paths.

  strncpy(_inoutName, _optString, FILENAME_MAX + 1);

  snprintf(_indexName, FILENAME_MAX, "%s/merylIndex", _optString);
  snprintf(_sqInfName, FILENAME_MAX, "%s/info",       _optString);
  snprintf(_sqRdsName, FILENAME_MAX, "%s/reads",      _optString);

  //  Now use a gigantic short-circuiting if test to actually parse the word.
  //  Each isFunction() will return true if it consumed the word we're trying
  //  to parse, and that will cause the whole test to fail.  If all of the
  //  isFunctions() return false, the word isn't recognized and an error is
  //  generated.
  //
  //  isInput() MUST come after isPrinter() so we can correctly handle
  //  'print some.meryl'.

  if ((isEmpty()            == false) &&   //  Consumes empty words.
      (isOption()           == false) &&   //  Consumes key=value options.
      (isAlias()            == false) &&   //  Consumes 'union', 'intersect', et cetera, aliases
      (isValueSelect()      == false) &&   //  Consumes 'value='
      (isLabelSelect()      == false) &&   //  Consumes 'label='
      (isValueFilter()      == false) &&   //  Consumes 'value:' filters
      (isLabelFilter()      == false) &&   //  Consumes 'label:' filters
      (isBasesFilter()      == false) &&   //  Consumes 'bases:' filters
      (isInputFilter()      == false) &&   //  Consumes 'input:' filters
      (isFilterConnective() == false) &&   //  Consumes 'and', 'or', 'not'
      (isCount()            == false) &&   //  Consumes 'count', 'count-forward', 'count-reverse'
      (isOutput()           == false) &&   //  Consumes 'output' and related database name
      (isPrinter()          == false) &&   //  Consumes 'print' and related output name
      (isInput()            == false))     //  Consumes inputs
    addError("Can't interpret '%s': not a meryl command, option, or recognized input file.", _optString);

  //  Process any operation close events (']') that were at the end of this
  //  word.  We discovered and stripped them out above.

  for (; terminate > 0; terminate--) {
    if (verbosity.showConstruction() == true)
      fprintf(stderr, "processWord()- Pop operation from top of stack.\n");
    if (_opStack.size() > 0)
      _opStack.pop();
    else
      addError("processWord()- Extra ']' encountered in command line.\n");
  }

  if (verbosity.showConstruction() == true)
    fprintf(stderr, "\n");


  if (verbosity.showConstruction() == true)
    for (uint32 rr=0; rr<numTrees(); rr++)
      printTree(getTree(rr), 0, 0);
}



