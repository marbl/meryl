
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


//  This is called by processWord() when an '[' is encountered in the input.
//  It will _always_ push a new operation onto the stack (unless the existing
//  operation is completely empty).  Depending on the stack contents, there
//  are two variants:
//    stack empty - add the new op to the list of root operations
//    stack not   - add the new op to the list of inputs for the current top op
//
void
merylCommandBuilder::addNewOperation(void) {
  merylOpTemplate *eop = getCurrent();

  if ((eop != nullptr) &&                    //  If an operation exists and it is
      (eop->_isSelector == false) &&         //  unused, use it instead of making
      (eop->_isCounting == false))           //  a new one.
    return;                                  //

  merylOpTemplate *nop = new merylOpTemplate(_opList.size() + 1);

  _opStack.push(nop);                        //  Add new node to the stack
  _opList.push_back(nop);                    //  and list of nodes.

  if (eop)                                   //  If there is an existing node,
    eop->addInput(new merylInput(nop));      //  add 'us' as input to it,
  else                                       //  otherwise, add 'us' as a new
    _opTree.push_back(_opList.size() - 1);   //  root node.

  if (globals.showConstruction() == true)    //  Now just some logging.
    fprintf(stderr, "addOperation()- Add new operation %s %ld.\n",
            (eop) ? "to stack at position" : "as root operation for tree",
            (eop) ? _opStack.size()        : _opTree.size());
}



//  Terminate up to n (default 1) operations.
//
//  Returns false if n is larger than the number of operations on the stack,
//  which is probably from too many ']'s on the command line.
//
bool
merylCommandBuilder::terminateOperations(uint32 nter, bool pop) {
  merylOpTemplate  *op = getCurrent();

  if (nter == 0)              return true;    //  Successfully terminated the requested number.
  if (_opStack.size() == 0)   return false;   //  Oops, too many terminates requested.

  if (globals.showConstruction() == true)
    fprintf(stderr, "terminateOperations()- Terminate action #%u at stack position %ld.\n",
            op->_ident, _opStack.size()-1);

 if (_curClass != opClass::clNone)
    sprintf(_errors, "terminateOperations()- Action #%u at stack position %ld did not process _curClass '%s'.\n",
            op->_ident, _opStack.size()-1, toString(_curClass));
  if (_curPname != opPname::pnNone)
    sprintf(_errors, "terminateOperations()- Action #%u at stack position %ld did not process _curPname '%s'.\n",
            op->_ident, _opStack.size()-1, toString(_curPname));

  if ((op->_isCounting == false) &&    //  If neither set, this action might be doing nothing.
      (op->_isSelector == false)) {    //  It could still be doing just histo/stats or output.
    fprintf(stderr, "terminateOperations()- Action #%u at stack position %ld possibly does nothing.\n",
            op->_ident, _opStack.size()-1);
    op->_isSelector = true;
  }

  //assert(_curClass == opClass::clNone);       //  If either of these are set, the action
  //assert(_curPname == opPname::pnNone);       //  was never fully parsed/added.


  if (op->_valueAssign == merylAssignValue::valueNOP)  op->_valueAssign = merylAssignValue::valueFirst;
  if (op->_labelAssign == merylAssignLabel::labelNOP)  op->_labelAssign = merylAssignLabel::labelFirst;

  //  Forget any files we've been remembering.

  _printACGTorder = false;

  //  Pop the action if told to, then recursively call to handle the
  //  remaining terminates.

  if (pop == true)
    _opStack.pop();

  return terminateOperations(nter-1, pop);
}



#if 0
bool
merylCommandBuilder::isCountingWord(void) {
  merylOpTemplate  *op = getCurrent();

  if ((strcmp(_optString, "count")         != 0) &&
      (strcmp(_optString, "count-forward") != 0) &&
      (strcmp(_optString, "count-reverse") != 0))
    return(false);

  if      (op->_isCounting)
    sprintf(_errors, "ERROR: operation is already a counting operation.\n");

  else if (op->_isSelector)
    sprintf(_errors, "ERROR: operation is a select operation, cannot also be a counting operation.\n");

  else {
    op->_isCounting = true;
    op->_counting   = new merylOpCounting(_optString);
  }

  return true;  //  Always consume the word, even if it triggered an error above.
}
#endif








//  We've buit the processing tree structure.  Walk through it looking for
//  errors or inconsistencies and opening inputs.
//
//   - terminateOperations() will close any still-open operations.  This can
//     be from (valid, legacy) usage or from (invalid but allowed) usage:
//       meryl act in1 in2
//       meryl [ act1 [ act2 in1 in2
//
//   - Tell each action that it has all of its inputs, and it is safe to open
//     them.  This will also check that the number of inputs is valid and add
//     an error if not.
//
//   - ensure that a kmer size is known
//
void
merylCommandBuilder::finalizeTrees(void) {

  terminateOperations(_opStack.size(), true);

  for (uint32 oo=0; oo<_opList.size(); oo++)
    _opList[oo]->finalizeTemplateInputs(oo, _errors);

  if (kmer::merSize() == 0)
    sprintf(_errors, "ERROR: mer size not set.\n");
}


//  Do counting.
void
merylCommandBuilder::performCounting(uint64 allowedMemory, uint32 allowedThreads) {
  for (uint32 oo=0; oo<_opList.size(); oo++)
    _opList[oo]->doCounting(allowedMemory, allowedThreads);
}


//  Tell all the actions to do their final initialization.
void
merylCommandBuilder::finalizeParameters(void) {
  for (uint32 oo=0; oo<_opList.size(); oo++)
    _opList[oo]->finalizeTemplateParameters();
}








bool
merylCommandBuilder::displayTreesAndErrors(void) {

  if (globals.showStandard() == true) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Found %u command tree%s.\n", numTrees(), (numTrees() == 1) ? "" : "s");

    for (uint32 ii=0; ii<numTrees(); ii++)
      printTree(ii);
  }

  if (numErrors() == 0)
    return false;

  fprintf(stderr, "Errors detected in the action tree%s:\n", (numTrees() == 1) ? "" : "s");

  for (uint32 ii=0; ii<numErrors(); ii++)
    if (getErrors()[ii] != nullptr)
      fprintf(stderr, "  %s\n", getErrors()[ii]);

  return true;
}






void
merylCommandBuilder::runThreads(uint32 allowedThreads) {

  setNumThreads(globals.allowedThreads());

  for (uint32 rr=0; rr<numTrees(); rr++) {
    merylOpTemplate *tpl = getTree(rr);

    if (tpl->_inputs.size() == 0)   //  A count operation (at the root node) but it's
      continue;                     //  all done now and there is nothing to do.

    if (globals.showStandard() == true) {
      fprintf(stderr, "\n");
      fprintf(stderr, "PROCESSING TREE #%u using %u thread%s.\n", rr+1, getMaxThreadsAllowed(), getMaxThreadsAllowed() == 1 ? "" : "s");
    }
    //printTree(tpl, 0, 11);

#pragma omp parallel for schedule(dynamic, 1)
    for (uint32 ff=0; ff<64; ff++) {
      merylOpCompute *cpu = getTree(rr, ff);

      while (cpu->nextMer() == true)
        ;
    }

    //  Signal that we're done processing.  This will (recursively) collect
    //  any statistics the user has requested be generated.

    tpl->finishAction();

    //  Simply deleting the root template node is enough to delete the entire
    //  tree, including the compute nodes.
    //
    //  NOTE that the pointers in _opList and _thList are left dangling!

    delete tpl;
  }

}
