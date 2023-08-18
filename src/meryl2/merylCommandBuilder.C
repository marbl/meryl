
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
//  It will _always_ push a new operation onto the stack.  Depending on the
//  stack contents, there are two variants:
//    stack empty - add the new op to the list of root operations
//    stack not   - add the new op to the list of inputs for the current top op
//
void
merylCommandBuilder::addNewOperation(void) {
  merylOpTemplate *nop = new merylOpTemplate(_opList.size() + 1);

  //  If the stack is empty, just add 'nop' to it as a root node.

  if (_opStack.size() == 0) {
    if (globals.showConstruction() == true)
      fprintf(stderr, "addOperation()- Add new operation as root operation.\n");

    _opStack.push(nop);
    _opList.push_back(nop);
    _opTree.push_back(_opList.size() - 1);

    return;
  }

  //  Otherwise, an operation exists on the stack.

  merylOpTemplate *eop = getCurrent();

  if (globals.showConstruction() == true) {
    fprintf(stderr, "addOperation()- Existing op: type '%s' ident %u value %s/%u label %s/%lu\n",
            toString(eop->_type),        eop->_ident,
            toString(eop->_valueAssign), eop->_valueConstant,
            toString(eop->_labelAssign), eop->_labelConstant);
  }

  //  If tht operation hasn't been used yet, do not add a new operation, just
  //  reuse what is there.

  if (eop->_type == merylOpType::opNothing) {
    fprintf(stderr, "addOperation()- Reuse existing empty operation at position %ld.\n",
            _opStack.size()-1);

    delete nop;

    return;
  }

  //  Finally, push 'nop' onto the stack as a child of the current operation.

  eop->addInput(new merylInput(nop));

  if (globals.showConstruction() == true)
    fprintf(stderr, "addOperation()- Add new operation to stack at position %ld.\n",
            _opStack.size());

  _opStack.push(nop);
  _opList.push_back(nop);
}



//  Terminate up to n (default 1) operations.
//
//  Returns false if n is larger than the number of operations on the stack,
//  which is probably from too many ']'s on the command line.
//
bool
merylCommandBuilder::terminateOperations(uint32 nter, bool pop) {

  if (nter == 0)              return true;    //  Successfully terminated the requested number.
  if (_opStack.size() == 0)   return false;   //  Oops, too many terminates requested.

  assert(_curClass == opClass::clNone);       //  If either of these are set, the action
  assert(_curPname == opPname::pnNone);       //  was never fully parsed/added.

  merylOpTemplate  *op = getCurrent();

  if (globals.showConstruction() == true)
    fprintf(stderr, "terminateOperations()- Terminate operation #%u at stack position %ld.\n",
            op->_ident, _opStack.size()-1);

  //  Set unset things in the action to defaults.

  if (op->_type        == merylOpType::opNothing)      op->_type        = merylOpType::opFilter;
  if (op->_valueAssign == merylAssignValue::valueNOP)  op->_valueAssign = merylAssignValue::valueFirst;
  if (op->_labelAssign == merylAssignLabel::labelNOP)  op->_labelAssign = merylAssignLabel::labelFirst;

  //  Finish file operations.

#if 0
  if (_isPrint)       op->addPrinter(nullptr, _printACGTorder, _errors);
  if (_isHistogram)   op->addHistogram("-", false, _errors);
  if (_isStatistics)  op->addHistogram("-", true, _errors);
#endif

  //  Forget any files we've been remembering.

  _isPrint        = false;
  _printACGTorder = false;
  _isHistogram    = false;
  _isStatistics   = false;

  //  Pop the action if told to, then recursively call to handle the
  //  remaining terminates.

  if (pop == true)
    _opStack.pop();

  return terminateOperations(nter-1, pop);
}



bool
merylCommandBuilder::isCount(void) {
  merylOpTemplate  *op = getCurrent();

  if ((strcmp(_optString, "count")         != 0) &&
      (strcmp(_optString, "count-forward") != 0) &&
      (strcmp(_optString, "count-reverse") != 0))
    return(false);

  if (op->_type == merylOpType::opCounting)
    sprintf(_errors, "ERROR: operation is already a counting operation.\n");

  op->_type     = merylOpType::opCounting;
  op->_counting = new merylOpCounting(_optString);

  return(true);
}









//  This doesn't really build the trees, since they're built as words are
//  added.
//
//  This does go through all the actions to make sure the inputs, outputs and
//  selectors are appropriate for that action, generating error messagess if
//  needed.
//
void
merylCommandBuilder::buildTrees(void) {

  //
  //  Clear the stack, we're done with it.  This handles, I think, only the
  //  case where the command line has no '[' or ']' symbols:
  //    meryl act in1 in2
  //  but does allow technically invalid commands that do not explicitly
  //  terminate any actions:
  //    meryl [ act1 [ act2 in1 in2
  //

  terminateOperations(_opStack.size(), true);

  //  Tell each action that it has all the inputs it is going to get,
  //  and add errors if there are too many or too few of them.

  for (uint32 oo=0; oo<_opList.size(); oo++) {
    merylOpTemplate  *ot = _opList[oo];

    ot->finalizeTemplateInputs(_errors);

    if (ot->_inputs.size() < ot->_inputsMin)
      sprintf(_errors, "ERROR: operation at position %u has %u inputs, but requires at least %u.\n",
               oo, ot->_inputs.size(), ot->_inputsMin);

    if (ot->_inputs.size() > ot->_inputsMax)
      sprintf(_errors, "ERROR: operation at position %u has %u inputs, but requires at most %u.\n",
               oo, ot->_inputs.size(), ot->_inputsMax);
  }

  //  Check COUNTING operations:
  //   - the mer size must be known
  //   - an output must be defined
  //   - inputs must be only 'sequence' types

  for (uint32 oo=0; oo<_opList.size(); oo++) {
    merylOpTemplate  *ot = _opList[oo];

    if (ot->_type != merylOpType::opCounting)
      continue;

    if (kmer::merSize() == 0)  //  ERROR: Kmer size not supplied with modifier k=<kmer-size>
      sprintf(_errors, "ERROR: counting operation at position %u doesn't have a mer size set.\n", oo);

    if (ot->_outDbseName == nullptr)
      sprintf(_errors, "ERROR: counting operation at position %u must have an output database.\n", oo);

    for (uint32 ii=0; ii<ot->_inputs.size(); ii++)
      if ((ot->_inputs[ii]->isFromSequence() == false) &&
          (ot->_inputs[ii]->isFromStore()    == false))
        sprintf(_errors, "ERROR: counting operation at position %u input %u must be from a sequence file or Canu seqStore.", oo, ii);
  }

  //  Check STATISTICS or HISTOGRAM operations:
  //   - print and output are allowed
  //   - inputs must not be 'sequence' type
  for (uint32 oo=0; oo<_opList.size(); oo++) {
    merylOpTemplate  *ot = _opList[oo];

    if ((ot->_type != merylOpType::opStatistics) &&
        (ot->_type != merylOpType::opHistogram))
      continue;

    for (uint32 ii=0; ii<ot->_inputs.size(); ii++)
      if ((ot->_inputs[ii]->isFromTemplate() == false) &&
          (ot->_inputs[ii]->isFromDatabase() == false))
        sprintf(_errors, "ERROR: stats/histo operation at position %u input %u must be from a meryl database or another operation.", oo, ii);
  }

  //  Check SELECTOR operations:
  //   - inputs must not be 'sequence' type

  for (uint32 oo=0; oo<_opList.size(); oo++) {
    merylOpTemplate  *ot = _opList[oo];

    if (ot->_type != merylOpType::opFilter)
      continue;

    for (uint32 ii=0; ii<ot->_inputs.size(); ii++)
      if ((ot->_inputs[ii]->isFromTemplate() == false) &&
          (ot->_inputs[ii]->isFromDatabase() == false))
        sprintf(_errors, "ERROR: select operation at position %u input %u must be from a meryl database or another operation.", oo, ii);
  }
}



void
merylCommandBuilder::performCounting(uint64 allowedMemory, uint32 allowedThreads) {
  for (uint32 oo=0; oo<numOperations(); oo++)
    getOperation(oo)->doCounting(allowedMemory, allowedThreads);
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

    //  Actions that were previously a count or a histo/stats on a database
    //  are all done and do not need to process their inputs.
    //
    //  Better - just remove the inputs from them?  But then we need to know
    //  if it supplies input to anything else.

#warning remove empty trees
    if (tpl->_type == merylOpType::opNothing)   //  Was previously a count, histo or stats
      continue;                                 //  operation, but it's all done now.

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
