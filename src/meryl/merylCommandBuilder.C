
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

#include <stdarg.h>




//  All the compute actions are finished and cleaned up.  The template action
//  need to harvest any summary statistics from them.
void
merylCommandBuilder::finishTree(uint32 rr) {
  uint32  tt = _opTree[rr];

  //  Harvest summary statistics from all the slices.

#warning collectStatistics()
  //for (uint32 ss=0; ss<64; ss++)
  //  _opList[tt]->collectStatistics(_thList[ss][tt]);

  //  Remove the slice computes.  This will close any open output database
  //  and delete the input objects which will propagate the delete/close to
  //  the rest of the objects in the tree.
  //
  //  Beware that the pointers in _thList[] are dangling now.
  //
  for (uint32 ss=0; ss<64; ss++)
    delete _thList[ss][tt];

  //  Emit statistics

  //if (_ot->_type == merylOpType::opHistogram)
  //  reportHistogram();
  //if (_ot->_type == merylOpType::opStatistics)
  //  reportStatistics();

  //_opList[tt]->emitStatistics();
  //_opList[tt]->emitHistogram();

  //  Delete the template.

  delete _opList[tt];
}




void
merylCommandBuilder::addError(char const *fmt, ...) {
  char    *err = new char [1024];
  va_list  ap;

  va_start(ap, fmt);
  vsnprintf(err, 1024, fmt, ap);
  va_end(ap);

  _errors.push_back(err);
}



//  This is called by processWord() when an '[' is encountered in the input.
//  It will _always_ push a new operation onto the stack.  Depending on the
//  stack contents, there are two variants:
//    empty - add the new op to the list of root operations
//    not   - add the new op to the list of inputs for the current top op
//
void
merylCommandBuilder::addNewOperation(void) {
  merylOpTemplate *op = new merylOpTemplate(_opList.size() + 1);

  if (_opStack.size() == 0) {
    if (verbosity.showConstruction() == true)
      fprintf(stderr, "addOperation()- Add new operation as root operation.\n");

    _opStack.push(op);
    _opList.push_back(op);
    _opTree.push_back(_opList.size() - 1);
  }

  else {
    if (verbosity.showConstruction() == true)
      fprintf(stderr, "addOperation()- Add new operation to stack at position %ld.\n",
              _opStack.size());

    getCurrent()->addInputFromOp(op, _errors);

    _opStack.push(op);
    _opList.push_back(op);
  }
}



bool
merylCommandBuilder::isCount(void) {
  merylOpTemplate  *op = getCurrent();

  if ((strcmp(_optString, "count")         != 0) &&
      (strcmp(_optString, "count-forward") != 0) &&
      (strcmp(_optString, "count-reverse") != 0))
    return(false);

  if (op->_type == merylOpType::opCounting)
    addError("operation is already a counting operation");

  op->_type     = merylOpType::opCounting;
  op->_counting = new merylOpCounting(_optString);

  return(true);
}



bool
merylCommandBuilder::isOutput(void) {

  //  If we see 'output' remember that the next arg is expected to be the
  //  database or filename.

  if (strcasecmp(_optString, "output") == 0) {
    _isOutput       = true;
    return(true);
  }

  //  If not expecting an output database in this word, return that we didn't
  //  consume the word.

  if (_isOutput == false)
    return(false);

  //  Otherwise, we're expecting an output database name.  Add it to the
  //  operation and return that we consumed the word.

  _isOutput = false;

  getCurrent()->addOutput(_inoutName, _errors);

  return(true);
}



bool
merylCommandBuilder::isPrinter(void) {
  merylOpTemplate  *op = getCurrent();

  //  If we see 'output' or 'print' or 'printacgt', remember that the next
  //  arg is expected to be the database or filename.

  if (strcasecmp(_optString, "print") == 0) {
    _printACGTorder = false;
    _isPrint        = true;
    return(true);
  }

  if (strcasecmp(_optString, "printACGT") == 0) {
    _printACGTorder = true;
    _isPrint        = true;
    return(true);
  }

  //  If not expecting a printer outpot name, return that we didn't
  //  consume the word.

  if (_isPrint == false)
    return(false);

  //  Otherwise, we're expecting a printer name.  This is complicated by
  //  wanting to allow 'print database.meryl' (to print the kmers to stdout),
  //  and 'print kmers.txt database.meryl' (to print to a file).
  //
  //  The former is handled by noticing that this word is a meryl database,
  //  adding the print option to the current operation and returning FALSE to
  //  let merylInput() catch the db and add it to the command as normal.
  //
  //  The latter simply adds the output name as a printer output.

  _printACGTorder = false;
  _isPrint        = false;

  if (fileExists(_indexName) == true) {
    op->addPrinter(nullptr, _printACGTorder, _errors);
    return(false);
  }

  op->addPrinter(_inoutName, _printACGTorder, _errors);

  return(true);
}



//  If this word is a meryl database, add it as an input to the current
//  operation.  If the current operation isn't defined yet, define it
//  to accept exactly one input.
//
//  This MUST come after isPrinter() so we can handle 'print some.meryl'
//  correctly.
bool
merylCommandBuilder::isInput(void) {
  merylOpTemplate  *op = getCurrent();

  bool   merylInput = ((fileExists(_indexName) == true));
  bool   canuInput  = ((fileExists(_sqInfName) == true) &&
                       (fileExists(_sqRdsName) == true));
  bool   seqInput   = ((fileExists(_inoutName) == true));

  //  If not an input, don't consume the word.
  if ((merylInput == false) &&
      (canuInput  == false) &&
      (seqInput   == false))
    return(false);

  //  Otherwise, it's an input.  If the op isn't set up yet, set it up to
  //  take exactly one input.  This _should_ be a print operation.
  //
  //if (op->_operation == merylOp::opNothing) {
  //}

  //  Now we can add the input to the operation.

  if (merylInput == true) {
    op->addInputFromDB(_inoutName, _errors);
  }

#ifdef CANU

  if (canuInput == true) {
    op->addInputFromCanu(_inoutName, _segment, _segmentMax, _errors);

    _segment    = 1;
    _segmentMax = 1;
  }

#else

  if (canuInput == true) {
    addError("ERROR: Detected Canu seqStore input '%s', but no Canu support is available.\n", _inoutName);
  }

#endif

  if (seqInput == true) {
    op->addInputFromSeq(_inoutName, _doCompression, _errors);
  }

  //  Reset state.

  return(true);
}


//  This doesn't really build the trees, since they're build as words are
//  added.
//
//  This does go through all the actions to make sure the inputs, outputs and
//  filters are appropriate for that action, generating error messagess if
//  needed.
//
void
merylCommandBuilder::buildTrees(void) {

  //  Clear the stack, we're done with it.
  while (_opStack.size() > 0)
    _opStack.pop();

  //  Tell each action that it has all the inputs it is going to get,
  //  and add errors if there are too many or too few of them.

  for (uint32 oo=0; oo<_opList.size(); oo++) {
    merylOpTemplate  *ot = _opList[oo];

    ot->finalizeTemplateInputs(_errors);

    if (ot->_inputs.size() < ot->_inputsMin)
      addError("operation at position %u has %u inputs, but requires at least %u\n",
               oo, ot->_inputs.size(), ot->_inputsMin);

    if (ot->_inputs.size() > ot->_inputsMax)
      addError("operation at position %u has %u inputs, but requires at most %u\n",
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
      addError("counting operation at position %u doesn't have a mer size set\n", oo);

    if (ot->_output == nullptr)
      addError("counting operation at position %u must have an output database\n", oo);

    for (uint32 ii=0; ii<ot->_inputs.size(); ii++)
      if ((ot->_inputs[ii]->isFromSequence() == false) &&
          (ot->_inputs[ii]->isFromStore()    == false))
        addError("counting operation at position %u input %u must be from a sequence file or Canu seqStore\n", oo, ii);
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
        addError("stats/histo operation at position %u input %u must be from a meryl database or another operation\n", oo, ii);
  }

  //  Check FILTER operations:
  //   - inputs must not be 'sequence' type

  for (uint32 oo=0; oo<_opList.size(); oo++) {
    merylOpTemplate  *ot = _opList[oo];

    if (ot->_type != merylOpType::opFilter)
      continue;

    for (uint32 ii=0; ii<ot->_inputs.size(); ii++)
      if ((ot->_inputs[ii]->isFromTemplate() == false) &&
          (ot->_inputs[ii]->isFromDatabase() == false))
        addError("filter operation at position %u input %u must be from a meryl database or another operation\n", oo, ii);
  }
}



void
merylCommandBuilder::performCounting(uint64 allowedMemory, uint32 allowedThreads) {
  for (uint32 oo=0; oo<numOperations(); oo++)
    getOperation(oo)->doCounting(allowedMemory, allowedThreads);
}





void
merylCommandBuilder::printTree(merylOpTemplate *op, uint32 inputNum, uint32 indent) {
  char   sA[1024];

  //  Scan the list of roots to decide if this is the root of a tree.

  bool   isTree  = false;
  uint32 rootNum = 0;

  for (uint32 rr=0; rr<numTrees(); rr++) {
    if (getTree(rr) == op) {
      isTree  = true;
      rootNum = rr;
    }
  }

  //  The action is either the root (and so the stream output from it is
  //  ignored) or is an input to some other action.

  if (isTree == true) {
    if (op->_type == merylOpType::opNothing)      fprintf(stderr, "%*s|- TREE %u: action #%u <empty>\n",           indent, "", rootNum, op->_ident);
    if (op->_type == merylOpType::opCounting)     fprintf(stderr, "%*s|- TREE %u: action #%u count kmers\n",       indent, "", rootNum, op->_ident);
    if (op->_type == merylOpType::opStatistics)   fprintf(stderr, "%*s|- TREE %u: action #%u report statistics\n", indent, "", rootNum, op->_ident);
    if (op->_type == merylOpType::opHistogram)    fprintf(stderr, "%*s|- TREE %u: action #%u report histogram\n",  indent, "", rootNum, op->_ident);
    if (op->_type == merylOpType::opFilter)       fprintf(stderr, "%*s|- TREE %u: action #%u filter kmers\n",      indent, "", rootNum, op->_ident);
  }
  else {
    if (op->_type == merylOpType::opNothing)      fprintf(stderr, "%*s<- INPUT @%u: action #%u <empty>\n",           indent, "", inputNum, op->_ident);
    if (op->_type == merylOpType::opCounting)     fprintf(stderr, "%*s<- INPUT @%u: action #%u count kmers\n",       indent, "", inputNum, op->_ident);
    if (op->_type == merylOpType::opStatistics)   fprintf(stderr, "%*s<- INPUT @%u: action #%u report statistics\n", indent, "", inputNum, op->_ident);
    if (op->_type == merylOpType::opHistogram)    fprintf(stderr, "%*s<- INPUT @%u: action #%u report histogram\n",  indent, "", inputNum, op->_ident);
    if (op->_type == merylOpType::opFilter)       fprintf(stderr, "%*s<- INPUT @%u: action #%u filter kmers\n",      indent, "", inputNum, op->_ident);
  }

  //  SELECTION
  //

  switch (op->_valueSelect) {
    case merylModifyValue::valueNOP:
      break;
    case merylModifyValue::valueSet:
      fprintf(stderr, "%*s|- SET value to constant %u\n", indent+3, "", op->_valueConstant);
      break;
    case merylModifyValue::valueSelected:
      fprintf(stderr, "%*s|- SET value to that of the kmer selected by label filter\n", indent+3, "");
      break;
    case merylModifyValue::valueFirst:
      fprintf(stderr, "%*s|- SET value to that of the kmer in the first input\n", indent+3, "");
      break;
    case merylModifyValue::valueMin:
      fprintf(stderr, "%*s|- SET value to the minimum of all kmers and constant %u\n", indent+3, "", op->_valueConstant);
      break;
    case merylModifyValue::valueMax:
      fprintf(stderr, "%*s|- SET value to the maximum of all kmers and constant %u\n", indent+3, "", op->_valueConstant);
      break;
    case merylModifyValue::valueAdd:
      fprintf(stderr, "%*s|- SET value to the sum of all kmers and constant %u\n", indent+3, "", op->_valueConstant);
      break;
    case merylModifyValue::valueSub:
      fprintf(stderr, "%*s|- SET value to the selected kmer minus all others and constant %u\n", indent+3, "", op->_valueConstant);
      break;
    case merylModifyValue::valueMul:
      fprintf(stderr, "%*s|- SET value to the product of all kmers and constant %u\n", indent+3, "", op->_valueConstant);
      break;
    case merylModifyValue::valueDiv:
      fprintf(stderr, "%*s|- SET value to the selected kmer divided by all others and constant %u\n", indent+3, "", op->_valueConstant);
      break;
    case merylModifyValue::valueDivZ:
      fprintf(stderr, "%*s|- SET value to the selected kmer divided by all others and constant %u\n", indent+3, "", op->_valueConstant);
      fprintf(stderr, "%*s    (rounding zero up to one)\n", indent+3, "");
      break;
    case merylModifyValue::valueMod:
      fprintf(stderr, "%*s|- SET value to the remainder of the selected kmer divided by all others and constant %u\n", indent+3, "", op->_valueConstant);
      break;
    case merylModifyValue::valueCount:
      fprintf(stderr, "%*s|- SET value to the number of input databases the kmer is present in\n", indent+3, "");
      break;
    default:
      break;
  }

  switch (op->_labelSelect) {
    case merylModifyLabel::labelNOP:
      break;
    case merylModifyLabel::labelSet:
      fprintf(stderr, "%*s|- SET label to constant %lu\n", indent+3, "", op->_labelConstant);
      break;
    case merylModifyLabel::labelSelected:
      fprintf(stderr, "%*s|- SET label to selected -- constant %lu\n", indent+3, "", op->_labelConstant);
      break;
    case merylModifyLabel::labelFirst:
      fprintf(stderr, "%*s|- SET label to first -- constant %lu\n", indent+3, "", op->_labelConstant);
      break;
    case merylModifyLabel::labelMin:
      fprintf(stderr, "%*s|- SET label to min -- constant %lu\n", indent+3, "", op->_labelConstant);
      break;
    case merylModifyLabel::labelMax:
      fprintf(stderr, "%*s|- SET label to max -- constant %lu\n", indent+3, "", op->_labelConstant);
      break;
    case merylModifyLabel::labelAnd:
      fprintf(stderr, "%*s|- SET label to and -- constant %lu\n", indent+3, "", op->_labelConstant);
      break;
    case merylModifyLabel::labelOr:
      fprintf(stderr, "%*s|- SET label to or -- constant %lu\n", indent+3, "", op->_labelConstant);
      break;
    case merylModifyLabel::labelXor:
      fprintf(stderr, "%*s|- SET label to xor -- constant %lu\n", indent+3, "", op->_labelConstant);
      break;
    case merylModifyLabel::labelDifference:
      fprintf(stderr, "%*s|- SET label to difference -- constant %lu\n", indent+3, "", op->_labelConstant);
      break;
    case merylModifyLabel::labelLightest:
      fprintf(stderr, "%*s|- SET label to lightest -- constant %lu\n", indent+3, "", op->_labelConstant);
      break;
    case merylModifyLabel::labelHeaviest:
      fprintf(stderr, "%*s|- SET label to heaviest -- constant %lu\n", indent+3, "", op->_labelConstant);
      break;
    case merylModifyLabel::labelInvert:
      fprintf(stderr, "%*s|- SET label to invert -- constant %lu\n", indent+3, "", op->_labelConstant);
      break;
    case merylModifyLabel::labelShiftLeft:
      fprintf(stderr, "%*s|- SET label to shift left -- constant %lu\n", indent+3, "", op->_labelConstant);
      break;
    case merylModifyLabel::labelShiftRight:
      fprintf(stderr, "%*s|- SET label to shift right -- constant %lu\n", indent+3, "", op->_labelConstant);
      break;
    case merylModifyLabel::labelRotateLeft:
      fprintf(stderr, "%*s|- SET label to rotate left -- constant %lu\n", indent+3, "", op->_labelConstant);
      break;
    case merylModifyLabel::labelRotateRight:
      fprintf(stderr, "%*s|- SET label to rotate right -- constant %lu\n", indent+3, "", op->_labelConstant);
      break;
    default:
      break;
  }



  //  Report filters

  for (uint32 ii=0; ii<op->_filter    .size(); ii++) {
    fprintf(stderr, "%*s|- FILTER %u\n", indent+3, "", ii+1);
    for (uint32 jj=0; jj<op->_filter[ii].size(); jj++) {
      fprintf(stderr, "%*s|- %s", indent+6, "", op->_filter[ii][jj].describe(sA));
    }

    //if (ii + 1 < op->_filter.size())
    //  fprintf(stderr, "%*sor\n", indent+4, "");
  }



#if 0
  if (op->_presentInFirst == true)
    fprintf(stderr, "%*s|-kmer must be present in first input\n", indent+3, "");

  if (op->_presentInAll == true)
    fprintf(stderr, "%*s|-kmer has no database presence requirements\n", indent+3, "");

  if (op->_presentInAll == true)
    fprintf(stderr, "%*s|-kmer must be present in all inputs\n", indent+3, "");

  for (uint32 pp=0; pp<op->_inputs.size(); pp++)
    if (op->_presentIn[pp] == true)
      fprintf(stderr, "%*s|-kmer must be present in %u inputs\n", indent+3, "", pp);
#endif

#if 0
  if (op->_mathConstant > 0)
    fprintf(stderr, "%*sconstant=%lu\n", indent+3, "", op->_mathConstant);

  if (op->_threshold != UINT64_MAX)
    fprintf(stderr, "%*sthreshold=%lu\n", indent+3, "", op->_threshold);

  if (op->_fracDist != DBL_MAX)
    fprintf(stderr, "%*sfraction-distinct=%f\n", indent+3, "", op->_fracDist);

  if (op->_wordFreq != DBL_MAX)
    fprintf(stderr, "%*sword-frequenct=%f\n", indent+3, "", op->_wordFreq);
#endif


  //  Report inputs

  for (uint32 ii=0; ii<op->_inputs.size(); ii++) {
    merylInput  *in = op->_inputs[ii];

    if (in->isFromTemplate() == true) {
      fprintf(stderr, "%*s<- INPUT @%u: action #%u\n", indent+3, "", ii+1, in->_template->_ident);
      printTree(in->_template, ii+1, indent+3);
    }

    if (in->isFromDatabase() == true) {
      fprintf(stderr, "%*s<- INPUT @%u: meryl database '%s'\n", indent+3, "", ii+1, in->name());
    }

    if (in->isFromSequence() == true) {
      fprintf(stderr, "%*s<- INPUT @%u: sequence file '%s'%s\n", indent+3, "", ii+1, in->name(), in->_homopolyCompress ? " (homopoly compressed)" : "");
    }

    if (in->isFromStore() == true) {
      fprintf(stderr, "%*s<- INPUT @%u: Canu sqStore '%s' (using reads %u through %u)\n", indent+3, "", ii+1, in->name(), in->_sqBgn, in->_sqEnd);
    }
  }

  //  Report outputs

  if (op->_output != nullptr) {
    fprintf(stderr, "%*s-> database '%s'\n", indent+3, "", op->_output->filename());
  }

  //  Report printing

  if (op->_printerName != nullptr) {
    fprintf(stderr, "%*s-> text file '%s'\n", indent+3, "", op->_printerName);
  }

  if (isTree == true) {
    fprintf(stderr, "TREE %u ends.\n", rootNum);
    fprintf(stderr, "\n");
  }
}



//  Clone the command tree(s) into thread-specific copies, one tree per thread.
//
//
void
merylCommandBuilder::spawnThreads(uint32 allowedThreads) {

  //  Tell all the actions to do their final initialization - this will
  //  create master input/output objects, query histograms to get any
  //  paramters based on those, and anything else that needs to be done once
  //  per action.

  for (uint32 rr=0; rr<numTrees(); rr++)
    getTree(rr)->finalizeTemplateParameters();

  //  Allocate compute objects for each of our 64 file slices, then copy the
  //  list of templates into each slice.  These need to exist before we start
  //  creating inputs.

  for (uint32 ss=0; ss<64; ss++) {
    _thList[ss] = new merylOpCompute * [_opList.size()];

    for (uint32 oo=0; oo<_opList.size(); oo++)
      _thList[ss][oo] = new merylOpCompute(_opList[oo], ss, _opList[oo]->_inputs.size());
  }

  //  Update all the input/output objects to be per-thread.

  for (uint32 ss=0; ss<64; ss++) {
    for (uint32 oo=0; oo<_opList.size(); oo++) {
      merylOpTemplate  *tpl = _opList[oo];       //  The template operation
      merylOpCompute   *cpu = _thList[ss][oo];   //  The per-thread operation we're creating.

      //  For each input to the template, add a new input to the compute.
      //
      for (uint32 ii=0; ii<tpl->_inputs.size(); ii++) {
        merylInput  *in = tpl->_inputs[ii];

        //  If the template input is from an action, we need to search for
        //  that action in the master list of templates, then set the
        //  per-thread input to be the corresponding object in the per-thread
        //  list.
        //

        if (in->isFromTemplate() == true) {
          uint32  inidx = UINT32_MAX;

          for (uint32 xx=0; xx<_opList.size(); xx++)
            if (in->_template == _opList[xx])
              inidx = xx;

          if (inidx == UINT32_MAX)
            fprintf(stderr, "Failed to find corresponding operation.\n"), exit(1);

          fprintf(stderr, "addInputFromOp oo %d ss %d inidx %d\n", oo, ss, inidx);
          cpu->addInputFromOp(_thList[ss][inidx]);
        }

        //  If the template input is from a database, make a new input for
        //  just the piece we're processing in this thread - this is done in
        //  addInputFromDB().
        //
        if (in->isFromDatabase() == true) {
          cpu->addInputFromDB(in->name(), ss);
        }

        //  If the template input is from anything else, it's an error.
        //
        assert(in->isFromSequence() == false);
        assert(in->isFromStore()    == false);
      }

      //  Add a writer or printer (if needed) for the slice we're operating on.
      //
      cpu->addOutput(tpl, ss);
      cpu->addPrinter(tpl, ss);
    }
  }
}

