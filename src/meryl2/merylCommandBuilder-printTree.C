
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
//  Print a full command tree with indent 'indent'.
//
void
merylCommandBuilder::printTree(uint32 tID, uint32 indent) {
  printTree(getTree(tID), tID+1, 0, indent);
}



//
//  Print a single tree node, recursing into child nodes.
//    op     - the tree node to print
//    tID    - integer index representing which tree is being printed
//    iID    - integer index representing which input this node is to some other node
//    indent - amount to indent the display
//
//  tID and iID are mutually exclusive.
//    if tID > 0 - this node is the root of a tree, and its output goes nowhere.
//    if iID > 0 - this node is the input to some other node.
//

void
merylCommandBuilder::printTree(merylOpTemplate *op, uint32 tID, uint32 iID, uint32 indent) {
  char   sA[1024];

  //
  //  Report where our kmers go.
  //

  {
    bool ist       = (tID > 0);       //  Display where our 'stream output' goes.
    char const *to = "|-+ TREE ";      //    For trees, it goes nowhere.
    char const *no = "|<- INPUT @";    //    For nodes, it is passed up to the parent node.

    switch (op->_type) {
      case merylOpType::opNothing:      fprintf(stderr, "%*s%s%u: action #%u <empty>\n",           indent, "", (ist) ? to : no, (ist) ? tID : iID, op->_ident);   break;
      case merylOpType::opCounting:     fprintf(stderr, "%*s%s%u: action #%u count kmers\n",       indent, "", (ist) ? to : no, (ist) ? tID : iID, op->_ident);   break;
      case merylOpType::opStatistics:   fprintf(stderr, "%*s%s%u: action #%u report statistics\n", indent, "", (ist) ? to : no, (ist) ? tID : iID, op->_ident);   break;
      case merylOpType::opHistogram:    fprintf(stderr, "%*s%s%u: action #%u report histogram\n",  indent, "", (ist) ? to : no, (ist) ? tID : iID, op->_ident);   break;
      case merylOpType::opPrint:        fprintf(stderr, "%*s%s%u: action #%u print kmers\n",       indent, "", (ist) ? to : no, (ist) ? tID : iID, op->_ident);   break;
      case merylOpType::opFilter:       fprintf(stderr, "%*s%s%u: action #%u select kmers\n",      indent, "", (ist) ? to : no, (ist) ? tID : iID, op->_ident);   break;
    }
  }

  //
  //  Report all OUTPUTS.
  //

  if (op->_outDbseName  != nullptr)   fprintf(stderr, "%*s|-> meryl database '%s'\n",    indent+3, "", op->_outDbseName);
  if (op->_outListName  != nullptr)   fprintf(stderr, "%*s|-> text list '%s'%s\n",       indent+3, "", op->_outListName, (op->_outList == nullptr) ? "" : " (parallel mode)");
  if (op->_outShowName  != nullptr)   fprintf(stderr, "%*s|-> text list (STDOUT)\n",     indent+3, "");
  if (op->_outPipeName  != nullptr)   fprintf(stderr, "%*s|-> meryl pipe '%s'\n",        indent+3, "", op->_outPipeName);
  if (op->_outStatsName != nullptr)   fprintf(stderr, "%*s|-> statistics report '%s'\n", indent+3, "", op->_outStatsName);
  if (op->_outHistoName != nullptr)   fprintf(stderr, "%*s|-> histogram report '%s'\n",  indent+3, "", op->_outHistoName);

  //
  //  Report kmer/value/label assignment.
  //

  switch (op->_valueAssign) {
    case merylAssignValue::valueNOP:
      break;
    case merylAssignValue::valueSet:
      fprintf(stderr, "%*s|-- SET value to constant %u\n", indent+3, "", op->_valueConstant);
      break;
    case merylAssignValue::valueSelected:
      fprintf(stderr, "%*s|-- SET value to that of the kmer selected by label selector\n", indent+3, "");
      break;
    case merylAssignValue::valueFirst:
      fprintf(stderr, "%*s|-- SET value to that of the kmer in the first input\n", indent+3, "");
      break;
    case merylAssignValue::valueMin:
      fprintf(stderr, "%*s|-- SET value to the minimum of all kmers and constant %u\n", indent+3, "", op->_valueConstant);
      break;
    case merylAssignValue::valueMax:
      fprintf(stderr, "%*s|-- SET value to the maximum of all kmers and constant %u\n", indent+3, "", op->_valueConstant);
      break;
    case merylAssignValue::valueAdd:
      fprintf(stderr, "%*s|-- SET value to the sum of all kmers and constant %u\n", indent+3, "", op->_valueConstant);
      break;
    case merylAssignValue::valueSub:
      fprintf(stderr, "%*s|-- SET value to the selected kmer minus all others and constant %u\n", indent+3, "", op->_valueConstant);
      break;
    case merylAssignValue::valueMul:
      fprintf(stderr, "%*s|-- SET value to the product of all kmers and constant %u\n", indent+3, "", op->_valueConstant);
      break;
    case merylAssignValue::valueDiv:
      fprintf(stderr, "%*s|-- SET value to the selected kmer divided by all others and constant %u\n", indent+3, "", op->_valueConstant);
      break;
    case merylAssignValue::valueDivZ:
      fprintf(stderr, "%*s|-- SET value to the selected kmer divided by all others and constant %u\n", indent+3, "", op->_valueConstant);
      fprintf(stderr, "%*s    (rounding zero up to one)\n", indent+3, "");
      break;
    case merylAssignValue::valueMod:
      fprintf(stderr, "%*s|-- SET value to the remainder of the selected kmer divided by all others and constant %u\n", indent+3, "", op->_valueConstant);
      break;
    case merylAssignValue::valueCount:
      fprintf(stderr, "%*s|-- SET value to the number of input databases the kmer is present in\n", indent+3, "");
      break;
    default:
      break;
  }

  switch (op->_labelAssign) {
    case merylAssignLabel::labelNOP:
      break;
    case merylAssignLabel::labelSet:
      fprintf(stderr, "%*s|-- SET label to constant %lu\n", indent+3, "", op->_labelConstant);
      break;
    case merylAssignLabel::labelSelected:
      fprintf(stderr, "%*s|-- SET label to selected -- constant %lu\n", indent+3, "", op->_labelConstant);
      break;
    case merylAssignLabel::labelFirst:
      fprintf(stderr, "%*s|-- SET label to first -- constant %lu\n", indent+3, "", op->_labelConstant);
      break;
    case merylAssignLabel::labelMin:
      fprintf(stderr, "%*s|-- SET label to min -- constant %lu\n", indent+3, "", op->_labelConstant);
      break;
    case merylAssignLabel::labelMax:
      fprintf(stderr, "%*s|-- SET label to max -- constant %lu\n", indent+3, "", op->_labelConstant);
      break;
    case merylAssignLabel::labelAnd:
      fprintf(stderr, "%*s|-- SET label to and -- constant %lu\n", indent+3, "", op->_labelConstant);
      break;
    case merylAssignLabel::labelOr:
      fprintf(stderr, "%*s|-- SET label to or -- constant %lu\n", indent+3, "", op->_labelConstant);
      break;
    case merylAssignLabel::labelXor:
      fprintf(stderr, "%*s|-- SET label to xor -- constant %lu\n", indent+3, "", op->_labelConstant);
      break;
    case merylAssignLabel::labelDifference:
      fprintf(stderr, "%*s|-- SET label to difference -- constant %lu\n", indent+3, "", op->_labelConstant);
      break;
    case merylAssignLabel::labelLightest:
      fprintf(stderr, "%*s|-- SET label to lightest -- constant %lu\n", indent+3, "", op->_labelConstant);
      break;
    case merylAssignLabel::labelHeaviest:
      fprintf(stderr, "%*s|-- SET label to heaviest -- constant %lu\n", indent+3, "", op->_labelConstant);
      break;
    case merylAssignLabel::labelInvert:
      fprintf(stderr, "%*s|-- SET label to invert -- constant %lu\n", indent+3, "", op->_labelConstant);
      break;
    case merylAssignLabel::labelShiftLeft:
      fprintf(stderr, "%*s|-- SET label to shift left -- constant %lu\n", indent+3, "", op->_labelConstant);
      break;
    case merylAssignLabel::labelShiftRight:
      fprintf(stderr, "%*s|-- SET label to shift right -- constant %lu\n", indent+3, "", op->_labelConstant);
      break;
    case merylAssignLabel::labelRotateLeft:
      fprintf(stderr, "%*s|-- SET label to rotate left -- constant %lu\n", indent+3, "", op->_labelConstant);
      break;
    case merylAssignLabel::labelRotateRight:
      fprintf(stderr, "%*s|-- SET label to rotate right -- constant %lu\n", indent+3, "", op->_labelConstant);
      break;
    default:
      break;
  }

  //
  //  Report selectors
  //

  for (uint32 ii=0; ii<op->_select.size(); ii++) {
    fprintf(stderr, "%*s|-- SELECTOR %u\n", indent+3, "", ii+1);
    for (uint32 jj=0; jj<op->_select[ii].size(); jj++) {
      fprintf(stderr, "%*s|-- %s", indent+6, "", op->_select[ii][jj].describe(sA));
    }

    //if (ii + 1 < op->_select.size())
    //  fprintf(stderr, "%*sor\n", indent+4, "");
  }

  //
  //  Report inputs
  //

  for (uint32 ii=0; ii<op->_inputs.size(); ii++) {
    merylInput  *in = op->_inputs[ii];

    if (in->isFromTemplate() == true) {
      //  The following line is printed by the recursive call to printTree().
      //fprintf(stderr, "%*s<- INPUT @%u: action #%u\n", indent+3, "", ii+1, in->_template->_ident);
      printTree(in->_template, 0, ii+1, indent+3);
    }

    if (in->isFromDatabase() == true) {
      fprintf(stderr, "%*s|<- INPUT @%u: meryl database '%s'\n", indent+3, "", ii+1, in->inputName());
    }

    if (in->isFromSequence() == true) {
      fprintf(stderr, "%*s|<- INPUT @%u: sequence file '%s'%s\n", indent+3, "", ii+1, in->inputName(), in->_squish ? " (homopoly compressed)" : "");
    }

    if (in->isFromStore() == true) {
      fprintf(stderr, "%*s|<- INPUT @%u: Canu sqStore '%s' (using reads %u through %u)\n", indent+3, "", ii+1, in->inputName(), in->_sqBgn, in->_sqEnd);
    }
  }


#if 0
  if (op->_presentInFirst == true)
    fprintf(stderr, "%*s|--kmer must be present in first input\n", indent+3, "");

  if (op->_presentInAll == true)
    fprintf(stderr, "%*s|--kmer has no database presence requirements\n", indent+3, "");

  if (op->_presentInAll == true)
    fprintf(stderr, "%*s|--kmer must be present in all inputs\n", indent+3, "");

  for (uint32 pp=0; pp<op->_inputs.size(); pp++)
    if (op->_presentIn[pp] == true)
      fprintf(stderr, "%*s|--kmer must be present in %u inputs\n", indent+3, "", pp);
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

  //
  //  If a tree, flag the end.
  //

  if (tID > 0) {
    fprintf(stderr, "%*s|-- TREE %u ends.\n", indent, "", tID);
    fprintf(stderr, "\n");
  }
}
