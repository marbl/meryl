
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
merylCommandBuilder::printTree(merylOpTemplate *op, uint32 inputNum, uint32 indent) {
  char   sA[1024];

  //
  //  Scan the list of roots to decide if this is the root of a tree.
  //
  //  If it is a root, the stream output is ignored;
  //  if it not, the stream output is the input to some other action.
  //

  bool   isTree  = false;
  uint32 rootNum = 0;

  for (uint32 rr=0; rr<numTrees(); rr++) {
    if (getTree(rr) == op) {
      isTree  = true;
      rootNum = rr;
    }
  }

  if (isTree == true) {
    switch (op->_type) {
      case merylOpType::opNothing:      fprintf(stderr, "\n%*s|- TREE %u: action #%u <empty>\n",           indent, "", rootNum, op->_ident);   break;
      case merylOpType::opCounting:     fprintf(stderr, "\n%*s|- TREE %u: action #%u count kmers\n",       indent, "", rootNum, op->_ident);   break;
      case merylOpType::opStatistics:   fprintf(stderr, "\n%*s|- TREE %u: action #%u report statistics\n", indent, "", rootNum, op->_ident);   break;
      case merylOpType::opHistogram:    fprintf(stderr, "\n%*s|- TREE %u: action #%u report histogram\n",  indent, "", rootNum, op->_ident);   break;
      case merylOpType::opPrint:        fprintf(stderr, "\n%*s|- TREE %u: action #%u print kmers\n",       indent, "", rootNum, op->_ident);   break;
      case merylOpType::opFilter:       fprintf(stderr, "\n%*s|- TREE %u: action #%u filter kmers\n",      indent, "", rootNum, op->_ident);   break;
    }
  }
  else {
    switch (op->_type) {
      case merylOpType::opNothing:      fprintf(stderr, "%*s^- INPUT @%u: action #%u <empty>\n",           indent, "", inputNum, op->_ident);   break;
      case merylOpType::opCounting:     fprintf(stderr, "%*s^- INPUT @%u: action #%u count kmers\n",       indent, "", inputNum, op->_ident);   break;
      case merylOpType::opStatistics:   fprintf(stderr, "%*s^- INPUT @%u: action #%u report statistics\n", indent, "", inputNum, op->_ident);   break;
      case merylOpType::opHistogram:    fprintf(stderr, "%*s^- INPUT @%u: action #%u report histogram\n",  indent, "", inputNum, op->_ident);   break;
      case merylOpType::opPrint:        fprintf(stderr, "%*s^- INPUT @%u: action #%u print kmers\n",       indent, "", inputNum, op->_ident);   break;
      case merylOpType::opFilter:       fprintf(stderr, "%*s^- INPUT @%u: action #%u filter kmers\n",      indent, "", inputNum, op->_ident);   break;
    }
  }

  //
  //  Report database, printing, histogram and statistics outputs.
  //

  if (op->_writer != nullptr) {
    fprintf(stderr, "%*s|> database '%s'\n", indent+3, "", op->_writer->filename());
  }

  if (op->_printerName != nullptr) {
    fprintf(stderr, "%*s|> text file '%s'\n", indent+3, "", op->_printerName);
  }

  if (op->_statsFile != nullptr) {
    fprintf(stderr, "%*s|> statistics to file '%s'\n", indent+3, "", op->_statsFile->filename());
  }

  if (op->_histoFile != nullptr) {
    fprintf(stderr, "%*s|> histogram to file '%s'\n", indent+3, "", op->_histoFile->filename());
  }

  //
  //  Report kmer/value/label selection.
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

  //
  //  Report filters
  //

  for (uint32 ii=0; ii<op->_filter    .size(); ii++) {
    fprintf(stderr, "%*s|- FILTER %u\n", indent+3, "", ii+1);
    for (uint32 jj=0; jj<op->_filter[ii].size(); jj++) {
      fprintf(stderr, "%*s|- %s", indent+6, "", op->_filter[ii][jj].describe(sA));
    }

    //if (ii + 1 < op->_filter.size())
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
      printTree(in->_template, ii+1, indent+3);
    }

    if (in->isFromDatabase() == true) {
      fprintf(stderr, "%*s^- INPUT @%u: meryl database '%s'\n", indent+3, "", ii+1, in->name());
    }

    if (in->isFromSequence() == true) {
      fprintf(stderr, "%*s^- INPUT @%u: sequence file '%s'%s\n", indent+3, "", ii+1, in->name(), in->_homopolyCompress ? " (homopoly compressed)" : "");
    }

    if (in->isFromStore() == true) {
      fprintf(stderr, "%*s^- INPUT @%u: Canu sqStore '%s' (using reads %u through %u)\n", indent+3, "", ii+1, in->name(), in->_sqBgn, in->_sqEnd);
    }
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

  //
  //  If a tree, flag the end.
  //

  if (isTree == true) {
    fprintf(stderr, "%*s|- TREE %u ends.\n", indent, "", rootNum);
    fprintf(stderr, "\n");
  }
}
