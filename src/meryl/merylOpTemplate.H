
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

#ifndef MERYLOPTEMPLATE_H
#define MERYLOPTEMPLATE_H

#include "meryl.H"

#include <vector>

//  This class is used for setting up the processing tree, checking args,
//  keeping track of parameters, etc.  Objects are made/modified while
//  the command line is being parsed.  Once parsing is done, they are
//  used to make merylOpCompute objects for each thread.
//
//  Counting operations are also done here.

//  The filters are a sum-of-products expression returning TRUE if the kmer should
//  be output and FALSE if the kmer should be ignored.
//
//  isKmerFilteredOut() computes the sum of those product filters applied
//  to the kmer we are thinking of outputting (_kmer).
//
//  Ideally, this should be captured in its own object, but it needs access
//  to the _act lists.
//
//  isTrue(k) returns true if the output kmer 'k' meets the conditions of
//  the filter and should be passed to the output.  It has access to all
//  the other input kmers through the _act inputs.
//
//  addNewFilterTerm() will add a new product term to the sum list.
//  addNewFilter()     will add a new filter to the most recent product term.

class merylOpTemplate {
public:
public:
  merylOpTemplate(uint32 ident);
  ~merylOpTemplate();

public:
  void    addInputFromOp  (merylOpTemplate *operation, std::vector<char const *> &err);
  void    addInputFromDB  (char const *dbName, std::vector<char const *> &err);
private:
  void    addInputFromDB  (char const *dbName);
public:
  void    addInputFromCanu(char const *stName, uint32 segment, uint32 segmentMax, std::vector<char const *> &err);
  void    addInputFromSeq (char const *sqName, bool doCompression, std::vector<char const *> &err);

  void    addOutput(char const *wrName, std::vector<char const *> &err);
  void    addPrinter(char const *prName, bool ACGTorder, std::vector<char const *> &err);


  void    finalizeTemplateInputs(std::vector<char const *> &err);
  void    finalizeTemplateParameters(void);

  //  Counting operations.
  void    doCounting(uint64 allowedMemory,
                     uint32 allowedThreads);

  merylOpCounting               *_counting = nullptr;

  //  Type of this operation.
  merylOpType                    _type      = merylOpType::opNothing;
  uint32                         _ident     =  0;
  //char                           _label[32] = {0};

  bool                           _onlyConfig = false;

  //  Inputs and a minimum/maximum number of inputs allowed.
  std::vector<merylInput *>      _inputs;

  uint32                         _inputsMin = 0;
  uint32                         _inputsMax = uint32max;

  //  Outputs.
  merylHistogram                *_stats   = nullptr;    //  For collecting stats on input kmers
  merylFileWriter               *_output  = nullptr;    //  Main output object

  compressedFileWriter          *_printer        = nullptr;
  char                          *_printerName    = nullptr;
  bool                           _printACGTorder = false;

  //  A description of how to generate output values/labels.
  merylModifyValue               _valueSelect   = merylModifyValue::valueNOP;
  kmvalu                         _valueConstant = 0;

  merylModifyLabel               _labelSelect   = merylModifyLabel::labelNOP;
  kmlabl                         _labelConstant = 0;

  //  Filters.
  void   addNewFilterProduct(void);
  void   addFilterToProduct(merylFilter const f);

  //  

private:
  std::vector< std::vector<merylFilter> >  _filter;

  //  Only so it can access _filter when displaying the command tree.
  friend class merylCommandBuilder;
};






inline
void
merylOpTemplate::addNewFilterProduct(void) {
  _filter.emplace_back();
}



inline
void
merylOpTemplate::addFilterToProduct(merylFilter const f) {

  _type = merylOpType::opFilter;   //  This op processes kmers.

  if (_filter.size() == 0)         //  If there are no filters yet,
    addNewFilterProduct();         //  add a new product term.

  _filter.back().emplace_back(f);  //  Add 'f' to the last product term.
}




#endif  //  MERYLOPTEMPLATE_H
