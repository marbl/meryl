
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

#ifndef MERYLFILTER_H
#define MERYLFILTER_H

#ifndef MERYLINCLUDE
#error "Do not use merylFilter.H, use meryl.H instead."
#endif

//  Relationships.  Tells how we should compare a kmer against a constant.
//
enum class merylFilterRelation {
  isNOP,         //  Nothing.  Generates an error.

  isEq,          //  True if the variable == constant
  isNeq,         //  True if the variable != constant
  isLeq,         //  True if the variable <= constant
  isGeq,         //  True if the variable >= constant
  isLt,          //  True if the variable <  constant
  isGt,          //  True if the variable >  constant
};


//  Quantities.  Tells what property of the kmer we're testing.
//
enum class merylFilterQuantity {
  isNOP,         //  Nothing.  Generates an error.

  isValue,       //  Compare kmer value
  isLabel,       //  Compare kmer label
  isBases,       //  Compare kmer bases
  isIndex,       //  Compare file index
};



//  A single filter; a term in the product expression.
//
//  The filter is true if 'kmer quantity' 'relation' 'constant',
//  for example "kmerValue == 4"
//
//  Filters can be combined in sum-of-product expressions.
//
//  The complicated part is countNonZeroBases().  This does some bit fiddling
//  to count the number of A's, C's, G's or T's in a kmer.
//
//  The basic idea is to count the number of bases that are represented by
//  to unset bits.  This is done by first squashing the two bits together,
//  masking out the unsquashed bit, then counting the number of set bits.
//
//    if the base is     A, '00'            will be squashed to '0'
//    if the base is not A, '01', '10', 11' are all squashed to '1'.
//
//  For other bases, we 'convert' that base into an A then count A's.
//    A = 00 - xor with 00  ->  A == 00  C == 01  T == 10  G == 11
//    C = 01 - xor with 01  ->  A == 01  C == 00  T == 11  G == 10
//    T = 10 - xor with 10  ->  A == 10  C == 11  T == 00  G == 01
//    G = 11 - xor with 11  ->  A == 11  C == 10  T == 01  G == 00
//
class merylFilter {
public:
  merylFilter(merylFilterQuantity type,
              merylFilterRelation rela,
              bool                invert,
              char const         *str);
  ~merylFilter();

private:
  template<typename X>
  bool
  compare(X x, X y) const {
    switch (_r) {
      case merylFilterRelation::isNOP:   assert(0);        break;
      case merylFilterRelation::isEq:    return(x == y);   break;
      case merylFilterRelation::isNeq:   return(x != y);   break;
      case merylFilterRelation::isLeq:   return(x <= y);   break;
      case merylFilterRelation::isGeq:   return(x >= y);   break;
      case merylFilterRelation::isLt:    return(x <  y);   break;
      case merylFilterRelation::isGt:    return(x >  y);   break;
      default:                           return(false);    break;
    }
    return(false);
  }

private:
  uint32 countNonZeroBases(kmdata k, uint64 cvt) const {   //  See above for docs.
    uint32   c = kmer::merSize();
    uint64  bh = k >> 64;                 //  Save the right-most 32 bases
    uint64  bl = k;                       //  Save the  left-most 32 bases

    bh ^= cvt;                            //  Apply any conversion from
    bl ^= cvt;                            //  base encoding to '00'

    bh |= (bh << 1);                      //  Combine two bits representing
    bl |= (bl << 1);                      //  a single base into one bit

    bh &= 0xaaaaaaaaaaaaaaaallu;          //  Strip out the non-combined bits.
    bl &= 0xaaaaaaaaaaaaaaaallu;          //

    return(countNumberOfSetBits64(bh) +   //  Sun the number of bases that are
           countNumberOfSetBits64(bl));   //  '01', '10' or '11'.
  }

private:
  uint32 countA(kmer k) const { return(kmer::merSize() - countNonZeroBases(k, 0x0000000000000000llu)); }
  uint32 countC(kmer k) const { return(kmer::merSize() - countNonZeroBases(k, 0x5555555555555555llu)); }
  uint32 countT(kmer k) const { return(kmer::merSize() - countNonZeroBases(k, 0xaaaaaaaaaaaaaaaallu)); }
  uint32 countG(kmer k) const { return(kmer::merSize() - countNonZeroBases(k, 0xffffffffffffffffllu)); }

  //  Evaluate the filter on kmer k, comparing against the constants saved in
  //  this filter object and the other kmer instances in the 'act' list.
public:
  bool
  isTrue(kmer k, uint32 actLen, kmer *act, uint32 *actIdx, uint32 *actRdx) const;


  //  finalizeFilterInputs() is used to let the filter decide if its
  //  configuration is sane, and to emit errors if needed.  This is also
  //  where we can figure out what 'all' means, or to fail if the filter
  //  specified a reference to file @4 but only 3 files exist.
  //
  //  finalizeFilterParameters() is used to set any run-time values (loaded
  //  from a meryl database, for example, to find the threshold that contains
  //  95% of the kmers).
  //
public:
  void   finalizeFilterInputs(merylOpTemplate *mot, std::vector<char const *> &err);
  void   finalizeFilterParameters(merylOpTemplate *mot);

  char *
  describe(char *str);


  //  What and how to compare:

  merylFilterRelation  _r = merylFilterRelation::isNOP;  //  Relation we need to satisfy.
  merylFilterQuantity  _q = merylFilterQuantity::isNOP;  //  Quantity we're comparing against.
  bool                 _t = true;                        //  Desired result for a 'true' filter result.

  //  A copy of the command line string, for reporting errors.
  char                           _str[FILENAME_MAX + 1] = { 0 };

  //  Source for left side of comparison:
  //
  //  vIndex:  if uint32max, compare against the appropriate
  //                         constant (as determined by _q).
  //           if 0,         compare against the output kmer
  //           otherwise,    compare against the input kmer from that db

  uint32               _vIndex1 = uint32max;   //  from kmer in specified database
  kmvalu               _vValue1 = 0;           //  from constant 'value' stored here
  kmlabl               _vLabel1 = 0;           //  from constant 'label' stored here
  uint32               _vBases1 = 0;           //  from constant 'bases' stored here ('bases count')

  //  Source for right side of comparison:

  uint32               _vIndex2 = uint32max;   //  from kmer in specified database
  kmvalu               _vValue2 = 0;           //  from constant 'value' stored here
  kmlabl               _vLabel2 = 0;           //  from constant 'label' stored here
  uint32               _vBases2 = 0;           //  from constant 'bases' stored here ('bases count')

  //  For filters on the value, the _vValue2 can be set to a constant or set
  //  based on the histogram of the input database.

  double               _vValue2Distinct = -1;
  double               _vValue2WordFreq = -1;

  //  Some flags on what to compare, so we don't explode the
  //  merylFilterQuantity enum with 15 more types.

  bool                 _countA = false;  //  These four tell us if we need to count the
  bool                 _countC = false;  //  number of A's etc in the kmer (_q == isBases).
  bool                 _countG = false;  //  If so we then compare this count against
  bool                 _countT = false;  //  vBases above.

  //  The 'input' filter needs to store a bunch of info until the whole
  //  command tree is constructed so that it can figure out what 'all inputs'
  //  means.
  //
  //  The 'input' filter tests one of two things:
  //    is the kmer present in a certin number of input databases?
  //    is the kmer present in a set of specific input databases?
  //
  //  The second one is tricky (and expensive).  It requires that the kmer be
  //  present in all the listed databases, but doesn't require that it be
  //  present ONLY in those databases.
  //
  //  Testing this is expensive, and a bit mask would be quicker (for up to
  //  64 inputs).

  //  While parsing the command line, we store the filter specification in two
  //  vectors:
  //    _input_num - a list of the number of inputs a kmer must occur in
  //               - 'n', 'n-m', 'all', 'any', 'n-all'
  //
  //    If _input_num_at_least if not uint32max then this specifies that the
  //    kmer must occur in at least this many inputs -- 'n-all' above.
  //
  //    If _input_num_all is true, then the kmer must occur in all inputs.
  //
  //    If _input_num_any is true, then there are no requirements for
  //    how many databases the kmer must occur in.
  //
  //    _input_idx - a list of specific input indexes a kmer must occur in
  //               - 'first', '@a', '@a-@b'
  //
  std::vector<uint32>            _input_num;
  uint32                         _input_num_at_least = uint32max;
  uint32                         _input_num_all      = false;
  uint32                         _input_num_any      = false;

  std::vector<uint32>            _input_idx;

  //  After the command line is parsed, finalizeFilterInputs() converts the
  //  above into convenient lookup structures.
  //
  //  The filter fails if:
  //       _presentInNum[_actLen]    == false
  //       _presentInIdx[_actIdx[i]] == false (for 0 <= i < _actLen)
  //
  //  The first condition tests if the kmer is present in the specified number of inputs.
  //
  //  The second condition tests if the kmer is present in all of the specified inputs.
  //
  bool                          *_presentInNum = nullptr;
  bool                          *_presentInIdx = nullptr;

  uint32                         _presentInLen  = 0;
  uint32                        *_presentInList = nullptr;
};


#endif  //  MERYLFILTER_H
