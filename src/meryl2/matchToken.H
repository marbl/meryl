
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


//  Expects inputs of
//    opt='<option>[:=]<value>' - the input word to test
//    pat='<pattrn>[:=]'        - the pattern to test against
//  and returns true if <option> is a prefix of <pattern> and the separator
//  matches.  optout is set to the letter after the separator.
//
//  Examples: 1) matchToken("out:", optout, "output:") == true
//            2) matchToken("out",  optout, "output:") == false
//
//  Example 1 is true because 'out' is a prefix of 'output' and both end
//  with a ':'; likewise, example 2 is false because the <opt> string
//  doesn't end with a ':'.
//
inline
bool
matchToken(char const *opt, char const *&optout, char const *pat, bool exact=false) {

  while ((*opt !=  0) && (*opt != ':') && (*opt != '=') &&   //  Advance opt and pat
         (*pat !=  0) &&                                     //  to the separator or
         (*opt == *pat))                                     //  first disagreement.
    opt++, pat++;

  if ((*opt != 0) && (*opt != ':') && (*opt != '='))         //  If opt isn't at its end,
    return false;                                            //  there is no match.

  if (exact == false)                                        //  Advance pat to its end
    while ((*pat != 0) && (*pat != ':') && (*pat != '='))    //  unless we're looking
      pat++;                                                 //  for an exact match.

  if (*opt != *pat)                                          //  If they disagree, the
    return false;                                            //  separators don't match.

  if (*opt)  optout = opt + 1;                               //  Return the next letter,
  else       optout = opt;                                   //  if one exists.

  return true;
}
