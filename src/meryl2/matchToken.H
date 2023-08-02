
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
//    opt='<option>[:=]'
//    str='<string>[:=]'
//  and returns true if <option> is a prefix of <string> and the separator
//  matches.  optout is set to the letter after the separator.
//
//  Of note, opt='word' and str='word:' are NOT matches.  We're expecting opt
//  to end in a ':'.
//
inline
bool
matchToken(char const *opt, char const *&optout, char const *str, bool exact=false) {

  while ((*opt !=  0) && (*opt != ':') && (*opt != '=') &&   //  Advance opt and str
         (*str !=  0) &&                                     //  to the separator or
         (*opt == *str))                                     //  first disagreement.
    opt++, str++;

  if ((*opt != 0) && (*opt != ':') && (*opt != '='))         //  If opt isn't at its end,
    return false;                                            //  there is no match.

  if (exact == false)                                        //  Advance str to its end
    while ((*str != 0) && (*str != ':') && (*str != '='))    //  unless we're looking
      str++;                                                 //  for an exact match.

  if (*opt != *str)                                          //  If they disagree, the
    return false;                                            //  separators don't match.

  if (*opt)  optout = opt + 1;                               //  Return the next letter,
  else       optout = opt;                                   //  if one exists.

  return true;
}
