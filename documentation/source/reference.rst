.. _reference:

Meryl is a tool for creating and working with DNA sequence k-mers.  K-mers
are typically created by **counting** how many times each k-mer sequence
occurs in some collection of sequences.  Meryl refers to this as the
**value** of the k-mer.  Each k-mer can also be annotated with up to 64-bits
of binary data, its **label**.  The label can be interpreted as a collection of
yes/no flags or binary-coded data, for example, 7 bits could be used to indicate
which of 7 samples the k-mer is present in, or a set of 10 bits could be used to
`unary encode <https://en.wikipedia.org/wiki/Unary_coding>`_ the 
`melting temperature <https://www.sigmaaldrich.com/US/en/technical-documents/protocol/genomics/pcr/oligos-melting-temp>`_
of the k-mer.

.. warning::
  NEED TO REFER TO ACTIONS ABOVE

Databases
=========

A set of k-mers, each k-mer with a value and a label, is stored in a
**database**.  The database is a set of 129 binary files -- 64 data files, 64
index files and one master index -- in a directory.
Meryl processes each of these files independently, using up to 64 compute
threads.

Databases also store the k-mer size (**-k** option) and any simple sequence
reductions (**--compress** and **--ssr** options) applied.  It is not
possible to combine databases with different parameters.

Counting K-mers
===============

The **count** action reads sequence from any number of input files and counts
the number of times each k-mer occurs.

By default, meryl uses **canonical** k-mers.  A canonical k-mer is the
lexicographically smaller of the forward-complement and revese-complement
k-mer.  Actions **count-forward** and **count-reverse** will instead count
k-mers as they are oriented in the input or the reverse-complement,
respectively.

Input sequences can be in either FASTA, FASTQ, raw bases, or if compiled with
Canu support, in a Canu seqStore database.  Sequence files can be gzip, bzip2
or xz compressed.

An output database must be supplied to all count actions.  K-mers are both
written to the output database and provided as input to destination actions.

Count actions, unless accompanied by an action that reads input from an
existing database, MUST specify the k-mer size on the command line with the
**-k** option.

Count actions can include a value or label selector, but cannot include any
filters.  A value selector could be used to assign each k-mer a constant value
instead of the count; a label selector could be used to assign each k-mer
a constant representing the input file.

Counting is resource intense.  Meryl will use memory and threads up to a
limit supplied by: the operating system (usually physical memory and the
number of CPUs), a grid manager (such as Slurm, PBS or SGE) or a command line
option (**-m** and **-t**).

Two algorithms are used for counting kmers.  The algorithm that is expected
to use the least memory is used.  The choice depends on the size of the input
sequences and the k-mer size.

Counting Small k-mers (k < 17)
------------------------------

For k at most 16, meryl counts k-mers directly, that is, by associating an
integer count with each possible k-mer.  This has the benefit of being simple
and uses a constant amount of memory regardless of the size of the input, but
quickly exhausts memory for even moderate kmer sizes.

.. warning::
  BUT THIS METHOD ISN'T used if there isn't a lot of input sequence

There are 4\ :sup:`k` kmers of size k; for k=16, there are 4,294,967,296
possible kmers.  Counting 16-mers with this method will use at least 8
GB of memory, independent of input size: counting 16-mers in an E.coli genome
will use 8 GB of memory, despite there being only 5 million or so k-mers.  Further,
memory usage can increase depending on the maximum count value.

This method uses only a single thread to read the input sequence and
increment counters in the array, but multiple threads can be used to generate
the output database.

.. warning::
  does count really only use one thread here?

Details: Each integer counter is initially a 16-bit value.  Once any count
exceeds 2\ :sup:`16` = 65,535 another bit is added to all value, resulting in
17-bit values for every kmer.  Once any count then exceeds 2\ :sup:`17` =
131,072, another bit is added, and so on.  Thus, memory usage is 512 MB *
log\ :sub:`2` maximum_count_value

Counting Large k-mers (k > 15)
------------------------------

For k larger than 15, or for small amounts of input sequence, meryl counts
k-mers by first converting the sequence to a list of k-mers, duplicates
included, then sorts the list, then scans the list to count the number of
times each k-mer is present.

If all k-mers in an input sequence do not fit in memory, a partial result is
written to disk.  After all input sequences have been processed, the partial
results are combined into a single output database.  In practice, this method
requires several additional gigabytes of memory to minimize the overhead of
writing and merging partial results.

This method can use multiple threads for every stage.

Details: Each kmer is split into a prefix and a suffix.  The prefix is used
to select a list to which the suffix is added.  A trade off is made between a
small prefix (resulting in few lists that store large suffixes) and a large
prefix (resulting in many lists where the overhead of each list could use
more space than the lists themselves).  When the (approximate) size of all
lists exceeds a user-supplied threshold, each list is sorted, the suffixes
are counted, and output to an intermediate database.  After all kmers are
processed, the intermediate databases are merged into one.

Actions
=======

Meryl processing is built around **actions**.  An action loads a k-mer from
one or multiple databases (or, for counting actions, computes the k-mer from
a sequence file) selects a new value and label for it, decides if it should
be output or discarded (e.g., "if the new value is more than 100, output the
k-mer"), and prints it to the screen, saves it to a new database, or
passes it on to another action for further processing.

K-mers are read "in order" from the inputs.  If an input does not contain
the "next" k-mer, it does not participate in the action processing.  For example,
suppose we have three input databases with the following 4-mers and their counts:

.. code-block:: none
  :caption: Sample databases.
  :linenos:

  input-1  input-2  input-3
  AAAA/1   AAAA/2   AAAA/3
  AAAC/1   CAAT/2   CCCC/3
  CAAT/1            GGGG/3
  GGGG/1

A 'union-sum' action with these three databases as input will output:

.. code-block:: none
  :caption: Sample output from union-sum action.
  :linenos:

  AAAA/6 (using the k-mer from input-1, input-2 and input-3)
  AAAC/1 (... from input-1)
  CAAT/3 (... from input-1 and input-2)
  CCCC/3 (... from input-3)
  GGGG/4 (... from input-1 and input-3)

A **selector** selects or computes the output value (label) for each k-mer
from among the input values (labels), or computes an output value (label)
from the input values (labels).  At most one selector can be supplied for the
value or label.

A **filter** decides if the k-mer should be output or discarded.  Filters can
use input values (labels), the new output value (label), the base composition
of the k-mer and how many and which specific inputs the k-mer was present in.
Any number of filters can be supplied, linked with **and**, **or** and
**not** operators.  See FILTERS.

Though it is possible to specify all those choices explicitly, **aliases** are
provided for most common operations.

.. code-block:: none
  :caption: Action aliases.
  :linenos:

  union A B ...         (output the kmer if it is in any input)
  union-min A B ...
  union-max A B ...
  union-sum A B ...

  intersect A B ...     (output the kmer if it is in all inputs)
  intersect-min A B ...
  intersect-max A B ...
  intersect-sum A B ...

  subtract A B ...    (value = A - B - ...)

  difference A B ...  (kmer occurs only in A)

  less-than X DB
  greater-than X DB
  at-least X DB
  at-most X DB
  equal-to X DB
  not-equal-to X DB

  increase X DB
  decrease X DB
  multiple X DB
  divide X DB
  divide-round X DB
  modulo X DB

Aliases exist to support common operations.  An alias sets the 'value',
'label' and 'input' options and so these are not allowed to be used with
aliases.  Examples of aliases and their explicit configuration:

.. warning::
  THIS is NOT CORRECT or COMPLETE!

.. warning::
  Entries with @ in them get rendered as email links

.. table:: Action Aliases
  :widths: 19 17 17 17 15 15

  +----------------+------------------------------------------------------------------------------+
  |                |                                 Action                                       |
  | Alias          +----------------+----------------+--------------+--------------+--------------+
  |                + Value Selector | Label Selector | Input Filter | Value Filter | Label Filter |
  +----------------+----------------+----------------+--------------+--------------+--------------+
  | union          | value=sum      | label=or       | input:any    | value:       | label:       |
  +----------------+----------------+----------------+--------------+--------------+--------------+
  | union-min      | value=min      | label=selected | input:any    | value:       | label:       |
  +----------------+----------------+----------------+--------------+--------------+--------------+
  | union-max      | value=max      | label=selected | input:any    | value:       | label:       |
  +----------------+----------------+----------------+--------------+--------------+--------------+
  | union-sum      | value=sum      | label=or       | input:any    | value:       | label:       |
  +----------------+----------------+----------------+--------------+--------------+--------------+
  +----------------+----------------+----------------+--------------+--------------+--------------+
  | intersect      | value=min      | label=and      | input:all    | value:       | label:       |
  +----------------+----------------+----------------+--------------+--------------+--------------+
  | intersect-min  | value=min      | label=selected | input:all    | value:       | label:       |
  +----------------+----------------+----------------+--------------+--------------+--------------+
  | intersect-max  | value=max      | label=selected | input:all    | value:       | label:       |
  +----------------+----------------+----------------+--------------+--------------+--------------+
  | intersext-sum  | value=sum      | label=and      | input:all    | value:       | label:       |
  +----------------+----------------+----------------+--------------+--------------+--------------+
  +----------------+----------------+----------------+--------------+--------------+--------------+
  | subtract       | value=sub      | label=first    | input:first  | value:       | label:       |
  +----------------+----------------+----------------+--------------+--------------+--------------+
  | difference     | value=sub      | label=first    | input:first  | value:       | label:       |
  +----------------+----------------+----------------+--------------+--------------+--------------+
  +----------------+----------------+----------------+--------------+--------------+--------------+
  | less-than X    | value=first    | label=first    | input:only   | value:<X     | label:       |
  +----------------+----------------+----------------+--------------+--------------+--------------+
  | greater-than X | value=first    | label=first    | input:only   | value:>X     | label:       |
  +----------------+----------------+----------------+--------------+--------------+--------------+
  | at-least X     | value=first    | label=first    | input:only   | value:>=X    | label:       |
  +----------------+----------------+----------------+--------------+--------------+--------------+
  | at-most X      | value=first    | label=first    | input:only   | value:<=X    | label:       |
  +----------------+----------------+----------------+--------------+--------------+--------------+
  | equal-to X     | value=first    | label=first    | input:only   | value:==X    | label:       |
  +----------------+----------------+----------------+--------------+--------------+--------------+
  | not-equal-to X | value=first    | label=first    | input:only   | value:!=X    | label:       |
  +----------------+----------------+----------------+--------------+--------------+--------------+
  +----------------+----------------+----------------+--------------+--------------+--------------+
  | increase X     | value=@1+X     | label=first    | input:only   | value:       | label:       |
  +----------------+----------------+----------------+--------------+--------------+--------------+
  | decrease X     | value=@1-X     | label=first    | input:only   | value:       | label:       |
  +----------------+----------------+----------------+--------------+--------------+--------------+
  | multiply X     | value=@1*X     | label=first    | input:only   | value:       | label:       |
  +----------------+----------------+----------------+--------------+--------------+--------------+
  | divide X       | value=@1/X     | label=first    | input:only   | value:       | label:       |
  +----------------+----------------+----------------+--------------+--------------+--------------+
  | divide-round X | value=@1/X ??  | label=first    | input:only   | value:       | label:       |
  +----------------+----------------+----------------+--------------+--------------+--------------+
  | modulo X       | value=@1%X     | label=first    | input:only   | value:       | label:       |
  +----------------+----------------+----------------+--------------+--------------+--------------+

('disjoint' is the former 'symmetric-difference' alias.)

Value Selectors
---------------

A **value selector** selects (or computes) the output value of the k-mer
based on the input values and possibly a single integer constant.

.. note::
  `(#X)` means to also include constant `X` in the computation.  Constants
  can be decimal integers (`123` or `123d`), hexadecimal (`abch`),
  octal (`147o`) or binary (`0101010b`).  SI suffixes can be used on plain decimal integers
  (`123k` == 123,000; `1mi` == 1,048,576).

.. warning::
  Things like value=@1-@2 are NOT supported.  Even the potentially useful
  value=@1 isn't supported (though it is listed below).

.. warning::
  value=selected isn't implemented.

.. table:: Value Selectors
  :widths: 20 80

  +--------------------+-------------------------------------------------+
  | Selector           | Set value to ...                                |
  +====================+=================================================+
  | value=#X           | ...constant X.                                  |
  +--------------------+-------------------------------------------------+
  | value=@X           | ...that of the k-mer in the Xth input           |
  +--------------------+-------------------------------------------------+
  | value=first        | ...that of the k-mer in the first input.        |
  +--------------------+-------------------------------------------------+
  | value=selected     | ...that of the k-mer selected by the label=     |
  |                    | selector.  When multiple k-mers are selected,   |
  |                    | the value of the first is used.                 |
  +--------------------+-------------------------------------------------+
  | value=min(#X)      | ...the minimum of all input values.             |
  +--------------------+-------------------------------------------------+
  | value=max(#X)      | ...the maximum of all input values.             |
  +--------------------+-------------------------------------------------+
  | value=add(#X)      | ...the sum of all input values.                 |
  +--------------------+-------------------------------------------------+
  | value=sum(#X)      | ...the sum of all input values.                 |
  +--------------------+-------------------------------------------------+
  | value=sub(#X)      | ...the value of the k-mer in the first input    |
  |                    | minus all other values.                         |
  +--------------------+-------------------------------------------------+
  | value=dif(#X)      | ...the value of the k-mer in the first input    |
  |                    | minus all other values.                         |
  +--------------------+-------------------------------------------------+
  | value=mul(#X)      | ...the product of all input values.             |
  +--------------------+-------------------------------------------------+
  | value=div(#X)      | ...the value of the k-mer in the first input    |
  |                    | divided by all other values.                    |
  +--------------------+-------------------------------------------------+
  | value=divzero(#X)  | ...the value of the k-mer in the first input    |
  |                    | divided by all other values, rounding zero up   |
  |                    | to one.                                         |
  +--------------------+-------------------------------------------------+
  | value=mod(#X)      | ...the remainder after the value of the k-mer in|
  |                    | the first input is divided by all other values. |
  +--------------------+-------------------------------------------------+
  | value=rem(#X)      | ...the remainder after the value of the k-mer in|
  |                    | the first input is divided by all other values. |
  +--------------------+-------------------------------------------------+
  | value=count        | ...the number of inputs the k-mer is present in.|
  +--------------------+-------------------------------------------------+


Label Selectors
---------------

A **label selector** selects (or computes) the output label of the k-mer
based on the input label and possibly a single 64-bit constant.

.. table:: Value Selectors
  :widths: 20 80

  +------------------------+-------------------------------------------------+
  | Selector               | Set label to ...                                |
  +========================+=================================================+
  | label=#X               | ...constant X.                                  |
  +------------------------+-------------------------------------------------+
  | label=@X               | ...that of the k-mer in the Xth input           |
  +------------------------+-------------------------------------------------+
  | label=first            | ...that of the k-mer in the first input.        |
  +------------------------+-------------------------------------------------+
  | label=selected         | ...that of the k-mer selected by the value=     |
  |                        | selector.  When multiple k-mers are selected,   |
  |                        | the label of the first is used.                 |
  +------------------------+-------------------------------------------------+
  | label=min(#X)          | ...the minimum of all input labels.             |
  +------------------------+-------------------------------------------------+
  | label=max(#X)          | ...the maximum of all input labels.             |
  +------------------------+-------------------------------------------------+
  | label=and(#X)          | ...the bitwise AND of all input labels.         |
  +------------------------+-------------------------------------------------+
  | label=or(#X)           | ...the bitwise OR of all input labels.          |
  +------------------------+-------------------------------------------------+
  | label=xor(#X)          | ...the bitwise XOR of all input labels.         |
  +------------------------+-------------------------------------------------+
  | label=difference(#X)   | ... ????                                        |
  +------------------------+-------------------------------------------------+
  | label=lightest(#X)     | ...the label with the fewest bit set.           |
  +------------------------+-------------------------------------------------+
  | label=heaviest(#X)     | ...the label with the most bits set.            |
  +------------------------+-------------------------------------------------+
  | label=invert(#X)       | ...the bitwise invert of the first input.       |
  +------------------------+-------------------------------------------------+
  | label=shift-left(#X)   | ...the first input shifted left by X places.    |
  +------------------------+-------------------------------------------------+
  | label=shift-right(#X)  | ...the first input shifted right by X places.   |
  +------------------------+-------------------------------------------------+
  | label=rotate-left(#X)  | ...the first input rotated left by X places.    |
  +------------------------+-------------------------------------------------+
  | label=rotate-right(#X) | ...the first input rotated right by X places.   |
  +------------------------+-------------------------------------------------+

Value Filters
-------------

Label Filters
-------------

Base Composition Filters
------------------------

Input Filters
-------------

Processing Trees
================

Meryl processes kmers using a tree of actions.  An action reads kmers from
multiple inputs, computes a function on the values and labels of all inputs
with the same kmer, and outputs a single kmer with a single value and a
single label.

(An action can also read sequence files and count the kmers.)

Each action in the tree is enclosed in square brackets.  Square brackets
around the top-level / outermost action are optional.

The input to an action can be either a meryl database on disk or the output of
a different action.

The 'union' action below reads input from meryl databases 'input-1.meryl' and
'input-2.meryl'.  All three forms below are equivalent.

.. code-block:: none
  :caption: A simple union action reading from two inputs.
  :linenos:

  [ union input-1.meryl input-2.meryl ]

.. code-block:: none
  :caption: A simple union action reading from two inputs, but formatted.
  :linenos:

    union
      input-1.meryl
      input-2.meryl

.. code-block:: none
  :caption: A simple union action reading from two inputs, as sub-actions.
  :linenos:

    union
      [ input-1.meryl ]   //  This form technically makes input-1 and input-2 into
      [ input-2.meryl ]   //  sub-actions instead of direct inputs to 'union'.

Sub-actions can pre-process inputs.  The 'intersect' action below reads input
from two counting actions, and the one after computes a `union` before the
`intersection`.

.. code-block:: none
  :caption: Sample databases.
  :linenos:

  intersect 
    [ count input-1.fasta ]
    [ count input-2.fasta ]

Each action will automatically pass its output kmers to the parent action,
and can optionally write them to an output database.

.. code-block:: none
  :caption: Sample databases.
  :linenos:

  intersect output=abINT12.meryl
    [ union input-a.meryl input-b.meryl output=ab.meryl ]
    [ union input-1.meryl input-2.meryl output=12.meryl ]

The original meryl allowed sub-actions to be supplied without surrounding
square brackets, but this led to great ambiguity in which action the output
modifier was associated with.  Without brackets, the following is ambiguous:

.. code-block:: none
  :caption: Sample databases.
  :linenos:

  meryl
    union
      intersect
        a.meryl
        b.meryl
      intersect
        c.meryl
        d.meryl

As written, the intent is clear, but meryl interprets the second 'intersect' action
as an input to the first:

.. code-block:: none
  :caption: Sample databases.
  :linenos:

  meryl
    union
      intersect
        a.meryl
        b.meryl
        intersect
          c.meryl
          d.meryl

Therefore, meryl2 **requires** actions (except the very first) to be
surrounded by square brackets.
