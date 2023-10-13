.. _grammar:

The command line is scanned for global options, and those are removed.  Other words
are appended to a single string separated by \n letters.  Internal spaces are
left as is.

For programs from files, a similar flow, except that
spaces between single and double quotes are preserved while
all other while space is replaced by newlines.

Aliases and alias options are handled by exact match.  This requires a flag
in command-builder to indicate if we're in an alias and if we're waiting for
an option.



GLOBAL OPTIONS
--------------

  ``-k <k>`` | ``k=<k>`` (legacy)
    Set k-mer size to ``<k>``.  If any databases are supplied this value is
    automatically discovered.  Any ``<k>`` supplied on the command line must
    agree with all databases.

    .. code-block:: none

        > meryl2 -k 23 o:stats basic-dh-1.meryl
        terminateOperations()- Action #1 at stack position 0 possibly does nothing.
        mer size mismatch, can't process this set of files.

  ``-c`` | ``--hpc`` | ``--homopoly-compress``
    Compress runs of the same base to a single base.

  ``--ssr <specification>``
    Apply a simple-sequence-reduction ``<specification>`` to input sequence before
    a counting operation.  [NOT IMPLEMENTED]

  ``-l <l>``
    Set label size to ``<l>`` bits.  If any databases are supplied this value
    is automatically discovered.  The maximum of all input atabases and any
    ``<l>`` parameter is used for the output label size.

  ``-m <m>`` | ``-memory <m>`` | ``--memory <m>`` | ``memory=<m>`` (legacy)
    Use up to ``<m>`` bytes memory.  This primarily affects counting operations.
    SI suffixes and common bases are allowed in the form ``<number><base><suffix>``:

    .. code-block:: none

      5g    ==   5 * 10^3 bytes  exactly 5 billion bytes
      101bm ==   5 * 10^2 bytes  exactly 5 million bytes
      5gi   ==   5 * 2^30 bytes  approximately 5.37 bilion bytes
      fhgi  == 0xf * 2^30 bytes  approximately 17.2 billion bytes

    If not specified, the limit is set to the smaller of the physical memory
    and any limit imposed by a grid scheduler.

  ``-t <t>`` | ``-threads <t>`` | ``--threads <t>`` | ``threads=<t>`` (legacy)
    Use ``<t>`` threads for processing, both during counting and normal processing.
    If not specified, the limit is set to either the number of physical cores present,
    or any limit imposed by a grid scheduler.

  ``-V``
    Increase verbosity of processing by the length of the option.  ``-V``
    increases verbosity by one, ``-VV`` by two, ``-Verbose`` by 7 and
    ``-verisimilitudinous`` will leave no detail unreported.  The maximum
    useful value is 5:

    .. code-block:: none

      0 - No messages (set by -Q)
      1 - Show compiled program before processing, the default
      2 - Show computation progress
      3 - Show compilation progress
      4 - Show details of compilation
      5 - Everything

  ``-Q``
    Display no logging or progress messages at all.

  ``-P``
    Enable progress reports.

  ``-C``
    Only compile the meryl processing tree.

  ``-h`` | ``-help`` | ``--help``
    Dislay command line usage.


ALIASES
-------

Aliases provide an easy way to specify common actions.  They configure both
the **assignment** and **selection** of an action.  They do not specify any
inputs or outputs.  Aliases are built-in to meryl or can be created by the user.

.. table:: Action Aliases
  :widths: 20 20 30 30

  +--------------------+--------------------+----------------------------------------------+----------------------------------------------+
  | Alias              | Output k-mer if... | Sets value to the...                         | Sets label to the...                         |
  +====================+====================+==============================================+==============================================+
  | union              | ...k-mer is in any | ...number of databases the k-mer is in.      | ...bitwise OR of all labels.                 |
  +--------------------+ input database.    +----------------------------------------------+----------------------------------------------+
  | union-min          |                    | ...smallest input value.                     | ...label of the selected k-mer.              |
  +--------------------+                    +----------------------------------------------+----------------------------------------------+
  | union-max          |                    | ...largest input value.                      | ...label of the selected k-mer.              |
  +--------------------+                    +----------------------------------------------+----------------------------------------------+
  | union-sum          |                    | ...sum of the input values.                  | ...bitwise OR of all labels.                 |
  +--------------------+--------------------+----------------------------------------------+----------------------------------------------+
  +--------------------+--------------------+----------------------------------------------+----------------------------------------------+
  | intersect          | ...k-mer is in all | ...value of the k-mer in the first database. | ...bitwise AND of all labels.                |
  +--------------------+ input databases.   +----------------------------------------------+----------------------------------------------+
  | intersect-min      |                    | ...smallest input value.                     | ...label of the selected k-mer.              |
  +--------------------+                    +----------------------------------------------+----------------------------------------------+
  | intersect-max      |                    | ...largest input value.                      | ...label of the selected k-mer.              |
  +--------------------+                    +----------------------------------------------+----------------------------------------------+
  | intersect-sum      |                    | ...sum of the input values.                  | ...bitwise AND of all labels.                |
  +--------------------+--------------------+----------------------------------------------+----------------------------------------------+
  +--------------------+--------------------+----------------------------------------------+----------------------------------------------+
  | subtract           | ...k-mer is in the | ...value of the k-mer in the first database  | ...bits set only in the k-mer in the         |
  |                    | first database.    | minus all other values.                      | first database.                              |
  +--------------------+                    +----------------------------------------------+----------------------------------------------+
  | difference         |                    | ...value of the k-mer in the first database. | ...label of the k-mer in the first database. |
  +--------------------+--------------------+----------------------------------------------+----------------------------------------------+
  +--------------------+--------------------+----------------------------------------------+----------------------------------------------+
  | less-than X        | ...k-mer is in the | ...value of the k-mer.                       | ...label of the k-mer.                       |
  +--------------------+ first and only     |                                              |                                              |
  | greater-than X     | database and the   |                                              |                                              |
  +--------------------+ value meets the    |                                              |                                              |
  | at-least X         | speficied          |                                              |                                              |
  +--------------------+ condition.         |                                              |                                              |
  | at-most X          |                    |                                              |                                              |
  +--------------------+                    |                                              |                                              |
  | equal-to X         |                    |                                              |                                              |
  +--------------------+                    |                                              |                                              |
  | not-equal-to X     |                    |                                              |                                              |
  +--------------------+--------------------+----------------------------------------------+----------------------------------------------+
  +--------------------+--------------------+----------------------------------------------+----------------------------------------------+
  | increase X         | ...k-mer is in the | ...value of the k-mer modified by            | ...label of the k-mer.                       |
  +--------------------+ first and only     | the specified operation.                     |                                              |
  | decrease X         | database.          |                                              |                                              |
  +--------------------+                    | (divide-round rounds 0 up to 1)              |                                              |
  | multiple X         |                    |                                              |                                              |
  +--------------------+                    |                                              |                                              |
  | divide X           |                    |                                              |                                              |
  +--------------------+                    |                                              |                                              |
  | divide-round X     |                    |                                              |                                              |
  +--------------------+                    |                                              |                                              |
  | modulo X           |                    |                                              |                                              |
  +--------------------+--------------------+----------------------------------------------+----------------------------------------------+

.. table:: Action Aliases
  :widths: 28 18 18 18 18

  +----------------------+-----------------------------------------------------------+
  |                      |                                  Action                   |
  | Alias                +---------------------------------+-------------------------+
  |                      + Assignment                      | Selector [#a]_          |
  +----------------------+-----------------+---------------+------------+------------+
  | union                | a:v=count       | a:l=or        | s:i:any    |            |
  +----------------------+-----------------+---------------+------------+------------+
  | union-min            | a:v=min         | a:l=selected  | s:i:any    |            |
  +----------------------+-----------------+---------------+------------+------------+
  | union-max            | a:v=max         | a:l=selected  | s:i:any    |            |
  +----------------------+-----------------+---------------+------------+------------+
  | union-sum            | a:v=sum         | a:l=or        | s:i:any    |            |
  +----------------------+-----------------+---------------+------------+------------+
  +----------------------+-----------------+---------------+------------+------------+
  | intersect            | a:v=first       | a:l=and       | s:i:all    |            |
  +----------------------+-----------------+---------------+------------+------------+
  | intersect-min        | a:v=min         | a:l=selected  | s:i:all    |            |
  +----------------------+-----------------+---------------+------------+------------+
  | intersect-max        | a:v=max         | a:l=selected  | s:i:all    |            |
  +----------------------+-----------------+---------------+------------+------------+
  | intersext-sum        | a:v=sum         | a:l=and       | s:i:all    |            |
  +----------------------+-----------------+---------------+------------+------------+
  +----------------------+-----------------+---------------+------------+------------+
  | subtract             | a:v=sub         | a:l=difference| s:i:first  |            |
  +----------------------+-----------------+---------------+------------+------------+
  | difference           | a:v=sub         | a:l=first     | s:i:first  |            |
  +----------------------+-----------------+---------------+------------+------------+
  +----------------------+-----------------+---------------+------------+------------+
  | less-than X [#b]_    | a:v=first       | a:l=first     | s:i:only   | s:v:<X     |
  +----------------------+-----------------+---------------+------------+------------+
  | greater-than X [#b]_ | a:v=first       | a:l=first     | s:i:only   | s:v:>X     |
  +----------------------+-----------------+---------------+------------+------------+
  | at-least X [#b]_     | a:v=first       | a:l=first     | s:i:only   | s:v:>=X    |
  +----------------------+-----------------+---------------+------------+------------+
  | at-most X [#b]_      | a:v=first       | a:l=first     | s:i:only   | s:v:<=X    |
  +----------------------+-----------------+---------------+------------+------------+
  | equal-to X [#b]_     | a:v=first       | a:l=first     | s:i:only   | s:v:==X    |
  +----------------------+-----------------+---------------+------------+------------+
  | not-equal-to X [#b]_ | a:v=first       | a:l=first     | s:i:only   | s:v:!=X    |
  +----------------------+-----------------+---------------+------------+------------+
  +----------------------+-----------------+---------------+------------+------------+
  | increase X           | a:v=\@1+X       | a:l=first     | s:i:only   |            |
  +----------------------+-----------------+---------------+------------+------------+
  | decrease X           | a:v=\@1-X       | a:l=first     | s:i:only   |            |
  +----------------------+-----------------+---------------+------------+------------+
  | multiply X           | a:v=\@1*X       | a:l=first     | s:i:only   |            |
  +----------------------+-----------------+---------------+------------+------------+
  | divide X             | a:v=\@1/X       | a:l=first     | s:i:only   |            |
  +----------------------+-----------------+---------------+------------+------------+
  | divide-round X       | a:v=\@1/X [#c]_ | a:l=first     | s:i:only   |            |
  +----------------------+-----------------+---------------+------------+------------+
  | modulo X             | a:v=\@1%X       | a:l=first     | s:i:only   |            |
  +----------------------+-----------------+---------------+------------+------------+

.. [#a] Selectors not listed are not used.
.. [#b] ``X`` can be an integer or ``distinct=<fraction>``, ``word-frequency=<fraction>``, or ``threshold=<integer>``.  See select:value for details.
.. [#c] The ``divide-round`` alias rounds values of 0 up to 1.

