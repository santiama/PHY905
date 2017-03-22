Flat profile:

Each sample counts as 0.01 seconds.
  %   cumulative   self              self     total           
 time   seconds   seconds    calls  ms/call  ms/call  name    
100.07      0.53     0.53        1   530.36   530.36  __linalg_MOD_jacobi_cyclic
  0.00      0.53     0.00        5     0.00     0.00  __linalg_MOD_offdiag_max
  0.00      0.53     0.00        1     0.00   530.36  MAIN__
  0.00      0.53     0.00        1     0.00     0.00  __input_MOD_dealloc
  0.00      0.53     0.00        1     0.00     0.00  __input_MOD_init
  0.00      0.53     0.00        1     0.00     0.00  __input_MOD_read_inputs
  0.00      0.53     0.00        1     0.00     0.00  __input_MOD_savedata
  0.00      0.53     0.00        1     0.00   530.36  __input_MOD_solve
  0.00      0.53     0.00        1     0.00     0.00  __linalg_MOD_jacobi_classical
  0.00      0.53     0.00        1     0.00     0.00  __unit_tests_MOD_gram_test
  0.00      0.53     0.00        1     0.00     0.00  __unit_tests_MOD_init_tests
  0.00      0.53     0.00        1     0.00     0.00  __unit_tests_MOD_max_test

 %         the percentage of the total running time of the
time       program used by this function.

cumulative a running sum of the number of seconds accounted
 seconds   for by this function and those listed above it.

 self      the number of seconds accounted for by this
seconds    function alone.  This is the major sort for this
           listing.

calls      the number of times this function was invoked, if
           this function is profiled, else blank.
 
 self      the average number of milliseconds spent in this
ms/call    function per call, if this function is profiled,
	   else blank.

 total     the average number of milliseconds spent in this
ms/call    function and its descendents per call, if this 
	   function is profiled, else blank.

name       the name of the function.  This is the minor sort
           for this listing. The index shows the location of
	   the function in the gprof listing. If the index is
	   in parenthesis it shows where it would appear in
	   the gprof listing if it were to be printed.

		     Call graph (explanation follows)


granularity: each sample hit covers 2 byte(s) for 1.89% of 0.53 seconds

index % time    self  children    called     name
                0.00    0.53       1/1           main [2]
[1]    100.0    0.00    0.53       1         MAIN__ [1]
                0.00    0.53       1/1           __input_MOD_solve [3]
                0.00    0.00       1/1           __unit_tests_MOD_init_tests [16]
                0.00    0.00       1/1           __unit_tests_MOD_max_test [17]
                0.00    0.00       1/1           __unit_tests_MOD_gram_test [15]
                0.00    0.00       1/1           __input_MOD_read_inputs [12]
                0.00    0.00       1/1           __input_MOD_init [11]
                0.00    0.00       1/1           __input_MOD_savedata [13]
                0.00    0.00       1/1           __input_MOD_dealloc [10]
-----------------------------------------------
                                                 <spontaneous>
[2]    100.0    0.00    0.53                 main [2]
                0.00    0.53       1/1           MAIN__ [1]
-----------------------------------------------
                0.00    0.53       1/1           MAIN__ [1]
[3]    100.0    0.00    0.53       1         __input_MOD_solve [3]
                0.53    0.00       1/1           __linalg_MOD_jacobi_cyclic [4]
-----------------------------------------------
                0.53    0.00       1/1           __input_MOD_solve [3]
[4]    100.0    0.53    0.00       1         __linalg_MOD_jacobi_cyclic [4]
-----------------------------------------------
                0.00    0.00       1/5           __unit_tests_MOD_max_test [17]
                0.00    0.00       4/5           __linalg_MOD_jacobi_classical [14]
[9]      0.0    0.00    0.00       5         __linalg_MOD_offdiag_max [9]
-----------------------------------------------
                0.00    0.00       1/1           MAIN__ [1]
[10]     0.0    0.00    0.00       1         __input_MOD_dealloc [10]
-----------------------------------------------
                0.00    0.00       1/1           MAIN__ [1]
[11]     0.0    0.00    0.00       1         __input_MOD_init [11]
-----------------------------------------------
                0.00    0.00       1/1           MAIN__ [1]
[12]     0.0    0.00    0.00       1         __input_MOD_read_inputs [12]
-----------------------------------------------
                0.00    0.00       1/1           MAIN__ [1]
[13]     0.0    0.00    0.00       1         __input_MOD_savedata [13]
-----------------------------------------------
                0.00    0.00       1/1           __unit_tests_MOD_gram_test [15]
[14]     0.0    0.00    0.00       1         __linalg_MOD_jacobi_classical [14]
                0.00    0.00       4/5           __linalg_MOD_offdiag_max [9]
-----------------------------------------------
                0.00    0.00       1/1           MAIN__ [1]
[15]     0.0    0.00    0.00       1         __unit_tests_MOD_gram_test [15]
                0.00    0.00       1/1           __linalg_MOD_jacobi_classical [14]
-----------------------------------------------
                0.00    0.00       1/1           MAIN__ [1]
[16]     0.0    0.00    0.00       1         __unit_tests_MOD_init_tests [16]
-----------------------------------------------
                0.00    0.00       1/1           MAIN__ [1]
[17]     0.0    0.00    0.00       1         __unit_tests_MOD_max_test [17]
                0.00    0.00       1/5           __linalg_MOD_offdiag_max [9]
-----------------------------------------------

 This table describes the call tree of the program, and was sorted by
 the total amount of time spent in each function and its children.

 Each entry in this table consists of several lines.  The line with the
 index number at the left hand margin lists the current function.
 The lines above it list the functions that called this function,
 and the lines below it list the functions this one called.
 This line lists:
     index	A unique number given to each element of the table.
		Index numbers are sorted numerically.
		The index number is printed next to every function name so
		it is easier to look up where the function in the table.

     % time	This is the percentage of the `total' time that was spent
		in this function and its children.  Note that due to
		different viewpoints, functions excluded by options, etc,
		these numbers will NOT add up to 100%.

     self	This is the total amount of time spent in this function.

     children	This is the total amount of time propagated into this
		function by its children.

     called	This is the number of times the function was called.
		If the function called itself recursively, the number
		only includes non-recursive calls, and is followed by
		a `+' and the number of recursive calls.

     name	The name of the current function.  The index number is
		printed after it.  If the function is a member of a
		cycle, the cycle number is printed between the
		function's name and the index number.


 For the function's parents, the fields have the following meanings:

     self	This is the amount of time that was propagated directly
		from the function into this parent.

     children	This is the amount of time that was propagated from
		the function's children into this parent.

     called	This is the number of times this parent called the
		function `/' the total number of times the function
		was called.  Recursive calls to the function are not
		included in the number after the `/'.

     name	This is the name of the parent.  The parent's index
		number is printed after it.  If the parent is a
		member of a cycle, the cycle number is printed between
		the name and the index number.

 If the parents of the function cannot be determined, the word
 `<spontaneous>' is printed in the `name' field, and all the other
 fields are blank.

 For the function's children, the fields have the following meanings:

     self	This is the amount of time that was propagated directly
		from the child into the function.

     children	This is the amount of time that was propagated from the
		child's children to the function.

     called	This is the number of times the function called
		this child `/' the total number of times the child
		was called.  Recursive calls by the child are not
		listed in the number after the `/'.

     name	This is the name of the child.  The child's index
		number is printed after it.  If the child is a
		member of a cycle, the cycle number is printed
		between the name and the index number.

 If there are any cycles (circles) in the call graph, there is an
 entry for the cycle-as-a-whole.  This entry shows who called the
 cycle (as parents) and the members of the cycle (as children.)
 The `+' recursive calls entry shows the number of function calls that
 were internal to the cycle, and the calls entry for each member shows,
 for that member, how many times it was called from other members of
 the cycle.


Index by function name

   [1] MAIN__                 [13] __input_MOD_savedata    [9] __linalg_MOD_offdiag_max
  [10] __input_MOD_dealloc     [3] __input_MOD_solve      [15] __unit_tests_MOD_gram_test
  [11] __input_MOD_init       [14] __linalg_MOD_jacobi_classical [16] __unit_tests_MOD_init_tests
  [12] __input_MOD_read_inputs [4] __linalg_MOD_jacobi_cyclic [17] __unit_tests_MOD_max_test
