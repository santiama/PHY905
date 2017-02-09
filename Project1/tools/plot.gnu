set terminal epslatex color

set key top right
set xlabel "$x_i$"
set ylabel "$v_i$"

set out "solutionN10.tex" 
plot 'general_solution-00000010.dat' u 1:2 w lp title "General GaussElim", 'general_solution-00000010.dat' u 1:3 w lp title "Exact", 'LU_solution-00000010.dat' u 1:2 w lp title "LU method", 'simple_solution-00000010.dat' u 1:2 w lp title "Simple GaussElim"
set out

set out "solutionN100.tex" 
plot 'general_solution-00000100.dat' u 1:2 w lp title "General GaussElim", 'general_solution-00000100.dat' u 1:3 w lp title "Exact", 'LU_solution-00000100.dat' u 1:2 w lp title "LU method", 'simple_solution-00000100.dat' u 1:2 w lp title "Simple GaussElim"
set out

set out "solutionN1000.tex" 
plot 'general_solution-00001000.dat' u 1:2 w lp title "General GaussElim", 'general_solution-00001000.dat' u 1:3 w lp title "Exact", 'LU_solution-00001000.dat' u 1:2 w lp title "LU method", 'simple_solution-00001000.dat' u 1:2 w lp title "Simple GaussElim"
set out

set out "solutionN10000.tex" 
plot 'general_solution-00010000.dat' u 1:2 w lp title "General GaussElim", 'general_solution-00010000.dat' u 1:3 w lp title "Exact", 'LU_solution-00010000.dat' u 1:2 w lp title "LU method", 'simple_solution-00010000.dat' u 1:2 w lp title "Simple GaussElim"
set out
