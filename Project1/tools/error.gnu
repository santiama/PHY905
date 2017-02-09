set terminal epslatex color
set out 'error.tex'

set key top right

set xlabel "$\log_{10}(h)$"
set ylabel "Max Error"

plot 'error.dat' u 1:2 w lp title "Max Relative Error"
set out
