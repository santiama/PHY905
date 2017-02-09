program solver 
    use input
    use HarmonicOscillator

    print *, "==============================================="
    print *, "ELLIPTIC 1D SOLVER : TRIDIAGONAL MATRIX ROUTINES"
    print *, "==============================================="
    print *, ""
    
    call read_inputs()
 
    if (solver_type .EQ.'PO') then
        do j=Nmin,Nmax
            N = 10.d0**j
            h = 1.d0/(1.d0+N)
            hsq = h*h       
            call solve_general() 
            call solve_simple()
		    call solve_lu()
        end do

    elseif (solver_type .EQ.'HO') then
        do j=Nmin,Nmax
            N = 10.d0**j
            h = 1.d0/(1.d0+N)
            hsq = h*h   
        end do
    end if            
end program solver
