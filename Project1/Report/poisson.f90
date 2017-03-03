program poisson
    use globals
    
    integer(kind=4) :: Nmax, j
    Nmax = 4
    
    print *, "==============================================="
    print *, "POISSON 1D SOLVER : TRIDIAGONAL MATRIX ROUTINES"
    print *, " PROJECT 1 - PHY480 - Marco Santia"
    print *, "==============================================="
    print *, ""
    print *, "Max k = ",Nmax
    print *, "where 10^k = N"
    print *, ""
    
    do j=1,Nmax
        N = 10.d0**j
        h = 1.d0/(1.d0+N)
        hsq = h*h       
        call solve_general() 
        call solve_simple()
		call solve_lu()
    end do
end program poisson
