module HarmonicOscillator 
    use math
    use linalg
    implicit none

    integer(kind=4) :: i,jj, N
	integer(kind=4), allocatable :: indx(:)
    real(kind=8) :: h, hsq, rho
    real(kind=8), allocatable :: xgrid(:), v_sol(:), rhs(:), vgrid(:)
	real(kind=8), allocatable :: Amat(:,:)
    real(kind=8), allocatable :: u_exact(:)
    real(kind=8), allocatable :: a(:),b(:),c(:)
    character(len=8) :: fmtstr
    character(32) x1
    character(len=128) :: filename

contains

    subroutine solve_oscillator()
    !{====================================
    ! = EIGENSOLVE HARMONICOSCILLATOR     =
    ! ====================================
        implicit none
        print *, "Begin Eigensolver N = ",N
        print *, ""
        call mem_init()
        print *, "~~MEMORY ALLOCATED~~"
        call initialize()
        print *, "~~VARIABLES INITIALIZED~~"
        
        call tridiag_general(a, b, c, v_sol, rhs, N)
        call mem_clear()
    !}
    end subroutine solve_general
    
    subroutine mem_init()
    !{====================================
    ! =   ALLOCATES MEMORY FOR VARIABLES =
    ! ====================================
        implicit none
        allocate(u_exact(N))
		allocate(Amat(N,N))
        allocate(xgrid(N))
        allocate(vgrid(N))
        allocate(v_sol(N))
        allocate(rhs(N))
        allocate(a(N))
        allocate(b(N))
        allocate(c(N))
		allocate(indx(N))
    !}
    end subroutine mem_init

    subroutine initialize()
    !{====================================
    ! =   INITIALIZES GRID AND VECTORS   =
    ! ====================================
        implicit none
		Amat(:,:) = 0.d0
		indx(:) = 0.d0
        do i=1,N        
            xgrid(i)= i*h
            rhs(i)= hsq*f(xgrid(i))
            u_exact(i) = exact(xgrid(i)) 
            a(i) = -1.d0
            b(i) = 2.d0
            c(i) = -1.d0
			Amat(i,i)= 2.d0
        end do

		! this is just to avoid if statements
		! for initializing Amat 
		Amat(2,1) = -1.d0
		Amat(N-1,N) = -1.d0
		do i=2,N-1
			Amat(i+1,i) = -1.d0
			Amat(i-1,i) = -1.d0
		end do
    !}
    end subroutine initialize

    subroutine mem_clear()
    !{====================================
    ! =   FREES MEMORY OF VARIABLES      =
    ! ====================================
        implicit none
        deallocate(u_exact)
		deallocate(Amat)
        deallocate(xgrid)
        deallocate(rhs)
        deallocate(vgrid)
        deallocate(v_sol)
        deallocate(a)
        deallocate(b)
        deallocate(c)
		deallocate(indx)
    !}
    end subroutine mem_clear
    
    subroutine potential()
    !{=======================================================================
    ! =  COMPUTES POTENTIAL - ALLOWS FOR SIMPLE MODIFICATION
    ! ==========================================================================
        implicit none
        do i=1,N
            vgrid(i) = rho**2
        end do
    end subroutine potential
    !}

end module HarmonicOscillator 
