module globals 
    use math
    use linalg
    implicit none

    integer(kind=4) :: i,jj, N
	integer(kind=4), allocatable :: indx(:)
    real(kind=8) :: h, hsq
    real(kind=8), allocatable :: xgrid(:), v_sol(:), rhs(:), error(:)
	real(kind=8), allocatable :: Amat(:,:)
    real(kind=8), allocatable :: u_exact(:)
    real(kind=8), allocatable :: a(:),b(:),c(:)
    character(len=8) :: fmtstr
    character(32) x1
    character(len=128) :: filename

contains
    subroutine solve_general()
    !{====================================
    ! = GENERAL TRIDIAG METHOD SOLVE     =
    ! ====================================
        implicit none
        real(kind=8) :: t1,t2
        print *, "Begin General Solver N = ",N
        print *, ""
        call mem_init()
        print *, "~~MEMORY ALLOCATED~~"
        call initialize()
        
        fmtstr = '(I8.8)'
        write(x1,fmtstr) N
        filename='general_solution-'//trim(x1)//'.dat'
        open(unit=73, file=filename, status="replace", action="write")

        print *, "~~VARIABLES INITIALIZED~~"
        
        call CPU_TIME(t1)
        call tridiag_general(a, b, c, v_sol, rhs, N)
        call CPU_TIME(t2)
        
        do i=1, N
            error(i) = log10(abs((v_sol(i)-u_exact(i))/u_exact(i)))
        end do

        print '("  SOLVE COMPLETED IN ",E16.5," SECONDS")',t2-t1
        print *, " writing to file... "

        write(73,'("# N = ",I8)') N
        write(73,'("# log10(h) = ",f14.10)') log10(h)
        write(73,'("# Max Error = ",f14.10)') maxval(error) 
        write(73,'("# Runtime: ",E16.5," s")') t2-t1
        write(73,'("#=====================================================")')
        write(73,'("#   x             v(x)             u(x)       Rel.Err")') 
        write(73,'("#=====================================================")')

        do i=1,N
            write(73,'(4f14.10)') xgrid(i), v_sol(i), u_exact(i), error(i)
        end do
        print *, " results written to file: ",filename
        print *, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        print *, ""
        close(unit=73)
        call mem_clear()
    !}
    end subroutine solve_general
    
    subroutine solve_simple()
    !{====================================
    ! = SIMPLIFIED TRIDIAG METHOD SOLVE  =
    ! ====================================
        implicit none
        real(kind=8) :: t1,t2
        print *, "Begin Simple Solver N = ",N
        print *, ""
        call mem_init()
        print *, "~~MEMORY ALLOCATED~~"
        call initialize()
        
        fmtstr = '(I8.8)'
        write(x1,fmtstr) N
        filename='simple_solution-'//trim(x1)//'.dat'
        open(unit=73, file=filename, status="replace", action="write")

        print *, "~~VARIABLES INITIALIZED~~"
        
        call CPU_TIME(t1)
        call tridiag_special(a, b, c, v_sol, rhs, N)
        call CPU_TIME(t2)
        
        do i=1, N
            error(i) = log10(abs((v_sol(i)-u_exact(i))/u_exact(i)))
        end do

        print '("  SOLVE COMPLETED IN ",E16.5," SECONDS")',t2-t1
        print *, " writing to file... "

        write(73,'("# N = ",I8)') N
        write(73,'("# log10(h) = ",f14.10)') log10(h)
        write(73,'("# Max Error = ",f14.10)') maxval(error) 
        write(73,'("# Runtime: ",E16.5," s")') t2-t1
        write(73,'("#=====================================================")')
        write(73,'("#   x             v(x)             u(x)       Rel.Err")') 
        write(73,'("#=====================================================")')

        do i=1,N
            write(73,'(4f14.10)') xgrid(i), v_sol(i), u_exact(i), error(i)
        end do
        print *, " results written to file: ",filename
        print *, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        print *, ""
        close(unit=73)
        call mem_clear()
    !}
    end subroutine solve_simple

	subroutine solve_LU()
    !{====================================
    ! = LU DECOMPOSITION METHOD SOLVE    =
    ! ====================================
		implicit none

		real(kind=8) :: det
        real(kind=8) :: t1,t2
        print *, "Begin LU Solver N = ",N
        print *, ""
        call mem_init()
        print *, "~~MEMORY ALLOCATED~~"
        call initialize()
        
        fmtstr = '(I8.8)'
        write(x1,fmtstr) N
        filename='LU_solution-'//trim(x1)//'.dat'
        open(unit=73, file=filename, status="replace", action="write")

        print *, "~~VARIABLES INITIALIZED~~"
        
        call CPU_TIME(t1)
		call lu_decompose(Amat,N,indx,det)
		call lu_linear_equation(Amat,N,indx,rhs)
        call CPU_TIME(t2)
        
        do i=1, N
            error(i) = log10(abs((rhs(i)-u_exact(i))/u_exact(i)))
        end do

        print '("  SOLVE COMPLETED IN ",E16.5," SECONDS")',t2-t1
        print *, " writing to file... "

        write(73,'("# N = ",I8)') N
        write(73,'("# log10(h) = ",f14.10)') log10(h)
        write(73,'("# Max Error = ",f14.10)') maxval(error) 
        write(73,'("# Runtime: ",E16.5," s")') t2-t1
        write(73,'("#=====================================================")')
        write(73,'("#   x             v(x)             u(x)       Rel.Err")') 
        write(73,'("#=====================================================")')

        do i=1,N
            write(73,'(4f14.10)') xgrid(i), rhs(i), u_exact(i), error(i)
        end do
        print *, " results written to file: ",filename
        print *, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        print *, ""
        close(unit=73)
        call mem_clear()

	!}
	end subroutine solve_LU

    subroutine mem_init()
    !{====================================
    ! =   ALLOCATES MEMORY FOR VARIABLES =
    ! ====================================
        implicit none
        allocate(u_exact(N))
		allocate(Amat(N,N))
        allocate(xgrid(N))
        allocate(error(N))
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
        deallocate(error)
        deallocate(v_sol)
        deallocate(a)
        deallocate(b)
        deallocate(c)
		deallocate(indx)
    !}
    end subroutine mem_clear

end module globals 
