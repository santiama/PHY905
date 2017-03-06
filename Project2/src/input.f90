module input
  use linalg
  implicit none    
	real(kind=8) :: rho_max, omega_r, h, hsq, pot 
  integer(kind=4) :: Ngrid, Nsort, Nbody, Niter, i, j
  real(kind=8), allocatable :: a(:,:), eigvecs(:,:), eigvals(:), rho(:)
  namelist /inputs/ Ngrid, rho_max, Nsort, Nbody, omega_r
        
contains
    subroutine read_inputs(infile)
    !{ set user input stuff 
      implicit none
      character infile*40
      open(10,file=infile)
      read(10,nml=inputs)
    !} 
    end subroutine read_inputs
    
    subroutine solve()
    !{ Solves via jacobi
      implicit none
      call jacobi(a, Ngrid, eigvals, eigvecs, Niter, Nsort)
    !}
    end subroutine solve

    subroutine savedata()
    !{ writes data to file
      implicit none

      write (*,*) ''
      write (*,'(A, I8)') ' Number of iterations = ', Niter
			write (*,202)
			write (*,201) (eigvals(i),i=1,Ngrid)
			write (*,203)
			do i = 1,Ngrid
			   write (*,201)  (eigvecs(i,j),j=1,Ngrid)
			end do
			
			201 format (6f12.6)
			202 format (/,' Eigenvalues')
			203 format (/,' Eigenvectors')

    !}
    end subroutine savedata

    subroutine init()
    !{ initialize variable
      implicit none
    
      allocate(a(Ngrid,Ngrid))
      allocate(rho(Ngrid))
      allocate(eigvals(Ngrid))
      allocate(eigvecs(Ngrid,Ngrid))
      
      a = 0.d0
      h = rho_max/(Ngrid+1)
      hsq = h*h
      do i = 1, Ngrid
        rho(i) = h*i
        pot = omega_r**2*rho(i)**2+(Nbody-1)/rho(i)
        a(i,i) = 2.d0/hsq+pot
      end do

		  a(2,1) = -1.d0/hsq
		  a(Ngrid-1,Ngrid) = -1.d0/hsq
		  do i=2,Ngrid-1
		  	a(i+1,i) = -1.d0/hsq
		  	a(i-1,i) = -1.d0/hsq
		  end do
      
			! print a header and the original matrix
      write (*,200)
      do i=1,Ngrid
         write (*,201) (a(i,j),j=1,Ngrid)
      end do

			200 format (' Quantum Harmonic Oscillator (Jacobi method) ',/, &
			            ' Matrix A')
			201 format (12f8.4)


    !}
    end subroutine init

    subroutine dealloc()
    !{ free memory
      implicit none
      deallocate(a)
      deallocate(rho)
      deallocate(eigvals)
      deallocate(eigvecs)
    !}
    end subroutine dealloc

end module input
