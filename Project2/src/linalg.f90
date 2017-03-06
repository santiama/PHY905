module linalg
    implicit none

contains

	subroutine jacobi(a,n,eigvals,x,Niter,Nsort)
	!{ ===========================================================
	implicit none
	integer :: nmax, i, j, k, sort_id
	integer(kind=4), intent(inout) ::  n, Nsort, Niter
  real(kind=8) :: sort_tmp(n)
	real(kind=8), intent(inout) ::  a(n,n),x(n,n),eigvals(n)
	real(kind=8) ::  tol, b2, bar
	real(kind=8) ::  beta, coeff, c, s, cs, sc, sort_min
	
	parameter(nmax=10000)
	parameter(tol=1.0d-10)
	
	Niter = 0
	x = 0.0
	do i=1,n
	  x(i,i) = 1.0
	end do

	b2 = 0.0
	do i=1,n
	  do j=1,n
	    if (i.ne.j) b2 = b2 + a(i,j)**2
	  end do
	end do
	
	if (b2 <= tol) return
	bar = 0.5*b2/float(n*n)
	
	do while (b2.gt.tol)
		Niter = Niter + 1
	  do i=1,n-1
	    do j=i+1,n
	      if (a(j,i)**2 <= bar) cycle 
	      b2 = b2 - 2.0*a(j,i)**2
	      bar = 0.5*b2/float(n*n)
	      
        beta = (a(j,j)-a(i,i))/(2.0*a(j,i))
	      coeff = 0.5*beta/sqrt(1.0+beta**2)
	      s = sqrt(max(0.5+coeff,0.0))
	      c = sqrt(max(0.5-coeff,0.0))
	      
        do k=1,n
	        cs =  c*a(i,k)+s*a(j,k)
	        sc = -s*a(i,k)+c*a(j,k)
	        a(i,k) = cs
	        a(j,k) = sc
	      end do
	      
        do k=1,n
	        cs =  c*a(k,i)+s*a(k,j)
	        sc = -s*a(k,i)+c*a(k,j)
	        a(k,i) = cs
	        a(k,j) = sc
	        cs =  c*x(k,i)+s*x(k,j)
	        sc = -s*x(k,i)+c*x(k,j)
	        x(k,i) = cs
	        x(k,j) = sc
	      end do
	    end do
	  end do
	end do
    
  do i = 1, n
  	eigvals(i) = a(i, i)
  end do

  ! sorting here
  do i = 1, Nsort
  	sort_min = eigvals(i)
  	sort_id = i
  	do j = i+1, n
  		if (eigvals(j) < sort_min) then
  			sort_min = eigvals(j)
  			sort_id = j
  		end if
  	end do

  	sort_tmp(1) = eigvals(i)
  	eigvals(i) = eigvals(sort_id)
  	eigvals(sort_id) = sort_tmp(1)

  	sort_tmp = x(:, i)
  	x(:, i) = x(:, sort_id)
  	x(:, sort_id) = sort_tmp
  end do   
	return
	!}
	end subroutine

end module linalg
