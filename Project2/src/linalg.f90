module linalg
    implicit none

contains
  subroutine jacobi_classical(a,n,eigvals,x,Niter,Nsort)
  !{ Classical Jacobi method. Finds pivots via max off-diagonal at each iteration. 
  !  Converges in fewer iterations but generally is slower than the cyclic method

  implicit none
  integer :: nmax, i, j, k, l, sort_id
  integer(kind=4), intent(inout) ::  n, Nsort, Niter
  real(kind=8), intent(inout) ::  a(n,n),x(n,n),eigvals(n)
  real(kind=8) ::  beta, tol, maxoffdiag
  real(kind=8) ::  c, s, t, sort_min, sort_tmp(n)
	real(kind=8) ::  a_ll, a_kk, a_ik, a_il, e_ik, e_il 
  parameter(nmax=50000)
  parameter(tol=1.0d-5)

	Niter = 0
  k = 0
  l = 0
	x = 0.0
	do i=1,n
	  x(i,i) = 1.0
	end do
  
  do
    call offdiag_max(a,n,k,l,maxoffdiag)
    if (Niter > nmax) then
      exit
    endif

    if (maxoffdiag**2 <= tol) then
      exit 
    end if

		if(a(k,l)/=0.d0) then
			maxoffdiag=0.d0
			beta=(a(l,l)-a(k,k))/(2.d0*a(k,l))
			if (beta>= 0.d0) then
				t = 1.d0/(beta+dsqrt(1.d0+beta**2))
			else
				t =-1.d0/(-beta+dsqrt(1.d0+beta**2))
			end if
			c=1.d0/dsqrt(1.d0+t**2) 
			s=c*t
		else
			c=1.d0
			s=0.d0
		end if

		a_kk=a(k,k)
		a_ll=a(l,l)
		a(k,k)=c**2*a_kk-2.d0*c*s*a(k,l)+s**2*a_ll
		a(l,l)=s**2*a_kk+2.d0*c*s*a(k,l)+c**2*a_ll
		a(k,l)=0.d0
		a(l,k)=0.d0

		do i=1,n-1
			if((i/=k).and.(i/=l))then
				a_ik=a(i,k)
        a_il=a(i,l)
				a(i,k)=c*a_ik-s*a_il
				a(k,i)=a(i,k)
				a(i,l)=s*a_ik+c*a_il
				a(l,i)=a(i,l)
			end if
			e_ik=x(k,i)
			e_il=x(l,i)
			x(k,i)=c*e_ik-s*e_il
			x(l,i)=c*e_il+s*e_ik
		end do
		Niter=Niter+1
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
  end subroutine jacobi_classical

	subroutine jacobi_cyclic(a,n,eigvals,x,Niter,Nsort)
  !{ Cyclic jacobi method. Instead of finding max off diagonal, 
  !  simply goes row by row and convergence criteria is based off
  !  sum of off diagonals instead

	implicit none
	integer :: nmax, i, j, k, sort_id
	integer(kind=4), intent(inout) ::  n, Nsort, Niter
  real(kind=8) :: sort_tmp(n)
	real(kind=8), intent(inout) ::  a(n,n),x(n,n),eigvals(n)
	real(kind=8) :: tol, b2, bar
	real(kind=8) :: beta, t, c, s, sort_min
  real(kind=8) :: a_ik, a_jk, a_ki, a_kj, e_ki, e_kj
	
	parameter(nmax=10000)
	parameter(tol=1.0d-10)
	
	Niter = 0
	x = 0.0
	do i=1,n
	  x(i,i) = 1.0
	end do
  
  ! cyclic method uses SUM of off diagonals not max. Goes row-by-row
	b2 = 0.0
	do i=1,n
	  do j=1,n
	    if (i.ne.j) b2 = b2 + a(i,j)**2
	  end do
	end do
	
	bar = 0.5*b2/float(n*n)
	do while (b2.gt.tol)
		Niter = Niter + 1
	  do i=1,n-1
	    do j=i+1,n
	      if (a(j,i)**2 <= bar) cycle 
	      b2 = b2 - 2.0*a(j,i)**2
	      bar = 0.5*b2/float(n*n)
	      
        beta = (a(j,j)-a(i,i))/(2.0*a(j,i))
	      t = 0.5*beta/sqrt(1.0+beta**2)
	      s = sqrt(max(0.5+t,0.0))
	      c = sqrt(max(0.5-t,0.0))
	      
        do k=1,n
          a_ik = a(i,k)
          a_jk = a(j,k) 
	        a(i,k) =  c*a_ik+s*a_jk
	        a(j,k) = -s*a_ik+c*a_jk
	      end do
	      
        do k=1,n
          a_ki = a(k,i)
          a_kj = a(k,j)
	        a(k,i) =  c*a_ki+s*a_kj
	        a(k,j) = -s*a_ki+c*a_kj

          e_ki = x(k,i)
          e_kj = x(k,j)
	        x(k,i) =  c*e_ki+s*e_kj
	        x(k,j) = -s*e_ki+c*e_kj
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

  subroutine offdiag_max(a,n,k,l,maxoff)
  !{ compute offdiagonal maximum
    implicit none
	  integer :: i, j
    integer(kind=4), intent(inout) :: k,l
    integer(kind=4), intent(in) :: n
    real(kind=8) :: tmp
    real(kind=8), intent(inout) :: a(n,n), maxoff
    maxoff = 0.d0
    do j=2,n
      do i=1,j-1
        tmp = dabs(a(i,j))
        if (tmp > maxoff) then
          l = j
          k = i
          maxoff = tmp
        end if
      end do
    end do
  !}
  end subroutine offdiag_max

end module linalg
