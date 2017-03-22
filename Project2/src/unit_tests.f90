module unit_tests
  use linalg 
  implicit none
	integer(kind=4) :: i,j,k,l,n,nit,nmax_test
  real(kind=8) :: maxoff, a(3,3), eigvals(3), tol_test
	character(len=20) :: fmt
	real(kind=8) :: x(3,3)
contains

subroutine init_tests()
  implicit none
  nmax_test = 1000000
  tol_test = 1.d-9
  a(1,1)=1.0
  a(1,2)=2.0
  a(1,3)=3.0
  a(2,1)=2.0
  a(2,2)=2.0
  a(2,3)=-2.0
  a(3,1)=3.0
  a(3,2)=-2.0
  a(3,3)=4.0
  n = 3
end subroutine init_tests

subroutine max_test()
  implicit none
  call offdiag_max(a,n,k,l,maxoff)
  print*, "Max off-diag test (answer=3) : Computed = ",maxoff
end subroutine max_test

subroutine gram_test()
  implicit none
  call jacobi_classical(a,n,eigvals,x,nit,n,nmax_test,tol_test)
  print *, "Orthogonality Testing : "
	do i = 1, n
		do j = i+1, n
			write(*,'(I4, A, I4, A, G15.6)')i,'  dot ', j,'  = ', dot_product(x(:, i), x(:, j))
		end do
	end do
end subroutine gram_test 
end module
