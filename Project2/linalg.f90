module linalg
    implicit none

contains

    subroutine tridiag_general(a, b, c, v_sol, rhs, N)
    !{============================================
    !  Solves the general tridiagonal setup where         
    !  a_i and c_i do not need to retain the same 
    !  values down their diagonals
    ! ============================================
        implicit none
        integer(kind=4) :: i
        integer(kind=4), intent(in) :: N
        real(kind=8), dimension(N), intent(inout) :: a, b, c, rhs
        real(kind=8), dimension(N), intent(inout) :: v_sol
        real(kind=8), dimension(N) :: rhs_new 
        real(kind=8) :: m 

        ! FORWARD SUBSTITUTION
        rhs_new(1) = rhs(1) 
        do i=2, N
            m = a(i)/b(i-1) 
            b(i) = b(i)-m*c(i-1)
            rhs_new(i) = rhs(i) - m*rhs_new(i-1)
        end do
        
        !BACKWARD SUB STEP
        v_sol(N) = rhs_new(N)/b(N)
        do i = N-1, 1, -1
           v_sol(i) = (rhs_new(i)-c(i)*v_sol(i+1))/b(i)
        end do
        return
    !}
    end subroutine tridiag_general

    subroutine tridiag_special(a, b, c, v_sol, rhs, N)
    !{=================================================
    !  Solves the specialized tridiagonal setup where         
    !  a_i and c_i have the same value down their entire 
    !  respective diagonals 
    ! ==================================================
        implicit none
        integer(kind=4) :: i
        integer(kind=4), intent(in) :: N
        real(kind=8), dimension(N), intent(inout) :: a, b, c, rhs
        real(kind=8), dimension(N), intent(inout) :: v_sol
        real(kind=8), dimension(N) :: rhs_new 

        ! FORWARD SUBSTITUTION
        rhs_new(1) = rhs(1) 
        do i=2, N
            b(i) = b(i)-1.d0/b(i-1)
            rhs_new(i) = rhs(i) + rhs_new(i-1)/b(i-1)
        end do
        
        !BACKWARD SUB STEP
        v_sol(N) = rhs_new(N)/b(N)
        do i = N-1, 1, -1
           v_sol(i) = (rhs_new(i)+v_sol(i+1))/b(i)
        end do
        return
    !}
    end subroutine tridiag_special
    
    SUBROUTINE lu_decompose(a,n,indx,d)
    !{    Given an NxN matrix A(N,N), this routine replaces it by the LU 
    !     decomposed one, where the matrix elements are stored in the same 
    !     matrix A. The array indx is  an output vector which records the row
    !     permutation effected by the partial pivoting. d is the determinant
    !
      IMPLICIT NONE
      INTEGER :: n, i, j, k, imax
      REAL(kind=8) :: sum , tiny, aamax, dum, d
      REAL(kind=8), DIMENSION(n,n) :: a
      INTEGER, DIMENSION(n) :: indx
      REAL(kind=8), ALLOCATABLE :: vv(:)
    
      tiny=1.0e-20
      ALLOCATE ( vv(n) )
      D=1.
      DO i=1,n
         aamax=0.
         DO j=1,n
            IF (ABS(a(i,j)) > aamax) aamax=ABS(a(i,j))
         ENDDO
         !     Zero is the largest element
         IF (aamax == 0.) STOP 'Singular matrix.'
         !     No nonzero largest element
         vv(i)=1./aamax
      ENDDO
      !     loop over columns
      DO j=1,n
         !     solves equation 2.3.12 except for i=j of Numerical Recipes
         IF (j > 1) THEN
            DO i=1,j-1
               sum=a(i,j)
               IF (i > 1)THEN
                  DO k=1,i-1
                     sum=sum-a(i,k)*a(k,j)
                  ENDDO
                  a(i,j)=sum
               ENDIF
            ENDDO
         ENDIF
         !    start searching for largest pivot element
         aamax=0.
         DO i=j,n
            sum=a(i,j)
            IF (j > 1)THEN
               DO k=1,j-1
                  sum=sum-a(i,k)*a(k,j)
               ENDDO
               a(i,j)=sum
            ENDIF
            dum=vv(i)*ABS(sum)
            IF (dum >= aamax) THEN
               imax=i
               aamax=dum
            ENDIF
         ENDDO
         !    interchange of rows
         IF (j /= imax)THEN
            DO k=1,n
               dum=a(imax,k)
               a(imax,k)=a(j,k)
               a(j,k)=dum
            ENDDO
            !    change of parity for determinant
            d=-d
            vv(imax)=vv(j)
         ENDIF
         indx(j)=imax
         IF(j /= n) THEN
            IF(a(j,j) == 0.) a(j,j)=tiny
            dum=1./a(j,j)
            DO i=j+1,n
               a(i,j)=a(i,j)*dum
            ENDDO
         ENDIF
         !    set up determinant
         d=d*a(j,j)
      ENDDO
      IF(a(n,n) == 0.)  a(n,n)=tiny
      DEALLOCATE ( vv)
    !}
    END SUBROUTINE lu_decompose
    
    SUBROUTINE lu_linear_equation(a,n,indx,b)
    !{    Solves set of linear equations Ax=b, A is input as an LU decompomsed
    !     matrix and indx keeps track of the permutations of the rows. b is input
    !     as the right-hand side vector b and returns the solution x. A, n and indx
    !     are not modified by this routine. This function takes into that b can contain
    !     many zeros and is therefore suitable for matrix inversion
      IMPLICIT NONE
      INTEGER :: n, ii, ll, i, j
      REAL(kind=8) :: sum 
      REAL(kind=8), DIMENSION(n,n) :: a
      REAL(kind=8), DIMENSION(n) :: b
      INTEGER, DIMENSION(n) :: indx
    
      ii=0
      !     First we solve equation 2.3.6 of numerical recipes 
      DO i=1,n
         ll=indx(i)
         sum=b(ll)
         b(ll)=b(i)
         IF (ii /= 0)THEN
            DO j=ii,i-1
               sum=sum-a(i,j)*b(j)
            ENDDO
         ELSEIF (sum /= 0.) THEN
            ii=i
         ENDIF
         b(i)=sum
      ENDDO
      !     then we solve equation 2.3.7
      DO i=n,1,-1
         sum=b(i)
         IF (i < n) THEN
            DO j=i+1,n
               sum=sum-a(i,j)*b(j)
            ENDDO
         ENDIF
         !     store a component of the solution x in the same place as b
         b(i)=sum/a(i,i)
      ENDDO
    !}
    END SUBROUTINE lu_linear_equation

end module linalg
