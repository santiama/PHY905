module math
implicit none

contains 
function f(x)
    implicit none
    real(kind=8) :: f
    real(kind=8), intent(in) :: x
    f = 100.d0*exp(-10.d0*x)
end function f

function exact(x)
    implicit none 
    real(kind=8) :: exact
    real(kind=8), intent(in) :: x
    exact = 1.d0-(1.d0-exp(-10.d0))*x-exp(-10.d0*x)
end function exact

end module math
