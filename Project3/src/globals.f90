module globals
  implicit none
  integer,   parameter :: DP = kind(0.0d0)           
  integer,   parameter :: LP = selected_int_kind(18) 
  character(len=100)   :: fmt
  character, parameter :: CR = char(13)              
  real(DP), parameter  :: PI = 4.0d0*atan(1.0d0)
  real(DP), parameter  :: G = 4.0d0*PI*PI
end module globals
