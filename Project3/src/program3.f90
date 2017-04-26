program program3
use globals
use bodies_class
implicit none

type(solar_system)             :: orbits
character(len=23), allocatable :: names(:)
character(len=23)              :: ini_arg
integer                        :: Nbody, Nsteps, Nargs, i, iostat, method
real(DP)                       :: runtime, toler
namelist /inputs/ method, toler, runtime, Nsteps, Nbody

open(10,file="param.in")
read(10,nml=inputs)

allocate(names(Nbody))
do i = 1, Nbody
    call get_command_argument(i, ini_arg, status=iostat)
    names(i) = ini_arg
end do

call orbits%init(names)
print *, " INITIALIZE SUCCESS"
call orbits%solve(runtime, toler, Nsteps, method)
print *, " SOLVE SUCCESS"

deallocate(names)
call orbits%deinit()

end program program3 
