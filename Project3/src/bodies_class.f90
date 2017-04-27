module bodies_class
use globals
implicit none
private

type :: particle 
    real(DP) :: pos(3), vel(3)
end type particle

type, extends(particle) :: body
    real(DP) :: mass
end type body

type, public :: solar_system 
!{
    type(body), allocatable :: bodies(:)
    integer :: Nbody
contains
    private
    procedure, public :: init
    procedure, public :: verlet_solver
    procedure, public :: euler_solver
    procedure, public :: rk4_solver
    procedure, public :: compute_accel
    procedure, public :: get_pos
    procedure, public :: getPE
    procedure, public :: getKE
    procedure, public :: getE
    procedure, public :: solve
    procedure, public :: deinit
!}
end type solar_system 

contains
subroutine init(this, names)
!{
    class(solar_system), intent(inout) :: this
    character(len=*),    intent(in)    :: names(:)
    character(len=23) :: body_in
    real(DP)          :: mass_in, pos_in(3), vel_in(3)
    integer           :: i, unit, iostat

    this%Nbody = size(names)
    allocate(this%bodies(this%Nbody))
    open(newunit=unit, file='initial.dat', status='old', action='read', iostat=iostat)

    do i = 1, this%Nbody
        do
            read(unit, '(A23, 7ES23.15)', iostat=iostat) body_in, pos_in, vel_in, mass_in
            if (trim(adjustl(body_in)) == trim(adjustl(names(i)))) exit
        end do
        this%bodies(i)%mass   = mass_in
        this%bodies(i)%pos(:) = pos_in
        this%bodies(i)%vel(:) = vel_in
        rewind(unit)
    end do
    close(unit)
!}
end subroutine init 

function getKE(this) result(kinE)
!{
    class(solar_system), intent(in) :: this
    real(DP)                  :: kinE(this%Nbody)

    kinE = 0.5*this%bodies(:)%mass*(this%bodies(:)%vel(1)*this%bodies(:)%vel(1) + &
                                 this%bodies(:)%vel(2)*this%bodies(:)%vel(2) + &
                                 this%bodies(:)%vel(3)*this%bodies(:)%vel(3))
!}
end function getKE

function getPE(this) result(potE)
!{
    class(solar_system), intent(inout) :: this
    real(DP) :: potE(this%Nbody), tempPE(this%Nbody)
    real(DP) :: disp(this%Nbody, 3)
    integer  :: i, n

    n = this%Nbody

    do i = 1, n-1

        disp(i+1:n, 1) = this%bodies(i+1:n)%pos(1) - this%bodies(i)%pos(1)
        disp(i+1:n, 2) = this%bodies(i+1:n)%pos(2) - this%bodies(i)%pos(2)
        disp(i+1:n, 3) = this%bodies(i+1:n)%pos(3) - this%bodies(i)%pos(3)

        tempPE(i+1:n) = (G*this%bodies(i)%mass)*this%bodies(i+1:n)%mass &
                                         /sqrt(disp(i+1:n,1)*disp(i+1:n,1) +&
                                               disp(i+1:n,2)*disp(i+1:n,2) +&
                                               disp(i+1:n,3)*disp(i+1:n,3))

        potE(i) = potE(i) - sum(tempPE(i+1:n))
        potE(i+1:n) = potE(i+1:n) - tempPE(i+1:n)
    end do
!}
end function getPE

function get_pos(this) result(positions)
!{
    class(solar_system), intent(inout) :: this
    real(DP) :: positions(this%Nbody, 3)
    integer :: i

    do i=1, this%Nbody
      positions(i,:) = this%bodies(i)%pos
    end do
!}
end function get_pos

subroutine compute_accel(this, pos, acc)
!{
    class(solar_system), intent(in) :: this
    real(DP),      intent(in)    :: pos(this%Nbody, 3)
    real(DP),      intent(out)   :: acc(this%Nbody, 3)
    real(DP) :: tmp(this%Nbody)
    real(DP) :: dists(this%Nbody, 3)
    integer  :: i, n
    n = this%Nbody
    accelerations = 0.0d0

    do i = 1, n-1
        dists(i+1:n, :) = pos(i, :) - pos(i+1:n, :)
        tmp(i+1:n) = dists(i+1:n, 1)**2 + dists(i+1:n,2)**2 + dists(i+1:n,3)**2
        tmp(i+1:n) = G/(tmp(i+1:n)*sqrt(tmp(i+1:n)))
        dists(i+1:n, :) = tmp(i+1:n)*dists(i+1:n, :)
        acc(i, :) = acc(i, :) - sum(this%bodies(i+1:n)%mass*dists(i+1:n, :))
        acc(i+1:n, :) = acc(i+1:n, :) + this%bodies(i)%mass*dists(i+1:n, :)
    end do
!}
end subroutine compute_accel

function getE(this) result(ene)
!{
    class(solar_system), intent(inout) :: this
    real(DP)                     :: ene(this%Nbody)

    ene = this%getKE() + this%getPE()
!}
end function getE

subroutine solve(this, total_time, tolerance, Nsteps, method)
!{
    class(solar_system), intent(inout) :: this
    real(DP)             :: dt
    real(DP), intent(in) :: tolerance, total_time
    integer, intent(in)  :: Nsteps, method !euler=0 verlet=1 rk4=2
    integer              :: unit, i

    dt = dble(total_time/Nsteps)
    open(newunit=unit, file='orbits.dat', status='replace', action='write')
    if (method==1) then
      call this%verlet_solver(dt, Nsteps, unit, fmt)
    else if (method==0) then
      call this%euler_solver(dt, Nsteps, unit, fmt)
    else
      call this%rk4_solver(dt, Nsteps, unit, fmt)
    end if
    close(unit)
!}
end subroutine solve

subroutine euler_solver(this, dt, Nsteps, unit, fmt)
!{
    class(solar_system),      intent(inout) :: this
    real(DP),           intent(in)    :: dt
    integer,            intent(in)    :: unit, Nsteps
    character(len=100), intent(in)    :: fmt
    real(DP) :: accelerations(this%Nbody, 3)
    integer  :: i, j, k

    write(unit, fmt) 0.d0, (this%bodies(j)%pos(1), this%bodies(j)%pos(2), this%bodies(j)%pos(3), j = 1, this%Nbody)

    do i=1,Nsteps

      call this%compute_accel(this%get_pos(), accelerations)
      do k=1,this%Nbody
        this%bodies(k)%pos(:) = this%bodies(k)%pos(:) + dt*this%bodies(k)%vel(:)
        this%bodies(k)%vel(:) = this%bodies(k)%vel(:) + dt*accelerations(k,:)
      end do

      write(unit, fmt) i*dt, (this%bodies(j)%pos(1), this%bodies(j)%pos(2), this%bodies(j)%pos(3), j = 1, this%Nbody)
    end do
!}
end subroutine euler_solver

subroutine rk4_solver(this, dt, Nsteps, unit, fmt)
!{
    class(solar_system), intent(inout) :: this
    real(DP),            intent(in)    :: dt
    integer,             intent(in)    :: unit, Nsteps
    character(len=100),  intent(in)    :: fmt
    integer                            :: i, j, k

    real(DP) :: k1v(this%Nbody, 3), k1p(this%Nbody, 3), k2v(this%Nbody, 3), k2p(this%Nbody, 3)
    real(DP) :: k3v(this%Nbody, 3), k3p(this%Nbody, 3), k4v(this%Nbody, 3), k4p(this%Nbody, 3)

    do i=1,Nsteps

      call this%compute_accel(this%get_pos(), k1v)
      k1p = this%vel
      call this%compute_accel(this%get_pos()+k1p*0.5*dt, k2v)
      k2p = this%vel+k1v*0.5*dt
      call this%compute_accel(this%get_pos()+k2p*0.5*dt, k3v)
      k3p = this%vel+k2v*0.5*dt
      call this%compute_accel(this%get_pos()+k3p*dt, k4v)
      k4p = this%vel+k3v*dt 

      this%bodies(k)%pos = this%bodies(k)%pos + (dt/6.d0)*(k1v+2*k2v+2*k3v+k4v)
      this%bodies(k)%vel = this%bodies(k)%vel + (dt/6.d0)*(k1p+2*k2p+2*k3p+k4p)

      write(unit, fmt) i*dt, (this%bodies(j)%pos(1), this%bodies(j)%pos(2), this%bodies(j)%pos(3), j = 1, this%Nbody)
    end do
!}
end subroutine rk4_solver

subroutine verlet_solver(this, dt, Nsteps, unit, fmt)
!{
    class(solar_system),      intent(inout) :: this
    real(DP),                 intent(in)    :: dt
    integer,                  intent(in)    :: unit, Nsteps
    character(len=100),       intent(in)    :: fmt
    real(DP) :: accelerations(this%Nbody, 3), pos(this%Nbody, 3)
    integer  :: i, j, k

    pos = this%get_pos()
    call this%compute_accel(pos, accelerations)

    write(unit, fmt) 0.d0, (this%bodies(j)%pos(1), this%bodies(j)%pos(2), this%bodies(j)%pos(3), j = 1, this%Nbody)

    do i=1,Nsteps
      do k=1,this%Nbody
        this%bodies(k)%vel(:) = this%bodies(k)%vel(:) + 0.5*dt*accelerations(k,:)
        pos(k,:) = this%bodies(k)%pos(:) + dt*this%bodies(k)%vel(:)
      end do

      call this%compute_accel(pos, accelerations)
      
      do k=1,this%Nbody
        this%bodies(k)%vel(:) = this%bodies(k)%vel(:) + 0.5*dt*accelerations(k,:)
        this%bodies(k)%pos(:) = pos(k,:)
      end do

      write(unit, fmt) i*dt, (this%bodies(j)%pos(1), this%bodies(j)%pos(2), this%bodies(j)%pos(3), j = 1, this%Nbody)
    end do
!}
end subroutine verlet_solver

subroutine deinit(this)
!{
    class(solar_system), intent(inout) :: this
    deallocate(this%bodies)
!}
end subroutine deinit

end module bodies_class
