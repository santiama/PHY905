program argon
  use plot
  use randomseed
  implicit none

  real(kind=8), parameter :: T_IC = 1.000D0                 
  real(kind=8), parameter :: rho = .88D0                     
  real(kind=8), parameter :: r_verlet = 4.5D0               
  real(kind=8), parameter :: r_verlet = 4.5D0                
  real(kind=8), parameter :: r_cut =  2.5D0                 
  real(kind=8), parameter :: dt = 0.004                      
  integer, parameter :: N_cell = 4                         
  integer, parameter :: N_bin = 5000                    
  integer, parameter :: time_end = 3000                 
  integer, parameter :: hist_time = 100                 
  integer, parameter :: eq_time = 500                   
  integer :: time, hist_runs = 0      
  integer, parameter :: N = 4*N_cell**3                 
  integer, parameter :: pairs_max = N*(N-1)/2           
  real(kind=8), parameter :: L_box = (N/rho)**(1.0D0/3.0D0)  
  real(kind=8), parameter :: L_cell = L_box/N_cell           
  real(kind=8), parameter :: pi = 4*atan(1.0D0)              
  real(kind=8), dimension(3,N) :: pos, vel, acc, tot_dis=0   
  real(kind=8), dimension(time_end) :: T, U, kin_energy,,disp_sq    
  real(kind=8) :: histogram(N_bin) = 0                       
  real(kind=8) :: T_scaling_factor                           
  real(kind=8) :: msq_vel                                    
  integer :: pairs(2,pairs_max)                         
  
  call initialize_system
  call plot_init(L_box)
  do time = 1, time_end
    if (modulo(time,50) == 0) call calculate_pairs      
    call update_pos_vel_acc
  end do
  call plot_close()
  call print_data
 contains

subroutine initialize_system
  time = 1
  call init_random_seed()                   
  call IC_vel                               
  call calculate_pairs                      
  call update_T                             
  call update_acc                           
end subroutine

subroutine IC_pos 
  integer :: x, y, z, pc 
  real(kind=8) :: fcc(3,4), origin(3)
  fcc(:,1) = [0.0D0, 0.0D0, 0.0D0]
  fcc(:,2) = [0.0D0, 0.5D0, 0.5D0]
  fcc(:,3) = [0.5D0, 0.0D0, 0.5D0]
  fcc(:,4) = [0.5D0, 0.5D0, 0.0D0]
  fcc = fcc*L_cell
  pc = 0
  do x = 1, N_cell
    do y = 1, N_cell
      do z = 1, N_cell
        origin = [x-1.D0, y-1.D0, z-1.D0]*L_cell
        pc = pc + 1
        pos (:,pc) = origin + fcc(:,1)
        pc = pc + 1
        pos (:,pc) = origin + fcc(:,2)
        pc = pc + 1
        pos (:,pc) = origin + fcc(:,3)
        pc = pc + 1
        pos (:,pc) = origin + fcc(:,4)
      end do
    end do
  end do
end subroutine IC_pos

subroutine random_gauss_routine(random_gauss)
  real(kind=8) :: v1, v2
  real(kind=8), intent(out) :: random_gauss
  call RANDOM_NUMBER(v1)
  call RANDOM_NUMBER(v2)
  random_gauss = sqrt(-2.0*log(v1))*sin(2.0D0*pi*v2)
end subroutine

subroutine IC_vel
  real(kind=8) :: avg_vel(3)
  integer :: i
  avg_vel = 0
  do i = 1,N        
    call random_gauss_routine(vel(1,i))
    call random_gauss_routine(vel(2,i))
    call random_gauss_routine(vel(3,i))
    avg_vel(:) = avg_vel(:) + vel(:,i)
  end do 
  avg_vel = avg_vel/N
  !subracting the average velocity so the center of mass does not move
  do i = 1, N
    vel(:,i) = vel(:,i) - avg_vel(:)
  end do  
end subroutine

subroutine calculate_pairs
  integer :: i, j, counter
  real(kind=8) :: r_sq, r(3)
  pairs = 0
  counter = 1
  do i = 1, N-1
    do j = i+1,N
      r = pos(:,i) - pos(:,j)
      r = r - Nint(r/L_box) * L_box
      r_sq = dot_product(r,r)
      if (r_sq < r_verlet**2) then
        pairs(1,counter) = i
        pairs(2,counter) = j
        counter = counter + 1
      end if
    end do
  end do
end subroutine

subroutine update_acc
  integer :: pc, i, j
  real(kind=8) :: r_sq, F(3), r(3)
  acc = 0
  U(time) = 0
  kin_energy(time) = 0.5*msq_vel
  do pc = 1, pairs_max
    i = pairs(1,pc)
    j = pairs(2,pc)
    if (i /= 0 .and. j /= 0) then 
      r = pos(:,i) - pos(:,j)
      r = r - Nint(r/L_box) * L_box      
      r_sq = dot_product(r,r)
      if (r_sq < r_cut**2 ) then
        F = 4*(12/r_sq**7-6/r_sq**4)*r
        acc(:,i) = acc(:,i) + F
        acc(:,j) = acc(:,j) - F
        P(time) = P(time) + dot_product(F,r)
      end if
      U(time) = U(time) + 4*(1/r_sq**6-1/r_sq**3)
      call make_histogram(N_bin, sqrt(r_sq), 0.D0, r_verlet, pc)
    end if
  end do
end subroutine

subroutine update_pos_vel_acc
   vel = vel + 0.5*acc*dt
   pos = modulo(pos + vel*dt, L_box)
   if (time >= eq_time) then 
    tot_dis = tot_dis + vel*dt
    disp_sq(time) = sum(tot_dis*tot_dis)/N
   end if
   call update_T
   call update_acc
   vel = vel + 0.5*acc*dt
end subroutine

subroutine update_T
  msq_vel = sum(vel*vel)
  T(time) = msq_vel/(3.D0*N)
  T_scaling_factor = sqrt(T_IC/T(time))
  if ((mod(time,20) == 0) .and. time < eq_time-100 ) vel = vel * T_scaling_factor
end subroutine

subroutine make_histogram(N_bin, x, xmin, xmax, pc)
  real(kind=8) :: bin_size, x, xmin, xmax
  integer :: histogram_entry, N_bin, pc
  bin_size = (xmax - xmin) / N_bin
  if (modulo(time,hist_time) == 0 .and. time > eq_time) then
    histogram_entry = ceiling(x/bin_size) 
    if (histogram_entry < N_bin) histogram(histogram_entry) = histogram(histogram_entry) + 1 
    if (pc == 1) hist_runs = hist_runs + 1 
  end if
end subroutine

subroutine print_data
    real(kind=8) :: array(5, size(T(eq_time:time_end)))
    histogram = histogram/hist_runs
    array(1,:) = T(eq_time:time_end)
    array(2,:) = kin_energy(eq_time:time_end)
    array(3,:) = U(eq_time:time_end)
    array(5,:) = disp_sq(eq_time:time_end)
    
    open (unit=10,file="array.dat")
    write (10,"(5F15.5)") array
    open (unit=12,file="histogram.dat")
    write (12,"(1F15.5)") histogram
    open (unit=18,file="data.dat")
    write (18,*) rho, T_IC, eq_time, time_end, dt, r_cut, r_verlet, N_bin, N
    print *, "avg T:        ", avg(T(eq_time:time_end))
    print *, "avg Kin:      ", avg(kin_energy(eq_time:time_end))
end subroutine

real(kind=8) function avg(x) result(average)
  real(kind=8), intent(in) :: x(:)
  average = sum(x) / size(x)
end function avg

real(kind=8) function var(x) result(variance)
  real(kind=8), intent(in) :: x(:)
  variance = sum(x*x) / size(x) - avg(x)**2
end function var
end program argon
