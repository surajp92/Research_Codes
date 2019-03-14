!-------------------------------------------------------------------------------
! This is a serial version of Advection-diffusion equation to validate 
! the MPI paralle version
!-------------------------------------------------------------------------------

program advection_diffusion

implicit none

! Computational parameters related to advection-diffusion equation
real*8, parameter				:: c = 1.0			! wave speed
real*8, parameter				:: alpha = 0.1		! diffusion coefficient 
real*8, parameter				:: cfl = 0.1		! CFL condition for stability
real*8, parameter				:: t_norm = 1.0		! Normalized time
real*8  					    :: kappa 			! wave number
real*8, dimension(:,:), allocatable	   	:: u_old
real*8, dimension(:,:), allocatable	   	:: u_new
real*8, dimension(:,:), allocatable	   	:: u_exact
real*8, dimension(:), allocatable	   	:: x, t
real*8, parameter				:: pi = 3.14159265358979323846264338D0
real*8, parameter				:: x_min = 0.0, x_max = 2.0*pi
integer, parameter				:: nx = 64
real*8						    :: dx, dt, x_loc_min
integer						    :: i, k 			! space and time index
integer						    :: time_steps 		! total number of time steps from 0 to T = 1.0
real*8						    :: phi				! phase angle
real*8						    :: start_time, stop_time, time 
real*8						    :: error_sum, avg_error
integer, parameter				:: angle_num = 10 	! number of randoom phase angle s for ensamble averaging
real*8, dimension(angle_num)	:: rand_angle, phis	! array for storing generated random numbers and phase angles
integer							:: angle_par	! index for phis

common /wavenumber/ kappa
kappa = 2.0
error_sum = 0.0

call Random_number(rand_angle)
! print *, rand_angle

do angle_par = 1, angle_num
	
	phis(angle_par) = 20.0*(rand_angle(angle_par) - 0.5)*2.0
	
end do

! grid size
dx 			= (x_max - x_min)/(nx)
dt 			= cfl*(dx**2)/alpha
time_steps 	= int(t_norm*2.0*pi/(abs(c)*dt))
! time_steps = 1


!Plot field
open(10,file='field_final.dat')


do angle_par = 1, angle_num
	phi 		= phis(angle_par)	! only one phase angle is considered. 

! allocate the local array for field variables and grid positions
	allocate(      u_old(0:nx, 0:time_steps))
	allocate(      u_new(0:nx, 0:time_steps))
	allocate(    u_exact(0:nx, 0:time_steps))
	allocate(          x(0:nx))
	allocate(    t(0:time_steps))

! allocate grid position
	do i = 0,nx
		x(i) 	  = x_min + dx*(i)		! position of x for each grid
	end do
	
	
	! print *, 'Rank = ', rank, 'X', x
	call initialize( x, nx, time_steps, u_old, u_exact, phi)
	
	call cpu_time(start_time)
	
	t(0) = 0.0
	
!------------------------------------------------------------------
! start time loop 
	do k = 1,time_steps
!------------------------------------------------------------------
! calculate exact solution at kth time step needed for left and right boundary
		time = k*dt
		t(k) = time
		do i = 0,nx
			u_exact(i,k) = exp(-time*alpha*kappa*kappa)*sin(kappa*(x(i)-c*time) + pi*phi/180.0)
		end do	
	
!----------------------------------------------------------------------	
! do calculation for all other points in the local domain
! skip points for i = 0 and i = nx+1 which are calculated by neigbouring domains
! and is exchanged by send/ receive	
	do i = 0,nx
			! update left boundary i = 0  
			if (i == 0) then
				u_new(i,k) = u_exact(i,k)
	        
	        ! update the right boundary i = nx
	        elseif (i == nx) then
	            u_new(nx,k) = u_exact(nx,k)
	        
	        ! internal points    
	        else
				u_new(i,k) = u_old(i,k-1) - c*dt*(u_old(i+1,k-1) - u_old(i-1,k-1))/(2*dx) &
						 + alpha*dt*(u_old(i+1,k-1) - 2.0*u_old(i,k-1) + u_old(i-1,k-1))/(dx*dx)
			end if
		end do
	 
!------------------------------------------------------------------    
! update old field to new field 
		do i = 0,nx
			u_old(i,k) = u_new(i,k)
	    end do
	end do	
!------------------------------------------------------------------
! calculate the local error between numerical solution and exact solution    
    do i = 0,nx
        error_sum = error_sum + abs(u_old(i,time_steps) - u_exact(i,time_steps)) 
    end do
    
    write(10,*) 'variables ="x","u","ue","error"'
    write(10,*)'zone T=',angle_par
	do i=0,nx
	! write(10,*) x(i),",",u_old(i, time_steps),",",u_exact(i, time_steps),",",u_old(i, time_steps)-u_exact(i, time_steps)
	write(10,*) x(i),u_old(i, time_steps),u_exact(i, time_steps),u_old(i, time_steps)-u_exact(i, time_steps)
	end do

    
	deallocate(u_old, u_new, u_exact)
	deallocate(x)
	deallocate(t)
end do

close(10)
call cpu_time(stop_time)
    
    avg_error = error_sum/(nx*angle_num)
    print *, 'Nx', nx
	print *,'Average error', avg_error
	print *, 'Time', stop_time-start_time

end program advection_diffusion

!------------------------------------------------------------------
! Subroutine to initialize the initial solution = exact solution 
!------------------------------------------------------------------
subroutine initialize(x, nx, time_steps, u_old, u_exact, phi)

implicit none

integer, intent(in)                         :: nx, time_steps
real*8, intent(in), dimension(0:nx)		    :: x
real*8, intent(out), dimension(0:nx, 0:time_steps)  :: u_old
real*8, intent(out), dimension(0:nx, 0:time_steps)  :: u_exact
real*8							    :: phi, kappa

integer							    :: i, k
real*8, parameter                   :: pi = 3.141592653589793238

common /wavenumber/ kappa

do i = 0,nx
	u_exact(i,0) = sin(kappa*x(i) + pi*phi/180)
	u_old(i,0)	= 0.0
	! u_old(i,0)   = u_exact(i,0)
end do

!Plot field
open(20,file='initial_field.dat')
write(20,*) 'variables ="x","u","ue","error"'
do i=0,nx
write(20,*) x(i),",",u_old(i, 0),",",u_exact(i, 0),",",u_old(i, 0)-u_exact(i, 0)
end do
close(20)

end subroutine initialize
