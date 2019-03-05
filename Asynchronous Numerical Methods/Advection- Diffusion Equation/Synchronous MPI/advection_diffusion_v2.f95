!-------------------------------------------------------------------------------
! The code solves the advection diffusion equation using Proxy equation
! approach given in the below reference paper using MPI. The code is
! still under development. 
!
!---------------------------------------------------------------------!
! Reference Paper: 
! Mittal, Ankita, and Sharath Girimaji. "Proxy-equation paradigm:
! A strategy for massively parallel asynchronous computations."
! Physical Review E 96.3 (2017): 033304.
!
!---------------------------------------------------------------------!
!	Author: Suraj Pawar
!-------------------------------------------------------------------------------

program advection_diffusion

implicit none

include 'mpif.h'

! MPI reated variables
integer 	                        :: ierr ! error signal variable, Standard value - 0
integer 	                        :: rank ! process ID (pid) / Number
integer		                        :: nprocs ! number of processors
integer, dimension(MPI_STATUS_SIZE)    	:: status ! status for receive
integer, parameter			:: id_left_to_right = 10 ! tag for send/ receive
integer, parameter			:: id_right_to_left = 20 ! tag for send/ receive
integer, dimension(4)                  	:: req = MPI_REQUEST_NULL
integer, dimension(MPI_STATUS_SIZE, 4)	:: status_array

! Computational parameters related to advection-diffusion equation
real*8, parameter				:: c = 1.0			! wave speed
real*8, parameter				:: alpha = 0.1		! diffusion coefficient 
real*8, parameter				:: cfl = 0.1		! CFL condition for stability
real*8, parameter				:: t_norm = 1.0		! Normalized time
real*8  					:: kappa 			! wave number
real*8, dimension(:,:), allocatable	   	:: u_old
real*8, dimension(:,:), allocatable	   	:: u_new
real*8, dimension(:,:), allocatable	   	:: u_exact
real*8, dimension(:), allocatable	   	:: x, t
real*8, parameter				:: pi = 3.14159265358979323846264338D0
real*8, parameter				:: x_min = 0.0, x_max = 2.0*pi
integer, parameter				:: nx_global = 512
integer						:: nx_local
real*8						:: dx, dt, x_loc_min
integer						:: i, k 
integer						:: time_steps 		! total number of time steps from 0 to T = 1.0
integer, parameter				:: istart = 1, kstart = 1
integer						:: i_global_low, i_global_high
real*8						:: phi				! phase angle
real*8						:: start_time, stop_time, time 
real*8						:: local_error_sum, global_error_sum, avg_global_error
integer, parameter				:: angle_num = 1 	! number of randoom phase angle s for ensamble averaging
real*8, dimension(angle_num)	:: rand_angle, phis	! array for storing generated random numbers and phase angles
integer							:: angle_par	! index for phis

common /wavenumber/ kappa
kappa = 2.0

!------------------------------------------------------------------
! MPI Starting


call MPI_INIT(ierr) ! Initialize MPI

call MPI_COMM_SIZE(MPI_COMM_WORLD, &
                           nprocs, &    ! number of processors
                           ierr     )   ! ierr = 0 if successful

call MPI_COMM_RANK(MPI_COMM_WORLD, &
                             rank, &    ! rank of the processor
                             ierr   )   ! ierr = 0 if successful

if (rank == 0) then
    write(*,*) 'No. of processor = ', nprocs, 'Rank = ', rank
    write(*,*) 'Nx_Global = ', nx_global
    write(*,*) 'Total size = ', nx_global
end if

! if (rank == 0) then

	call Random_number(rand_angle)
	
	do angle_par = 1, angle_num
		phis(angle_par) = 20.0*(rand_angle(angle_par) - 0.5)*2.0
	end do

! end if
!-------------------------------------------------------------------------------
! Allocate x-position to each grid point
! The domain is partitioned in x-direction and there are two global ghost nodes
! i = 0 and i = nx_global + 1
! Each local domain has two ghost nodes i = 0 and i = nx_local + 1
! (i=0)_rank = (i = nx_local)_(rank-1)
! (i = nx+1)_rank = (i = 0)_(rank+1)


! local size
nx_local 			= nx_global/nprocs
local_error_sum 	= 0.0
global_error_sum 	= 0.0
! print *, rank, nx_global, nx_local

do angle_par = 1, angle_num

! grid size
dx 			= (x_max - x_min)/(nx_global)
dt 			= cfl*(dx**2)/alpha

time_steps 	= int(t_norm*2.0*pi/(abs(c)*dt))

phi = phis(angle_par)
! print *, rank, phi

! allocate the local array for field variables and grid positions
allocate(      u_old(0:nx_local+1, 0:time_steps))
allocate(      u_new(0:nx_local+1, 0:time_steps))
allocate(    u_exact(0:nx_local+1, 0:time_steps))
allocate(          x(0:nx_local+1))
allocate(          t(0:time_steps))

! allocate grid position
do i = 0,nx_local+1
	x_loc_min = rank*((x_max - x_min)/nprocs) ! dimesnion in y-direction is 1.0
	x(i) 	  = x_loc_min + dx*(i-1)		! position of x for each grid
end do


! print *, 'Rank = ', rank, 'X', x
call initialize(rank, nprocs, x, nx_local, time_steps, u_old, u_exact, phi)

!------------------------------------------------------------------
! start recording actual execution time

call cpu_time(start_time)


t(0) = 0.0
do k = 1,time_steps

! calculate exact solution at kth time step
	time = k*dt
	t(k) = time
	do i = 0,nx_local+1
		
		u_exact(i,k) = exp(-time*alpha*kappa*kappa)*sin(kappa*(x(i)-c*time) + pi*phi/180.0)
	
	end do	
	
! do calculation for all other points in the local domain
! skip points for i = 0 and i = nx+1 which are calculated by neigbouring domains
! and is exchanged by send/ receive
	
	do i = 1,nx_local
		! update left boundary i = 1 for rank = 0 
		if (rank == 0 .and. i == 1) then
			u_new(i,k) = u_exact(i,k)
		else
			u_new(i,k) = u_old(i,k-1) - c*dt*(u_old(i+1,k-1) - u_old(i-1,k-1))/(2*dx) &
					 + alpha*dt*(u_old(i+1,k-1) - 2.0*u_old(i,k-1) + u_old(i-1,k-1))/(dx*dx)
		end if
	end do
	
	! update the right boundary i = nx_local+1 and rank = nprocs-1
	if (rank == nprocs-1) then
		u_new(nx_local+1,k) = u_exact(nx_local+1,k) 
	end if

!------------------------------------------------------------------
! communicate the solution between ghost nodes	
! left and right boundary of internal domains gets updated for u_old

	! send solution of last grid point of local domain (rank 0,1..) (i = nx) to
    ! left ghost nodes of neighbouring right domain ((rank+1) 1,2..) (i = 0)
    if (rank /= nprocs-1) then
        call MPI_Isend(u_new(nx_local,k), & ! starting address of data to be send
                                       1, & ! number of bytes of data to be send
                    MPI_DOUBLE_PRECISION, & ! datatype of data to be send
                                  rank+1, & ! target processor rank
                        id_left_to_right, & ! tag of the send signal
                          MPI_COMM_WORLD, & !
                                  req(1), & !
                                    ierr)   ! ierr = 0 if successful
    end if
	
	! send solution of first grid point of local domain (rank 1,2..) (i = 1) to
    ! right ghost nodes of neighbouring left domain ((rank-1) 0,1..) (i = nx+1)
    if (rank /= 0) then
        call MPI_Isend(      u_new(1,k), & ! starting address of data to be send
                                      1, & ! number of bytes of data to be send
                   MPI_DOUBLE_PRECISION, & ! datatype of data to be send
                                 rank-1, & ! target processor rank
                       id_right_to_left, & ! tag of the send signal
                         MPI_COMM_WORLD, & !
                                 req(2), & !
                                   ierr)   ! ierr = 0 if successful
    end if
	
	! receive data from left domain (rank-1 0,1) (i = nx)  and
    ! store it in left ghost nodes (i = 0)
    if (rank /= 0) then
        call MPI_Irecv(      u_old(0,k), & ! starting address of data to be receive
                                      1, & ! number of bytes of data to be receive
                   MPI_DOUBLE_PRECISION, & ! datatype of data to be receive
                                 rank-1, & ! source processor rank
                       id_left_to_right, & ! tag of the send signal
                         MPI_COMM_WORLD, & !
                                 req(3), & !
                                   ierr)   ! ierr = 0 if successful
    end if
    
    ! receive data from right domain (rank+1 1,2) (i = 1)  and
    ! store it in right ghost nodes (i = nx+1)
    if (rank /= nprocs-1) then
        call MPI_Irecv(u_old(nx_local+1,k), & ! starting address of data to be send
                                         1, & ! number of bytes of data to be send
                      MPI_DOUBLE_PRECISION, & ! datatype of data to be send
                                   rank+1 , & ! source processor rank
                          id_right_to_left, & ! tag of the send signal
                            MPI_COMM_WORLD, & !
                                    req(4), & !
                                      ierr)   ! ierr = 0 if successful
    end if
    
! wait till Isend and Irecv above is completed
    call MPI_Waitall(4, req, status_array, ierr)
    
! update old field to new field 
    do i = 1,nx_local
		u_old(i,k) = u_new(i,k)
    end do
    
! update old field for left boundary
    if (rank == 0) then
		u_old(1,k) = u_new(1,k)
    end if
    
! update old field for right boundary
    if (rank == nprocs-1) then
		u_old(nx_local+1,k) = u_new(nx_local+1,k)
    end if 
end do

! calculate error using u_old and u_exact 
	do i = 1,nx_local
		if (rank == 0 .and. i ==0) cycle
		local_error_sum = local_error_sum + abs(u_old(i,time_steps) - u_exact(i,time_steps)) 
	end do

! take sum of error from each processor  
	call MPI_Reduce(local_error_sum, & ! variable to be collected from all processors
                   global_error_sum, & ! variable in which result is to be stored
                                  1, & ! number of values
               MPI_DOUBLE_PRECISION, & ! datatype of variable
                            MPI_SUM, & ! type of operation
                                  0, & ! rank of processor in which variable is to be stored
                     MPI_COMM_WORLD, & !
                                ierr)  ! ierr = 0 if successful
                                
	
	deallocate(u_old, u_new, u_exact)
	deallocate(x)
	deallocate(t)
	
! wait till all processors come to this points 
! this is for each processor will start with same randomly generated phase angle
! this also ensures accurate timing and clean output
call MPI_Barrier(MPI_COMM_WORLD, & !
                            ierr) ! ierr = 0 if successful
                            

end do

call cpu_time(stop_time)

! ensamble average of error
if (rank == 0) then
	avg_global_error = global_error_sum/(nx_global*angle_num)
	print *, 'average global', avg_global_error
	print *, 'Time', stop_time-start_time
end if

! call subroutine to write the results for plotting
! call write_tecplot_file(x, t, u_new, u_exact, nx_local, time_steps, rank, nprocs)



!------------------------------------------------------------------
! MPI final calls
call MPI_FINALIZE(ierr) ! ierr = 0 if successful

!deallocate(u_old)
!deallocate(u_new)
!deallocate(u_exact)

end program advection_diffusion


!----------------------------- Initialize subroutine ---------------------------
! This subroutine initializes the exact solution at 0th time step and initial 
! condition at 0th time step 
!-------------------------------------------------------------------------------
subroutine initialize(rank, nprocs, x, nx_local, time_steps, u_old, u_exact, phi)

implicit none

integer, intent(in)                             	    :: nx_local, time_steps, rank, nprocs
real*8, intent(in), dimension(0:nx_local+1)		    	:: x
real*8, intent(out), dimension(0:nx_local+1, 0:time_steps)  :: u_old
real*8, intent(out), dimension(0:nx_local+1, 0:time_steps)  :: u_exact
real*8							    :: phi, kappa

integer							    :: i, k
real*8, parameter                   :: pi = 3.141592653589793238

common /wavenumber/ kappa
! print *, 'rank = ', rank, 'u_exact = ', u_exact
do i = 0,nx_local+1
	u_exact(i,0) = sin(kappa*x(i) + pi*phi/180)
	u_old(i,0)   = u_exact(i,0)
end do

end subroutine initialize


!------------------------------- Export results subroutine ---------------------
! This subroutine writes the results in .dat file which can be plotted using
! tecplot
!-------------------------------------------------------------------------------
subroutine write_tecplot_file(x, t, u_new, u_exact, nx_local, time_steps, rank, nprocs)

implicit none

integer, intent(in)	:: nx_local, time_steps, rank, nprocs
real*8, intent(in), dimension(0:nx_local+1, 0:time_steps)	:: u_new, u_exact
real*8, intent(in), dimension(0:nx_local+1)			:: x
real*8, intent(in), dimension(0:time_steps)			:: t

integer								:: i, k
character(80)						:: char_temp, filename


!-------------------------------------------------------------------------------
! store the value of rank as a character in the the character variable char_temp
write(char_temp, '(i5)') rank

!-------------------------------------------------------------------------------
! define the file name
filename = "tecplot_" // trim(adjustl(char_temp)) // '.dat'

!-------------------------------------------------------------------------------
! Open the file and start writing the data
open(unit = 1, file = filename, status = "unknown")

!-------------------------------------------------------------------------------
! header
write(1,*) 'title =', '"Partition_'// trim(adjustl(char_temp)), '"'
write(1,*) 'variables = "x", "U", "U_Exact", "U-U_Exact"'

!if (rank == 0) then
!	write(1,*) 'zone T =', 'Partition_'// trim(adjustl(char_temp)),' i', nx_local+1, 'j=', time_steps+1
!else
!	write(1,*) 'zone T =', 'Partition_'// trim(adjustl(char_temp)),' i', nx_local+2, 'j=', time_steps+1
!end if

! do k = 0, time_steps
	do i = 0, nx_local+1
		if (rank == 0 .and. i ==0) cycle
		
		write(1,*) x(i), u_new(i,time_steps), u_exact(i,time_steps), u_new(i,time_steps) - u_exact(i,time_steps)
	
	end do
!end do

! close the file
close(1)
end subroutine write_tecplot_file
