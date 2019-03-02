!-------------------------------------------------------------------------------
! The code solves the advection diffusion equation using Proxy equation 
! approach given in the below reference paper using MPI. The code is 
! still under development
!
!---------------------------------------------------------------------!
!
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
integer 	                           	:: ierr ! error signal variable, Standard value - 0
integer 	                           	:: rank ! process ID (pid) / Number
integer		                           	:: nprocs ! number of processors
integer, dimension(MPI_STATUS_SIZE)    	:: status ! status for receive
integer, parameter					   	:: id_left_to_right = 10 ! tag for send/ receive
integer, parameter					   	:: id_right_to_left = 20 ! tag for send/ receive
integer, dimension(4)                  	:: req = MPI_REQUEST_NULL
integer, dimension(MPI_STATUS_SIZE, 4)	:: status_array

! Computational parameters related to advection-diffusion equation
real*8, parameter					   	:: T = 1.0, c = 1.0, alpha = 0.1, cfl = 0.1
real*8, parameter					   	:: kappa = 2.0
real*8, dimension(:), allocatable	   	:: u_old 
real*8, dimension(:), allocatable	   	:: u_new
real*8, dimension(:), allocatable	   	:: u_exact
real*8, dimension(:), allocatable	   	:: x
real*8, parameter					   	:: pi = 3.14159265358979323846264338D0
real*8, parameter					   	:: x_min = 0.0, x_max = 2.0*pi
integer, parameter					   	:: nx_global = 64
integer									:: nx_local
real*8									:: dx, dt, x_loc_min
integer									:: i, k
integer, parameter						:: istart = 1, kstart = 1
integer									:: i_global_low, i_global_high

!-------------------------------------------------------------------------------
! MPI Starting
!-------------------------------------------------------------------------------

call MPI_INIT(ierr) ! Initialize MPI

call MPI_COMM_SIZE(MPI_COMM_WORLD, &
                           nprocs, &    ! number of processors
                           ierr     )   ! ierr = 0 if successful

call MPI_COMM_RANK(MPI_COMM_WORLD, &
                             rank, &    ! rank of the processor
                             ierr   )   ! ierr = 0 if successful

if (rank == 0) then
    write(*,*) 'No. of processor = ', nprocs, 'Rank = ', rank
    write(*,*)
    write(*,*) 'Nx_Global = ', nx_global
    write(*,*) 'Total size = ', nx_global
end if

!-------------------------------------------------------------------------------
! Allocate x-position to each grid point
! The domain is partitioned in x-direction and there are two global ghost nodes
! i = 0 and i = nx_global + 1
! Each local domain has two ghost nodes i = 0 and i = nx_local + 1
! (i=0)_rank = (i = nx_local)_(rank-1)
! (i = nx+1)_rank = (i = 0)_(rank+1)
!-------------------------------------------------------------------------------

! local size
nx_local = nx_glovbal/nprocs

! grid size
dx = (x_max - x_min)/(nx_global)

! allocate the local array for field variables and grid positions
allocate(      u_old(0:nx_local+1))
allocate(      u_new(0:nx_local+1))
allocate(    u_exact(0:nx_local+1))
allocate(          x(0:nx_local+1))

! allocate grid position
do i = 0,nx+1
	x_loc_min = rank*((x_max - x_min)/nprocs) ! dimesnion in y-direction is 1.0
	x(i) = x_loc_min + dx*(i-1)		! position of x for each grid
end do




















