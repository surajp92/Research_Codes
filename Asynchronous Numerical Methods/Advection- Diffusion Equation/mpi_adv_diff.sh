clear
ls
# Remove old executables
rm *.exe
ls

# COmpile and Build the executable
/usr/bin/mpif90 -fbounds-check -o mpi_adv_diff.exe advection_diffusion_v2.f95
ls
# Run the executable
/usr/bin/mpirun -n 2 ./mpi_adv_diff.exe
