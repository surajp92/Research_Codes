clear
ls
rm *.o *,mod *.exe
ls
gfortran -c *.f95
ls
gfortran *.o -o serial_adv_diff.exe
ls
./serial_adv_diff.exe
