### One-dimensional diffusion equation

The LAPACK package (lapack-3.10.0.tar.gz - download via http://www.netlib.org/lapack/#_lapack_version_3_10_0_2) has been added to the folder, 
as well as the compiled executable (diffusion.out). So, to compile the programm just unzip the LAPACK package (lapack-3.10.0.tar.gz) and write
for "~" the directory path where the LAPACK folder is stored: 

- gfortran diffusion.f90 -llapack -lblas -L/~/lapack-3.10.0 -O3 -W

For further details see the report.
