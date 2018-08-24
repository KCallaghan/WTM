all: equilibrium_w_sw.exe priority.exe #hadcm3.exe

equilibrium_w_sw.exe: Equilibrium_version_with_SW.f90          
	mpif90 --version
	mpif90 -f90 gfortran-7  -O3 -g             Equilibrium_version_with_SW.f90           -o equilibrium_w_sw.exe  -I/usr/include -L/usr/lib/x86_64-linux-gnu/ -lnetcdff -ffast-math -march=native -mtune=native #-fsanitize=address

priority.exe:         priority-SW-coupled-GW-and-SW_WORKING.f90
	mpif90 --version
	mpif90 -f90 gfortran-7 -O3 -g -std=f2008 priority-SW-coupled-GW-and-SW_WORKING.f90 -o priority.exe          -I/usr/include -L/usr/lib/x86_64-linux-gnu/ -lnetcdff -ffast-math -march=native -mtune=native #-fsanitize=address

#hadcm3.exe:           Equilibrium_version_HADCM3.f90           
#	mpif90 -O3 -g Equilibrium_version_HADCM3.f90            -o hadcm3.exe            -I/usr/include -L/usr/lib/x86_64-linux-gnu/ -lnetcdff #-fsanitize=address

clean:
	rm -f *.exe