all:
	mpif90 -O3 -g Equilibrium_version_with_SW.f90  -I/usr/include -L/usr/lib/x86_64-linux-gnu/ -lnetcdff -fsanitize=address
	#mpif90 -O3 -g priority-SW-coupled-GW-and-SW_WORKING.f90 -I/usr/include -L/usr/lib/x86_64-linux-gnu/ -lnetcdff
	#mpif90 -O3 -g Equilibrium_version_HADCM3.f90 -I/usr/include -L/usr/lib/x86_64-linux-gnu/ -lnetcdff
	
