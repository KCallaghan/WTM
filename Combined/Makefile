# Default to optimized compile
opt?=-O3
fsanitize?=
$(info Optimization flag: $(opt))
$(info Defaults to "-O3"; use "opt=X" to change it.)

PETSC_DIR=/home/kerry/petsc
PETSC_ARCH=arch-linux-c-debug

CXXFLAGS="-g -O3 -Wall -Wno-unknown-pragmas -pedantic -I../common/richdem/include -I/usr/include/gdal"
CXXFLAGS="$(CXXFLAGS) -DRICHDEM_LOGGING"

export GDAL_LIBS=`gdal-config --libs` -DUSEGDAL -DUSE_GDAL
export GDAL_CFLAGS=`gdal-config --cflags` -DUSE_GDAL -DUSEGDAL

RICHDEM_GIT_HASH="-NA-"
RICHDEM_COMPILE_TIME=`date -u +'%Y-%m-%d %H:%M:%S UTC'`
export RD_CXX_FLAGS=-I../common/richdem/include -DUSE_GDAL\
	-DRICHDEM_GIT_HASH="\"$(RICHDEM_GIT_HASH)\"" \
	-DRICHDEM_COMPILE_TIME="\"$(RICHDEM_COMPILE_TIME)\""
export CXXFLAGS=$(GDAL_CFLAGS) --std=c++17 ${opt} -g -Wall -Wno-unknown-pragmas
# use -fsanitize=address (above) to check for memory errors
export LIBS= $(GDAL_LIBS) -fopenmp -lpetsc

objectsTWSM = twsm.o parameters.o richdem.o ArrayPack.o transient_groundwater.o
objectsDH = run_dephier.o parameters.o richdem.o ArrayPack.o

all : twsm.x dephier.x

dephier.x : $(objectsDH)
	mpicxx $(CXXFLAGS) $(RD_CXX_FLAGS) $(objectsDH) $(LIBS) -o dephier.x

run_dephier.o : run_dephier.cpp
	mpicxx -c -fPIC $(CXXFLAGS) $(RD_CXX_FLAGS) run_dephier.cpp $(LIBS) \
		-o run_dephier.o

twsm.x : $(objectsTWSM)
	mpicxx $(CXXFLAGS) $(RD_CXX_FLAGS) $(objectsTWSM) $(LIBS) $(GDAL_LIBS)  \
			-o twsm.x

twsm.o : TWSM.cpp
	mpicxx -c -fPIC $(CXXFLAGS) $(RD_CXX_FLAGS) TWSM.cpp $(LIBS) $(GDAL_LIBS)  -o twsm.o

richdem.o : ../common/richdem/src/richdem.cpp
	mpicxx -fPIC -c $(CXXFLAGS) $(RD_CXX_FLAGS) $(GDAL_LIBS) ../common/richdem/src/richdem.cpp $(LIBS) -o richdem.o

transient_groundwater.o : transient_groundwater.cpp
	mpicxx -c -fPIC $(CXXFLAGS) $(RD_CXX_FLAGS) transient_groundwater.cpp $(LIBS)  -o transient_groundwater.o


ArrayPack.o : ArrayPack.cpp  ArrayPack.hpp
	mpicxx -fPIC -c $(CXXFLAGS) $(RD_CXX_FLAGS) $(GDAL_LIBS) ArrayPack.cpp $(LIBS) -o ArrayPack.o

parameters.o : parameters.cpp parameters.hpp
	mpicxx -fPIC -c $(CXXFLAGS) $(RD_CXX_FLAGS) $(GDAL_LIBS) parameters.cpp $(LIBS) -o parameters.o

clean:
	rm twsm.x
	rm dephier.x

distclean:
	rm $(objectsTWSM)
	rm run_dephier.o

