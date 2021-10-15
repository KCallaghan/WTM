# Default to optimized compile
opt?=-O3
fsanitize?=
$(info Optimization flag: $(opt))
$(info Defaults to "-O3"; use "opt=X" to change it.)

PETSC_DIR=/home/kerry/petsc
PETSC_ARCH=arch-linux-c-debug


GDAL_LIBS="`gdal-config --libs`"
GDAL_CFLAGS="`gdal-config --cflags` "
LIBS="$GDAL_LIBS"
CXXFLAGS="$CXXFLAGS $GDAL_CFLAGS"

RICHDEM_GIT_HASH="-NA-"
RICHDEM_COMPILE_TIME=`date -u +'%Y-%m-%d %H:%M:%S UTC'`
export RD_CXX_FLAGS=-I../common/richdem/include -DUSE_GDAL\
	-DRICHDEM_GIT_HASH="\"$(RICHDEM_GIT_HASH)\"" \
	-DRICHDEM_COMPILE_TIME="\"$(RICHDEM_COMPILE_TIME)\""
export CXXFLAGS=--std=c++17 ${opt} -g -Wall -Wno-unknown-pragmas
# use -fsanitize=address (above) to check for memory errors
export LIBS= -fopenmp -lpetsc

objectsTWSM = twsm.o parameters.o richdem.o ArrayPack.o transient_groundwater.o
objectsDH = run_dephier.o parameters.o richdem.o ArrayPack.o

all : twsm.x dephier.x

dephier.x : $(objectsDH)
	mpicxx $(CXXFLAGS) $(RD_CXX_FLAGS) $(objectsDH) $(LIBS) -o dephier.x

run_dephier.o : run_dephier.cpp
	mpicxx -c -fPIC $(CXXFLAGS) $(RD_CXX_FLAGS) run_dephier.cpp $(LIBS) \
		-o run_dephier.o

twsm.x : $(objectsTWSM)
	mpicxx $(CXXFLAGS) $(RD_CXX_FLAGS) $(objectsTWSM) $(LIBS) ${PETSC_KSP_LIB} -I$(PETSC_DIR)/${PETSC_ARCH}/include  -I ${PETSC_DIR}/include \
			-o twsm.x

twsm.o : TWSM.cpp
	mpicxx -c -fPIC $(CXXFLAGS) $(RD_CXX_FLAGS) TWSM.cpp $(LIBS) ${PETSC_KSP_LIB} -I$(PETSC_DIR)/${PETSC_ARCH}/include  -I ${PETSC_DIR}/include  -o twsm.o

transient_groundwater.o : transient_groundwater.cpp
	mpicxx -c -fPIC $(CXXFLAGS) $(RD_CXX_FLAGS) transient_groundwater.cpp $(LIBS) ${PETSC_KSP_LIB} -I$(PETSC_DIR)/${PETSC_ARCH}/include  -I ${PETSC_DIR}/include  -o transient_groundwater.o

richdem.o : ../common/richdem/src/richdem.cpp
	mpicxx -fPIC -c $(CXXFLAGS) $(RD_CXX_FLAGS) ../common/richdem/src/richdem.cpp $(LIBS) -o richdem.o

ArrayPack.o : ArrayPack.cpp  ArrayPack.hpp
	mpicxx -fPIC -c $(CXXFLAGS) $(RD_CXX_FLAGS) ArrayPack.cpp $(LIBS) -o ArrayPack.o

parameters.o : parameters.cpp parameters.hpp
	mpicxx -fPIC -c $(CXXFLAGS) $(RD_CXX_FLAGS) parameters.cpp $(LIBS) -o parameters.o

clean:
	rm twsm.x
	rm dephier.x

distclean:
	rm $(objectsTWSM)
	rm run_dephier.o


