RICHDEM_GIT_HASH="-NA-"
RICHDEM_COMPILE_TIME=`date -u +'%Y-%m-%d %H:%M:%S UTC'`
export RD_CXX_FLAGS=-I../common/richdem/include -DRICHDEM_GIT_HASH="\"$(RICHDEM_GIT_HASH)\"" -DRICHDEM_COMPILE_TIME="\"$(RICHDEM_COMPILE_TIME)\""
export CXXFLAGS=--std=c++17 -O3 -g -Wall -Wno-unknown-pragmas #-fsanitize=address
export LIBS=-lnetcdf

twsm: DisjointDenseIntSet.hpp ArrayPack.cpp  ArrayPack.hpp  evaporation.hpp transient_groundwater.hpp TWSM.cpp irf.cpp parameters.cpp parameters.hpp Makefile ../common/netcdf.hpp
	g++-7 $(CXXFLAGS) $(RD_CXX_FLAGS) TWSM.cpp parameters.cpp ArrayPack.cpp ../common/richdem/include/richdem/richdem.cpp $(LIBS) -o twsm

clean:
	rm a.out
