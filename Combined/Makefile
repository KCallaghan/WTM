RICHDEM_GIT_HASH="-NA-"
RICHDEM_COMPILE_TIME=`date -u +'%Y-%m-%d %H:%M:%S UTC'`
export RD_CXX_FLAGS=-I../common/richdem/include -DRICHDEM_GIT_HASH="\"$(RICHDEM_GIT_HASH)\"" -DRICHDEM_COMPILE_TIME="\"$(RICHDEM_COMPILE_TIME)\""
export CXXFLAGS=--std=c++17 -g -Wall -Wno-unknown-pragmas #-fsanitize=address
export LIBS=-lnetcdf


a.out: DisjointDenseIntSet.hpp ArrayPack.cpp  ArrayPack.hpp  evaporation.hpp transient_groundwater.hpp run_transient_like_equilibrium_coupled.cpp parameters.cpp parameters.hpp Makefile ../common/netcdf.hpp
	g++-7 $(CXXFLAGS) $(RD_CXX_FLAGS) run_transient_like_equilibrium_coupled.cpp parameters.cpp ArrayPack.cpp ../common/richdem/include/richdem/richdem.cpp $(LIBS)

clean:
	rm a.out
