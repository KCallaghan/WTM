RICHDEM_GIT_HASH="-NA-"
RICHDEM_COMPILE_TIME=`date -u +'%Y-%m-%d %H:%M:%S UTC'`
export RD_CXX_FLAGS=-I../richdem/include -DRICHDEM_GIT_HASH="\"$(RICHDEM_GIT_HASH)\"" -DRICHDEM_COMPILE_TIME="\"$(RICHDEM_COMPILE_TIME)\""
export CXXFLAGS=--std=c++17 -g -O3 -Wall -Wno-unknown-pragmas #-fsanitize=address
export LIBS=-lnetcdf


a.out: ArrayPack.cpp  ArrayPack.hpp  groundwater.cpp parameters.cpp  parameters.hpp Makefile
	g++-7 $(CXXFLAGS) $(RD_CXX_FLAGS) groundwater.cpp parameters.cpp ArrayPack.cpp ../richdem/include/richdem/richdem.cpp $(LIBS)

clean:
	rm a.out
