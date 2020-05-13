RICHDEM_GIT_HASH="-NA-"
RICHDEM_COMPILE_TIME=`date -u +'%Y-%m-%d %H:%M:%S UTC'`
export RD_CXX_FLAGS=-I../common/richdem/include -DRICHDEM_GIT_HASH="\"$(RICHDEM_GIT_HASH)\"" -DRICHDEM_COMPILE_TIME="\"$(RICHDEM_COMPILE_TIME)\""
export CXXFLAGS=--std=c++17 -g -Wall -Wno-unknown-pragmas #-fsanitize=address
export LIBS=-lnetcdf

objects = test.o parameters.o richdem.o ArrayPack.o transient_groundwater.o
test-aw : $(objects)
	g++ $(CXXFLAGS) $(RD_CXX_FLAGS) $(objects) $(LIBS) \
			-o test-aw.x

test.o : test.cpp
	g++ -c -fPIC $(CXXFLAGS) $(RD_CXX_FLAGS) $(PYBINDFLAGS) test.cpp $(LIBS) -o test.o

twsm.o : TWSM.cpp
	g++ -c -fPIC $(CXXFLAGS) $(RD_CXX_FLAGS) TWSM.cpp $(LIBS) -o twsm.o

transient_groundwater.o : transient_groundwater.cpp
	g++ -c -fPIC $(CXXFLAGS) $(RD_CXX_FLAGS) transient_groundwater.cpp $(LIBS) -o transient_groundwater.o

richdem.o : ../common/richdem/include/richdem/richdem.cpp
	g++ -fPIC -c $(CXXFLAGS) $(RD_CXX_FLAGS) ../common/richdem/include/richdem/richdem.cpp $(LIBS) -o richdem.o

ArrayPack.o : ArrayPack.cpp  ArrayPack.hpp
	g++ -fPIC -c $(CXXFLAGS) $(RD_CXX_FLAGS) ArrayPack.cpp $(LIBS) -o ArrayPack.o

parameters.o : parameters.cpp parameters.hpp
	g++ -fPIC -c $(CXXFLAGS) $(RD_CXX_FLAGS) parameters.cpp $(LIBS) -o parameters.o

clean:
	rm test-aw.x

distclean:
	rm test-aw.x
	rm $(objects)

