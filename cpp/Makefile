export GDAL_LIBS=`gdal-config --libs`
export GDAL_CFLAGS=`gdal-config --cflags` -DUSEGDAL
RICHDEM_GIT_HASH=`git rev-parse HEAD`
RICHDEM_COMPILE_TIME=`date -u +'%Y-%m-%d %H:%M:%S UTC'`
export CXXFLAGS=$(GDAL_CFLAGS) --std=c++17 -g -O3 -Wall -Wno-unknown-pragmas -I../include -DRICHDEM_GIT_HASH="\"$(RICHDEM_GIT_HASH)\"" -DRICHDEM_COMPILE_TIME="\"$(RICHDEM_COMPILE_TIME)\""
export LIBS=$(GDAL_LIBS)

a.out: dep-hier.cpp
	g++-7 $(CXXFLAGS) dep-hier.cpp -lnetcdf -I/home/rick/projects/watershed/richdem/include -DUSEGDAL $(LIBS) -Wall

clean:
	rm a.out