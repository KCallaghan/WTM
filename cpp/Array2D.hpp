#ifndef _array2d_hpp_
#define _array2d_hpp_

#include <string>
#include <iostream>
#include <netcdf.h>
#include <stdexcept>
#include <memory>

template<class T>
class Array2D {
 public:
  T *data;
 private:
  int mywidth  = -1;
  int myheight = -1;
  void getDimLength(const int ncid, const int dimnum);
 public:
  Array2D() = default;
  Array2D(std::string filename, std::string datavar);
  Array2D(const int width0, const int height0, const T val0);
  ~Array2D();
  T&   operator()(const int x, const int y);
  T    operator()(const int x, const int y) const;
  T&   operator()(const int i);
  T    operator()(const int i) const;
  int  width()  const;
  int  height() const;
  int  size()   const;
  bool isEdgeCell(const int x, const int y) const;
  bool inGrid(const int x, const int y) const;
};

template<class T>
void Array2D<T>::getDimLength(const int ncid, const int dimnum){
  char dimnamebuf[100];
  size_t dimlen;
  int retval;
  if((retval= nc_inq_dimname(ncid, dimnum, dimnamebuf)))
    throw std::runtime_error("Couldn't get name of dimension!");

  std::string dimname(dimnamebuf);

  if((retval=nc_inq_dimlen(ncid, dimnum, &dimlen)))
    throw std::runtime_error("Couldn't get length for dimension '"+dimname+"'!");

  if(dimname=="lat")
    myheight = dimlen;
  else if(dimname=="lon")
    mywidth = dimlen;
}

template<class T>
Array2D<T>::Array2D(std::string filename, std::string datavar){
  /* This will be the netCDF ID for the file and data variable. */
  int ncid, varid, retval, dim_count;

  /* Open the file. NC_NOWRITE tells netCDF we want read-only access
   * to the file.*/
  if ((retval = nc_open(filename.c_str(), NC_NOWRITE, &ncid)))
    throw std::runtime_error("Failed to open file '" + filename + "'!");

  if ((retval = nc_inq_ndims(ncid, &dim_count)))
    throw std::runtime_error("Failed to get number of dimensions from file '" + filename + "'!");

  if(dim_count!=3)
    throw std::runtime_error("File '" + filename + "' did not have 2 dimensions!");    

  getDimLength(ncid, 0);
  getDimLength(ncid, 1);
  getDimLength(ncid, 2);

  if(mywidth==-1 || myheight==-1)
    throw std::runtime_error("File '" + filename + "' did not have a lat or lon dimension!");    

  data = new float[mywidth*myheight];

  /* Get the varid of the data variable, based on its name. */
  if ((retval = nc_inq_varid(ncid, datavar.c_str(), &varid)))
    throw std::runtime_error("Failed to get dataset '"+datavar+"' from file '" + filename + "'!");

  /* Read the data. */
  //TODO: Check data type
  if ((retval = nc_get_var_float(ncid, varid, data)))
    throw std::runtime_error("Failed to read data from '"+datavar+" from file '" + filename + "'! Error: " + nc_strerror(retval));

  /* Close the file, freeing all resources. */
  if ((retval = nc_close(ncid)))
    throw std::runtime_error("Failed to close file '" + filename + "'!");
}

template<class T>
Array2D<T>::Array2D(const int width0, const int height0, const T val0) {
  data = new T[width0*height0];
  mywidth  = width0;
  myheight = height0;
  for(int i=0;i<mywidth*myheight;i++)
    data[i] = val0;
}

template<class T>
Array2D<T>::~Array2D() {
  delete[] data;
}


template<class T>
T& Array2D<T>::operator()(const int x, const int y) {
  return data[y*mywidth+(x%mywidth)];
}

template<class T>
T  Array2D<T>::operator()(const int x, const int y) const {
  return data[y*mywidth+(x%mywidth)];
}

template<class T>
T& Array2D<T>::operator()(const int i) {
  return data[i];
}

template<class T>
T  Array2D<T>::operator()(const int i) const {
  return data[i];
}

template<class T>
int Array2D<T>::width() const {
  return mywidth;
}

template<class T>
int Array2D<T>::height() const {
  return myheight;
}

template<class T>
int Array2D<T>::size() const {
  return mywidth*myheight;
}

template<class T>
bool Array2D<T>::isEdgeCell(const int x, const int y) const {
  return y==0 || y==myheight-1;
}

template<class T>
bool Array2D<T>::inGrid(const int x, const int y) const {
  return y>=0 || y<=myheight-1;
}

#endif
