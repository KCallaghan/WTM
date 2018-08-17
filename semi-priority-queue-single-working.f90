subroutine Merge(A,NA,B,NB,C,NC)
 
   integer, intent(in) :: NA,NB,NC         ! Normal usage: NA+NB = NC
   integer, intent(in out) :: A(NA)        ! B overlays C(NA+1:NC)
   integer, intent(in)     :: B(NB)
   integer, intent(in out) :: C(NC)
 
   integer :: I,J,K
 
   I = 1; J = 1; K = 1;
   do while(I <= NA .and. J <= NB)
      if (A(I) <= B(J)) then
         C(K) = A(I)
         I = I+1
      else
         C(K) = B(J)
         J = J+1
      endif
      K = K + 1
   enddo
   do while (I <= NA)
      C(K) = A(I)
      I = I + 1
      K = K + 1
   enddo
   return
 
end subroutine merge
 
recursive subroutine MergeSort(A,N,T,indices)
 
   integer, intent(in) :: N
   integer, dimension(N), intent(in out) :: A
   real,dimension(2,N) :: indices
   integer, dimension((N+1)/2), intent (out) :: T
 
   integer :: NA,NB,V
 
   if (N < 2) return
   if (N == 2) then
      if (A(1) > A(2)) then
         V = A(1)
         A(1) = A(2)
         A(2) = V

         V_row = indices(1,1)
         v_col = indices(2,1)
         indices(1,1) = indices(1,2)
         indices(1,2) = V_row

         indices(2,1) = indices(2,2)
         indices(2,1) = V_col



      endif
      return
   endif      
   NA=(N+1)/2
   NB=N-NA
 
   call MergeSort(A,NA,T,indices)
   call MergeSort(A(NA+1),NB,T,indices)
 
   if (A(NA) > A(NA+1)) then
      T(1:NA)=A(1:NA)
      call Merge(T,NA,A(NA+1),NB,A,N)
   endif
   return
 
end subroutine MergeSort
 





program surface_water

use netcdf
use mpi

implicit none


integer :: n2,n3,ncid,varid,iret,i,j



character*20 :: surfdatadir,time_start
character*100 :: filetopo_start

real,allocatable,dimension(:,:) :: topo
real,allocatable,dimension(:) :: hz_1D

  integer, parameter :: Merge_size = 4000
   real, dimension(Merge_size) :: result_topo
   integer, dimension ((Merge_size+1)/2) :: T



integer :: col,row,converged,counter


real :: upvalue,downvalue,leftvalue,rightvalue,water

real maxdiff

real,allocatable,dimension(:,:) :: h,hold,topo_read,h_values,&
hz_read,truth_array,diff


real,allocatable,dimension(:,:) :: arr

REAL,PARAMETER :: UNDEF = -1.0E+7



!setup data ***************************************************************************************************

print *,'initialise'


n2 = 2000            !number of columns in the rotated version of madagascar. Fortran thinks these are ROWS
n3 = 1000            !number of rows in the rotated version of madagascar. Fortran thinks these are COLUMNS.

time_start = 'Mad_000000'
surfdatadir = 'surfdata/'

upvalue = 0.
downvalue = 0.
leftvalue = 0.
rightvalue = 0.
counter = 0
converged = 0


filetopo_start = trim(surfdatadir)//trim(time_start)//'_topo_rotated.nc'



allocate(h(n2,n3))
allocate(hold(n2,n3))
h(:,:) = 1.
hold = h



allocate(hz_read(n2,n3))
allocate(h_values(n2,n3))
h_values = h

allocate(truth_array(n2,n3))
truth_array = 1

print *,'line 166'


allocate(topo(n2,n3))
     
 iret = nf90_open(filetopo_start,0,ncid) !reading in the topo
 call check_err(iret)

 iret = nf90_inq_varid(ncid,'value',varid)
 call check_err(iret)

 iret = nf90_get_var(ncid,varid,topo)
 call check_err(iret)

 iret = nf90_close(ncid)
 call check_err(iret)


where(topo .le. UNDEF) topo = 0. !change undefined cells to 0


hz_read = topo + h_values   

allocate(hz_1D(n2*n3))


do i=1,n2
    hz_1D(((i-1)*n3)+1:i*n3) = hz_read(i,:)

end do

allocate(arr(2,n2*n3))

do row=1,n2
 arr(1,((row-1)*n3)+1:row*n3)=row
end do

do row=1,n2
  do col = 1,n3
    arr(2,(row-1)*n3+col)=col
end do
end do



MAIN: do while(converged .eq.0)


hold=h_values
counter = counter + 1
do i=1,n2
  do j=1,n3
    row = arr(1,n2*n3-(j-1)-(i-1)*n3)
    col = arr(2,n2*n3-(j-1)-(i-1)*n3)



           upvalue = hz_read(row,col) - hz_read(row-1,col)
            downvalue = hz_read(row,col) - hz_read(row+1,col)
            leftvalue = hz_read(row,col) - hz_read(row,col-1)
            rightvalue = hz_read(row,col) - hz_read(row,col+1)

            if(max(upvalue,downvalue,leftvalue,rightvalue) .le. 0.) then 
                 truth_array(row,col) = 0

                CYCLE   

            elseif(max(upvalue,downvalue,leftvalue,rightvalue) .eq. upvalue) then
                water = min(h_values(row,col),upvalue/2.)
                h_values(row,col) = h_values(row,col) - water
                h_values(row-1,col) = h_values(row-1,col) + water


            elseif(max(upvalue,downvalue,leftvalue,rightvalue) .eq. downvalue) then
                water = min(h_values(row,col),downvalue/2.)
                h_values(row,col) = h_values(row,col) - water
                h_values(row+1,col) = h_values(row+1,col) + water


            elseif(max(upvalue,downvalue,leftvalue,rightvalue) .eq. rightvalue) then
                water = min(h_values(row,col),rightvalue/2.)
                h_values(row,col) = h_values(row,col) - water
                h_values(row,col+1) = h_values(row,col+1) + water


            elseif(max(upvalue,downvalue,leftvalue,rightvalue) .eq. leftvalue) then
                water = min(h_values(row,col),leftvalue/2.)
                h_values(row,col) = h_values(row,col) - water
                h_values(row,col-1) = h_values(row,col-1) + water


            endif


  end do
end do


diff = abs(hold-h_values)
maxdiff = 0
maxdiff = maxval(diff)
print *,'maxdiff',maxdiff

hz_read = topo + h_values  

do i=1,n2
    hz_1D(((i-1)*n3)+1:i*n3) = hz_read(i,:)

end do


do row=1,n2
 arr(1,((row-1)*n3)+1:row*n3)=row
end do

do row=1,n2
  do col = 1,n3
    arr(2,(row-1)*n3+col)=col
end do
end do



   call MergeSort(hz_1D,Merge_size,T,arr)

print *,counter
if(counter.eq.1000)then
converged=1
endif


end do MAIN


 print *,'done'

  
 
    open(23,file = 'test_priority_version.dat',form='unformatted',access='stream')!access='direct',recl=n2*n3)!access='stream')!do the final write - create the file


    print *,'the file has been opened'

    write(23,rec=1)((h_values(i,j),i=1,n2),j=1,n3) !and write it to file

    print *,'the file has been written'

    close(23)

    print *,'written to file'





end program surface_water


!********************************************************************************************************

subroutine check_err(statusnc)

  use netcdf
  integer statusnc

  if(statusnc.ne. nf90_noerr) then
      stop 'Stopped due to a catch in the check_err subroutine'
  endif

end subroutine check_err




