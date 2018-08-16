!use rotated versions of data for pacman flow
!This version of the code sorts cells by priority before moving water, for the surface water component.
!It is significantly faster than just running through the array start to finish. 
!We use a finite difference method for groundwater flow
!groundwater and surface water are on the same array now 


!*************************************************************************************************************

subroutine Merge(A,NA,B,NB,C,NC)
 
   integer, intent(in)     :: NA,NB,NC     ! Normal usage: NA+NB = NC
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
 
recursive subroutine MergeSort(A,N,T,indices)           !A is the hz array, N is the length of the array (or a shorter length), T is an empty array, indices is the indices array to also be sorted. 
 
   integer, intent(in) :: N
   integer, dimension(N), intent(in out) :: A
   real,dimension(2,N) :: indices
   integer, dimension((N+1)/2), intent (out) :: T
 
   integer :: NA,NB,V
 
   if (N < 2) return
   if (N == 2) then
      if (A(1) > A(2)) then                 !Both the hz array with topo + h info and an array recording the indices are sorted
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



!********************************************************************************************************


subroutine dividedomain(n2,n3,numtasks,nini,filemask,ntotal)
  use netcdf
  implicit none

 
  integer :: n2,n3,numtasks
  integer :: nini(1:numtasks-1)
  real,allocatable,dimension(:,:) :: varread
  integer,allocatable,dimension(:) :: ncells
  integer :: iret,ncid,varid,ncount,n,j
  integer(KIND=8) :: ntotal
  character(len=100) :: filemask

  allocate(varread(n2,n3))

  write(6,*)'reading in the mask to divide the domain'
  print *,filemask
  iret = nf90_open(filemask,0,ncid)  !open the mask file
  call check_err(iret, "Unable to load file '" // filemask // "'!")
  write(6,*)'first call'

  iret = nf90_inq_varid(ncid,'value',varid) !get the ID of the value layer
  call check_err(iret, "Problem getting ID of mask layer!")
  write(6,*) 'second call'

  iret = nf90_get_var(ncid,varid,varread) !read the actual values into the array called varread
  call check_err(iret, "Problem reading the values from the mask layer!")
  write(6,*) 'third call'

  iret = nf90_close(ncid) !close the mask file
  call check_err(iret, "Problem closing the mask file!")
  write(6,*)'fourth call'

  ntotal = count(varread>-100) !count the number of land cells. I changed this slightly since I am using mask here rather than topo; all cells with a value of 1 should be included. 

  print *,'ntotal',ntotal,numtasks,'numtasks'

  allocate(ncells(n3))

  ncells = count(varread>-100,1) !The number of cells which are defined in the 1 dimension

  ncount=0

  nini(1) = 2

  n=2 !counter
  
  do j=1,n3
      ncount=ncount+ncells(j) !add the number of cells defined in the current column
      if (ncount .ge. ntotal/(numtasks-1)) then !>= total number of defined cells/number of threads
        nini(n) = j-1 !Telling it which row it needs to start in, to divide up the work equally, so each task handles the same number of cells.
        ncount  = ncells(j) !reset the ncount for the next task
        n       = n+1
      endif
      if(n .eq. numtasks) exit
  end do

  print *,'down here in the subroutine',nini

  deallocate(ncells)
  deallocate(varread)

  return

end subroutine dividedomain


!********************************************************************************************************

subroutine check_err(statusnc, msg)
  use netcdf
  integer statusnc
  character(len=*) msg

  if(statusnc.ne.nf90_noerr) then
    IF (LEN(msg)>0) then
      ERROR STOP msg
    ELSE
      ERROR STOP 'Stopped due to a catch in the check_err subroutine'
    ENDIF
  endif

end subroutine check_err

!********************************************************************************************************



program GWSW

use netcdf
use mpi

implicit none

! Define all variables: ***********************************************************************


integer :: i,j,n2,n3,ncid,varid,error,iret,col,row,converged,counter,nmax,n,iter

integer(KIND=8) :: ntotal

integer ierr, pid,numtasks,tasktype,status(MPI_STATUS_SIZE),rc,columntype,columntypeint,iterations

real :: upvalue,downvalue,leftvalue,rightvalue,water,water1,water2,delta_xy

real diff_total,maxdiff,dx,dy,xn,xs,deltat,placeholder

character(len=20) :: surfdatadir,time_start,initdatadir,time_end

character(len=100) :: filetopo_start,filemask,filetopo_end,filerech_start,filerech_end, &
file_fslope_start, file_fslope_end, filetemp_start, filetemp_end, fileksat,filewtd

real,allocatable,dimension(:,:) :: topo,h,hold,diff,hz_read,hold_read,mask,mask_read

real, allocatable,dimension(:,:) :: topo_start_read, topo_end_read, rech_start_read,&
 rech_end_read,fslope_start_read, fslope_end_read,temp_start_read, &
 temp_end_read,ksat_read,wtd_read,wtd, wtdnew,ksat,fdepth,wtdglob,&
topo_now,rech_now,fslope_now,temp_now,fdepth_now,fdepth_start,&
topo_sent,fdepth_sent,rech_sent,head,kcell,qlat_north,qlat_south,&
qlat_east,qlat_west,landmask


real,allocatable,dimension(:) :: hz_1D,area

integer :: Merge_size
integer,allocatable,dimension(:)  :: T

real,allocatable,dimension(:,:) :: arr


real,allocatable :: xlat(:),alpha(:)

REAL,PARAMETER :: UNDEF = -1.0E+7

REAL(KIND=8) :: SEDGE

REAL(KIND=8),PARAMETER :: pi=3.141592653589793D0+0

integer,allocatable :: domblock(:),domblocksmall(:),nini(:),nend(:)


call MPI_INIT(ierr)      !Initialise MPI
if (ierr .ne. MPI_SUCCESS) then     !Error catching - this part should never run.
    print *,'Error starting MPI program. Terminating.'
    call MPI_ABORT(MPI_COMM_WORLD,rc,ierr)
end if


!*******************************************************************************************

!Data setup:

!Values which may need to change from run to run:

iterations = 500*365  !500 years * 12 months = 6000 iterations for 500-year runs
n2         = 2000     !number of columns in the rotated version of madagascar. Fortran thinks these are ROWS
n3         = 1000     !number of rows in the rotated version of madagascar. Fortran thinks these are COLUMNS.
time_start = 'Mad_021000'
time_end   = 'Mad_020500'

SEDGE = -27. !latitude of the southern edge of Madagascar. 


deltat = (365.*24.*3600./365.) !Seconds in a monthly timestep


!Files and folders:

surfdatadir = 'surfdata/'
initdatadir = 'initdata/'

filetopo_start = trim(surfdatadir)//trim(time_start)//'_topo_rotated.nc'
filetopo_end   = trim(surfdatadir)//trim(time_end)//'_topo_rotated.nc'
  
filerech_start = trim(surfdatadir)//trim(time_start)//'_rech_rotated.nc'
filerech_end   = trim(surfdatadir)//trim(time_end)//'_rech_rotated.nc'

file_fslope_start = trim(surfdatadir)//trim(time_start)//'_fslope_rotated.nc'
file_fslope_end   = trim(surfdatadir)//trim(time_end)//'_fslope_rotated.nc'

filetemp_start = trim(surfdatadir)//trim(time_start)//'_temp_rotated.nc'
filetemp_end   = trim(surfdatadir)//trim(time_end)//'_temp_rotated.nc'

fileksat = trim(surfdatadir)//'Mad_ksat_rotated.nc' !same ksat is always used
filemask = trim(surfdatadir)//trim(time_end)//'_mask_rotated.nc' !ask Andy if he agrees with using the end time for mask
filewtd  = trim(surfdatadir)//trim(time_start)//'_wtd_rotated.nc' !water table determined from the previous time step


!Hard-coded values:

delta_xy = 120 !there are 120 30 arc-second pieces in one degree

converged = 0 !flag for the surface water loop

!initialise a few values:

upvalue    = 0.
downvalue  = 0.
leftvalue  = 0.
rightvalue = 0.
counter    = 0

iter       = 0
diff_total = 0


open(15,file=trim(time_end)//'output_GW_SW.txt') !text file to save iteration data - rather than printing on the screen


call MPI_COMM_RANK(MPI_COMM_WORLD,pid,ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numtasks,ierr)
print *,'Number of tasks = ',numtasks,'My rank = ',pid

!allocate array space, the size of the number of tasks available

allocate(nini         (1:numtasks-1))
allocate(nend         (1:numtasks-1))
allocate(domblock     (1:numtasks-1))
allocate(domblocksmall(1:numtasks-1))
               

!creating different data types and things we need to do the calculation:

call MPI_TYPE_CONTIGUOUS(numtasks-1,MPI_INTEGER,tasktype,ierr) !creates a contiguous datatype. The new datatype is called tasktype
call MPI_Type_commit(tasktype,ierr)


!divide the domain among tasks:

if(pid .eq. 0) then
    write(6,*) 'PID = 0'
    call dividedomain(n2,n3,numtasks,nini,filemask,ntotal) !Divides up the work equally among all of the ranks, by number of defined land cells.
    write (6,*) 'Dividedomain done'

    nend(numtasks-1) = n3+1 !define where each task must finish
    do n=2,numtasks-1
        nend(n-1) = nini(n) -1 !moving everything in the list along one space - the column to end in for each task.
    end do

    do n=1,numtasks-1
        call MPI_send(nini(1),1,tasktype,n,1,MPI_COMM_WORLD,ierr) !because only PID=0 has these values right now, so we have to send them out. 
        call MPI_send(nend(1),1,tasktype,n,20,MPI_COMM_WORLD,ierr) 
        call MPI_send(ntotal,1,MPI_INTEGER,n,2001,MPI_COMM_WORLD,ierr)
    end do
else
    call MPI_recv(nini(1),1,tasktype,0,1,MPI_COMM_WORLD,status,ierr) !receive what was sent above.
    call MPI_recv(nend(1),1,tasktype,0,20,MPI_COMM_WORLD,status,ierr)
    call MPI_recv(ntotal,1,MPI_INTEGER,0,2001,MPI_COMM_WORLD,status,ierr)
endif    

   

do n=1,numtasks-1
    nmax = nend(n) - nini(n) + 4 !max number of columns we have in this group

    call MPI_TYPE_CONTIGUOUS(n2*nmax,MPI_REAL,domblock(n),ierr) !creating the domblock and columntype data types
    call MPI_type_commit(domblock(n),ierr)

    nmax = nmax-3
    call MPI_TYPE_CONTIGUOUS(n2*nmax,MPI_REAL,domblocksmall(n),ierr)
    call MPI_type_commit(domblocksmall(n),ierr)
end do


call MPI_TYPE_CONTIGUOUS(n2,MPI_REAL,columntype,ierr)
call MPI_type_commit(columntype,ierr)

call MPI_TYPE_CONTIGUOUS(n2,MPI_INTEGER1,columntypeint,ierr)
call MPI_type_commit(columntypeint,ierr)


!read in all of the data:

if(pid .eq. 0) then

    allocate(h(n2,n3))
    allocate(hold(n2,n3))
    hold(:,:) = 0.

    allocate(wtdglob(n2,n3))
    wtdglob = 0.

    allocate(topo_start_read(n2,n3)) !starting time topography
     
      iret = nf90_open(filetopo_start,0,ncid) !reading in the topo
      call check_err(iret, "Unable to load file '" // filetopo_start // "'!")

      iret = nf90_inq_varid(ncid,'value',varid)
      call check_err(iret, "")

      iret = nf90_get_var(ncid,varid,topo_start_read)
      call check_err(iret, "")

      iret = nf90_close(ncid)
      call check_err(iret, "Unable to close file!")

    where(topo_start_read .le. UNDEF) topo_start_read = 0. !change undefined cells to 0

    allocate(topo_end_read(n2,n3)) !ending time topography
     
      iret = nf90_open(filetopo_end,0,ncid) !reading in the topo
      call check_err(iret, "Unable to load file '" // filetopo_end // "'!")

      iret = nf90_inq_varid(ncid,'value',varid)
      call check_err(iret, "")

      iret = nf90_get_var(ncid,varid,topo_end_read)
      call check_err(iret, "")

      iret = nf90_close(ncid)
      call check_err(iret, "Unable to close file!")

    where(topo_end_read .le. UNDEF) topo_end_read = 0. !change undefined cells to 0

    allocate(mask_read(n2,n3))  !mask to mask out ocean cells
     
      iret = nf90_open(filemask,0,ncid) !reading in the mask
      call check_err(iret, "Unable to load file '" // filemask // "'!")

      iret = nf90_inq_varid(ncid,'value',varid)
      call check_err(iret, "")

      iret = nf90_get_var(ncid,varid,mask_read)
      call check_err(iret, "")

      iret = nf90_close(ncid)
      call check_err(iret, "Unable to close file!")
      
  !    where(mask.eq.0) h = 0.
  
    allocate(wtd_read(n2,n3))

      wtd_read(:,:) = -5*mask_read !for now - in the future we should read in the wtd from the previous time step

   !   iret = nf90_open(filewtd,0,ncid) !reading in the mask
    !  call check_err(iret)
!print *,'259'
 !     iret = nf90_inq_varid(ncid,'value',varid)
  !    call check_err(iret)
!print *,'262'
 !     iret = nf90_get_var(ncid,varid,mask_read)
  !    call check_err(iret)
!print *,'265'
  !    iret = nf90_close(ncid)
 !     call check_err(iret)

    !  wtd_read = -wtd_read !it is read as a positive value, but in the model it's negative
!false - it is read in as a negative value from now on! 
          
    allocate(ksat_read(n2,n3))
        
      iret = nf90_open(fileksat,0,ncid) !reading in the ksat
      call check_err(iret, "Unable to load file '" // fileksat // "'!")

      iret = nf90_inq_varid(ncid,'value',varid)
      call check_err(iret,"")

      iret = nf90_get_var(ncid,varid,ksat_read)
      call check_err(iret,"")

      iret = nf90_close(ncid)
      call check_err(iret, "Unable to close file!")


    allocate(rech_start_read(n2,n3)) !recharge from the starting time
      
      iret = nf90_open(filerech_start,0,ncid) !reading in the recharge file
      call check_err(iret, "Unable to load file '" // filerech_start // "'!")

      iret = nf90_inq_varid(ncid,'value',varid)
      call check_err(iret, "")

      iret = nf90_get_var(ncid,varid,rech_start_read)
      call check_err(iret, "")

      iret = nf90_close(ncid)
      call check_err(iret, "Unable to close file!")

    rech_start_read = max(rech_start_read,0.)
    rech_start_read = (rech_start_read/365.) !converting to monthly - CHECK what units have I done for the new ones? I think it's per month already! 

    print *,'311'

    allocate(rech_end_read(n2,n3)) !recharge from the end time
      
      iret = nf90_open(filerech_end,0,ncid) !reading in the recharge file
      call check_err(iret, "Unable to load file '" // filerech_end // "'!")

      iret = nf90_inq_varid(ncid,'value',varid)
      call check_err(iret, "")

      iret = nf90_get_var(ncid,varid,rech_end_read)
      call check_err(iret, "")

      iret = nf90_close(ncid)
      call check_err(iret, "Unable to close file!")

    rech_end_read = max(rech_end_read,0.)
    rech_end_read = (rech_end_read/365.)


    allocate(fslope_start_read(n2,n3)) !the idea is to load in the f values for these, i.e. f = 100/(1+150*LGM_slope) so that it doesn't still have to be calculated in here. 
      
      iret = nf90_open(file_fslope_start,0,ncid) !reading in the fslope file
      call check_err(iret, "Unable to load file '" // file_fslope_start // "'!")

      iret = nf90_inq_varid(ncid,'value',varid)
      call check_err(iret, "")

      iret = nf90_get_var(ncid,varid,fslope_start_read)
      call check_err(iret, "")

      iret = nf90_close(ncid)
      call check_err(iret, "Unable to close file!")

    print *,'347'

    allocate(fslope_end_read(n2,n3))
      
      iret = nf90_open(file_fslope_end,0,ncid) !reading in the fslope file
      call check_err(iret, "Unable to load file '" // file_fslope_end // "'!")

      iret = nf90_inq_varid(ncid,'value',varid)
      call check_err(iret, "")

      iret = nf90_get_var(ncid,varid,fslope_end_read)
      call check_err(iret, "")

      iret = nf90_close(ncid)
      call check_err(iret, "Unable to close file!")

    print *,'362'

    allocate(temp_start_read(n2,n3)) !temp is used to calculate the final e-folding depth (along with fslope)
      
      iret = nf90_open(filetemp_start,0,ncid) !reading in the starting temperature file
      call check_err(iret, "Unable to load file '" // filetemp_start // "'!")

      iret = nf90_inq_varid(ncid,'value',varid)
      call check_err(iret, "")

      iret = nf90_get_var(ncid,varid,temp_start_read)
      call check_err(iret, "")

      iret = nf90_close(ncid)
      call check_err(iret, "Unable to close file!")

    allocate(temp_end_read(n2,n3))
      
      iret = nf90_open(filetemp_end,0,ncid) !reading in the ending temperature file
      call check_err(iret, "Unable to load file '" // filetemp_end // "'!")

      iret = nf90_inq_varid(ncid,'value',varid)
      call check_err(iret, "")

      iret = nf90_get_var(ncid,varid,temp_end_read)
      call check_err(iret, "")

      iret = nf90_close(ncid)
      call check_err(iret, "Unable to close file!")


!send and receive the data from above

    do n=1,numtasks-2
        call MPI_send(mask_read(1,nini(n)-1),1,domblock(n),n, 1,MPI_COMM_WORLD,ierr)
        call MPI_send(wtd_read (1,nini(n)-1),1,domblock(n),n, 2,MPI_COMM_WORLD,ierr)       
        call MPI_send(ksat_read(1,nini(n)-1),1,domblock(n),n, 3,MPI_COMM_WORLD,ierr)
        call MPI_send(wtd_read (1,nini(n)-1),1,domblock(n),n, 4,MPI_COMM_WORLD,ierr)  
        call MPI_send(hold     (1,nini(n)-1),1,domblock(n),n,10,MPI_COMM_WORLD,ierr)
    end do

    call MPI_send(mask_read(1,nini(numtasks-1)-1),1,domblocksmall(numtasks-1),numtasks-1, 1,MPI_COMM_WORLD,ierr)
    call MPI_send(wtd_read (1,nini(numtasks-1)-1),1,domblocksmall(numtasks-1),numtasks-1, 2,MPI_COMM_WORLD,ierr)
    call MPI_send(ksat_read(1,nini(numtasks-1)-1),1,domblocksmall(numtasks-1),numtasks-1, 3,MPI_COMM_WORLD,ierr)
    call MPI_send(wtd_read (1,nini(numtasks-1)-1),1,domblocksmall(numtasks-1),numtasks-1, 4,MPI_COMM_WORLD,ierr)
    call MPI_send(hold     (1,nini(numtasks-1)-1),1,domblocksmall(numtasks-1),numtasks-1,10,MPI_COMM_WORLD,ierr)
  
    call MPI_send(mask_read(1,1),1,columntype,numtasks-1, 1,MPI_COMM_WORLD,ierr)
    call MPI_send(wtd_read (1,1),1,columntype,numtasks-1, 2,MPI_COMM_WORLD,ierr)
    call MPI_send(ksat_read(1,1),1,columntype,numtasks-1, 3,MPI_COMM_WORLD,ierr)
    call MPI_send(wtd_read (1,1),1,columntype,numtasks-1, 4,MPI_COMM_WORLD,ierr)
    call MPI_send(hold     (1,1),1,columntype,numtasks-1,10,MPI_COMM_WORLD,ierr)

    call MPI_send(mask_read(1,2),1,columntype,numtasks-1, 1,MPI_COMM_WORLD,ierr)
    call MPI_send(wtd_read (1,2),1,columntype,numtasks-1, 2,MPI_COMM_WORLD,ierr)
    call MPI_send(ksat_read(1,2),1,columntype,numtasks-1, 3,MPI_COMM_WORLD,ierr)
    call MPI_send(wtd_read (1,2),1,columntype,numtasks-1, 4,MPI_COMM_WORLD,ierr)
    call MPI_send(hold     (1,2),1,columntype,numtasks-1,10,MPI_COMM_WORLD,ierr)


    call MPI_send(mask_read(1,3),1,columntype,numtasks-1, 1,MPI_COMM_WORLD,ierr)
    call MPI_send(wtd_read (1,3),1,columntype,numtasks-1, 2,MPI_COMM_WORLD,ierr)
    call MPI_send(ksat_read(1,3),1,columntype,numtasks-1, 3,MPI_COMM_WORLD,ierr)
    call MPI_send(wtd_read (1,3),1,columntype,numtasks-1, 4,MPI_COMM_WORLD,ierr)
    call MPI_send(hold     (1,3),1,columntype,numtasks-1,10,MPI_COMM_WORLD,ierr)

    deallocate(mask_read)
    deallocate(wtd_read)
    deallocate(ksat_read)

else
  
    nmax = nend(pid) - nini(pid) +4
    allocate(wtd       (n2,nmax))
    allocate(wtdnew    (n2,nmax))
    allocate(ksat      (n2,nmax))
    allocate(landmask  (n2,nmax))
    allocate(head      (n2,nmax))
    allocate(kcell     (n2,nmax))
    allocate(qlat_north(n2,nmax))
    allocate(qlat_south(n2,nmax))
    allocate(qlat_east (n2,nmax))
    allocate(qlat_west (n2,nmax))
    allocate(hold_read (n2,nmax))
 
    allocate(topo_sent  (n2,nmax))
    allocate(rech_sent  (n2,nmax))
    allocate(fdepth_sent(n2,nmax))

    write(6,*)'allocated all'

      if(pid.lt.numtasks-1)then

        call MPI_recv(landmask (1,1),1,domblock(pid),0,1,MPI_COMM_WORLD,status,ierr)
        call MPI_recv(wtd      (1,1),1,domblock(pid),0,2,MPI_COMM_WORLD,status,ierr)
        call MPI_recv(ksat     (1,1),1,domblock(pid),0,3,MPI_COMM_WORLD,status,ierr)
        call MPI_recv(wtdnew   (1,1),1,domblock(pid),0,4,MPI_COMM_WORLD,status,ierr)
        call MPI_recv(hold_read(1,1),1,domblock(pid),0,10,MPI_COMM_WORLD,status,ierr)

   else

        call MPI_recv(landmask (1,1),1,domblocksmall(pid),0, 1,MPI_COMM_WORLD,status,ierr)
        call MPI_recv(wtd      (1,1),1,domblocksmall(pid),0, 2,MPI_COMM_WORLD,status,ierr)
        call MPI_recv(ksat     (1,1),1,domblocksmall(pid),0, 3,MPI_COMM_WORLD,status,ierr)
        call MPI_recv(wtdnew   (1,1),1,domblocksmall(pid),0, 4,MPI_COMM_WORLD,status,ierr)
        call MPI_recv(hold_read(1,1),1,domblocksmall(pid),0,10,MPI_COMM_WORLD,status,ierr)

        call MPI_recv(landmask (1,nmax-2),1,columntype,0, 1,MPI_COMM_WORLD,status,ierr)
        call MPI_recv(wtd      (1,nmax-2),1,columntype,0, 2,MPI_COMM_WORLD,status,ierr)
        call MPI_recv(ksat     (1,nmax-2),1,columntype,0, 3,MPI_COMM_WORLD,status,ierr)
        call MPI_recv(wtdnew   (1,nmax-2),1,columntype,0, 4,MPI_COMM_WORLD,status,ierr)
        call MPI_recv(hold_read(1,nmax-2),1,columntype,0,10,MPI_COMM_WORLD,status,ierr)

        call MPI_recv(landmask (1,nmax-1),1,columntype,0, 1,MPI_COMM_WORLD,status,ierr)
        call MPI_recv(wtd      (1,nmax-1),1,columntype,0, 2,MPI_COMM_WORLD,status,ierr)
        call MPI_recv(ksat     (1,nmax-1),1,columntype,0, 3,MPI_COMM_WORLD,status,ierr)
        call MPI_recv(wtdnew   (1,nmax-1),1,columntype,0, 4,MPI_COMM_WORLD,status,ierr)
        call MPI_recv(hold_read(1,nmax-1),1,columntype,0,10,MPI_COMM_WORLD,status,ierr)

        call MPI_recv(landmask (1,nmax),1,columntype,0, 1,MPI_COMM_WORLD,status,ierr)
        call MPI_recv(wtd      (1,nmax),1,columntype,0, 2,MPI_COMM_WORLD,status,ierr)
        call MPI_recv(ksat     (1,nmax),1,columntype,0, 3,MPI_COMM_WORLD,status,ierr)
        call MPI_recv(wtdnew   (1,nmax),1,columntype,0, 4,MPI_COMM_WORLD,status,ierr)
        call MPI_recv(hold_read(1,nmax),1,columntype,0,10,MPI_COMM_WORLD,status,ierr)
    endif
   
    write(6,*) 'variables received'
endif


!a little pre-loop math to get latitudes, total time in a time-step, and cell areas:


dy = 6370000.*pi/(180.*delta_xy) !radius of the earth * pi / number of possible cells in the y-direction. This should equal the height of each cell in the N-S direction.
dx=dy

if(pid .gt. 0) then
    allocate(xlat(n2))
    allocate(area(n2))
    allocate(alpha(n2))

    do j=1,n2 !changing area of cell depending on its latitude
        xlat(j)     = (float(j-2)/delta_xy+SEDGE)*pi/180. !latitude in radians
        xs          = (float(2*(j-1))/(delta_xy*2.)+SEDGE)*pi/180. !Latitude on southern cell edge in radians
        xn          = (float(2*(j+1))/(delta_xy*2.)+SEDGE)*pi/180. !latitude on northern cell edge in radians
        placeholder = dy*6370000.*(sin(xn)-sin(xs))/2.
        area(j)     = dy*6370000.*(sin(xn)-sin(xs))/2. !final cell area for that latitude: trapezoid dy * dx

        alpha(j) = 0.5*deltat/placeholder 
    end do
end if


!calculate the e-folding depth, and send this plus topo and rech

if (pid.eq.0) then

    allocate(fdepth_start(n2,n3))
     
    do j=1,n3
        do i=1,n2
            if (temp_start_read(i,j) .gt. -5) then
                fdepth_start(i,j) = fslope_start_read(i,j)
            elseif (temp_start_read(i,j) .lt. -14) then
                fdepth_start(i,j) = fslope_start_read(i,j) * (0.17+0.005*temp_start_read(i,j))
            else
                fdepth_start(i,j) = fslope_start_read(i,j) * (1.5 + 0.1*temp_start_read(i,j))
            endif
        end do
    end do

    do n=1,numtasks-2
        call MPI_send(topo_start_read(1,nini(n)-1),1,domblock(n),n,5,MPI_COMM_WORLD,ierr)
        call MPI_send(rech_start_read(1,nini(n)-1),1,domblock(n),n,6,MPI_COMM_WORLD,ierr)
        call MPI_send(fdepth_start   (1,nini(n)-1),1,domblock(n),n,7,MPI_COMM_WORLD,ierr)
    end do

    call MPI_send(topo_start_read(1,nini(numtasks-1)-1),1,domblocksmall(numtasks-1),numtasks-1,5,MPI_COMM_WORLD,ierr)
    call MPI_send(rech_start_read(1,nini(numtasks-1)-1),1,domblocksmall(numtasks-1),numtasks-1,6,MPI_COMM_WORLD,ierr)
    call MPI_send(fdepth_start   (1,nini(numtasks-1)-1),1,domblocksmall(numtasks-1),numtasks-1,7,MPI_COMM_WORLD,ierr)

    call MPI_send(topo_start_read(1,1),1,columntype,numtasks-1,5,MPI_COMM_WORLD,ierr)
    call MPI_send(rech_start_read(1,1),1,columntype,numtasks-1,6,MPI_COMM_WORLD,ierr)
    call MPI_send(fdepth_start   (1,1),1,columntype,numtasks-1,7,MPI_COMM_WORLD,ierr)

    call MPI_send(topo_start_read(1,2),1,columntype,numtasks-1,5,MPI_COMM_WORLD,ierr)
    call MPI_send(rech_start_read(1,2),1,columntype,numtasks-1,6,MPI_COMM_WORLD,ierr)
    call MPI_send(fdepth_start   (1,2),1,columntype,numtasks-1,7,MPI_COMM_WORLD,ierr)

    call MPI_send(topo_start_read(1,3),1,columntype,numtasks-1,5,MPI_COMM_WORLD,ierr)
    call MPI_send(rech_start_read(1,3),1,columntype,numtasks-1,6,MPI_COMM_WORLD,ierr)
    call MPI_send(fdepth_start   (1,3),1,columntype,numtasks-1,7,MPI_COMM_WORLD,ierr)


 
    write(6,*)'sent initial'

    deallocate(fdepth_start)

elseif(pid.lt.numtasks-1)then
  
    call MPI_recv(topo_sent  (1,1),1,domblock(pid),0,5,MPI_COMM_WORLD,status,ierr) !receiving everthing that was sent above     
    call MPI_recv(rech_sent  (1,1),1,domblock(pid),0,6,MPI_COMM_WORLD,status,ierr)
    call MPI_recv(fdepth_sent(1,1),1,domblock(pid),0,7,MPI_COMM_WORLD,status,ierr)

else

    call MPI_recv(topo_sent  (1,1),1,domblocksmall(pid),0,5,MPI_COMM_WORLD,status,ierr)
    call MPI_recv(rech_sent  (1,1),1,domblocksmall(pid),0,6,MPI_COMM_WORLD,status,ierr)
    call MPI_recv(fdepth_sent(1,1),1,domblocksmall(pid),0,7,MPI_COMM_WORLD,status,ierr)

    call MPI_recv(topo_sent  (1,nmax-2),1,columntype,0,5,MPI_COMM_WORLD,status,ierr)
    call MPI_recv(rech_sent  (1,nmax-2),1,columntype,0,6,MPI_COMM_WORLD,status,ierr)
    call MPI_recv(fdepth_sent(1,nmax-2),1,columntype,0,7,MPI_COMM_WORLD,status,ierr)

    call MPI_recv(topo_sent  (1,nmax-1),1,columntype,0,5,MPI_COMM_WORLD,status,ierr)
    call MPI_recv(rech_sent  (1,nmax-1),1,columntype,0,6,MPI_COMM_WORLD,status,ierr)
    call MPI_recv(fdepth_sent(1,nmax-1),1,columntype,0,7,MPI_COMM_WORLD,status,ierr)

    call MPI_recv(topo_sent  (1,nmax),1,columntype,0,5,MPI_COMM_WORLD,status,ierr)
    call MPI_recv(rech_sent  (1,nmax),1,columntype,0,6,MPI_COMM_WORLD,status,ierr)
    call MPI_recv(fdepth_sent(1,nmax),1,columntype,0,7,MPI_COMM_WORLD,status,ierr)


    write(6,*)'received initial'
endif


!allocate the arrays for the mergesort done later
if(pid.ne.0)then
    nmax = nend(pid) - nini(pid) +3 
    Merge_size = n2*nmax           !Note that this may be faster with a smaller value. Largely because the actual sort takes longer, although the number of iterations may sometimes also be higher despite that not making much sense
    !Also note if unsure - a smaller value may only partially sort the array, but a too-large value may fill with zeros! 
    allocate(T((Merge_size+1)/2))
endif

!ready to start the loops ************************************************************************************

print *,'about to start the groundwater loop'

print *,'(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((',nini(pid),nend(pid),nmax,pid

GROUNDWATER: DO while(iter<iterations) 

    iter = iter + 1

    if (pid .eq. 1) then !This is not completely necessary, but it's nice to see how far the run is. 
        write (6,*) 'iter number',iter
    endif


    if(mod(float(iter),3650.).eq.0.) then !update the changing topography and recharge, here being done every 10 years
        write(6,*)'pid',pid
        if (pid.eq.0) then !here I am going to do the adjustments to the topo,rech,fslope,temp. 
            allocate(topo_now(n2,n3))
            topo_now = (topo_start_read * (1- iter/iterations)) + (topo_end_read * (iter/iterations))
            write (6,*) 'allocated topo'

            allocate(rech_now(n2,n3))
            rech_now = (rech_start_read * (1- iter/iterations)) + (rech_end_read * (iter/iterations))
            write (6,*) 'allocated rech'

            allocate(fslope_now(n2,n3))
            fslope_now = (fslope_start_read * (1- iter/iterations)) + (fslope_end_read * (iter/iterations))
            write (6,*) 'allocated fslope'

            allocate(temp_now(n2,n3))
            temp_now = (temp_start_read * (1- iter/iterations)) + (temp_end_read * (iter/iterations))
            write (6,*) 'allocated temp'

            allocate(fdepth_now(n2,n3))
            do j=1,n3
                do i=1,n2
                    if (temp_now(i,j) .gt. -5) then
                        fdepth_now(i,j) = fslope_now(i,j)
                    elseif (temp_now(i,j) .lt. -14) then
                        fdepth_now(i,j) = fslope_now(i,j) * (0.17+0.005*temp_now(i,j))
                    else
                        fdepth_now(i,j) = fslope_now(i,j) * (1.5 + 0.1*temp_now(i,j))
                    endif
                end do
            end do


            do n=1,numtasks-2
                call MPI_send(topo_now  (1,nini(n)-1),1,domblock(n),n, 9,MPI_COMM_WORLD,ierr)
                call MPI_send(rech_now  (1,nini(n)-1),1,domblock(n),n,10,MPI_COMM_WORLD,ierr)
                call MPI_send(fdepth_now(1,nini(n)-1),1,domblock(n),n,11,MPI_COMM_WORLD,ierr)
            end do

            call MPI_send(topo_now  (1,nini(numtasks-1)-1),1,domblocksmall(numtasks-1),numtasks-1, 9,MPI_COMM_WORLD,ierr)
            call MPI_send(rech_now  (1,nini(numtasks-1)-1),1,domblocksmall(numtasks-1),numtasks-1,10,MPI_COMM_WORLD,ierr)
            call MPI_send(fdepth_now(1,nini(numtasks-1)-1),1,domblocksmall(numtasks-1),numtasks-1,11,MPI_COMM_WORLD,ierr)

            call MPI_send(topo_now  (1,1),1,columntype,numtasks-1, 9,MPI_COMM_WORLD,ierr)
            call MPI_send(rech_now  (1,1),1,columntype,numtasks-1,10,MPI_COMM_WORLD,ierr)
            call MPI_send(fdepth_now(1,1),1,columntype,numtasks-1,11,MPI_COMM_WORLD,ierr)
 
            call MPI_send(topo_now  (1,2),1,columntype,numtasks-1, 9,MPI_COMM_WORLD,ierr)
            call MPI_send(rech_now  (1,2),1,columntype,numtasks-1,10,MPI_COMM_WORLD,ierr)
            call MPI_send(fdepth_now(1,2),1,columntype,numtasks-1,11,MPI_COMM_WORLD,ierr)

            call MPI_send(topo_now  (1,3),1,columntype,numtasks-1, 9,MPI_COMM_WORLD,ierr)
            call MPI_send(rech_now  (1,3),1,columntype,numtasks-1,10,MPI_COMM_WORLD,ierr)
            call MPI_send(fdepth_now(1,3),1,columntype,numtasks-1,11,MPI_COMM_WORLD,ierr)

 
            deallocate(topo_now)
            deallocate(rech_now)
            deallocate(fslope_now)
            deallocate(temp_now)
            deallocate(fdepth_now)


        elseif(pid.lt.numtasks-1)then
            call MPI_recv(topo_sent  (1,1),1,domblock(pid),0, 9,MPI_COMM_WORLD,status,ierr) !receiving everthing that was sent above     
            call MPI_recv(rech_sent  (1,1),1,domblock(pid),0,10,MPI_COMM_WORLD,status,ierr)
            call MPI_recv(fdepth_sent(1,1),1,domblock(pid),0,11,MPI_COMM_WORLD,status,ierr)

        else
            call MPI_recv(topo_sent  (1,1),1,domblocksmall(pid),0, 9,MPI_COMM_WORLD,status,ierr) 
            call MPI_recv(rech_sent  (1,1),1,domblocksmall(pid),0,10,MPI_COMM_WORLD,status,ierr)
            call MPI_recv(fdepth_sent(1,1),1,domblocksmall(pid),0,11,MPI_COMM_WORLD,status,ierr)
 
            call MPI_recv(topo_sent  (1,nmax-1),1,columntype,0, 9,MPI_COMM_WORLD,status,ierr) 
            call MPI_recv(rech_sent  (1,nmax-1),1,columntype,0,10,MPI_COMM_WORLD,status,ierr)
            call MPI_recv(fdepth_sent(1,nmax-1),1,columntype,0,11,MPI_COMM_WORLD,status,ierr)
 
     
            call MPI_recv(topo_sent  (1,nmax),1,columntype,0, 9,MPI_COMM_WORLD,status,ierr) 
            call MPI_recv(rech_sent  (1,nmax),1,columntype,0,10,MPI_COMM_WORLD,status,ierr)
            call MPI_recv(fdepth_sent(1,nmax),1,columntype,0,11,MPI_COMM_WORLD,status,ierr)
 
            call MPI_recv(topo_sent  (1,nmax+1),1,columntype,0, 9,MPI_COMM_WORLD,status,ierr) 
            call MPI_recv(rech_sent  (1,nmax+1),1,columntype,0,10,MPI_COMM_WORLD,status,ierr)
            call MPI_recv(fdepth_sent(1,nmax+1),1,columntype,0,11,MPI_COMM_WORLD,status,ierr)
        endif
    endif

!what about issues with using the first and last cell but still checking surrounding cells here?

    IF (pid .gt. 0) then     
        nmax = nend(pid) - nini(pid) + 3
        do j=1,nmax  +1
            do i=1,n2
                if(fdepth_sent(i,j) .gt. 0. ) then 
                    wtd(i,j) = wtd(i,j) + rech_sent(i,j) !add the month's recharge to the cell
                    head(i,j) = topo_sent(i,j) + wtd(i,j) !gives the water table height (starts off the same as land surface if a wtd isn't loaded in), which = head

                    if(wtd(i,j) .lt. -1.5) then !work out hydraulic conductivity for each cell
                        kcell(i,j) = fdepth_sent(i,j) *ksat(i,j)*exp((wtd(i,j)+1.5)/fdepth_sent(i,j)) !This is equation S6 from the paper
                    elseif(wtd(i,j) .le. 0) then
                        kcell(i,j) = ksat(i,j)*(wtd(i,j)+1.5+fdepth_sent(i,j)) !equation S4 from the paper 
                    else
                        kcell(i,j) = ksat(i,j)*(0+1.5+fdepth_sent(i,j)) !maxes out when water is at the surface, to avoid instabilities in surface water movement.
                    endif
                endif
            end do
        end do

        write(6,*) 'kcell done',nmax,pid,nini(pid)

        do j=2,nmax!-1 !-2 !I think there may be some issues with the very last processor here overshooting
            do i=2,n2-1
                if(landmask(i,j) .gt. 0 ) then
                    ! alpha = time step / area
                 ! Approx. small distance between N and S edges, so just use cos (xlat(j)) at midpoint
                 ! N-S, then W-E
                 !finite difference method
                 wtdnew(i,j) = wtd(i,j) + kcell(i,j) * alpha(j) / cos(xlat(j)) &
                                    * (head(i,j-1) - 2*head(i,j) + head(i,j+1))/dx &
                               + kcell(i,j) * alpha(j) * cos(xlat(j)) &
                                    * (head(i-1,j) - 2*head(i,j) + head(i+1,j))/dx 
 
                endif
            end do
        end do

 
        wtd = wtdnew !update the whole wtd array

!sending and receiving the lines on either side of each section to allow for flow across those lines:

        if(pid.eq.1)then
            call MPI_recv(wtd(1,nmax),  1,columntype,pid+1,     8,MPI_COMM_WORLD,status,ierr)
            call MPI_recv(wtd(1,1),     1,columntype,pid+1,     9,MPI_COMM_WORLD,status,ierr)
            call MPI_send(wtd(1,2),     1,columntype,numtasks-1,8,MPI_COMM_WORLD,ierr)
            call MPI_send(wtd(1,nmax-1),1,columntype,numtasks-1,9,MPI_COMM_WORLD,ierr)

        elseif(pid.eq.numtasks-1)then
            if(mod(pid,2).eq.0)then
                call MPI_recv(wtd(1,nmax),  1,columntype,1,    8,MPI_COMM_WORLD,status,ierr)
                call MPI_recv(wtd(1,1),     1,columntype,1,    9,MPI_COMM_WORLD,status,ierr)
                call MPI_send(wtd(1,2),     1,columntype,pid-1,8,MPI_COMM_WORLD,ierr)   
                call MPI_send(wtd(1,nmax-1),1,columntype,pid-1,9,MPI_COMM_WORLD,ierr)

            else
                call MPI_send(wtd(1,2),     1,columntype,pid-1,8,MPI_COMM_WORLD,ierr)   
                call MPI_send(wtd(1,nmax-1),1,columntype,pid-1,9,MPI_COMM_WORLD,ierr)
                call MPI_recv(wtd(1,nmax),  1,columntype,1,    8,MPI_COMM_WORLD,status,ierr)
                call MPI_recv(wtd(1,1),     1,columntype,1,    9,MPI_COMM_WORLD,status,ierr)
            endif

        else
            if(mod(pid,2).eq.0)then
                call MPI_recv(wtd(1,nmax),  1,columntype,pid+1,8,MPI_COMM_WORLD,status,ierr)
                call MPI_recv(wtd(1,1),     1,columntype,pid+1,9,MPI_COMM_WORLD,status,ierr)
                call MPI_send(wtd(1,2),     1,columntype,pid-1,8,MPI_COMM_WORLD,ierr)
                call MPI_send(wtd(1,nmax-1),1,columntype,pid-1,9,MPI_COMM_WORLD,ierr)
            else
                call MPI_send(wtd(1,2),     1,columntype,pid-1,8,MPI_COMM_WORLD,ierr)
                call MPI_send(wtd(1,nmax-1),1,columntype,pid-1,9,MPI_COMM_WORLD,ierr)
                call MPI_recv(wtd(1,nmax),  1,columntype,pid+1,8,MPI_COMM_WORLD,status,ierr)
                call MPI_recv(wtd(1,1),     1,columntype,pid+1,9,MPI_COMM_WORLD,status,ierr)
            endif
        endif
    ENDIF 



    if(mod((iter),365).eq.0.) then
        print *,'starting with the surface water loop!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

!reset variables to 0
        converged  = 0
        counter    = 0
        diff_total = 0
      
        SURFACE: do while(converged .eq. 0) 
           
            counter = counter + 1


            if(counter .eq.5000)then    !select threshold here. Consider a better way to threshold or a max number of iterations if it isn't reaching that. BUT note it'll sometimes level out for a while and then still be able to take a big jump to improvement, so it's not so simple as just looking for when it levels out! 
                print *,'reached max iterations',diff_total
                converged = 1
            endif

            if(counter .gt.5 .and. diff_total .lt.1)then
                print *,'success',diff_total
                converged = 1
            endif


            if (pid .eq. 0) then

                print *,'Surface water counter',counter
                print *,'max',diff_total

 
            else !Any PID that is not 0

                hold_read = wtd
                nmax      = nend(pid) - nini(pid) +3  


                allocate(hz_read(n2,nmax+1))
                allocate(diff(n2,nmax+1))

                hz_read = max(topo_sent,topo_sent+wtd) !water moves on the surface, so it must be AT LEAST the height of the topography, or if there is surface water then it's topo + wtd

                allocate(hz_1D(n2*nmax))  !Convert the array to a 1D for the mergesort

                do i=1,n2
                    hz_1D(((i-1)*nmax)+1:i*nmax) = hz_read(i,:) !unpack hz into a 1D array
                end do

                allocate(arr(2,n2*nmax))    !Create an array of the indices, for picking which cells to process first
                do row=1,n2
                    arr(1,((row-1)*nmax)+1:row*nmax)=row
                end do

                do row=1,n2
                    do col = 1,nmax
                        arr(2,(row-1)*nmax+col)=col
                    end do
                end do

                call MergeSort(hz_1D,Merge_size,T,arr)  !Sort to obtain order for processing the array - water flows from the highest cells first


                COLS1: do i=2,nmax
                ROWS1: do j=2,n2
                    row = arr(1,n2*nmax-(j-1)-(i-1)*n2) !get the next item in the sorted list to be processed, from the sorted index array
                    col = arr(2,n2*nmax-(j-1)-(i-1)*n2)

                    if(col.ge.nmax-1)then !Doing the end two columns separately, so skip them here
                      CYCLE
                    endif

                    if(col.le.2)then !Doing the end two columns separately, so skip them here
                      CYCLE
                    endif

                    if(row.eq.1 .or.row.eq.n2)then !Can't correctly process the very edge rows
                      CYCLE
                    endif

                    if(wtd(row,col).le.0)then  !there is no surface water to move
                      CYCLE
                    endif

                    if(landmask(row,col) .eq. 0) then !we are in the ocean
                      wtd(row,col) = 0
                    endif

                    if(hz_read(row,col) .le. hz_read(row+1,col) .and. hz_read(row,col) &
                        .le. hz_read(row-1,col) .and. hz_read(row,col) .le. hz_read(row,col+1) .and. &
                        hz_read(row,col) .le. hz_read(row,col-1)) then !this is the lowest cell in the area, no water can leave
                      CYCLE
                    else
                      upvalue    = hz_read(row,col) - hz_read(row-1,col) !find the steepest direction, in which water would move
                      downvalue  = hz_read(row,col) - hz_read(row+1,col)
                      leftvalue  = hz_read(row,col) - hz_read(row,col-1)
                      rightvalue = hz_read(row,col) - hz_read(row,col+1)
      
                      if(max(upvalue,downvalue,leftvalue,rightvalue) .le. 0.) then     !should have already eliminated this option, but in case, it's the lowest local cell    
                        CYCLE   
                      elseif(max(upvalue,downvalue,leftvalue,rightvalue) .eq. upvalue) then  !choose direction which is steepest
                          water  = min(wtd(row,col)*area(row),upvalue*area(row)/2.)          !water is the minumum of the total water available, or of half the difference between 2 cells
                          water1 = water/area(row-1)                                         !get height of water change in the giving and receiving cells
                          water2 = water/area(row)

                          wtd(row,col)   = wtd(row,col) - water2                             !update water table in the cells
                          wtd(row-1,col) = wtd(row-1,col) + water1

                          hz_read(row,col)   = hz_read(row,col) - water2                     !update hz for the rest of the calculations
                          hz_read(row-1,col) = hz_read(row-1,col) + water1

                      elseif(max(upvalue,downvalue,leftvalue,rightvalue) .eq. downvalue) then
                          water          = min(wtd(row,col)*area(row),downvalue*area(row)/2.)    
                          water1         = water/area(row+1)  
                          water2         = water/area(row)
                          wtd(row,col)   = wtd(row,col) - water2
                          wtd(row+1,col) = wtd(row+1,col) + water1

                          hz_read(row,col) = hz_read(row,col) - water2
                          hz_read(row+1,col) = hz_read(row+1,col) + water1

                      elseif(max(upvalue,downvalue,leftvalue,rightvalue) .eq. rightvalue) then
                          water          = min(wtd(row,col)*area(row),rightvalue*area(row)/2.)    
                          water1         = water/area(row)  
                          water2         = water/area(row)
                          wtd(row,col)   = wtd(row,col) - water2
                          wtd(row,col+1) = wtd(row,col+1) + water1

                          hz_read(row,col) = hz_read(row,col) - water2
                          hz_read(row,col+1) = hz_read(row,col+1) + water1

                      elseif(max(upvalue,downvalue,leftvalue,rightvalue) .eq. leftvalue) then
                          water          = min(wtd(row,col)*area(row),leftvalue*area(row)/2.)    
                          water1         = water/area(row)  
                          water2         = water/area(row)
                          wtd(row,col)   = wtd(row,col) - water2
                          wtd(row,col-1) = wtd(row,col-1) + water1

                          hz_read(row,col) = hz_read(row,col) - water2
                          hz_read(row,col-1) = hz_read(row,col-1) + water1

                      endif
                         
                    endif

                end do ROWS1
                end do COLS1


                if(pid.eq.1)then        !send & receive the edge columns
                    call MPI_recv(wtd(1,nmax+1),1,columntype,pid+1,     8,MPI_COMM_WORLD,status,ierr)
                    call MPI_recv(wtd(1,nmax),  1,columntype,pid+1,     9,MPI_COMM_WORLD,status,ierr)
                    call MPI_send(wtd(1,3),     1,columntype,numtasks-1,8,MPI_COMM_WORLD,ierr)
                    call MPI_send(wtd(1,2),     1,columntype,numtasks-1,9,MPI_COMM_WORLD,ierr)

                    call MPI_recv(hz_read(1,nmax+1),1,columntype,pid+1,     8,MPI_COMM_WORLD,status,ierr)
                    call MPI_recv(hz_read(1,nmax),  1,columntype,pid+1,     9,MPI_COMM_WORLD,status,ierr)
                    call MPI_send(hz_read(1,3),     1,columntype,numtasks-1,8,MPI_COMM_WORLD,ierr)
                    call MPI_send(hz_read(1,2),     1,columntype,numtasks-1,9,MPI_COMM_WORLD,ierr)

                elseif(pid.eq.numtasks-1)then
                    if(mod(pid,2).eq.0)then
                      call MPI_recv(wtd(1,nmax+1),1,columntype,1,    8,MPI_COMM_WORLD,status,ierr)
                      call MPI_recv(wtd(1,nmax),  1,columntype,1,    9,MPI_COMM_WORLD,status,ierr)
                      call MPI_send(wtd(1,3),     1,columntype,pid-1,8,MPI_COMM_WORLD,ierr)   
                      call MPI_send(wtd(1,2),     1,columntype,pid-1,9,MPI_COMM_WORLD,ierr)
           
                      call MPI_recv(hz_read(1,nmax+1),1,columntype,1,    8,MPI_COMM_WORLD,status,ierr)
                      call MPI_recv(hz_read(1,nmax),  1,columntype,1,    9,MPI_COMM_WORLD,status,ierr)
                      call MPI_send(hz_read(1,3),     1,columntype,pid-1,8,MPI_COMM_WORLD,ierr)   
                      call MPI_send(hz_read(1,2),     1,columntype,pid-1,9,MPI_COMM_WORLD,ierr)
                    else
                      call MPI_send(wtd(1,3),     1,columntype,pid-1,8,MPI_COMM_WORLD,ierr)   
                      call MPI_send(wtd(1,2),     1,columntype,pid-1,9,MPI_COMM_WORLD,ierr)
                      call MPI_recv(wtd(1,nmax+1),1,columntype,1,    8,MPI_COMM_WORLD,status,ierr)
                      call MPI_recv(wtd(1,nmax),  1,columntype,1,    9,MPI_COMM_WORLD,status,ierr)
           
                      call MPI_send(hz_read(1,3),     1,columntype,pid-1,8,MPI_COMM_WORLD,ierr)   
                      call MPI_send(hz_read(1,2),     1,columntype,pid-1,9,MPI_COMM_WORLD,ierr)
                      call MPI_recv(hz_read(1,nmax+1),1,columntype,1,    8,MPI_COMM_WORLD,status,ierr)
                      call MPI_recv(hz_read(1,nmax),  1,columntype,1,    9,MPI_COMM_WORLD,status,ierr)
                    endif
                else
                    if(mod(pid,2).eq.0)then
                        call MPI_recv(wtd(1,nmax+1),1,columntype,pid+1,8,MPI_COMM_WORLD,status,ierr)
                        call MPI_recv(wtd(1,nmax),1,columntype,pid+1,9,MPI_COMM_WORLD,status,ierr)
                        call MPI_send(wtd(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,ierr)
                        call MPI_send(wtd(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,ierr)

                        call MPI_recv(hz_read(1,nmax+1),1,columntype,pid+1,8,MPI_COMM_WORLD,status,ierr)
                        call MPI_recv(hz_read(1,nmax),1,columntype,pid+1,9,MPI_COMM_WORLD,status,ierr)
                        call MPI_send(hz_read(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,ierr)
                        call MPI_send(hz_read(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,ierr)

                    else
                        call MPI_send(wtd(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,ierr)
                        call MPI_send(wtd(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,ierr)
                        call MPI_recv(wtd(1,nmax+1),1,columntype,pid+1,8,MPI_COMM_WORLD,status,ierr)
                        call MPI_recv(wtd(1,nmax),1,columntype,pid+1,9,MPI_COMM_WORLD,status,ierr)

                        call MPI_send(hz_read(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,ierr)
                        call MPI_send(hz_read(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,ierr)
                        call MPI_recv(hz_read(1,nmax+1),1,columntype,pid+1,8,MPI_COMM_WORLD,status,ierr)
                        call MPI_recv(hz_read(1,nmax),1,columntype,pid+1,9,MPI_COMM_WORLD,status,ierr)
                    endif
                endif


                COLS3: do col=nmax-1,nmax
                ROWS3: do row=2,n2
                    if(row.eq.1 .or.row.eq.n2)then 
                        CYCLE
                    endif


                    if(wtd(row,col).le.0)then
                      CYCLE
                    endif

                    if(landmask(row,col) .eq. 0) then
                      wtd(row,col) = 0
                    endif

                          
                    if(hz_read(row,col) .le. hz_read(row+1,col) .and. hz_read(row,col) &
                    .le. hz_read(row-1,col) .and. hz_read(row,col) .le. hz_read(row,col+1) .and. &
                    hz_read(row,col) .le. hz_read(row,col-1)) then
                        CYCLE
                    else

                        upvalue    = hz_read(row,col) - hz_read(row-1,col)
                        downvalue  = hz_read(row,col) - hz_read(row+1,col)
                        leftvalue  = hz_read(row,col) - hz_read(row,col-1)
                        rightvalue = hz_read(row,col) - hz_read(row,col+1)

                        if(max(upvalue,downvalue,leftvalue,rightvalue) .le. 0.) then         
                            CYCLE   
                        elseif(max(upvalue,downvalue,leftvalue,rightvalue) .eq. upvalue) then
                            water          = min(wtd(row,col)*area(row),upvalue*area(row)/2.)   
                            water1         = water/area(row-1)  
                            water2         = water/area(row)
                            wtd(row,col)   = wtd(row,col) - water2
                            wtd(row-1,col) = wtd(row-1,col) + water1

                            hz_read(row,col)   = hz_read(row,col) - water2
                            hz_read(row-1,col) = hz_read(row-1,col) + water1


                        elseif(max(upvalue,downvalue,leftvalue,rightvalue) .eq. downvalue) then
                            water          = min(wtd(row,col)*area(row),downvalue*area(row)/2.)   
                            water1         = water/area(row+1)  
                            water2         = water/area(row)
                            wtd(row,col)   = wtd(row,col) - water2
                            wtd(row+1,col) = wtd(row+1,col) + water1

                            hz_read(row,col)   = hz_read(row,col) - water2
                            hz_read(row+1,col) = hz_read(row+1,col) + water1

                        elseif(max(upvalue,downvalue,leftvalue,rightvalue) .eq. rightvalue) then
                            water          = min(wtd(row,col)*area(row),rightvalue*area(row)/2.)   
                            water1         = water/area(row)  
                            water2         = water/area(row)
                            wtd(row,col)   = wtd(row,col) - water2
                            wtd(row,col+1) = wtd(row,col+1) + water1

                            hz_read(row,col)   = hz_read(row,col) - water2
                            hz_read(row,col+1) = hz_read(row,col+1) + water1

                        elseif(max(upvalue,downvalue,leftvalue,rightvalue) .eq. leftvalue) then
                            water          = min(wtd(row,col)*area(row),leftvalue*area(row)/2.)   
                            water1         = water/area(row)  
                            water2         = water/area(row)
                            wtd(row,col)   = wtd(row,col) - water2
                            wtd(row,col-1) = wtd(row,col-1) + water1

                            hz_read(row,col) = hz_read(row,col) - water2
                            hz_read(row,col-1) = hz_read(row,col-1) + water1

                        endif

                                    
                    endif

                end do ROWS3
                end do COLS3
                     

               if(pid.eq.1)then  !Once again send & receive the edge columns
                    call MPI_send(wtd(1,nmax+1),1,columntype,pid+1,     8,MPI_COMM_WORLD,ierr)
                    call MPI_send(wtd(1,nmax),  1,columntype,pid+1,     9,MPI_COMM_WORLD,ierr)
                    call MPI_recv(wtd(1,3),     1,columntype,numtasks-1,8,MPI_COMM_WORLD,status,ierr)
                    call MPI_recv(wtd(1,2),     1,columntype,numtasks-1,9,MPI_COMM_WORLD,status,ierr)

                    call MPI_send(hz_read(1,nmax+1),1,columntype,pid+1,     8,MPI_COMM_WORLD,ierr)
                    call MPI_send(hz_read(1,nmax),  1,columntype,pid+1,     9,MPI_COMM_WORLD,ierr)
                    call MPI_recv(hz_read(1,3),     1,columntype,numtasks-1,8,MPI_COMM_WORLD,status,ierr)
                    call MPI_recv(hz_read(1,2),     1,columntype,numtasks-1,9,MPI_COMM_WORLD,status,ierr)

                elseif(pid.eq.numtasks-1)then
                    if(mod(pid,2).eq.0)then


                    call MPI_send(wtd(1,nmax+1),1,columntype,1,    8,MPI_COMM_WORLD,ierr)
                    call MPI_send(wtd(1,nmax),  1,columntype,1,    9,MPI_COMM_WORLD,ierr)
                    call MPI_recv(wtd(1,3),     1,columntype,pid-1,8,MPI_COMM_WORLD,status,ierr)   
                    call MPI_recv(wtd(1,2),     1,columntype,pid-1,9,MPI_COMM_WORLD,status,ierr)

                    call MPI_send(hz_read(1,nmax+1),1,columntype,1,    8,MPI_COMM_WORLD,ierr)
                    call MPI_send(hz_read(1,nmax),  1,columntype,1,    9,MPI_COMM_WORLD,ierr)
                    call MPI_recv(hz_read(1,3),     1,columntype,pid-1,8,MPI_COMM_WORLD,status,ierr)   
                    call MPI_recv(hz_read(1,2),     1,columntype,pid-1,9,MPI_COMM_WORLD,status,ierr)


                    else

                    call MPI_recv(wtd(1,3),     1,columntype,pid-1,8,MPI_COMM_WORLD,status,ierr)   
                    call MPI_recv(wtd(1,2),     1,columntype,pid-1,9,MPI_COMM_WORLD,status,ierr)
                    call MPI_send(wtd(1,nmax+1),1,columntype,1,    8,MPI_COMM_WORLD,ierr)
                    call MPI_send(wtd(1,nmax),  1,columntype,1,    9,MPI_COMM_WORLD,ierr)

                    call MPI_recv(hz_read(1,3),     1,columntype,pid-1,8,MPI_COMM_WORLD,status,ierr)   
                    call MPI_recv(hz_read(1,2),     1,columntype,pid-1,9,MPI_COMM_WORLD,status,ierr)
                    call MPI_send(hz_read(1,nmax+1),1,columntype,1,    8,MPI_COMM_WORLD,ierr)
                    call MPI_send(hz_read(1,nmax),  1,columntype,1,    9,MPI_COMM_WORLD,ierr)

                    endif

                else
                    if(mod(pid,2).eq.0)then
                        call MPI_recv(wtd(1,3),     1,columntype,pid-1,8,MPI_COMM_WORLD,status,ierr)
                        call MPI_recv(wtd(1,2),     1,columntype,pid-1,9,MPI_COMM_WORLD,status,ierr)
                        call MPI_send(wtd(1,nmax+1),1,columntype,pid+1,8,MPI_COMM_WORLD,ierr)
                        call MPI_send(wtd(1,nmax),  1,columntype,pid+1,9,MPI_COMM_WORLD,ierr)

                        call MPI_recv(hz_read(1,3),     1,columntype,pid-1,8,MPI_COMM_WORLD,status,ierr)
                        call MPI_recv(hz_read(1,2),     1,columntype,pid-1,9,MPI_COMM_WORLD,status,ierr)
                        call MPI_send(hz_read(1,nmax+1),1,columntype,pid+1,8,MPI_COMM_WORLD,ierr)
                        call MPI_send(hz_read(1,nmax),  1,columntype,pid+1,9,MPI_COMM_WORLD,ierr)

                    else
                        call MPI_send(wtd(1,nmax+1),1,columntype,pid+1,8,MPI_COMM_WORLD,ierr)
                        call MPI_send(wtd(1,nmax),  1,columntype,pid+1,9,MPI_COMM_WORLD,ierr)
                        call MPI_recv(wtd(1,3),     1,columntype,pid-1,8,MPI_COMM_WORLD,status,ierr)
                        call MPI_recv(wtd(1,2),     1,columntype,pid-1,9,MPI_COMM_WORLD,status,ierr)

                        call MPI_send(hz_read(1,nmax+1),1,columntype,pid+1,8,MPI_COMM_WORLD,ierr)
                        call MPI_send(hz_read(1,nmax),  1,columntype,pid+1,9,MPI_COMM_WORLD,ierr)
                        call MPI_recv(hz_read(1,3),     1,columntype,pid-1,8,MPI_COMM_WORLD,status,ierr)
                        call MPI_recv(hz_read(1,2),     1,columntype,pid-1,9,MPI_COMM_WORLD,status,ierr)

                    endif
                endif

                diff = abs(hold_read-wtd)
                maxdiff = 0
                maxdiff = maxval(diff)
                deallocate(hz_read)
                deallocate(diff)
                deallocate(hz_1D)
                deallocate(arr)

            endif

            call MPI_ALLREDUCE(maxdiff,diff_total,1,MPI_REAL,mpi_max,MPI_COMM_WORLD,ierr)
   
        end do SURFACE




        if (pid .eq. 0) then
            print *,'done yay yay yay'
 
            do n=1,numtasks-1
                call MPI_recv(wtdglob(1,nini(n)),1,domblocksmall(n),n,4,MPI_COMM_WORLD,status,ierr)!receive the final wtd data from everyone
            end do
     
 
            open(25,file = trim(time_end)//'_27June.dat',form='unformatted',access='stream')!access='direct',recl=n2*n3)!access='stream')!do the final write - create the file
            print *,'the file has been opened'
            write(25)((wtdglob(i,j),i=1,n2),j=1,n3) !and write it to file
            print *,'the file has been written'
            close(25)
            print *,'written to file',pid

        else
            call MPI_send(wtd(1,2),1,domblocksmall(pid),0,4,MPI_COMM_WORLD,ierr) !everyone sends out the final result to 
        endif

    endif

END DO GROUNDWATER



if(pid.eq.0) then !do the final water table file write
    write(6,*) 'done; iterations = ',iter   
    wtdglob = 0.

    open(23,file = trim(time_end)//'_water_test.dat',form='unformatted',access='stream')!'direct',recl=n2*n3) !do the final write - create the file
    do n=1,numtasks-1
        call MPI_recv(wtdglob(1,nini(n)),1,domblocksmall(n),n,8,MPI_COMM_WORLD,status,ierr) !receive the final wtd data from everyone
    end do

    write(23,rec=1)((wtdglob(i,j),i=1,n2),j=1,n3) !and write it to file
    close(23)
    print *,'wrote the groundwater file'
else

    call MPI_send(wtd(1,2),1,domblocksmall(pid),0,8,MPI_COMM_WORLD,ierr) !everyone sends out the final result to pid 0! 
endif



print *,'about to try deallocating',pid


 deallocate(topo_start_read,stat=error)
  if (error.ne.0) then
      print *, 'topo_start_read error'
  endif
 deallocate(topo_end_read,stat=error)
  if (error.ne.0) then
      print *, 'topo_end_read error'
  endif
  deallocate(rech_start_read,stat=error)
  if (error.ne.0) then
      print *, 'rech_start_read error'
  endif
 deallocate(rech_end_read,stat=error)
  if (error.ne.0) then
      print *, 'rech_end_read error'
  endif
 deallocate(fslope_start_read,stat=error)
  if (error.ne.0) then
      print *, 'fslope_start_read error'
  endif
 deallocate(fslope_end_read,stat=error)
  if (error.ne.0) then
      print *, 'fslope_end_read error'
  endif
 deallocate(temp_start_read,stat=error)
  if (error.ne.0) then
      print *, 'temp_start_read error'
  endif 
 deallocate(temp_end_read,stat=error)
  if (error.ne.0) then
      print *, 'temp_end_read error'
  endif
 deallocate(topo_sent,stat=error)
  if (error.ne.0) then
      print *, 'topo_sent error'
  endif
 deallocate(rech_sent,stat=error)
  if (error.ne.0) then
      print *, 'rech_sent error'
  endif
 deallocate(fdepth_sent,stat=error)
  if (error.ne.0) then
      print *, 'fdepth_sent error'
  endif
 deallocate(qlat_north,stat=error)
  if (error.ne.0) then
      print *, 'qlat_north error'
  endif
 deallocate(qlat_south,stat=error)
  if (error.ne.0) then
      print *, 'qlat_south error'
  endif
 deallocate(qlat_east,stat=error)
  if (error.ne.0) then
      print *, 'qlat_east error'
  endif
 deallocate(qlat_west,stat=error)
  if (error.ne.0) then
      print *, 'qlat_west error'
  endif
 

 ! deallocate(landmask,stat=error)
 ! if (error.ne.0)then
 !     print *,'landmask error'
 ! endif
!  deallocate(wtd,stat=error)
!  if (error.ne.0)then
!      print *,'wtd error'
!  endif
 ! deallocate(wtdnew,stat=error)
 ! if (error.ne.0) then
 !     print *, 'wtdnew error'
 ! endif
 
  deallocate(fdepth,stat=error)
  if (error.ne.0) then
      print *,'fdepth error'
  endif
  deallocate(ksat,stat=error)
  if (error.ne.0) then
      print *,'ksat error'
  endif
 
  deallocate(kcell,stat=error)
  if (error.ne.0) then
      print *,'kcell error'
  endif
  deallocate(head,stat=error)
  if (error.ne.0) then
      print *,'head error'
  endif
   deallocate(nini,stat=error)
  if (error.ne.0) then
      print *,'nini error'
  endif
      deallocate(nend,stat=error)
  if (error.ne.0) then
      print *,'nend error'
  endif
  deallocate(domblock,stat=error)
  if (error.ne.0) then
      print *,'domblock error'
  endif
 
  deallocate(domblocksmall,stat=error)
  if (error.ne.0) then
      print *,'domblocksmall error'
  end if
  deallocate(xlat,stat=error)
  if (error.ne.0) then
      print *,'xlat error'
  endif
  deallocate(alpha,stat=error)
  if (error.ne.0) then
      print *, 'alpha error'
  endif


  deallocate(topo,stat=error)
  if (error.ne.0) then
      print *, 'topo error'
  endif

  deallocate(mask,stat=error)
  if (error.ne.0) then
      print *,'mask error'
  endif

  deallocate(nini,stat=error)
  if (error.ne.0) then
      print *,'nini error'
  endif

  deallocate(nend,stat=error)
  if (error.ne.0) then
      print *,'nend error'
  endif

  deallocate(domblock,stat=error)
  if (error.ne.0) then
      print *,'domblock error'
  endif


  deallocate(domblocksmall,stat=error)
  if (error.ne.0) then
      print *,'domblocksmall error'
  end if




  print *,'done',pid


  call MPI_FINALIZE(ierr)

  print *,'finished',pid

end program GWSW
