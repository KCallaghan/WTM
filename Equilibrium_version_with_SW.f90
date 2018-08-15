!this version tries to obtain equilibrium for both groundwater and surface water
!this should be run at LGM to obtain the starting water tables for future transient runs
!one model run per setup - ie one for TRACE, one for HADCM3
!We will try running the surface water on a coarser grid then switch to finer, to help with processing times




!*************************************************************************************************************

subroutine Merge(A,NA,B,NB,C,NC) !used to partially sort the surface water array to speed processing
 
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


subroutine dividedomain(n2,n3,numtasks,nini,filemask,ntotal) !divide among processors
  use netcdf
  implicit none

 
  integer :: n2,n3,numtasks
  integer :: nini(1:numtasks-1)
  real,allocatable,dimension(:,:) :: varread
  integer,allocatable,dimension(:) :: ncells
  integer :: iret,ncid,varid,ntotal,ncount,n,j
  character*100 :: filemask

  allocate(varread(n2,n3))

  write(6,*)'reading in the mask to divide the domain'

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

  ntotal = count(varread>0.5) !count the number of land cells. I changed this slightly since I am using mask here rather than topo; all cells with a value of 1 should be included. 


  allocate(ncells(n3))
 
  ncells = count(varread>0.5,1) !The number of cells which are defined in the 1 dimension

  ncount=0

  nini(1) = 2

  n=2 !counter
  
  do j=1,n3!,5
      ncount=ncount+ncells(j) !add the number of cells defined in the current column
      if (ncount .ge. ntotal/(numtasks-1)) then !>= total number of defined cells/number of threads
          nini(n) = j-1 !Telling it which row it needs to start in, to divide up the work equally, so each task handles the same number of cells.
          ncount = ncells(j) !reset the ncount for the next task
          n = n+1
      endif
      if(n .eq. numtasks) exit
  end do

  deallocate(varread,ncells)

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


program equilibrium

  use netcdf
  use mpi

  implicit none !ensures you don't create any variables without explicitly defining their type

!define all variables:

  REAL,PARAMETER :: UNDEF = -1.0E+7
  REAL(KIND=8),PARAMETER :: pi=3.141592653589793D0+0

  integer :: n2,n3,ntotal,iret,ncid,varid,n,j,nmax,i,numbertotal,numbercount,iter,error,numberold
  integer ierr,pid,numtasks,rc,tasktype,status(MPI_STATUS_SIZE),columntype,columntypeint,&
fine,n_use,n2_use

  real :: dltxy = 120 !there are 120 30 arc-second pieces in one degree
  real thres,d0,d1,d2,d3,deltat,dx,dy,xn,xs,wtdmax,q,total

  character*60 :: surfdatadir,initdatadir
  character*8 :: region
  character*100 :: filetopo,filewtd,filefdepth,fileksat,filerech,filemask,filetemp

  REAL(KIND=8) :: SEDGE !kind defines precision, number of bytes

  integer,allocatable :: domblock(:),domblocksmall(:),nini(:),nend(:)
  integer*1,allocatable,dimension(:) :: maskline
  integer*1, allocatable,dimension(:,:) :: landmask,maskold,watermask

  real, allocatable :: alpha(:),xlat(:),alphamonth(:)
  real, allocatable, dimension(:,:) :: topo,wtd,kcell,ksat,head,qlat,rechmean, fdepth,&
varread, wtd_read, fdepth_read, ksat_read,rech_read,rech_month_read, &
rech_month, wtdglob,fdepth_start, temp_read,fslope,temp_sent, &
h,hold,diff,hz_read,hold_read,mask_read,mask,topo_coarse,wtd_use, &
diff_coarse

  integer :: col,row,converged,counter

  real :: upvalue,downvalue,leftvalue,rightvalue,water,water1,water2

  real diff_total,maxdiff,placeholder


  real,allocatable,dimension(:) :: hz_1D,area,area_coarse

  integer :: Merge_size
  integer,allocatable,dimension(:)  :: T

  real,allocatable,dimension(:,:) :: arr



!*******************************************************************************************

 
  call MPI_INIT(ierr)       !Initialise MPI
  if (ierr .ne. MPI_SUCCESS) then       !Error catching - this part should never run.
      print *,'Error starting MPI program. Terminating.'
      call MPI_ABORT(MPI_COMM_WORLD,rc,ierr)
  end if



!*******************************************************************************************

!Set global parameters:

  n2          = 2000 !number of cells in x and y directions
  n3          = 1000
  region      = 'Mad_' !folder name for the 4 files. Limit 8 characters.
  surfdatadir = 'surfdata/' !folder names for data, limit 60 characters. Topo, fdepth, and rech go here.
  initdatadir = 'initdata/' !mask and wtd go here
!ksat is not within a folder since the same one is being used for all times.


  Thres = 0.001 !When to stop calculating - can try with different values and see how it goes. This is the coarsest option.
  SEDGE = -27 !southern most latitude

!input files:

  filetopo   = trim(surfdatadir)//trim(region)//'021000_topo_rotated.nc'
  filefdepth = trim(surfdatadir)//trim(region)//'021000_fslope_rotated.nc'
  fileksat   = trim(surfdatadir)//'Mad_ksat_rotated.nc'
  filerech   = trim(surfdatadir)//trim(region)//'021000_rech_rotated.nc'
  filemask   = trim(surfdatadir)//trim(region)//'020500_mask_rotated.nc'
  filewtd    = trim(surfdatadir)//trim(region)//'021000_wtd_rotated.nc'
  filetemp   = trim(surfdatadir)//trim(region)//'021000_temp_rotated.nc'

!print *,filetopo,filefdepth,fileksat,filerech,filemask,filewtd,filetemp
!output text file

open (15,file=trim(region)//'_output_equilibrium_SW.txt') !creating a text file to store the values of iterations etc

!different increments to use to move water table depth:

  d0 = 0.01 !use /12 for monthly iterations
  d1 = 0.05 !no need to use /12 any more, we now switch to monthly within the loop. 
  d2 = 0.1
  d3 = 0.25

  IF (Thres .eq. 0.01) THEN
      WRITE(6,*) 'Doing the coarse adjustment based on the 10 mm threshold'
      WRITE(15,*) 'Doing the coarse adjustment based on the 10 mm threshold'
      d0=d0
  ELSEIF (Thres .eq. 0.005) THEN
      WRITE(6,*) 'Medium adjustment for the 5 mm threshold'
      WRITE(15,*) 'Medium adjustment for the 5 mm threshold'
      d0=d0/5.
  ELSEIF (Thres .eq. 0.001) THEN
      WRITE(6,*) 'Fine adjustment for the 1 mm threshold'
      WRITE(15,*) 'Medium adjustment for the 5 mm threshold'
      d0=d0/20.
  ELSE
      WRITE(6,*) 'The threshold input does not exist, the program will now quit.'
      WRITE(15,*) 'The threshold input does not exist, the program will now quit.'
      STOP
  ENDIF

 
!get the number of tasks available
  call MPI_COMM_RANK(MPI_COMM_WORLD,pid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,numtasks,ierr)
  print *,'Number of tasks=',numtasks,'My rank=',pid
 
!allocate array space, the size of the number of tasks available

  allocate(nini(1:numtasks-1))
  allocate(nend(1:numtasks-1))
  allocate(domblock(1:numtasks-1))
  allocate(domblocksmall(1:numtasks-1))
  allocate(maskline(n2))               

  call MPI_TYPE_CONTIGUOUS(numtasks-1,MPI_INTEGER,tasktype,ierr) !creates a contiguous datatype. The new datatype is called tasktype
  call MPI_Type_commit(tasktype,ierr)


!divide the domain among tasks:

  if(pid .eq. 0) then
      write(6,*) 'PID = 0'
      call dividedomain(n2,n3,numtasks,nini,filemask,ntotal) !Divides up the work equally among all of the ranks, by number of defined land cells.
      write (6,*) 'Dividedomain done'
      write (15,*) 'Dividedomain done' 

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
print *,nend(pid)

  endif    


 
 do n=1,numtasks-1
      nmax = nend(n) - nini(n) + 4 !max number of columns we have in this group
print *,'here',nmax,nini(n),nend(n),pid
      call MPI_TYPE_CONTIGUOUS(n2*nmax,MPI_REAL,domblock(n),ierr)
      call MPI_type_commit(domblock(n),ierr)

      nmax = nmax-3
      call MPI_TYPE_CONTIGUOUS(n2*nmax,MPI_REAL,domblocksmall(n),ierr)
      call MPI_type_commit(domblocksmall(n),ierr)
  end do


  call MPI_TYPE_CONTIGUOUS(n2,MPI_REAL,columntype,ierr)
  call MPI_type_commit(columntype,ierr)


  call MPI_TYPE_CONTIGUOUS(n2,MPI_INTEGER1,columntypeint,ierr)
  call MPI_type_commit(columntypeint,ierr)


  deltat = (365.*24.*3600.) !Seconds in an annual timestep

  dy = 6370000.*pi/(180.*dltxy) !radius of the earth * pi / number of possible cells in the y-direction. This should equal the height of each cell in the N-S direction.
  dx=dy
  
  if(pid .gt. 0) then
    !  nmax = nend(pid) - nini(pid) +1
      allocate(xlat(n2))
      allocate(area(n2))
      allocate(alpha(n2))
      allocate(alphamonth(n2))


      do j=1,n2  !changing area of cell depending on its latitude. Not totally sure what each of the different variables here represents...
          xlat(j)     = (float(j)/dltxy+SEDGE)*pi/180.
          xs          = (float(2*(j)-1)/(dltxy*2.)+SEDGE)*pi/180. !Latitude number 1
          xn          = (float(2*(j)+1)/(dltxy*2.)+SEDGE)*pi/180. !latitude number 2
          placeholder = dy*6370000.*(sin(xn)-sin(xs))/2.
          area        = dy*6370000.*(sin(xn)-sin(xs))/2. !final cell area for that latitude

          alpha(j) = deltat/placeholder 
          alphamonth(j) = (deltat/12)/placeholder 
      end do
  end if

!############################################################################# thing to change with lakes
  wtdmax = 0 !water table depth is 0. This should probably be changed when I start bringing in the lakes. 

  if(pid .eq. 0) then



    allocate(wtdglob(n2,n3))
    wtdglob = 0.




      allocate(varread(n2,n3))
     
      iret = nf90_open(filetopo,0,ncid) !reading in the topo
      call check_err(iret, "Unable to load file '" // filemask // "'!")

      iret = nf90_inq_varid(ncid,'value',varid)
      call check_err(iret, "")

      iret = nf90_get_var(ncid,varid,varread)
      call check_err(iret, "")

      iret = nf90_close(ncid)
      call check_err(iret, "Unable to close file!")


      where(varread .le. UNDEF) varread = 0. !change undefined cells to 0

 
      allocate(watermask(n2,n3)) 
     
      iret = nf90_open(filemask,0,ncid) !reading in the mask
      call check_err(iret, "Unable to load file '" // filemask // "'!")

      iret = nf90_inq_varid(ncid,'value',varid)
      call check_err(iret, "")

      iret = nf90_get_var(ncid,varid,watermask)
      call check_err(iret, "")

      iret = nf90_close(ncid)
      call check_err(iret, "Unable to close file!")
      


      allocate(wtd_read(n2,n3))

      wtd_read = 0           !no longer using an input water table! 

     
    
      allocate(ksat_read(n2,n3))
   
     
      iret = nf90_open(fileksat,0,ncid) !reading in the ksat
      call check_err(iret, "Unable to load file '" // filemask // "'!")

      iret = nf90_inq_varid(ncid,'value',varid)
      call check_err(iret,"")

      iret = nf90_get_var(ncid,varid,ksat_read)
      call check_err(iret,"")

      iret = nf90_close(ncid)
      call check_err(iret, "Unable to close file!")



      allocate(rech_read(n2,n3))
      allocate(rech_month_read(n2,n3))
     
      iret = nf90_open(filerech,0,ncid) !reading in the recharge file
      call check_err(iret, "Unable to load file '" // filemask // "'!")

      iret = nf90_inq_varid(ncid,'value',varid)
      call check_err(iret, "")

      iret = nf90_get_var(ncid,varid,rech_read)
      call check_err(iret, "")

      iret = nf90_close(ncid)
      call check_err(iret, "Unable to close file!")

!####################################################################################### potentially change with lakes
      rech_read = max(rech_read,0.) !setting it so recharge can only be positive

    rech_month_read = (rech_read/12.)  



      allocate(fdepth_read(n2,n3))
    
      iret = nf90_open(filefdepth,0,ncid) !reading in the fdepth
      call check_err(iret, "Unable to load file '" // filemask // "'!")

      iret = nf90_inq_varid(ncid,'value',varid)
      call check_err(iret, "")

      iret = nf90_get_var(ncid,varid,fdepth_read)
      call check_err(iret, "")

      iret = nf90_close(ncid)
      call check_err(iret, "")



allocate(temp_read(n2,n3))
    
      iret = nf90_open(filetemp,0,ncid) !reading in the fdepth
      call check_err(iret, "Unable to load file '" // filemask // "'!")

      iret = nf90_inq_varid(ncid,'value',varid)
      call check_err(iret, "")

      iret = nf90_get_var(ncid,varid,temp_read)
      call check_err(iret, "")

      iret = nf90_close(ncid)
      call check_err(iret, "")


!now send everything we have opened:

print *,'line 498'

      do n=1,numtasks-2
          call MPI_send(varread        (1,nini(n)-1),1,domblock(n),n,1,MPI_COMM_WORLD,ierr)
print *,'line 502'
          call MPI_send(watermask      (1,nini(n)-1),1,domblock(n),n,2,MPI_COMM_WORLD,ierr)
          call MPI_send(wtd_read       (1,nini(n)-1),1,domblock(n),n,3,MPI_COMM_WORLD,ierr)
          call MPI_send(ksat_read      (1,nini(n)-1),1,domblock(n),n,5,MPI_COMM_WORLD,ierr)
          call MPI_send(rech_read      (1,nini(n)-1),1,domblock(n),n,6,MPI_COMM_WORLD,ierr)
          call MPI_send(rech_month_read(1,nini(n)-1),1,domblock(n),n,7,MPI_COMM_WORLD,ierr)
  print *,'line 507'      
      end do


    call MPI_send(varread        (1,nini(numtasks-1)-1),1,domblocksmall(numtasks-1),numtasks-1,1,MPI_COMM_WORLD,ierr)
    call MPI_send(watermask      (1,nini(numtasks-1)-1),1,domblocksmall(numtasks-1),numtasks-1,2,MPI_COMM_WORLD,ierr)
    call MPI_send(wtd_read       (1,nini(numtasks-1)-1),1,domblocksmall(numtasks-1),numtasks-1,3,MPI_COMM_WORLD,ierr)
    call MPI_send(ksat_read      (1,nini(numtasks-1)-1),1,domblocksmall(numtasks-1),numtasks-1,5,MPI_COMM_WORLD,ierr)
    call MPI_send(rech_read      (1,nini(numtasks-1)-1),1,domblocksmall(numtasks-1),numtasks-1,6,MPI_COMM_WORLD,ierr)
    call MPI_send(rech_month_read(1,nini(numtasks-1)-1),1,domblocksmall(numtasks-1),numtasks-1,7,MPI_COMM_WORLD,ierr)
   

    call MPI_send(varread        (1,1),1,columntype,numtasks-1,1,MPI_COMM_WORLD,ierr)
    call MPI_send(watermask      (1,1),1,columntype,numtasks-1,2,MPI_COMM_WORLD,ierr)
    call MPI_send(wtd_read       (1,1),1,columntype,numtasks-1,3,MPI_COMM_WORLD,ierr)
    call MPI_send(ksat_read      (1,1),1,columntype,numtasks-1,5,MPI_COMM_WORLD,ierr)
    call MPI_send(rech_read      (1,1),1,columntype,numtasks-1,6,MPI_COMM_WORLD,ierr)
    call MPI_send(rech_month_read(1,1),1,columntype,numtasks-1,7,MPI_COMM_WORLD,ierr)
   

    call MPI_send(varread        (1,2),1,columntype,numtasks-1,1,MPI_COMM_WORLD,ierr)
    call MPI_send(watermask      (1,2),1,columntype,numtasks-1,2,MPI_COMM_WORLD,ierr)
    call MPI_send(wtd_read       (1,2),1,columntype,numtasks-1,3,MPI_COMM_WORLD,ierr)
    call MPI_send(ksat_read      (1,2),1,columntype,numtasks-1,5,MPI_COMM_WORLD,ierr)
    call MPI_send(rech_read      (1,2),1,columntype,numtasks-1,6,MPI_COMM_WORLD,ierr)
    call MPI_send(rech_month_read(1,2),1,columntype,numtasks-1,7,MPI_COMM_WORLD,ierr)
   

    call MPI_send(varread        (1,3),1,columntype,numtasks-1,1,MPI_COMM_WORLD,ierr)
    call MPI_send(watermask      (1,3),1,columntype,numtasks-1,2,MPI_COMM_WORLD,ierr)
    call MPI_send(wtd_read       (1,3),1,columntype,numtasks-1,3,MPI_COMM_WORLD,ierr)
    call MPI_send(ksat_read      (1,3),1,columntype,numtasks-1,5,MPI_COMM_WORLD,ierr)
    call MPI_send(rech_read      (1,3),1,columntype,numtasks-1,6,MPI_COMM_WORLD,ierr)
    call MPI_send(rech_month_read(1,3),1,columntype,numtasks-1,7,MPI_COMM_WORLD,ierr)
  

      deallocate(varread)
      deallocate(ksat_read)
      deallocate(wtd_read)
      deallocate(rech_read)
      deallocate(watermask)
      deallocate(rech_month_read)


  else
  
      nmax = nend(pid) - nini(pid) +4
      allocate(topo      (n2,nmax))
      allocate(landmask  (n2,nmax))
      allocate(wtd       (n2,nmax))
      allocate(ksat      (n2,nmax))
      allocate(rechmean  (n2,nmax))
      allocate(rech_month(n2,nmax))
      allocate(kcell     (n2,nmax))
      allocate(head      (n2,nmax))
      allocate(qlat      (n2,nmax))
      allocate(mask      (n2,nmax),maskold(n2,nmax))
      write(6,*) 'allocated all'

    if(pid.lt.numtasks-1)then

      call MPI_recv(topo      (1,1),1,domblock(pid),0,1,MPI_COMM_WORLD,status,ierr) !receiving everthing that was sent above     
      call MPI_recv(landmask  (1,1),1,domblock(pid),0,2,MPI_COMM_WORLD,status,ierr)
      call MPI_recv(wtd       (1,1),1,domblock(pid),0,3,MPI_COMM_WORLD,status,ierr)
      call MPI_recv(ksat      (1,1),1,domblock(pid),0,5,MPI_COMM_WORLD,status,ierr)
      call MPI_recv(rechmean  (1,1),1,domblock(pid),0,6,MPI_COMM_WORLD,status,ierr)
      call MPI_recv(rech_month(1,1),1,domblock(pid),0,7,MPI_COMM_WORLD,status,ierr)

   else

      call MPI_recv(topo      (1,1), 1, domblocksmall(pid),0,1,MPI_COMM_WORLD,status,ierr)
      call MPI_recv(landmask  (1,1), 1, domblocksmall(pid),0,2,MPI_COMM_WORLD,status,ierr)
      call MPI_recv(wtd       (1,1), 1, domblocksmall(pid),0,3,MPI_COMM_WORLD,status,ierr)
      call MPI_recv(ksat      (1,1), 1, domblocksmall(pid),0,5,MPI_COMM_WORLD,status,ierr)
      call MPI_recv(rechmean  (1,1), 1, domblocksmall(pid),0,6,MPI_COMM_WORLD,status,ierr)
      call MPI_recv(rech_month(1,1), 1, domblocksmall(pid),0,7,MPI_COMM_WORLD,status,ierr)
    
      call MPI_recv(topo      (1,nmax-2),1,columntype,0,1,MPI_COMM_WORLD,status,ierr)
      call MPI_recv(landmask  (1,nmax-2),1,columntype,0,2,MPI_COMM_WORLD,status,ierr)
      call MPI_recv(wtd       (1,nmax-2),1,columntype,0,3,MPI_COMM_WORLD,status,ierr)
      call MPI_recv(ksat      (1,nmax-2),1,columntype,0,5,MPI_COMM_WORLD,status,ierr)
      call MPI_recv(rechmean  (1,nmax-2),1,columntype,0,6,MPI_COMM_WORLD,status,ierr)
      call MPI_recv(rech_month(1,nmax-2),1,columntype,0,7,MPI_COMM_WORLD,status,ierr)
  
      call MPI_recv(topo      (1,nmax-1),1,columntype,0,1,MPI_COMM_WORLD,status,ierr)
      call MPI_recv(landmask  (1,nmax-1),1,columntype,0,2,MPI_COMM_WORLD,status,ierr)
      call MPI_recv(wtd       (1,nmax-1),1,columntype,0,3,MPI_COMM_WORLD,status,ierr)
      call MPI_recv(ksat      (1,nmax-1),1,columntype,0,5,MPI_COMM_WORLD,status,ierr)
      call MPI_recv(rechmean  (1,nmax-1),1,columntype,0,6,MPI_COMM_WORLD,status,ierr)
      call MPI_recv(rech_month(1,nmax-1),1,columntype,0,7,MPI_COMM_WORLD,status,ierr)

      call MPI_recv(topo      (1,nmax),1,columntype,0,1,MPI_COMM_WORLD,status,ierr)
      call MPI_recv(landmask  (1,nmax),1,columntype,0,2,MPI_COMM_WORLD,status,ierr)
      call MPI_recv(wtd       (1,nmax),1,columntype,0,3,MPI_COMM_WORLD,status,ierr)
      call MPI_recv(ksat      (1,nmax),1,columntype,0,5,MPI_COMM_WORLD,status,ierr)
      call MPI_recv(rechmean  (1,nmax),1,columntype,0,6,MPI_COMM_WORLD,status,ierr)
      call MPI_recv(rech_month(1,nmax),1,columntype,0,7,MPI_COMM_WORLD,status,ierr)
   
    endif
   
      write(6, *) 'variables received'
      write(15,*) 'variables received'

      maskold = 1 !This is just being done ahead of the big loop, maskold is used to show which cells still need to be processed. 

  endif



if (pid.eq.0) then


    allocate(fdepth_start(n2,n3))


     write(6,*)'fdepth_Start allocated'
 
          do j=1,n3
              do i=1,n2
                  if (temp_read(i,j) .gt. -5) then
                      fdepth_start(i,j) = fdepth_read(i,j)
                  elseif (temp_read(i,j) .lt. -14) then
                      fdepth_start(i,j) = fdepth_read(i,j) * (0.17+0.005*temp_read(i,j))
                  else
                      fdepth_start(i,j) = fdepth_read(i,j) * (1.5 + 0.1*temp_read(i,j))
                  endif
              end do
          end do
write(6,*) 'fdepth_start created'
  

where(fdepth_start .le. 0.0001) fdepth_start = 0.0001 !change undefined cells to 0



          do n=1,numtasks-2
write (6,*) 'trying to send'
             
              call MPI_send(fdepth_start(1,nini(n)-1),1,domblock(n),n,15,MPI_COMM_WORLD,ierr)
             
          end do 


   call MPI_send(fdepth_start(1,nini(numtasks-1)-1),1,domblocksmall(numtasks-1),numtasks-1,15,MPI_COMM_WORLD,ierr)
   call MPI_send(fdepth_start(1,1),1,columntype,numtasks-1,15,MPI_COMM_WORLD,ierr)
   call MPI_send(fdepth_start(1,2),1,columntype,numtasks-1,15,MPI_COMM_WORLD,ierr)
   call MPI_send(fdepth_start(1,3),1,columntype,numtasks-1,15,MPI_COMM_WORLD,ierr)

          write(6,*)'sent initial'

          deallocate(fdepth_start)

      elseif(pid.lt.numtasks-1)then
  

       allocate(fdepth(n2,nmax))
       
          call MPI_recv(fdepth(1,1),1,domblock(pid),0,15,MPI_COMM_WORLD,status,ierr)
 
          write(6,*)'received initial'

      else

    allocate(fdepth(n2,nmax))

   call MPI_recv(fdepth(1,1),1,domblocksmall(pid),0,15,MPI_COMM_WORLD,status,ierr)
   call MPI_recv(fdepth(1,nmax-2),1,columntype,0,15,MPI_COMM_WORLD,status,ierr)
   call MPI_recv(fdepth(1,nmax-1),1,columntype,0,15,MPI_COMM_WORLD,status,ierr)
   call MPI_recv(fdepth(1,nmax  ),1,columntype,0,15,MPI_COMM_WORLD,status,ierr)


       endif



!allocate the arrays for the mergesort done later



  iter = 0              !count the number of iterations - this is used to know when to switch from annual to monthly cycles. 
  numbercount = 0
  numbertotal = ntotal !we will use this to get numbertotal less than x % land cells as equilibrium condition
 

 
  GROUNDWATER : DO while(numbertotal > ntotal/200 .and. iter<10 ) !start of the main loop with condition - either number of iterations or % equilibrium. I am going for 98% equilibrium. 
      
 

      if (iter .eq. 30000) then           !I have switched to 30000 iterations before switching to monthly processing, because with 50000 it always seemed to be doing nothing for a long time. 
!here we automatically switch to monthly processing, 
          write(6,*) '30000 iterations,adjusting the values for monthly processing'
          write(15,*) '30000 iterations,adjusting the values for monthly processing'
          thres       = thres/12.
          d0          = d0/12.
          d1          = d1/12.
          d2          = d2/12.
          d3          = d3/12.
          deltat      = deltat/12.
          alpha       = alphamonth
          rechmean    = rech_month
          maskold     = 1          !because we changed the threshold, so we want to recheck all the cells. 
          numbertotal = ntotal
      endif

      iter = iter + 1 !count the iterations

      if (pid .eq. 0) then !This is not completely necessary, but it's nice to see how far the run is. 
          write (6,*) 'PID = 0. Numbertotal:',numbertotal,'ntotal',ntotal
          write (6,*) 'iter number',iter
          write (15,*) 'PID = 0. Numbertotal:',numbertotal,'ntotal',ntotal
          write (15,*) 'iter number',iter
      endif

      numbercount = 0


      if(mod(float(iter),10000.).eq.0.)then       !save the data
          if(pid.eq.0) then
              allocate(wtdglob(n2,n3))
              open(23,file=trim(region)//'equilibrium_SW.dat',form='unformatted',access='direct',recl=n2*n3)

              do n=1,numtasks-1
                  write(6,*)'receiving wtd',n
                  write(15,*)'receiving wtd',n
                  call MPI_recv(wtdglob(1,nini(n)),1,domblocksmall(n),n,8,MPI_COMM_WORLD,status,ierr)

              end do
              
              write(23,rec=1)((wtdglob(i,j),i=1,n2),j=1,n3)

              deallocate(wtdglob)
          else
              call MPI_send(wtd(1,2),1,domblocksmall(pid),0,8,MPI_COMM_WORLD,ierr)
          endif
      endif




!######################################################################################################## change for lakes
      if (pid .gt. 0) wtd = min(wtd,0.) !any water table above ground gets reset to 0 - this needs to be changed for lakes


      IF (pid .gt. 0) then
          mask = 0
          nmax = nend(pid) - nini(pid)
          do j=1,nmax+1
              do i=1,n2
                  if(fdepth(i,j) .gt. 0. .and. (maskold(i,j).gt.0.or.j.eq.1.or.j.eq.nmax+1)) then !maskold is used to stop us from repeating all of the steps on cells which have already reached equilibrium. 
                      head(i,j) = topo(i,j) + wtd(i,j) !gives the water table height (same as land surface if a wtd isn't loaded in), which = head
                      if(wtd(i,j) .lt. -1.5) then
!work out hydraulic conductivity for each cell
                          kcell(i,j) = fdepth(i,j) *ksat(i,j)*exp((wtd(i,j)+1.5)/fdepth(i,j)) !This is equation S6 from the paper
                      elseif(wtd(i,j).le.0)then
                          kcell(i,j) = ksat(i,j)*(wtd(i,j)+1.5+fdepth(i,j)) !equation S4 from the paper 
                    else
                        kcell(i,j) = ksat(i,j)*(0+1.5+fdepth(i,j)) !maxes out when water is at the surface, to avoid instabilities in surface water movement.
                      endif
                  endif
              end do
          end do


!empty spots in a matrix are zero. 


!head for ocean cells is empty ie zero, which is correct. I wonder if it is correct for kcell also to be 0?


          do j=2,nmax
              do i=2,n2-1

                  if(landmask(i,j) .gt. 0 .and. maskold(i,j).gt.0) then
                      q=0.
                      !north
                      q  = q + (kcell(i,j+1)+kcell(i,j))*(head(i,j+1)-head(i,j)) * cos(xlat(j)+pi/(180.*dltxy*2.))   !soo... we're just adding to the total q each time? we're getting a total discharge but not actually moving it in these directions?
!it seems like we are getting the total which will be discharged from each cell
                      !south
                      q  = q + (kcell(i,j-1)+kcell(i,j))*(head(i,j-1)-head(i,j)) * cos(xlat(j)-pi/(180.*dltxy*2.))
                      !west
                      q  = q + (kcell(i-1,j)+kcell(i,j))*(head(i-1,j)-head(i,j)) / cos(xlat(j))
                      !east
                      q  = q + (kcell(i+1,j)+kcell(i,j))*(head(i+1,j)-head(i,j)) / cos(xlat(j))
                      qlat(i,j)=alpha(j)*q  !and we multiply it with alpha, which somehow brings in the timestep?
!I think multiplying it with alpha gets the total that will be discharged as it builds up over that whole time. 

!mm -> m
                      total=rechmean(i,j)*1.e-3 + qlat(i,j)  
                     ! total = rechmean(i,j) + qlat(i,j)
                      numberold = numbercount

                    

              !       As recharge is fixed, the following applies:
                     !       (a) if total <0, meaning too much lateral flows i.e., water table is too high.
                     !       (b) if total >0, meaning too little lateral flow, i.e., water table is too low.

                      if(total .lt. -1.) then   !adjustment size depending on how far off it is
                          wtd(i,j) = wtd(i,j) -d3 !we use d0-d3 as different size increments of adjustment
                          numbercount = numbercount + 1 !and count the cell as not yet being in equilibrium. 
                      elseif (total .lt. -0.25) then
                          wtd(i,j) = wtd(i,j) -d2
                          numbercount = numbercount + 1
                      elseif (total .lt. -0.1) then
                          wtd(i,j) = wtd(i,j) -d1 
                          numbercount = numbercount + 1
                          
                      elseif (total .lt. -thres) then
                          wtd(i,j) = wtd(i,j) -d0
                          numbercount = numbercount + 1

!##############################################################################################################change for lakes
                      elseif(total .gt. 1. .and.wtd(i,j).lt.wtdmax) then                   !for now, we have to prevent more from pooling on the surface or we never get equilibrium. BUT this clearly has to change with lakes!
                          wtd(i,j) = wtd(i,j) +d3
                          numbercount = numbercount + 1
                      elseif (total .gt. 0.25 .and.wtd(i,j).lt.wtdmax) then
                          wtd(i,j) = wtd(i,j) +d2
                          numbercount = numbercount + 1
                      elseif (total .gt. 0.1 .and.wtd(i,j).lt.wtdmax) then
                          wtd(i,j) = wtd(i,j) +d1 
                          numbercount = numbercount + 1
                      elseif (total .gt. thres .and.wtd(i,j).lt.wtdmax) then
                          wtd(i,j) = wtd(i,j) +d0
                          numbercount = numbercount + 1
                      !things which are between -thres and thres do not change; they are in equilibrium.
                      endif


                     if (numberold .ne. numbercount) then
                         mask(i+1,j) = 1 !tag cells to show they are still not in equilibrium.
                         mask(i-1,j) = 1
                         mask(i,j+1) = 1
                         mask(i,j-1) = 1
                         mask(i,j)   = 1
                     endif
 

                  endif
              end do
          end do


      
         maskold = mask !to record which cells still need to be processed.

         numbertotal = numbercount !total number of cells still out of equilibrium this iteration
                             


!sending and receiving the lines on either side of each section to allow for flow across those lines. 

          if(pid .eq. 1) then
              call MPI_send(wtd(1,nmax-1),1,columntype,2,0,MPI_COMM_WORLD,ierr)
              call MPI_send(maskold(1,nmax),1,columntypeint,2,1,MPI_COMM_WORLD,ierr) !so we send out the border part of our wtd and our mask

              call MPI_recv(wtd(1,nmax),1,columntype,2,0,MPI_COMM_WORLD,status,ierr)
              call MPI_recv(maskline(1),1,columntypeint,2,1,MPI_COMM_WORLD,status,ierr) !and we receive it from the neighbour

              do i=1,n2
                  maskold(i,nend(1)-1) = maskold(i,nend(1)-1)+maskline(i) !and then we add that end line of the mask to what we already have, so that any cells which were out of equilibrium on EITHER side of the line get included in the next round of processing. 
              end do


          elseif (pid .eq. numtasks-1) then     !and, continue to do this for all of the different tasks. We separate odd and even to allow our sends and receives to work without creating a blockage. 
              if (mod(pid,2) .eq.numtasks - 1) then
                  call MPI_recv(wtd(1,1),1,columntype,pid-1,0,MPI_COMM_WORLD,status,ierr)
                  call MPI_recv(maskline(1),1,columntypeint,pid-1,1,MPI_COMM_WORLD,status,ierr)
              
                  call MPI_send(wtd(1,2),1,columntype,pid-1,0,MPI_COMM_WORLD,ierr)
                  call MPI_send(maskold(1,1),1,columntypeint,pid-1,1,MPI_COMM_WORLD,ierr)

                  do i=1,n2
                      maskold(i,2)=maskold(i,2)+maskline(i)
                  end do



              else
                  call MPI_send(wtd(1,2),1,columntype,pid-1,0,MPI_COMM_WORLD,ierr)
                  call MPI_send(maskold(1,1),1,columntypeint,pid-1,1,MPI_COMM_WORLD,ierr)

                  call MPI_recv(wtd(1,1),1,columntype,pid-1,0,MPI_COMM_WORLD,status,ierr)
                  call MPI_recv(maskline(1),1,columntypeint,pid-1,1,MPI_COMM_WORLD,status,ierr)

                  do i=1,n2
                      maskold(i,2)=maskold(i,2)+maskline(i)
                  end do



              endif
          else
              nmax = nend(pid) - nini(pid)+1
              if(mod(pid,2).eq.0) then
                  call MPI_recv(wtd(1,nmax),1,columntype,pid+1,0,MPI_COMM_WORLD,status,ierr)
                  call MPI_recv(wtd(1,1),1,columntype,pid-1,0,MPI_COMM_WORLD,status,ierr)
                  call MPI_recv(maskline(1),1,columntypeint,pid+1,1,MPI_COMM_WORLD,status,ierr)

                  do i=1,n2
                      maskold(i,nmax-1) = maskold(i,nmax-1)+maskline(i)
                  end do

                  call MPI_recv(maskline(1),1,columntypeint,pid-1,1,MPI_COMM_WORLD,status,ierr)

                  do i=1,n2
                      maskold(i,2) = maskold(i,2) + maskline(i)
                  end do

                  call MPI_send(wtd(1,nmax-1),1,columntype,pid+1,0,MPI_COMM_WORLD,ierr)
                  call MPI_send(wtd(1,2),1,columntype,pid-1,0,MPI_COMM_WORLD,ierr)
                  call MPI_send(maskold(1,nmax),1,columntypeint,pid+1,1,MPI_COMM_WORLD,ierr)
                  call MPI_send(maskold(1,1),1,columntypeint,pid-1,1,MPI_COMM_WORLD,ierr)
              else
                  call MPI_send(wtd(1,nmax-1),1,columntype,pid+1,0,MPI_COMM_WORLD,ierr)
                  call MPI_send(wtd(1,2),1,columntype,pid-1,0,MPI_COMM_WORLD,ierr)
                  call MPI_send(maskold(1,nmax),1,columntypeint,pid+1,1,MPI_COMM_WORLD,ierr)
                  call MPI_send(maskold(1,1),1,columntypeint,pid-1,1,MPI_COMM_WORLD,ierr)

                  call MPI_recv(wtd(1,nmax),1,columntype,pid+1,0,MPI_COMM_WORLD,status,ierr)
                  call MPI_recv(wtd(1,1),1,columntype,pid-1,0,MPI_COMM_WORLD,status,ierr)
                  call MPI_recv(maskline(1),1,columntypeint,pid+1,1,MPI_COMM_WORLD,status,ierr)


                  do i=1,n2
                      maskold(i,nmax-1)=maskold(i,nmax-1)+maskline(i)
                  end do



                  call MPI_recv(maskline(1),1,columntypeint,pid-1,1,MPI_COMM_WORLD,status,ierr)

                  do i=1,n2
                      maskold(i,2) = maskold(i,2)+maskline(i)
                  end do


              endif
        endif


      ENDIF 



      call MPI_ALLREDUCE(numbercount,numbertotal,1,MPI_INTEGER,mpi_sum,MPI_COMM_WORLD,ierr) !adds together results from multiple threads
!updates numbertotal





    if(mod((iter),10).eq.0.) then


        print *,'starting with the surface water loop!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

!reset variables to 0
        converged=0
        counter=0
        diff_total = 0
        fine=0

if(iter.gt.180000)then
        fine=1
endif
      
        SURFACE: do while(converged .eq. 0) 
           
            counter = counter + 1


            if(counter .eq.50)then    !select threshold here. Consider a better way to threshold or a max number of iterations if it isn't reaching that. BUT note it'll sometimes level out for a while and then still be able to take a big jump to improvement, so it's not so simple as just looking for when it levels out! 
                print *,'reached max iterations',diff_total
                converged = 1
            endif

           if(counter .gt.5 .and. diff_total .lt.1.and.fine.eq.0)then
                print *,'switching to fine grid',diff_total
                fine=1
           
            endif



            if(counter .gt.5 .and. diff_total .lt.1)then
                print *,'success',diff_total
                converged = 1
            endif


            if (pid .eq. 0) then

                print *,'Surface water counter',counter
                print *,'max',diff_total

 
            else !Any PID that is not 0

                nmax = nend(pid) - nini(pid) +3  
                

             


                if(fine.eq.0)then
                    n_use = nmax/5
                    n2_use = n2/5
                    allocate(topo_coarse(n2/5,(nmax+1)/5))
                    allocate(hz_read(n2/5,(nmax+1)/5))
                    allocate(diff(n2/5,(nmax+1)/5))
                    allocate(wtd_use(n2/5,(nmax+1)/5))
                    allocate(area_coarse(n2/5))
                    allocate(hold_read(n2/5,(nmax+1)/5))

print *,'line 1026'                    

                    row = 0
                     col = 0


                    do j=1,n2,5
                        col = 0
                        row = row + 1
                         do i=1,nmax+1,5
                             col = col +1


                             topo_coarse(row,col)=sum(topo(j:j+4,i:i+4))/25.
                            ! diff_coarse(row,col)=sum(diff(j:j+4,i:i+4))/25.
                             wtd_use(row,col)=sum(wtd(j:j+4,i:i+4))/25.
                             area_coarse = sum(area(j:j+4))/5

                        end do
                   end do
  hold_read=wtd_use
                hz_read = max(topo_coarse,topo_coarse+wtd_use) !water moves on the surface, so it must be AT LEAST the height of the topography, or if there is surface water then it's topo + wtd

  deallocate(topo_coarse)

print *,'line 1051'

                elseif(fine.eq.1)then
                 n_use=nmax
                 n2_use=n2

                allocate(hz_read(n2,nmax+1))
                allocate(diff(n2,nmax+1))
                allocate(wtd_use(n2,nmax+1))
allocate(hold_read(n2,nmax+1))
                wtd_use=wtd
                hz_read = max(topo,topo+wtd_use) !water moves on the surface, so it must be AT LEAST the height of the topography, or if there is surface water then it's topo + wtd
             hold_read=wtd_use

                endif

print *,'line 1067'
if(pid.ne.0)then
    Merge_size = n2_use*n_use           !Note that this may be faster with a smaller value. Largely because the actual sort takes longer, although the number of iterations may sometimes also be higher despite that not making much sense
    !Also note if unsure - a smaller value may only partially sort the array, but a too-large value may fill with zeros! 
    allocate(T((Merge_size+1)/2))
endif
             

                allocate(hz_1D(n2_use*n_use))  !Convert the array to a 1D for the mergesort

                do i=1,n2_use
                    hz_1D(((i-1)*n_use)+1:i*n_use) = hz_read(i,:) !unpack hz into a 1D array
                end do

                allocate(arr(2,n2_use*n_use))    !Create an array of the indices, for picking which cells to process first
                do row=1,n2_use
                    arr(1,((row-1)*n_use)+1:row*n_use)=row
                end do

                do row=1,n2_use
                    do col = 1,n_use
                        arr(2,(row-1)*n_use+col)=col
                    end do
                end do

                call MergeSort(hz_1D,Merge_size,T,arr)  !Sort to obtain order for processing the array - water flows from the highest cells first

print *,'starting cols1'
                COLS1: do i=2,n_use
                    ROWS1:do j=2,n2_use

                        row = arr(1,n2_use*n_use-(j-1)-(i-1)*n2_use) !get the next item in the sorted list to be processed, from the sorted index array
                        col = arr(2,n2_use*n_use-(j-1)-(i-1)*n2_use)


                 if(col.ge.n_use-1)then !Doing the end two columns separately, so skip them here
                    CYCLE
                endif

                if(col.le.2)then !Doing the end two columns separately, so skip them here
                    CYCLE
                endif

                if(row.eq.1 .or.row.eq.n2_use)then !Can't correctly process the very edge rows
                    CYCLE
                endif


                        if(wtd_use(row,col).le.0)then  !there is no surface water to move
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
                            upvalue = hz_read(row,col) - hz_read(row-1,col) !find the steepest direction, in which water would move
                            downvalue = hz_read(row,col) - hz_read(row+1,col)
                            leftvalue = hz_read(row,col) - hz_read(row,col-1)
                            rightvalue = hz_read(row,col) - hz_read(row,col+1)
            
                            if(max(upvalue,downvalue,leftvalue,rightvalue) .le. 0.) then     !should have already eliminated this option, but in case, it's the lowest local cell    
                                CYCLE   
                elseif(max(upvalue,downvalue,leftvalue,rightvalue) .eq. upvalue) then  !choose direction which is steepest
                    water = min(wtd_use(row,col)*area(row),upvalue*area(row)/2.)         !water is the minumum of the total water available, or of half the difference between 2 cells
                    water1 = water/area(row-1)                                       !get height of water change in the giving and receiving cells
                    water2 = water/area(row)

                    wtd_use(row,col) = wtd_use(row,col) - water2                              !update water table in the cells
                    wtd_use(row-1,col) = wtd_use(row-1,col) + water1

                    hz_read(row,col) = hz_read(row,col) - water2                      !update hz for the rest of the calculations
                    hz_read(row-1,col) = hz_read(row-1,col) + water1



                elseif(max(upvalue,downvalue,leftvalue,rightvalue) .eq. downvalue) then
                    water              = min(wtd_use(row,col)*area(row),downvalue*area(row)/2.)    
                    water1             = water/area(row+1)  
                    water2             = water/area(row)
                    wtd_use(row,col)   = wtd_use(row,col) - water2
                    wtd_use(row+1,col) = wtd_use(row+1,col) + water1

                    hz_read(row,col) = hz_read(row,col) - water2
                    hz_read(row+1,col) = hz_read(row+1,col) + water1

                elseif(max(upvalue,downvalue,leftvalue,rightvalue) .eq. rightvalue) then
                    water              = min(wtd_use(row,col)*area(row),rightvalue*area(row)/2.)    
                    water1             = water/area(row)  
                    water2             = water/area(row)
                    wtd_use(row,col)   = wtd_use(row,col) - water2
                    wtd_use(row,col+1) = wtd_use(row,col+1) + water1

                    hz_read(row,col) = hz_read(row,col) - water2
                    hz_read(row,col+1) = hz_read(row,col+1) + water1

                elseif(max(upvalue,downvalue,leftvalue,rightvalue) .eq. leftvalue) then
                    water              = min(wtd_use(row,col)*area(row),leftvalue*area(row)/2.)    
                    water1             = water/area(row)  
                    water2             = water/area(row)
                    wtd_use(row,col)   = wtd_use(row,col) - water2
                    wtd_use(row,col-1) = wtd_use(row,col-1) + water1

                    hz_read(row,col) = hz_read(row,col) - water2
                    hz_read(row,col-1) = hz_read(row,col-1) + water1

                endif
                         
                        endif

                    end do ROWS1
                end do COLS1
print *,'done with cols1'

       if(pid.eq.1)then        !send & receive the edge columns
            call MPI_recv(wtd_use(1,n_use+1),1,columntype,pid+1,8,MPI_COMM_WORLD,status,ierr)
            call MPI_recv(wtd_use(1,n_use  ),1,columntype,pid+1,9,MPI_COMM_WORLD,status,ierr)
            call MPI_send(wtd_use(1,3),1,columntype,numtasks-1,8,MPI_COMM_WORLD,ierr)
            call MPI_send(wtd_use(1,2),1,columntype,numtasks-1,9,MPI_COMM_WORLD,ierr)

            call MPI_recv(hz_read(1,n_use+1),1,columntype,pid+1,8,MPI_COMM_WORLD,status,ierr)
            call MPI_recv(hz_read(1,n_use  ),1,columntype,pid+1,9,MPI_COMM_WORLD,status,ierr)
            call MPI_send(hz_read(1,3),1,columntype,numtasks-1,8,MPI_COMM_WORLD,ierr)
            call MPI_send(hz_read(1,2),1,columntype,numtasks-1,9,MPI_COMM_WORLD,ierr)

        elseif(pid.eq.numtasks-1)then
            if(mod(pid,2).eq.0)then

            call MPI_recv(wtd_use(1,n_use+1),1,columntype,1,8,MPI_COMM_WORLD,status,ierr)
            call MPI_recv(wtd_use(1,n_use  ),1,columntype,1,9,MPI_COMM_WORLD,status,ierr)
            call MPI_send(wtd_use(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,ierr)   
            call MPI_send(wtd_use(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,ierr)
 
            call MPI_recv(hz_read(1,n_use+1),1,columntype,1,8,MPI_COMM_WORLD,status,ierr)
            call MPI_recv(hz_read(1,n_use  ),1,columntype,1,9,MPI_COMM_WORLD,status,ierr)
            call MPI_send(hz_read(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,ierr)   
            call MPI_send(hz_read(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,ierr)
 
            else

            call MPI_send(wtd_use(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,ierr)   
            call MPI_send(wtd_use(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,ierr)
            call MPI_recv(wtd_use(1,n_use+1),1,columntype,1,8,MPI_COMM_WORLD,status,ierr)
            call MPI_recv(wtd_use(1,n_use  ),1,columntype,1,9,MPI_COMM_WORLD,status,ierr)
 
            call MPI_send(hz_read(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,ierr)   
            call MPI_send(hz_read(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,ierr)
            call MPI_recv(hz_read(1,n_use+1),1,columntype,1,8,MPI_COMM_WORLD,status,ierr)
            call MPI_recv(hz_read(1,n_use  ),1,columntype,1,9,MPI_COMM_WORLD,status,ierr)


            endif

        else
            if(mod(pid,2).eq.0)then
                call MPI_recv(wtd_use(1,n_use+1),1,columntype,pid+1,8,MPI_COMM_WORLD,status,ierr)
                call MPI_recv(wtd_use(1,n_use  ),1,columntype,pid+1,9,MPI_COMM_WORLD,status,ierr)
                call MPI_send(wtd_use(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,ierr)
                call MPI_send(wtd_use(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,ierr)

                call MPI_recv(hz_read(1,n_use+1),1,columntype,pid+1,8,MPI_COMM_WORLD,status,ierr)
                call MPI_recv(hz_read(1,n_use  ),1,columntype,pid+1,9,MPI_COMM_WORLD,status,ierr)
                call MPI_send(hz_read(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,ierr)
                call MPI_send(hz_read(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,ierr)

            else
                call MPI_send(wtd_use(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,ierr)
                call MPI_send(wtd_use(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,ierr)
                call MPI_recv(wtd_use(1,n_use+1),1,columntype,pid+1,8,MPI_COMM_WORLD,status,ierr)
                call MPI_recv(wtd_use(1,n_use  ),1,columntype,pid+1,9,MPI_COMM_WORLD,status,ierr)

                call MPI_send(hz_read(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,ierr)
                call MPI_send(hz_read(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,ierr)
                call MPI_recv(hz_read(1,n_use+1),1,columntype,pid+1,8,MPI_COMM_WORLD,status,ierr)
                call MPI_recv(hz_read(1,n_use  ),1,columntype,pid+1,9,MPI_COMM_WORLD,status,ierr)

            endif
        endif

print *,'starting cols3'
                    COLS3: do col=n_use-1,n_use
                        ROWS3:do row=2,n2_use


                if(row.eq.1 .or.row.eq.n2_use)then 
                    CYCLE
                endif


                            if(wtd_use(row,col).le.0)then
                              CYCLE
                            endif
 
                            if(landmask(row,col) .eq. 0) then
                                wtd_use(row,col) = 0
                            endif

                  
                            if(hz_read(row,col) .le. hz_read(row+1,col) .and. hz_read(row,col) &
                            .le. hz_read(row-1,col) .and. hz_read(row,col) .le. hz_read(row,col+1) .and. &
                            hz_read(row,col) .le. hz_read(row,col-1)) then
                                CYCLE
                            else

                                upvalue = hz_read(row,col) - hz_read(row-1,col)
                                downvalue = hz_read(row,col) - hz_read(row+1,col)
                                leftvalue = hz_read(row,col) - hz_read(row,col-1)
                                rightvalue = hz_read(row,col) - hz_read(row,col+1)

                                if(max(upvalue,downvalue,leftvalue,rightvalue) .le. 0.) then         
                                    CYCLE   
                         elseif(max(upvalue,downvalue,leftvalue,rightvalue) .eq. upvalue) then
                    water = min(wtd_use(row,col)*area(row),upvalue*area(row)/2.)   
                    water1 = water/area(row-1)  
                    water2 = water/area(row)
                           wtd_use(row,col) = wtd_use(row,col) - water2
                            wtd_use(row-1,col) = wtd_use(row-1,col) + water1

                    hz_read(row,col) = hz_read(row,col) - water2
                    hz_read(row-1,col) = hz_read(row-1,col) + water1


                        elseif(max(upvalue,downvalue,leftvalue,rightvalue) .eq. downvalue) then
                    water = min(wtd_use(row,col)*area(row),downvalue*area(row)/2.)   
                    water1 = water/area(row+1)  
                    water2 = water/area(row)
                           wtd_use(row,col) = wtd_use(row,col) - water2
                            wtd_use(row+1,col) = wtd_use(row+1,col) + water1

                    hz_read(row,col) = hz_read(row,col) - water2
                    hz_read(row+1,col) = hz_read(row+1,col) + water1

                        elseif(max(upvalue,downvalue,leftvalue,rightvalue) .eq. rightvalue) then
                    water = min(wtd_use(row,col)*area(row),rightvalue*area(row)/2.)   
                    water1 = water/area(row)  
                    water2 = water/area(row)
                           wtd_use(row,col) = wtd_use(row,col) - water2
                            wtd_use(row,col+1) = wtd_use(row,col+1) + water1

                    hz_read(row,col) = hz_read(row,col) - water2
                    hz_read(row,col+1) = hz_read(row,col+1) + water1

                        elseif(max(upvalue,downvalue,leftvalue,rightvalue) .eq. leftvalue) then
                    water = min(wtd_use(row,col)*area(row),leftvalue*area(row)/2.)   
                    water1 = water/area(row)  
                    water2 = water/area(row)
                           wtd_use(row,col) = wtd_use(row,col) - water2
                            wtd_use(row,col-1) = wtd_use(row,col-1) + water1

                    hz_read(row,col) = hz_read(row,col) - water2
                    hz_read(row,col-1) = hz_read(row,col-1) + water1

                        endif

                            
                            endif

                        end do ROWS3
                    end do COLS3
             
print *,'done with cols3'
       if(pid.eq.1)then  !Once again send & receive the edge columns
            call MPI_send(wtd_use(1,n_use+1),1,columntype,pid+1,8,MPI_COMM_WORLD,ierr)
            call MPI_send(wtd_use(1,n_use  ),1,columntype,pid+1,9,MPI_COMM_WORLD,ierr)
            call MPI_recv(wtd_use(1,3),1,columntype,numtasks-1,8,MPI_COMM_WORLD,status,ierr)
            call MPI_recv(wtd_use(1,2),1,columntype,numtasks-1,9,MPI_COMM_WORLD,status,ierr)

            call MPI_send(hz_read(1,n_use+1),1,columntype,pid+1,8,MPI_COMM_WORLD,ierr)
            call MPI_send(hz_read(1,n_use  ),1,columntype,pid+1,9,MPI_COMM_WORLD,ierr)
            call MPI_recv(hz_read(1,3),1,columntype,numtasks-1,8,MPI_COMM_WORLD,status,ierr)
            call MPI_recv(hz_read(1,2),1,columntype,numtasks-1,9,MPI_COMM_WORLD,status,ierr)

        elseif(pid.eq.numtasks-1)then
            if(mod(pid,2).eq.0)then


            call MPI_send(wtd_use(1,n_use+1),1,columntype,1,8,MPI_COMM_WORLD,ierr)
            call MPI_send(wtd_use(1,n_use  ),1,columntype,1,9,MPI_COMM_WORLD,ierr)
            call MPI_recv(wtd_use(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,status,ierr)   
            call MPI_recv(wtd_use(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,status,ierr)

            call MPI_send(hz_read(1,n_use+1),1,columntype,1,8,MPI_COMM_WORLD,ierr)
            call MPI_send(hz_read(1,n_use  ),1,columntype,1,9,MPI_COMM_WORLD,ierr)
            call MPI_recv(hz_read(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,status,ierr)   
            call MPI_recv(hz_read(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,status,ierr)


            else

            call MPI_recv(wtd_use(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,status,ierr)   
            call MPI_recv(wtd_use(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,status,ierr)
            call MPI_send(wtd_use(1,n_use+1),1,columntype,1,8,MPI_COMM_WORLD,ierr)
            call MPI_send(wtd_use(1,n_use  ),1,columntype,1,9,MPI_COMM_WORLD,ierr)

            call MPI_recv(hz_read(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,status,ierr)   
            call MPI_recv(hz_read(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,status,ierr)
            call MPI_send(hz_read(1,n_use+1),1,columntype,1,8,MPI_COMM_WORLD,ierr)
            call MPI_send(hz_read(1,n_use  ),1,columntype,1,9,MPI_COMM_WORLD,ierr)

            endif

        else
            if(mod(pid,2).eq.0)then
                call MPI_recv(wtd_use(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,status,ierr)
                call MPI_recv(wtd_use(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,status,ierr)
                call MPI_send(wtd_use(1,n_use+1),1,columntype,pid+1,8,MPI_COMM_WORLD,ierr)
                call MPI_send(wtd_use(1,n_use  ),1,columntype,pid+1,9,MPI_COMM_WORLD,ierr)

                call MPI_recv(hz_read(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,status,ierr)
                call MPI_recv(hz_read(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,status,ierr)
                call MPI_send(hz_read(1,n_use+1),1,columntype,pid+1,8,MPI_COMM_WORLD,ierr)
                call MPI_send(hz_read(1,n_use  ),1,columntype,pid+1,9,MPI_COMM_WORLD,ierr)

            else
                call MPI_send(wtd_use(1,n_use+1),1,columntype,pid+1,8,MPI_COMM_WORLD,ierr)
                call MPI_send(wtd_use(1,n_use  ),1,columntype,pid+1,9,MPI_COMM_WORLD,ierr)
                call MPI_recv(wtd_use(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,status,ierr)
                call MPI_recv(wtd_use(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,status,ierr)

                call MPI_send(hz_read(1,n_use+1),1,columntype,pid+1,8,MPI_COMM_WORLD,ierr)
                call MPI_send(hz_read(1,n_use  ),1,columntype,pid+1,9,MPI_COMM_WORLD,ierr)
                call MPI_recv(hz_read(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,status,ierr)
                call MPI_recv(hz_read(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,status,ierr)

            endif
        endif


print *,'finish up the lasts'
row=0
                   do j=1,n2,5
                        col = 0
                        row = row + 1
                         do i=1,nmax+1,5
                             col = col +1
                             
                             wtd(j:j+4,i:i+4) = wtd_use(row,col)

                        end do
                   end do
 

print *,'wtd assigned'

            diff = abs(hold_read-wtd_use)
print *,'diff'
            maxdiff = 0
            maxdiff = maxval(diff)
print *,'deallocate'
   !         deallocate(hz_read)
print *,'hz_read'
    !        deallocate(diff)
print *,'diff del'
     !       deallocate(hz_1D)
print *,'hz_1D'
      !      deallocate(arr)
print *,'arr'
       !     deallocate(wtd_use)
print *,'wtd use'


        endif
print *,'do the allreduce'
        call MPI_ALLREDUCE(maxdiff,diff_total,1,MPI_REAL,mpi_max,MPI_COMM_WORLD,ierr)
   
        end do SURFACE


endif





  END DO GROUNDWATER !done with the big piece of processing! We now have 98% equilibrium! 


  if(pid.eq.0) then
      write(6,*) 'done; numbertotal = ',numbertotal,'iterations = ',iter
      write(15,*) 'done; numbertotal = ',numbertotal,'iterations = ',iter

      allocate(wtdglob(n2,n3))
      wtdglob = 0.

      open(23,file = trim(region)//'_equilibrium_SW.dat',form='unformatted',access='direct',recl=n2*n3) !do the final write - create the file

      do n=1,numtasks-1
          call MPI_recv(wtdglob(1,nini(n)+1),1,domblocksmall(n),n,8,MPI_COMM_WORLD,status,ierr) !receive the final wtd data from everyone
      end do

      write(23,rec=1)((wtdglob(i,j),i=1,n2),j=1,n3) !and write it to file

      close(23)

  else
      call MPI_send(wtd(1,2),1,domblocksmall(pid),0,8,MPI_COMM_WORLD,ierr) !everyone sends out the final result to pid 0! 

  endif


  write (6,*)'about to try deallocating'
  write (15,*)'about to try deallocating'
  deallocate(topo,stat=error)
  if (error.ne.0) then
      print *, 'topo error'
  endif
  deallocate(landmask,stat=error)
  if (error.ne.0)then
      print *,'landmask error'
  endif
  deallocate(wtd,stat=error)
  if (error.ne.0)then
      print *,'wtd error'
  endif
  deallocate(fdepth_start,stat=error)
  if (error.ne.0) then
      print *,'fdepth error'
  endif 
deallocate(fslope,stat=error)
  if (error.ne.0) then
      print *,'fslope error'
  endif
deallocate(fdepth,stat=error)
  if (error.ne.0) then
      print *,'fslope error'
  endif

  deallocate(ksat,stat=error)
  if (error.ne.0) then
      print *,'ksat error'
  endif
  deallocate(rechmean,stat=error)
  if (error.ne.0) then
      print *,'rechmean error'
  endif
  deallocate(kcell,stat=error)
  if (error.ne.0) then
      print *,'kcell error'
  endif
  deallocate(head,stat=error)
  if (error.ne.0) then
      print *,'head error'
  endif
  deallocate(qlat,stat=error)
  if (error.ne.0) then
      print *,'qlat error'
  endif
  deallocate(mask,maskold,stat=error)
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
  deallocate(xlat,stat=error)
  if (error.ne.0) then
      print *,'xlat error'
  endif
  deallocate(alpha,stat=error)
  if (error.ne.0) then
      print *, 'alpha error'
  endif

 deallocate(temp_sent,stat=error)
  if (error.ne.0) then
      print *, 'temp error'
  endif



  write(6,*)'done'
  write(15,*)'done'




  call MPI_FINALIZE(ierr)
end program equilibrium

