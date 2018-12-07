!this version tries to obtain equilibrium for both groundwater and surface water
!I am retyping/copying/reordering to try to make sure formatting looks good as well as to try iron out errors
!I am also trying to get rid of temperature and mask files, we don't need these, can do a completed fdepth and just use topo as mask. 
!this should be run at LGM to obtain the starting water tables for future transient runs
!one model run per setup - ie one for TRACE, one for HADCM3, although it can be done any time you want equilibrium
!For surface water, we run on a coarser grid then switch to finer, to help with processing times


!SUBROUTINES:

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


subroutine dividedomain(n2,n3,numtasks,nini,filetopo,ntotal) !divide among processors
  use netcdf
  implicit none

 
  integer :: n2,n3,numtasks
  integer :: nini(1:numtasks-1)
  real,allocatable,dimension(:,:) :: topo_read
  integer,allocatable,dimension(:) :: ncells
  integer :: iret,ncid,varid,ntotal,ncount,n,j
  character*100 :: filetopo
  real, parameter :: UNDEF = -1.0E+7

  allocate(topo_read(n2,n3))

  write(6,*)'reading in the mask to divide the domain'

  iret = nf90_open(filetopo,0,ncid)  !open the topography file
  call check_err(iret)
  write(6,*)'first call'

  iret = nf90_inq_varid(ncid,'value',varid) !get the ID of the value layer
  call check_err(iret)
  write(6,*) 'second call'

  iret = nf90_get_var(ncid,varid,topo_read) !read the actual values into the array called topo_read
  call check_err(iret)
  write(6,*) 'third call'

  iret = nf90_close(ncid) !close the mask file
  call check_err(iret)
  write(6,*)'fourth call'

  ntotal = count(topo_read>UNDEF) !count the number of land cells. I changed this slightly since I am using mask here rather than topo; all cells with a value of 1 should be included. 


  allocate(ncells(n3))
 
  ncells = count(topo_read>UNDEF,1) !The number of cells which are defined in the 1 dimension

  ncount=0

  nini(1) = 2

  n=2 !counter
  
  do j=1,n3
      ncount=ncount+ncells(j) !add the number of cells defined in the current column
      if (ncount .ge. ntotal/(numtasks-1)) then !>= total number of defined cells/number of threads
          nini(n) = j-1 !Telling it which row it needs to start in, to divide up the work equally, so each task handles the same number of cells.
          ncount = ncells(j) !reset the ncount for the next task
          n = n+1
      endif
      if(n .eq. numtasks) exit
  end do

  deallocate(topo_read,ncells)

  return

end subroutine dividedomain

  
!********************************************************************************************************

subroutine check_err(statusnc)

  use netcdf
  integer statusnc

  if(statusnc.ne. nf90_noerr) then
      stop 'Stopped due to a catch in the check_err subroutine'
  endif

end subroutine check_err

!********************************************************************************************************


program coupled_GW_SW

  use netcdf
  use mpi

  implicit none !ensures you don't create any variables without explicitly defining their type

!define all variables:

!character:
  character*8 :: region
character*60 :: surfdatadir
character*100 :: filetopo,filefdepth,fileksat,filerech


!integer:
integer n2, n3, ntotal, iret, ncid, varid, n, j, i, nmax, numbertotal, numbercount, iter, error, &
numberold, ierr, pid, numtasks, rc, tasktype, status(MPI_STATUS_SIZE), columntype, columntypeint, &
fine, n_use, n2_use, col, row, converged, counter, Merge_size

integer,allocatable :: domblock(:), domblocksmall(:), nini(:), nend(:)
integer*1, allocatable, dimension(:) :: maskline, T
integer*1, allocatable, dimension (:,:) :: landmask, maskold, watermask,mask


!real:
real thres, d0, d1, d2, d3, deltat, dx, dy, xn, xs, q, total, upvalue, downvalue, &
leftvalue, rightvalue, water, water1, water2, diff_total, maxdiff, placeholder, delta_xy, &
SEDGE,div_count

real, allocatable, dimension(:) :: alpha, xlat, alphamonth, hz_1D, area, area_use

real, allocatable, dimension(:,:) :: topo, wtd, kcell, ksat, head, qlat, rechmean, fdepth, &
topo_read, wtd_read, fdepth_read, ksat_read, rech_read, rech_month_read, rech_month, wtd_save, &
h, hold, diff, hz_read, hold_read, topo_use, wtd_use, diff_coarse, arr,wtdglob

real, parameter :: UNDEF = -1.0E+7
REAL(KIND=8),PARAMETER :: pi=3.141592653589793D0+0

!*******************************************************************************************
 
  call MPI_INIT(ierr)       !Initialise MPI
  if (ierr .ne. MPI_SUCCESS) then       !Error catching - this part should never run.
      print *,'Error starting MPI program. Terminating.'
      call MPI_ABORT(MPI_COMM_WORLD,rc,ierr)
  end if


!*******************************************************************************************

!set global parameters:
  n2 = 2000 !number of cells in x and y directions
  n3 = 1000
  region = 'Mad_' !start of file name for the files. Limit 8 characters. 
  surfdatadir = 'surfdata/' !folder names for the data, limit 60 characters. Topo, fdepth, and rech go here. 
!ksat is not within a folder since the same one is being used for all times.

  Thres = 0.001 !when to stop calculating - 0.01, 0.005, or 0.001. This is the finest option. 
  SEDGE = -27 !southern most latitude of your domain

  delta_xy = 120
  deltat = (365.*24.*3600.) !seconds in an annual timestep

!input files:

  filetopo = trim(surfdatadir)//trim(region)//'021000_topo_rotated.nc'
  filefdepth = trim(surfdatadir)//trim(region)//'021000_fslope_rotated.nc' !CHANGE TO FDEPTH
  filerech = trim(surfdatadir)//trim(region)//'021000_rech_rotated.nc'
  fileksat = trim(surfdatadir)//trim(region)//'ksat_rotated.nc'

!output text file:

  open (15,file=trim(region)//'output_equilibrium_coupled.txt') 

!different increments to use to move water table depth:

  d0 = Thres/2
  d1 = 0.05
  d2 = 0.1
  d3 = 0.25
 
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
      write(15,*) 'PID = 0'
      call dividedomain(n2,n3,numtasks,nini,filetopo,ntotal) !Divides up the work equally among all of the ranks, by number of defined land cells.
      write (15,*) 'Dividedomain done' 
      print *,'dividedomain done'

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


 
  do n=1, numtasks-1
      nmax = nend(n) - nini(n) + 4 !max number of columns we have in this group

      call MPI_TYPE_CONTIGUOUS(n2*nmax,MPI_REAL,domblock(n),ierr)
      call MPI_type_commit(domblock(n),ierr)


      call MPI_TYPE_CONTIGUOUS(n2*(nmax-3),MPI_REAL,domblocksmall(n),ierr)
      call MPI_type_commit(domblocksmall(n),ierr)
  end do

  call MPI_TYPE_CONTIGUOUS(n2,MPI_REAL,columntype,ierr)
  call MPI_type_commit(columntype,ierr)

  call MPI_TYPE_CONTIGUOUS(n2,MPI_INTEGER1,columntypeint,ierr)
  call MPI_type_commit(columntypeint,ierr)




!Get the changing cell areas across different latitudes:

  dy = 6370000.*pi/(180.*delta_xy) !radius of the earth * pi/number of cells in a half-circle of the globe. This should equal the height of each cell in the N-S direction. 

  dx = dy
  
  if(pid .gt. 0) then
      allocate(xlat(n2))
      allocate(area(n2))
      allocate(alpha(n2))
      allocate(alphamonth(n2))

      do j=1, n2 !changing area of cell depending on its latitude
          xlat(j) = (float(j)/delta_xy + SEDGE)*pi/180.
          xs = (float(2*(j)-1)/(delta_xy*2.)+SEDGE)*pi/180. !Latitude on southern cell edge
          xn = (float(2*(j)+1)/(delta_xy*2.)+SEDGE)*pi/180. !latitude on northern cell edge

          placeholder = dy*6370000.*(sin(xn)-sin(xs))/2.
          area(j) = placeholder !final cell area for that latitude: trapezoid dy * dx

          alpha(j) = deltat/placeholder 
          alphamonth(j) = (deltat/12.)/placeholder
      end do
  end if



!*******************************************************************************************


!import all of the data:

  if(pid .eq. 0) then

      allocate(wtd_read(n2,n3))
      wtd_read = 0. 

      allocate(topo_read(n2,n3))
     
      iret = nf90_open(filetopo,0,ncid) !reading in the topo
      call check_err(iret)

      iret = nf90_inq_varid(ncid,'value',varid)
      call check_err(iret)

      iret = nf90_get_var(ncid,varid,topo_read)
      call check_err(iret)

      iret = nf90_close(ncid)
      call check_err(iret)


      where(topo_read .le. UNDEF) topo_read = 0. !change undefined cells to 0
     
    
      allocate(ksat_read(n2,n3))
     
      iret = nf90_open(fileksat,0,ncid) !reading in the ksat
      call check_err(iret)

      iret = nf90_inq_varid(ncid,'value',varid)
      call check_err(iret)

      iret = nf90_get_var(ncid,varid,ksat_read)
      call check_err(iret)

      iret = nf90_close(ncid)
      call check_err(iret)


      allocate(rech_read(n2,n3))
      allocate(rech_month_read(n2,n3))
     
      iret = nf90_open(filerech,0,ncid) !reading in the recharge file
      call check_err(iret)

      iret = nf90_inq_varid(ncid,'value',varid)
      call check_err(iret)

      iret = nf90_get_var(ncid,varid,rech_read)
      call check_err(iret)

      iret = nf90_close(ncid)
      call check_err(iret)

!####################################################################################### potentially change with lakes
      rech_read = max(rech_read,0.) !setting it so recharge can only be positive
    rech_month_read = (rech_read/12.)  

      allocate(fdepth_read(n2,n3))
    
      iret = nf90_open(filefdepth,0,ncid) !reading in the fdepth
      call check_err(iret)

      iret = nf90_inq_varid(ncid,'value',varid)
      call check_err(iret)

      iret = nf90_get_var(ncid,varid,fdepth_read)
      call check_err(iret)

      iret = nf90_close(ncid)
      call check_err(iret)

      where(fdepth_read .le. 0.0001) fdepth_read = 0.0001

!now send everything we have opened:

      do n=1,numtasks-2
          call MPI_send(topo_read(1,nini(n)-1),1,domblock(n),n,1,MPI_COMM_WORLD,ierr)
          call MPI_send(wtd_read(1,nini(n)-1),1,domblock(n),n,3,MPI_COMM_WORLD,ierr)
          call MPI_send(ksat_read(1,nini(n)-1),1,domblock(n),n,5,MPI_COMM_WORLD,ierr)
          call MPI_send(rech_read(1,nini(n)-1),1,domblock(n),n,6,MPI_COMM_WORLD,ierr)
          call MPI_send(rech_month_read(1,nini(n)-1),1,domblock(n),n,7,MPI_COMM_WORLD,ierr)
          call MPI_send(fdepth_read(1,nini(n)-1),1,domblock(n),n,8,MPI_COMM_WORLD,ierr)
      end do

      call MPI_send(topo_read(1,nini(numtasks-1)-1),1,domblocksmall(numtasks-1),numtasks-1,1,MPI_COMM_WORLD,ierr)
      call MPI_send(wtd_read(1,nini(numtasks-1)-1),1,domblocksmall(numtasks-1),numtasks-1,3,MPI_COMM_WORLD,ierr)
      call MPI_send(ksat_read(1,nini(numtasks-1)-1),1,domblocksmall(numtasks-1),numtasks-1,5,MPI_COMM_WORLD,ierr)

      call MPI_send(rech_read(1,nini(numtasks-1)-1),1,domblocksmall(numtasks-1),numtasks-1,6,MPI_COMM_WORLD,ierr)
    call MPI_send(rech_month_read(1,nini(numtasks-1)-1),1,domblocksmall(numtasks-1),numtasks-1,7,MPI_COMM_WORLD,ierr)
      call MPI_send(fdepth_read(1,nini(numtasks-1)-1),1,domblocksmall(numtasks-1),numtasks-1,8,MPI_COMM_WORLD,ierr)
   

      call MPI_send(topo_read(1,1),1,columntype,numtasks-1,1,MPI_COMM_WORLD,ierr)
      call MPI_send(wtd_read(1,1),1,columntype,numtasks-1,3,MPI_COMM_WORLD,ierr)
      call MPI_send(ksat_read(1,1),1,columntype,numtasks-1,5,MPI_COMM_WORLD,ierr)
      call MPI_send(rech_read(1,1),1,columntype,numtasks-1,6,MPI_COMM_WORLD,ierr)
    call MPI_send(rech_month_read(1,1),1,columntype,numtasks-1,7,MPI_COMM_WORLD,ierr)
      call MPI_send(fdepth_read(1,1),1,columntype,numtasks-1,8,MPI_COMM_WORLD,ierr)
   
      call MPI_send(topo_read(1,2),1,columntype,numtasks-1,1,MPI_COMM_WORLD,ierr)
      call MPI_send(wtd_read(1,2),1,columntype,numtasks-1,3,MPI_COMM_WORLD,ierr)
      call MPI_send(ksat_read(1,2),1,columntype,numtasks-1,5,MPI_COMM_WORLD,ierr)
      call MPI_send(rech_read(1,2),1,columntype,numtasks-1,6,MPI_COMM_WORLD,ierr)
    call MPI_send(rech_month_read(1,2),1,columntype,numtasks-1,7,MPI_COMM_WORLD,ierr)
      call MPI_send(fdepth_read(1,2),1,columntype,numtasks-1,8,MPI_COMM_WORLD,ierr)
   

      call MPI_send(topo_read(1,3),1,columntype,numtasks-1,1,MPI_COMM_WORLD,ierr)
      call MPI_send(wtd_read(1,3),1,columntype,numtasks-1,3,MPI_COMM_WORLD,ierr)
      call MPI_send(ksat_read(1,3),1,columntype,numtasks-1,5,MPI_COMM_WORLD,ierr)
      call MPI_send(rech_read(1,3),1,columntype,numtasks-1,6,MPI_COMM_WORLD,ierr)
    call MPI_send(rech_month_read(1,3),1,columntype,numtasks-1,7,MPI_COMM_WORLD,ierr)
      call MPI_send(fdepth_read(1,3),1,columntype,numtasks-1,8,MPI_COMM_WORLD,ierr)
  

      deallocate(topo_read)
      deallocate(ksat_read)
      deallocate(wtd_read)
      deallocate(rech_read)
      deallocate(rech_month_read)
      deallocate(fdepth_read)

  else
  
      allocate(topo(n2,nmax))
      allocate(wtd(n2,nmax))
      allocate(ksat(n2,nmax))
      allocate(rechmean(n2,nmax))
      allocate(rech_month(n2,nmax))
      allocate(kcell(n2,nmax))
      allocate(head(n2,nmax))
      allocate(qlat(n2,nmax))
      allocate(maskold(n2,nmax))
      allocate(mask(n2,nmax))
      allocate(fdepth(n2,nmax))
      allocate(landmask(n2,nmax))

    if(pid.lt.numtasks-1)then

          call MPI_recv(topo(1,1),1,domblock(pid),0,1,MPI_COMM_WORLD,status,ierr) !receiving everthing that was sent above     
          call MPI_recv(wtd(1,1),1,domblock(pid),0,3,MPI_COMM_WORLD,status,ierr)
          call MPI_recv(ksat(1,1),1,domblock(pid),0,5,MPI_COMM_WORLD,status,ierr)
          call MPI_recv(rechmean(1,1),1,domblock(pid),0,6,MPI_COMM_WORLD,status,ierr)
      call MPI_recv(rech_month(1,1),1,domblock(pid),0,7,MPI_COMM_WORLD,status,ierr)
          call MPI_recv(fdepth(1,1),1,domblock(pid),0,8,MPI_COMM_WORLD,status,ierr)


   else

          call MPI_recv(topo(1,1),1,domblocksmall(pid),0,1,MPI_COMM_WORLD,status,ierr)
          call MPI_recv(wtd(1,1),1,domblocksmall(pid),0,3,MPI_COMM_WORLD,status,ierr)
          call MPI_recv(ksat(1,1),1,domblocksmall(pid),0,5,MPI_COMM_WORLD,status,ierr)

          call MPI_recv(rechmean(1,1),1,domblocksmall(pid),0,6,MPI_COMM_WORLD,status,ierr)
          call MPI_recv(rech_month(1,1),1,domblocksmall(pid),0,7,MPI_COMM_WORLD,status,ierr)
          call MPI_recv(fdepth(1,1),1,domblocksmall(pid),0,8,MPI_COMM_WORLD,status,ierr)
          call MPI_recv(topo(1,nmax-2),1,columntype,0,1,MPI_COMM_WORLD,status,ierr)
          call MPI_recv(wtd(1,nmax-2),1,columntype,0,3,MPI_COMM_WORLD,status,ierr)
          call MPI_recv(ksat(1,nmax-2),1,columntype,0,5,MPI_COMM_WORLD,status,ierr)
          call MPI_recv(rechmean(1,nmax-2),1,columntype,0,6,MPI_COMM_WORLD,status,ierr)
       call MPI_recv(rech_month(1,nmax-2),1,columntype,0,7,MPI_COMM_WORLD,status,ierr)
          call MPI_recv(fdepth(1,nmax-2),1,columntype,0,8,MPI_COMM_WORLD,status,ierr)
  
          call MPI_recv(topo(1,nmax-1),1,columntype,0,1,MPI_COMM_WORLD,status,ierr)
          call MPI_recv(wtd(1,nmax-1),1,columntype,0,3,MPI_COMM_WORLD,status,ierr)
          call MPI_recv(ksat(1,nmax-1),1,columntype,0,5,MPI_COMM_WORLD,status,ierr)
          call MPI_recv(rechmean(1,nmax-1),1,columntype,0,6,MPI_COMM_WORLD,status,ierr)
       call MPI_recv(rech_month(1,nmax-1),1,columntype,0,7,MPI_COMM_WORLD,status,ierr)
          call MPI_recv(fdepth(1,nmax-1),1,columntype,0,8,MPI_COMM_WORLD,status,ierr)

          call MPI_recv(topo(1,nmax),1,columntype,0,1,MPI_COMM_WORLD,status,ierr)
          call MPI_recv(wtd(1,nmax),1,columntype,0,3,MPI_COMM_WORLD,status,ierr)
          call MPI_recv(ksat(1,nmax),1,columntype,0,5,MPI_COMM_WORLD,status,ierr)
          call MPI_recv(rechmean(1,nmax),1,columntype,0,6,MPI_COMM_WORLD,status,ierr)
       call MPI_recv(rech_month(1,nmax),1,columntype,0,7,MPI_COMM_WORLD,status,ierr)
          call MPI_recv(fdepth(1,nmax),1,columntype,0,8,MPI_COMM_WORLD,status,ierr)   
    endif
   

      maskold = 1 !This is just being done ahead of the big loop, maskold is used to show which cells still need to be processed. 
      mask = 0
      landmask = topo

  endif


!Get ready to start the groundwater loop:

  iter = 0
  numbercount = 0
  numbertotal = ntotal
 

 
GROUNDWATER : DO while(numbertotal > ntotal/100 .and. iter < 40) !start of the main loop with conditions - number of iterations or % equilibrium. Currently set for 99% equilibrium or 400000 iterations; can be changed if this is too much
      
      if(iter .eq. 50000) then !here we switch from annual to monthly processing. 
          write(15,*) '50000 iterations, adjusting the values for monthly processing'
          thres = thres/12.  !these need to be adjusted since the total water moving for a month will be 12x< than that for a year, so the threshold must also be smaller because of the way we decide if it's at equilibrium. 
          d0 = d0/12.
          d1 = d1/12.
          d2 = d2/12.
          d3 = d3/12.
          deltat = deltat/12.
          alpha = alphamonth
          rechmean = rech_month
          maskold = 1          !because we changed the threshold, so we want to recheck all the cells. 
          numbertotal = ntotal
      endif

      iter = iter + 1 ! count the iterations


      if (pid.eq.0)then !just for checking how far the run is
          write(15,*) 'PID = 0. Numbertotal:', numbertotal, 'ntotal:', ntotal, 'iter number:', iter
      endif

      numbercount = 0

      if(mod(float(iter),10000.) .eq. 0.) then  !save the data in case it crashes or you need multiple runs to reach equilibrium
          if(pid.eq.0)then
              allocate(wtd_save(n2,n3))
              open(23,file=trim(region)//'equilibrium_GW_SW.dat',form='unformatted',access='stream')

              do n=1,numtasks-1
                  write(15,*) 'receiving wtd',n   
                  call MPI_recv(wtd_save(1,nini(n)),1,domblocksmall(n),n,8,MPI_COMM_WORLD,status,ierr)
              end do
              
              write (23,rec=1)((wtd_save(i,j),i=1,n2),j=1,n3)

              deallocate(wtd_save)

          else
              call MPI_send(wtd(1,2),1,domblocksmall(pid),0,8,MPI_COMM_WORLD,ierr)
          endif
      endif



     IF (pid.gt.0) then
          mask = 0
          nmax = nend(pid) - nini(pid)
          do j=1, nmax+1
              do i=1, n2
                  if(fdepth(i,j) .gt. 0 .and.(maskold(i,j).gt.0.or.j.eq.1.or.j.eq.nmax+1)) then !maskold is used to stop us from repeating all of the steps on cells which have already reached equilibrium. 
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


          do j=2,nmax
              do i=2,n2-1

                  if(landmask(i,j) .gt. 0 .and. maskold(i,j).gt.0) then
                      q=0.
                      !north
                      q  = q + (kcell(i,j+1)+kcell(i,j))*(head(i,j+1)-head(i,j)) * cos(xlat(j)+pi/(180.*delta_xy*2.))   
!it seems like we are getting the total which will be discharged from each cell
                      !south
                      q  = q + (kcell(i,j-1)+kcell(i,j))*(head(i,j-1)-head(i,j)) * cos(xlat(j)-pi/(180.*delta_xy*2.))
                      !west
                      q  = q + (kcell(i-1,j)+kcell(i,j))*(head(i-1,j)-head(i,j)) / cos(xlat(j))
                      !east
                      q  = q + (kcell(i+1,j)+kcell(i,j))*(head(i+1,j)-head(i,j)) / cos(xlat(j))
                      qlat(i,j)=alpha(j)*q  
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
                      elseif(total .gt. 1.) then                  
                          wtd(i,j) = wtd(i,j) +d3
                          numbercount = numbercount + 1
                      elseif (total .gt. 0.25 ) then
                          wtd(i,j) = wtd(i,j) +d2
                          numbercount = numbercount + 1
                      elseif (total .gt. 0.1 ) then
                          wtd(i,j) = wtd(i,j) +d1 
                          numbercount = numbercount + 1
                      elseif (total .gt. thres ) then
                          wtd(i,j) = wtd(i,j) +d0
                          numbercount = numbercount + 1
                      !things which are between -thres and thres do not change; they are in equilibrium.
                      endif


                     if (numberold .ne. numbercount) then
                         mask(i+1,j) =1 !tag cells to show they are still not in equilibrium.
                         mask(i-1,j) = 1
                         mask(i,j+1) = 1
                         mask(i,j-1) = 1
                         mask(i,j) = 1
                     endif
 

                  endif
              end do
          end do

          maskold = mask  !to record which cells still need to be processed
         numbertotal = numbercount !total number of cells still out of equilibrium this iteration
                             
!send and receive the lines on either side of each section to allow for flow across these lines:
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










      if(mod(iter,10).eq.0)then

          write (15,*) 'starting with the surface water loop'
!reset variables to 0
              converged = 0
              counter = 0
        diff_total = 0
              fine = 0
if(iter.gt.180000)then
        fine=1
endif


      
        SURFACE: do while(converged .eq. 0) 
           
            counter = counter + 1


            if(counter .eq.50)then    !select threshold here. Consider a better way to threshold or a max number of iterations if it isn't reaching that. BUT note it'll sometimes level out for a while and then still be able to take a big jump to improvement, so it's not so simple as just looking for when it levels out! 
                print *,'reached max iterations',diff_total
                converged = 1
            endif

            if(counter .gt.5 .and. diff_total .lt.1)then
                  print *,'success'!,diff_total
                converged = 1
              endif

              if(pid.ne.0)then
                  Merge_size = n2_use*n_use           !Note that this may be faster with a smaller value. Largely because the actual sort takes longer, although the number of iterations may sometimes also be higher despite that not making much sense
    !Also note if unsure - a smaller value may only partially sort the array, but a too-large value may fill with zeros! 
                  allocate(T((Merge_size+1)/2))
            endif


            if (pid .eq. 0) then

                print *,'Surface water counter',counter
                print *,'max',diff_total

 
            else !Any PID that is not 0

                  hold_read = wtd
                nmax = nend(pid) - nini(pid) +3  
                

                allocate(hz_read(n2,nmax+1))
                allocate(diff(n2,nmax+1))

                  hz_read = max(topo,topo+wtd) !water moves on the surface, so it must be AT LEAST the height of the topography, or if there is surface water then it's topo + wtd

                  allocate(hz_1D(n2*nmax))  !Convert the array to a 1D for the mergesort
             
                  do i=1,n2
!print *,'wth',i
                      hz_1D(((i-1)*nmax)+1:i*nmax) = hz_read(i,1:nmax) !unpack hz into a 1D array
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
                      ROWS1:do j=2,n2

                          row = arr(1,n2*nmax-(j-1)-(i-1)*n2) !get the next item in the sorted list to be processed, from the sorted index array
                          col = arr(2,n2*nmax-(j-1)-(i-1)*n2)


                          if(col.ge.nmax-1)then !Doing the end two columns separately, so skip them here
                    CYCLE
                endif

                if(col.le.2)then !Doing the end two columns separately, so skip them here
                    CYCLE
                endif

                    CYCLE
                endif


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
                                  water = min(wtd(row,col)*area(row),upvalue*area(row)/2.)         !water is the minumum of the total water available, or of half the difference between 2 cells
                    water1 = water/area(row-1)                                       !get height of water change in the giving and receiving cells
                    water2 = water/area(row)

                                  wtd(row,col) = wtd(row,col) - water2                              !update water table in the cells
                                  wtd(row-1,col) = wtd(row-1,col) + water1

                    hz_read(row,col) = hz_read(row,col) - water2                      !update hz for the rest of the calculations
                    hz_read(row-1,col) = hz_read(row-1,col) + water1



                elseif(max(upvalue,downvalue,leftvalue,rightvalue) .eq. downvalue) then
                                  water = min(wtd(row,col)*area(row),downvalue*area(row)/2.)    
                                  water1 = water/area(row+1)  
                                  water2 = water/area(row)
                                  wtd(row,col) = wtd(row,col) - water2
                                  wtd(row+1,col) = wtd(row+1,col) + water1

                    hz_read(row,col) = hz_read(row,col) - water2
                    hz_read(row+1,col) = hz_read(row+1,col) + water1

                elseif(max(upvalue,downvalue,leftvalue,rightvalue) .eq. rightvalue) then
                                  water = min(wtd(row,col)*area(row),rightvalue*area(row)/2.)    
                                  water1 = water/area(row)  
                                  water2 = water/area(row)
                                  wtd(row,col) = wtd(row,col) - water2
                                  wtd(row,col+1) = wtd(row,col+1) + water1

                    hz_read(row,col) = hz_read(row,col) - water2
                    hz_read(row,col+1) = hz_read(row,col+1) + water1

                elseif(max(upvalue,downvalue,leftvalue,rightvalue) .eq. leftvalue) then
                                  water = min(wtd(row,col)*area(row),leftvalue*area(row)/2.)    
                                  water1 = water/area(row)  
                                  water2 = water/area(row)
                                  wtd(row,col) = wtd(row,col) - water2
                                  wtd(row,col-1) = wtd(row,col-1) + water1

                    hz_read(row,col) = hz_read(row,col) - water2
                    hz_read(row,col-1) = hz_read(row,col-1) + water1

                endif
                        endif

  
                    end do ROWS1
                end do COLS1



       if(pid.eq.1)then        !send & receive the edge columns
                      call MPI_recv(wtd(1,nmax+1),1,columntype,pid+1,8,MPI_COMM_WORLD,status,ierr)
                      call MPI_recv(wtd(1,nmax),1,columntype,pid+1,9,MPI_COMM_WORLD,status,ierr)
                      call MPI_send(wtd(1,3),1,columntype,numtasks-1,8,MPI_COMM_WORLD,ierr)
                      call MPI_send(wtd(1,2),1,columntype,numtasks-1,9,MPI_COMM_WORLD,ierr)

                      call MPI_recv(hz_read(1,nmax+1),1,columntype,pid+1,8,MPI_COMM_WORLD,status,ierr)
                      call MPI_recv(hz_read(1,nmax),1,columntype,pid+1,9,MPI_COMM_WORLD,status,ierr)
            call MPI_send(hz_read(1,3),1,columntype,numtasks-1,8,MPI_COMM_WORLD,ierr)
            call MPI_send(hz_read(1,2),1,columntype,numtasks-1,9,MPI_COMM_WORLD,ierr)

        elseif(pid.eq.numtasks-1)then
            if(mod(pid,2).eq.0)then

                          call MPI_recv(wtd(1,nmax+1),1,columntype,1,8,MPI_COMM_WORLD,status,ierr)
                          call MPI_recv(wtd(1,nmax),1,columntype,1,9,MPI_COMM_WORLD,status,ierr)
                          call MPI_send(wtd(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,ierr)   
                          call MPI_send(wtd(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,ierr)
 
                          call MPI_recv(hz_read(1,nmax+1),1,columntype,1,8,MPI_COMM_WORLD,status,ierr)
                          call MPI_recv(hz_read(1,nmax),1,columntype,1,9,MPI_COMM_WORLD,status,ierr)
            call MPI_send(hz_read(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,ierr)   
            call MPI_send(hz_read(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,ierr)
 
            else

                          call MPI_send(wtd(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,ierr)   
                          call MPI_send(wtd(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,ierr)
                          call MPI_recv(wtd(1,nmax+1),1,columntype,1,8,MPI_COMM_WORLD,status,ierr)
                          call MPI_recv(wtd(1,nmax),1,columntype,1,9,MPI_COMM_WORLD,status,ierr)
 
            call MPI_send(hz_read(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,ierr)   
            call MPI_send(hz_read(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,ierr)
                          call MPI_recv(hz_read(1,nmax+1),1,columntype,1,8,MPI_COMM_WORLD,status,ierr)
                          call MPI_recv(hz_read(1,nmax),1,columntype,1,9,MPI_COMM_WORLD,status,ierr)


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
                      ROWS3:do row=2,n2


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

                                upvalue = hz_read(row,col) - hz_read(row-1,col)
                                downvalue = hz_read(row,col) - hz_read(row+1,col)
                                leftvalue = hz_read(row,col) - hz_read(row,col-1)
                                rightvalue = hz_read(row,col) - hz_read(row,col+1)

                                if(max(upvalue,downvalue,leftvalue,rightvalue) .le. 0.) then         
                                    CYCLE   
                         elseif(max(upvalue,downvalue,leftvalue,rightvalue) .eq. upvalue) then
                                  water = min(wtd(row,col)*area(row),upvalue*area(row)/2.)   
                    water1 = water/area(row-1)  
                    water2 = water/area(row)
                                  wtd(row,col) = wtd(row,col) - water2
                                  wtd(row-1,col) = wtd(row-1,col) + water1

                    hz_read(row,col) = hz_read(row,col) - water2
                    hz_read(row-1,col) = hz_read(row-1,col) + water1


                        elseif(max(upvalue,downvalue,leftvalue,rightvalue) .eq. downvalue) then
                                  water = min(wtd(row,col)*area(row),downvalue*area(row)/2.)   
                    water1 = water/area(row+1)  
                    water2 = water/area(row)
                                  wtd(row,col) = wtd(row,col) - water2
                                  wtd(row+1,col) = wtd(row+1,col) + water1

                    hz_read(row,col) = hz_read(row,col) - water2
                    hz_read(row+1,col) = hz_read(row+1,col) + water1

                        elseif(max(upvalue,downvalue,leftvalue,rightvalue) .eq. rightvalue) then
                                  water = min(wtd(row,col)*area(row),rightvalue*area(row)/2.)   
                    water1 = water/area(row)  
                    water2 = water/area(row)
                                  wtd(row,col) = wtd(row,col) - water2
                                  wtd(row,col+1) = wtd(row,col+1) + water1

                    hz_read(row,col) = hz_read(row,col) - water2
                    hz_read(row,col+1) = hz_read(row,col+1) + water1

                        elseif(max(upvalue,downvalue,leftvalue,rightvalue) .eq. leftvalue) then
                                  water = min(wtd(row,col)*area(row),leftvalue*area(row)/2.)   
                    water1 = water/area(row)  
                    water2 = water/area(row)
                                  wtd(row,col) = wtd(row,col) - water2
                                  wtd(row,col-1) = wtd(row,col-1) + water1

                    hz_read(row,col) = hz_read(row,col) - water2
                    hz_read(row,col-1) = hz_read(row,col-1) + water1

                        endif

                            
                            endif

                        end do ROWS3
                    end do COLS3
             

       if(pid.eq.1)then  !Once again send & receive the edge columns
                      call MPI_send(wtd(1,nmax+1),1,columntype,pid+1,8,MPI_COMM_WORLD,ierr)
                      call MPI_send(wtd(1,nmax),1,columntype,pid+1,9,MPI_COMM_WORLD,ierr)
                      call MPI_recv(wtd(1,3),1,columntype,numtasks-1,8,MPI_COMM_WORLD,status,ierr)
                      call MPI_recv(wtd(1,2),1,columntype,numtasks-1,9,MPI_COMM_WORLD,status,ierr)

                      call MPI_send(hz_read(1,nmax+1),1,columntype,pid+1,8,MPI_COMM_WORLD,ierr)
                      call MPI_send(hz_read(1,nmax),1,columntype,pid+1,9,MPI_COMM_WORLD,ierr)
            call MPI_recv(hz_read(1,3),1,columntype,numtasks-1,8,MPI_COMM_WORLD,status,ierr)
            call MPI_recv(hz_read(1,2),1,columntype,numtasks-1,9,MPI_COMM_WORLD,status,ierr)

        elseif(pid.eq.numtasks-1)then
            if(mod(pid,2).eq.0)then


                          call MPI_send(wtd(1,nmax+1),1,columntype,1,8,MPI_COMM_WORLD,ierr)
                          call MPI_send(wtd(1,nmax),1,columntype,1,9,MPI_COMM_WORLD,ierr)
                          call MPI_recv(wtd(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,status,ierr)   
                          call MPI_recv(wtd(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,status,ierr)

                          call MPI_send(hz_read(1,nmax+1),1,columntype,1,8,MPI_COMM_WORLD,ierr)
                          call MPI_send(hz_read(1,nmax),1,columntype,1,9,MPI_COMM_WORLD,ierr)
            call MPI_recv(hz_read(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,status,ierr)   
            call MPI_recv(hz_read(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,status,ierr)


            else

                          call MPI_recv(wtd(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,status,ierr)   
                          call MPI_recv(wtd(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,status,ierr)
                          call MPI_send(wtd(1,nmax+1),1,columntype,1,8,MPI_COMM_WORLD,ierr)
                          call MPI_send(wtd(1,nmax),1,columntype,1,9,MPI_COMM_WORLD,ierr)

            call MPI_recv(hz_read(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,status,ierr)   
            call MPI_recv(hz_read(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,status,ierr)
                          call MPI_send(hz_read(1,nmax+1),1,columntype,1,8,MPI_COMM_WORLD,ierr)
                          call MPI_send(hz_read(1,nmax),1,columntype,1,9,MPI_COMM_WORLD,ierr)

            endif

        else
            if(mod(pid,2).eq.0)then
                          call MPI_recv(wtd(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,status,ierr)
                          call MPI_recv(wtd(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,status,ierr)
                          call MPI_send(wtd(1,nmax+1),1,columntype,pid+1,8,MPI_COMM_WORLD,ierr)
                          call MPI_send(wtd(1,nmax),1,columntype,pid+1,9,MPI_COMM_WORLD,ierr)

                call MPI_recv(hz_read(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,status,ierr)
                call MPI_recv(hz_read(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,status,ierr)
                          call MPI_send(hz_read(1,nmax+1),1,columntype,pid+1,8,MPI_COMM_WORLD,ierr)
                          call MPI_send(hz_read(1,nmax),1,columntype,pid+1,9,MPI_COMM_WORLD,ierr)

            else
                          call MPI_send(wtd(1,nmax+1),1,columntype,pid+1,8,MPI_COMM_WORLD,ierr)
                          call MPI_send(wtd(1,nmax),1,columntype,pid+1,9,MPI_COMM_WORLD,ierr)
                          call MPI_recv(wtd(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,status,ierr)
                          call MPI_recv(wtd(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,status,ierr)

                          call MPI_send(hz_read(1,nmax+1),1,columntype,pid+1,8,MPI_COMM_WORLD,ierr)
                          call MPI_send(hz_read(1,nmax),1,columntype,pid+1,9,MPI_COMM_WORLD,ierr)
                call MPI_recv(hz_read(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,status,ierr)
                call MPI_recv(hz_read(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,status,ierr)

            endif
        endif


                             

 
                  diff = abs(hold_read-wtd)
            maxdiff = 0
            maxdiff = maxval(diff)


                  deallocate(hz_read)
                  deallocate(diff)
                  deallocate(hz_1D)
                  deallocate(arr)
                  deallocate(T)

        endif

        call MPI_ALLREDUCE(maxdiff,diff_total,1,MPI_REAL,mpi_max,MPI_COMM_WORLD,ierr)
   
        end do SURFACE
endif




end do GROUNDWATER




if(pid.eq.0) then !do the final water table file write
    write(6,*) 'done; iterations = ',iter   
      allocate(wtdglob(n2,n3))

    open(23,file = trim(region)//'_water_test.dat',form='unformatted',access='stream')!'direct',recl=n2*n3) !do the final write - create the file
      do n=1,numtasks-1
        call MPI_recv(wtdglob(1,nini(n)),1,domblocksmall(n),n,8,MPI_COMM_WORLD,status,ierr) !receive the final wtd data from everyone
      end do

      write(23,rec=1)((wtdglob(i,j),i=1,n2),j=1,n3) !and write it to file
      close(23)
    print *,'wrote the groundwater file'
    deallocate(wtdglob)
else

      call MPI_send(wtd(1,2),1,domblocksmall(pid),0,8,MPI_COMM_WORLD,ierr) !everyone sends out the final result to pid 0! 
  endif



print *,'about to try deallocating',pid

 
  deallocate(topo,stat=error)
  if (error.ne.0) then
      print *, 'topo error'
  endif
  deallocate(wtd,stat=error)
  if (error.ne.0) then
      print *, 'wtd error'
  endif 
  deallocate(ksat,stat=error)
  if (error.ne.0) then
      print *, 'ksat error'
  endif
 deallocate(rechmean,stat=error)
  if (error.ne.0) then
      print *, 'rechmean error'
  endif
 deallocate(rech_month,stat=error)
  if (error.ne.0) then
      print *, 'rech_month error'
  endif
  deallocate(kcell,stat=error)
  if (error.ne.0) then
      print *, 'kcell error'
  endif
  deallocate(head,stat=error)
  if (error.ne.0) then
      print *, 'head error'
  endif
  deallocate(qlat,stat=error)
  if (error.ne.0) then
      print *, 'qlat error'
  endif
 deallocate(maskold,stat=error)
  if (error.ne.0) then
      print *, 'maskold error'
  endif
 deallocate(mask,stat=error)
  if (error.ne.0) then
      print *, 'mask error'
  endif
 deallocate(fdepth,stat=error)
  if (error.ne.0) then
      print *, 'fdepth error'
  endif
deallocate(landmask,stat=error)
  if (error.ne.0) then
      print *, 'landmask error'
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










print *,'done'

  call MPI_FINALIZE(ierr)
end program coupled_GW_SW



