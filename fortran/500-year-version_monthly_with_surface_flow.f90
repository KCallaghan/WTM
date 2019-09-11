!In this version, I am going to adjust the model to do forward-differencing from an initial equilibrium state. 
!We will run to equilibrium at 22000 years ago, then do iterations in 500-year timesteps (6000 iterations for each run) to get actual, non-equilibrated water table. 
!This version is trying to allow lake formation for the first time, with surface flow. 
!In this version, I am treating surface flow just as another layer of groundwater. 
!Not ideal, but let's see what happens! 
!I have commented out any evaporation code for now to test one thing at a time. 

program ewtd

  use netcdf
  use mpi

  implicit none

!define all variables:

  REAL,PARAMETER :: UNDEF = -1.0E+7
  REAL(KIND=8),PARAMETER :: pi=3.141592653589793D0+0

  integer :: n2,n3,ntotal,ncid,iret,varid,n,nmax,iter,j,i,error
  integer ierr,pid,numtasks,tasktype,rc,status(MPI_STATUS_SIZE),columntype,columntypeint,iterations

  integer, allocatable :: domblock(:),domblocksmall(:),domblockint(:),nini(:),nend(:)
 
  integer*1,allocatable,dimension(:,:) :: landmask

  real deltat,dy,wtdmax,xn,xs,area,qnorth,qsouth,qeast,qwest,LH_vap
  real :: delta_xy
  real(KIND=8) :: SEDGE
  real, allocatable :: alpha(:),xlat(:)
  real, allocatable,dimension(:,:) :: topo_start_read, topo_end_read, rech_start_read, rech_end_read,&
 fslope_start_read,fslope_end_read,temp_start_read, temp_end_read,ksat_read, wtd_read,mask_read,wtd, wtdnew,&
 ksat,fdepth, wtdglob, topo_now,rech_now,fslope_now,temp_now,fdepth_now,fdepth_start,topo_sent,fdepth_sent,&
rech_sent, head,kcell, qlat_north,qlat_south,qlat_east,qlat_west!, latent_heat_start_read, latent_heat_end_read, LE_now,Evaporation, Evaporation_sent, Evaporation_start

  character*20 :: surfdatadir,initdatadir,time_start,time_end
  character*100 :: filetopo_start,filetopo_end,filerech_start,filerech_end,&
 file_fslope_start, file_fslope_end, filetemp_start, filetemp_end, fileksat,filemask,filewtd!,file_LE_start,file_LE_end


  call MPI_INIT(ierr)
  if (ierr .ne. MPI_SUCCESS) then
      print *,'Error starting MPI program. Terminating.'
      call MPI_ABORT(MPI_COMM_WORLD,rc,ierr)
  end if







iterations = 500.  !500 years * 12 months = 6000 iterations for 500-year runs
n2 = 2000  !number of columns in the rotated version of madagascar. Fortran thinks these are ROWS
n3 = 1000  !number of rows in the rotated version of madagascar. Fortran thinks these are COLUMNS.
time_start = 'Mad_021000'
time_end = 'Mad_020500'
  
surfdatadir = 'surfdata/'
initdatadir = 'initdata/'
delta_xy = 120 !there are 120 30 arc-second pieces in one degree

SEDGE = -27. !latitude of the southern edge of Madagascar. 
delta_xy = 120


!input files:

filetopo_start = trim(surfdatadir)//trim(time_start)//'_topo_rotated.nc'
filetopo_end = trim(surfdatadir)//trim(time_end)//'_topo_rotated.nc'
  
filerech_start = trim(surfdatadir)//trim(time_start)//'_rech_rotated.nc'
filerech_end = trim(surfdatadir)//trim(time_end)//'_rech_rotated.nc'

file_fslope_start = trim(surfdatadir)//trim(time_start)//'_fslope_rotated.nc'
file_fslope_end = trim(surfdatadir)//trim(time_end)//'_fslope_rotated.nc'

filetemp_start = trim(surfdatadir)//trim(time_start)//'_temp_rotated.nc'
filetemp_end = trim(surfdatadir)//trim(time_end)//'_temp_rotated.nc'

fileksat = trim(surfdatadir)//'Mad_ksat_rotated.nc' !same ksat is always used
filemask = trim(surfdatadir)//trim(time_end)//'_mask_rotated.nc' !ask Andy if he agrees with using the end time for mask
filewtd = trim(surfdatadir)//trim(time_start)//'_wtd_rotated.nc' !water table determined from the previous time step


!Set global parameters:


!input files:

  

!output text file

  open(15,file=trim(time_end)//'_output_forward_differenced.txt')


!get the number of tasks available

  call MPI_COMM_RANK(MPI_COMM_WORLD,pid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,numtasks,ierr)
  print *,'Number of tasks = ',numtasks,'My rank = ',pid

!allocate array space, the size of the number of tasks available

  allocate(nini(1:numtasks-1))
  allocate(nend(1:numtasks-1))
  allocate(domblock(1:numtasks-1))
  allocate(domblockint(1:numtasks-1))
  allocate(domblocksmall(1:numtasks-1))
               

!creating different data types and things we need to do the calculation:

  call MPI_TYPE_CONTIGUOUS(numtasks-1,MPI_INTEGER,tasktype,ierr) !creates a contiguous datatype. The new datatype is called tasktype
  call MPI_Type_commit(tasktype,ierr)


!divide the domain among tasks:

  if(pid .eq. 0) then
      write(6,*) 'PID = 0'
      call dividedomain(n2,n3,numtasks,nini,filemask,ntotal) !Divides up the work equally among all of the ranks, by number of defined land cells.
      write (6,*) 'Dividedomain done'
      write (15,*) 'Dividedomain done'
      
      do n=1,numtasks-1
          call MPI_send(nini(1),1,tasktype,n,1,MPI_COMM_WORLD,ierr) !because only PID=0 has these values right now, so we have to send them out. 
          call MPI_send(ntotal,1,MPI_INTEGER,n,2001,MPI_COMM_WORLD,ierr)
      end do
  else
      call MPI_recv(nini(1),1,tasktype,0,1,MPI_COMM_WORLD,status,ierr) !receive what was sent above.
      call MPI_recv(ntotal,1,MPI_INTEGER,0,2001,MPI_COMM_WORLD,status,ierr)

  endif    



  nend(numtasks-1) = n3 !define where each task must finish
  do n=2,numtasks-1
      nend(n-1) = nini(n) +1 !moving everything in the list along one space - the column to end in for each task.
  end do



  do n=1,numtasks-1
      nmax = nend(n) - nini(n) + 1 !max number of columns we have in this group

      call MPI_TYPE_CONTIGUOUS(n2*nmax,MPI_REAL,domblock(n),ierr)
      call MPI_type_commit(domblock(n),ierr)


      call MPI_TYPE_CONTIGUOUS(n2*nmax,MPI_INTEGER1,domblockint(n),ierr)
      call MPI_type_commit(domblockint(n),ierr)

      nmax = nmax-2
      call MPI_TYPE_CONTIGUOUS(n2*nmax,MPI_REAL,domblocksmall(n),ierr)
      call MPI_type_commit(domblocksmall(n),ierr)
  end do


  call MPI_TYPE_CONTIGUOUS(n2,MPI_REAL,columntype,ierr)
  call MPI_type_commit(columntype,ierr)


  call MPI_TYPE_CONTIGUOUS(n2,MPI_INTEGER1,columntypeint,ierr)
  call MPI_type_commit(columntypeint,ierr)



!read in all of the data:


  if(pid .eq. 0) then
      allocate(topo_start_read(n2,n3))
     
      iret = nf90_open(filetopo_start,0,ncid) !reading in the topo
      call check_err(iret)

      iret = nf90_inq_varid(ncid,'value',varid)
      call check_err(iret)

      iret = nf90_get_var(ncid,varid,topo_start_read)
      call check_err(iret)

      iret = nf90_close(ncid)
      call check_err(iret)


      where(topo_start_read .le. UNDEF) topo_start_read = 0. !change undefined cells to 0


      allocate(topo_end_read(n2,n3))
     
      iret = nf90_open(filetopo_end,0,ncid) !reading in the topo
      call check_err(iret)

      iret = nf90_inq_varid(ncid,'value',varid)
      call check_err(iret)

      iret = nf90_get_var(ncid,varid,topo_end_read)
      call check_err(iret)

      iret = nf90_close(ncid)
      call check_err(iret)


      where(topo_end_read .le. UNDEF) topo_end_read = 0. !change undefined cells to 0






      allocate(mask_read(n2,n3)) 
     
      iret = nf90_open(filemask,0,ncid) !reading in the mask
      call check_err(iret)

      iret = nf90_inq_varid(ncid,'value',varid)
      call check_err(iret)

      iret = nf90_get_var(ncid,varid,mask_read)
      call check_err(iret)

      iret = nf90_close(ncid)
      call check_err(iret)
      


      allocate(wtd_read(n2,n3))

 
      iret = nf90_open(filewtd,0,ncid) !reading in the mask
      call check_err(iret)

      iret = nf90_inq_varid(ncid,'value',varid)
      call check_err(iret)

      iret = nf90_get_var(ncid,varid,wtd_read)
      call check_err(iret)

      iret = nf90_close(ncid)
      call check_err(iret)

    !  wtd_read = -wtd_read !it is read as a positive value, but in the model it's negative
!false - it is read in as a negative value from now on! 
      
  
            
    
      allocate(ksat_read(n2,n3))
   
     
      iret = nf90_open(fileksat,0,ncid) !reading in the ksat
      call check_err(iret)

      iret = nf90_inq_varid(ncid,'value',varid)
      call check_err(iret)

      iret = nf90_get_var(ncid,varid,ksat_read)
      call check_err(iret)

      iret = nf90_close(ncid)
      call check_err(iret)





      allocate(rech_start_read(n2,n3))
      
      iret = nf90_open(filerech_start,0,ncid) !reading in the recharge file
      call check_err(iret)

      iret = nf90_inq_varid(ncid,'value',varid)
      call check_err(iret)

      iret = nf90_get_var(ncid,varid,rech_start_read)
      call check_err(iret)

      iret = nf90_close(ncid)
      call check_err(iret)

      rech_start_read = max(rech_start_read,0.)
      rech_start_read = (rech_start_read/12.) !converting to monthly



      allocate(rech_end_read(n2,n3))
      
      iret = nf90_open(filerech_end,0,ncid) !reading in the recharge file
      call check_err(iret)

      iret = nf90_inq_varid(ncid,'value',varid)
      call check_err(iret)

      iret = nf90_get_var(ncid,varid,rech_end_read)
      call check_err(iret)

      iret = nf90_close(ncid)
      call check_err(iret)

      rech_end_read = max(rech_end_read,0.)
      rech_end_read = (rech_end_read/12.)




      allocate(fslope_start_read(n2,n3)) !the idea is to load in the f values for these, i.e. f = 100/(1+150*LGM_slope) so that it doesn't still have to be calculated in here. 
      
      iret = nf90_open(file_fslope_start,0,ncid) !reading in the recharge file
      call check_err(iret)

      iret = nf90_inq_varid(ncid,'value',varid)
      call check_err(iret)

      iret = nf90_get_var(ncid,varid,fslope_start_read)
      call check_err(iret)

      iret = nf90_close(ncid)
      call check_err(iret)


      allocate(fslope_end_read(n2,n3))
      
      iret = nf90_open(file_fslope_end,0,ncid) !reading in the recharge file
      call check_err(iret)

      iret = nf90_inq_varid(ncid,'value',varid)
      call check_err(iret)

      iret = nf90_get_var(ncid,varid,fslope_end_read)
      call check_err(iret)

      iret = nf90_close(ncid)
      call check_err(iret)



      allocate(temp_start_read(n2,n3))
      
      iret = nf90_open(filetemp_start,0,ncid) !reading in the recharge file
      call check_err(iret)

      iret = nf90_inq_varid(ncid,'value',varid)
      call check_err(iret)

      iret = nf90_get_var(ncid,varid,temp_start_read)
      call check_err(iret)

      iret = nf90_close(ncid)
      call check_err(iret)


      allocate(temp_end_read(n2,n3))
      
      iret = nf90_open(filetemp_end,0,ncid) !reading in the recharge file
      call check_err(iret)

      iret = nf90_inq_varid(ncid,'value',varid)
      call check_err(iret)

      iret = nf90_get_var(ncid,varid,temp_end_read)
      call check_err(iret)

      iret = nf90_close(ncid)
      call check_err(iret)


 !     allocate(latent_heat_start_read(n2,n3))
      
 !     iret = nf90_open(file_LE_start,0,ncid) !reading in the recharge file
 !     call check_err(iret)

 !     iret = nf90_inq_varid(ncid,'value',varid)
 !     call check_err(iret)

 !     iret = nf90_get_var(ncid,varid,latent_heat_start_read)
 !     call check_err(iret)

 !     iret = nf90_close(ncid)
 !     call check_err(iret)


 !     allocate(latent_heat_end_read(n2,n3))
      
 !     iret = nf90_open(file_LE_end,0,ncid) !reading in the recharge file
 !     call check_err(iret)

 !     iret = nf90_inq_varid(ncid,'value',varid)
 !     call check_err(iret)

 !     iret = nf90_get_var(ncid,varid,latent_heat_end_read)
 !     call check_err(iret)

 !     iret = nf90_close(ncid)
 !     call check_err(iret)



!allocate(Evaporation_start(n2,n3))
 !       Evaporation_start = (latent_heat_start_read/LH_vap)*0.001*deltat !I'm getting confused with my units here but I think this is evaporation in metres per month...
!write(6,*) 'allocated Evaporation'




!now send everything we have opened:



      do n=1,numtasks-1
          call MPI_send(mask_read(1,nini(n)),1,domblockint(n),n,1,MPI_COMM_WORLD,ierr)
          call MPI_send(wtd_read(1,nini(n)),1,domblock(n),n,2,MPI_COMM_WORLD,ierr)       
          call MPI_send(ksat_read(1,nini(n)),1,domblock(n),n,3,MPI_COMM_WORLD,ierr)
          call MPI_send(wtd_read(1,nini(n)),1,domblock(n),n,4,MPI_COMM_WORLD,ierr)       
      end do

      deallocate(mask_read)
      deallocate(wtd_read)
      deallocate(ksat_read)

  else
  
      nmax = nend(pid) - nini(pid) +1
      allocate(wtd(n2,nmax))
      allocate(wtdnew(n2,nmax))
      allocate(ksat(n2,nmax))
      allocate(landmask(n2,nmax))
      allocate(head(n2,nmax))
      allocate(kcell(n2,nmax))
      allocate(qlat_north(n2,nmax))
      allocate(qlat_south(n2,nmax))
      allocate(qlat_east(n2,nmax))
      allocate(qlat_west(n2,nmax))
  

          allocate(topo_sent(n2,nmax))
          allocate(rech_sent(n2,nmax))
          allocate(fdepth_sent(n2,nmax))
!allocate(Evaporation_sent(n2,nmax))



      write(6,*)'allocated all'
      write(15,*) 'allocated all'




      call MPI_recv(landmask(1,1),1,domblockint(pid),0,1,MPI_COMM_WORLD,status,ierr)
      call MPI_recv(wtd(1,1),1,domblock(pid),0,2,MPI_COMM_WORLD,status,ierr)
      call MPI_recv(ksat(1,1),1,domblock(pid),0,3,MPI_COMM_WORLD,status,ierr)
      call MPI_recv(wtdnew(1,1),1,domblock(pid),0,4,MPI_COMM_WORLD,status,ierr)

     
      write(6,*) 'variables received'
      write(15,*) 'variables received'

  endif

!a little pre-loop math to get latitudes, total time in a time-step, and cell areas:


  deltat = (365.*24.*3600./12.) !Seconds in a monthly timestep


  dy = 6370000.*pi/(180.*delta_xy) !radius of the earth * pi / number of possible cells in the y-direction. This should equal the height of each cell in the N-S direction.

  if(pid .gt. 0) then
      nmax = nend(pid) - nini(pid) +1
      allocate(xlat(nmax))
      allocate(alpha(nmax))

      do j=1,nmax !changing area of cell depending on its latitude
          xlat(j) = (float(j+nini(pid)-2)/delta_xy+SEDGE)*pi/180. !latitude in radians
          xs = (float(2*(j+nini(pid)-2)-1)/(delta_xy*2.)+SEDGE)*pi/180. !latitude at one side of cell
          xn = (float(2*(j+nini(pid)-2)+1)/(delta_xy*2.)+SEDGE)*pi/180. !and at the other side
          area = dy*6370000.*(sin(xn)-sin(xs)) !final cell area for that latitude

          alpha(j) = 0.5*deltat/area !this is a weird way of doing it. I think I'll move this to later. 
      end do
  end if


 ! wtdmax = 0
  iter = 0


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



          do n=1,numtasks-1
              call MPI_send(topo_start_read(1,nini(n)),1,domblock(n),n,5,MPI_COMM_WORLD,ierr)
              call MPI_send(rech_start_read(1,nini(n)),1,domblock(n),n,6,MPI_COMM_WORLD,ierr)
              call MPI_send(fdepth_start(1,nini(n)),1,domblock(n),n,7,MPI_COMM_WORLD,ierr)
         !     call MPI_send(Evaporation_start(1,nini(n)),1,domblock(n),n,13,MPI_COMM_WORLD,ierr)
          end do
 
          write(6,*)'sent initial'



          deallocate(fdepth_start)

!deallocate(Evaporation_start)


      else
  
      

    !      write(6,*)'hmm'


          call MPI_recv(topo_sent(1,1),1,domblock(pid),0,5,MPI_COMM_WORLD,status,ierr) !receiving everthing that was sent above     
          call MPI_recv(rech_sent(1,1),1,domblock(pid),0,6,MPI_COMM_WORLD,status,ierr)
          call MPI_recv(fdepth_sent(1,1),1,domblock(pid),0,7,MPI_COMM_WORLD,status,ierr)
 !call MPI_recv(Evaporation_sent(1,1),1,domblock(pid),0,13,MPI_COMM_WORLD,status,ierr)

          write(6,*)'received initial'

       endif






  EWTLOOP : DO while(iter < iterations) !start of the main loop with condition - 500 years x 12 months = 6000 total iterations; 12000 for the 1000 year step
!write(6,*) 'start loop',pid

      iter = iter + 1

!if(mod(float(iter),100.).eq.0.)then
!write(6,*) 'pid',iter,pid
!endif


      if (pid .eq. 1) then !This is not completely necessary, but it's nice to see how far the run is. 
         
          write (6,*) 'iter number',iter
         
          write (15,*) 'iter number',iter
      endif


      if(mod(float(iter),100.).eq.0.)then
          if(pid.eq.0) then
              allocate(wtdglob(n2,n3))
              open(23,file=trim(time_end)//'_wtd_forward.dat',form='unformatted',access='stream')!'direct',recl=n2*n3)

              do n=1,numtasks-1
                  write(6,*)'receiving wtd',n
                  write(15,*)'receiving wtd',n
                  call MPI_recv(wtdglob(1,nini(n)+1),1,domblocksmall(n),n,8,MPI_COMM_WORLD,status,ierr)

              end do
              
              write(23,rec=1)((wtdglob(i,j),i=1,n2),j=1,n3)

              deallocate(wtdglob)
          else
              write (6,*) 'sending wtd',pid
              call MPI_send(wtd(1,2),1,domblocksmall(pid),0,8,MPI_COMM_WORLD,ierr)
              
          endif
      endif


      if(mod(float(iter),120.).eq.0.) then
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


!allocate(LE_now(n2,n3))
!          LE_now = (latent_heat_start_read * (1- iter/iterations)) + (latent_heat_end_read * (iter/iterations))
!write (6,*) 'allocated LE'






     
 
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

write (6,*) 'allocated fdepth'


!allocate(Evaporation(n2,n3))
!        Evaporation = (LE_now/LH_vap)*0.001*deltat !I'm getting confused with my units here but I think this is evaporation in metres per month...
!write(6,*) 'allocated Evaporation'



          do n=1,numtasks-1
write(6,*) 'ready to send'
              call MPI_send(topo_now(1,nini(n)),1,domblock(n),n,9,MPI_COMM_WORLD,ierr)
write (6,*) 'sent 1'
              call MPI_send(rech_now(1,nini(n)),1,domblock(n),n,10,MPI_COMM_WORLD,ierr)
write (6,*) 'sent 2'
              call MPI_send(fdepth_now(1,nini(n)),1,domblock(n),n,11,MPI_COMM_WORLD,ierr)
write (6,*) 'sent 3'
            !  call MPI_send(Evaporation(1,nini(n)),1,domblock(n),n,12,MPI_COMM_WORLD,ierr)
!write (6,*) 'sent 4'
          end do
 
          write(6,*)'sent'


          deallocate(topo_now)
          deallocate(rech_now)
          deallocate(fslope_now)
          deallocate(temp_now)
          deallocate(fdepth_now)
    !      deallocate(Evaporation)
    !      deallocate(LE_now)


      else
  
      

    !      write(6,*)'hmm'
write(6,*) 'ready to receive'


          call MPI_recv(topo_sent(1,1),1,domblock(pid),0,9,MPI_COMM_WORLD,status,ierr) !receiving everthing that was sent above     
          call MPI_recv(rech_sent(1,1),1,domblock(pid),0,10,MPI_COMM_WORLD,status,ierr)
          call MPI_recv(fdepth_sent(1,1),1,domblock(pid),0,11,MPI_COMM_WORLD,status,ierr)
 !call MPI_recv(Evaporation_sent(1,1),1,domblock(pid),0,12,MPI_COMM_WORLD,status,ierr)

          write(6,*)'received'

       endif
endif



  !   if (pid .gt. 0) wtd = min(wtd,0.) !any water table above ground gets reset to 0 - this needs to be changed for lakes





      IF (pid .gt. 0) then
          
          nmax = nend(pid) - nini(pid)
          do j=1,nmax+1
              do i=1,n2
                  if(fdepth_sent(i,j) .gt. 0. ) then !maskold is used to stop us from repeating all of the steps on cells which have already reached equilibrium. 
         !             wtd(i,j) = wtd(i,j) + rech_sent(i,j)  !adding in any recharge for that time step. 

                      head(i,j) = topo_sent(i,j) + wtd(i,j) + rech_sent(i,j) !gives the water table height (same as land surface if a wtd isn't loaded in), which = head

                      if(wtd(i,j) .lt. -1.5) then
!work out hydraulic conductivity for each cell
                          kcell(i,j) = fdepth_sent(i,j) *ksat(i,j)*exp((wtd(i,j)+1.5)/fdepth_sent(i,j)) !This is equation S6 from the paper
                      else
                          kcell(i,j) = ksat(i,j)*(wtd(i,j)+1.5+fdepth_sent(i,j)) !equation S4 from the paper 
!*******************************************************************************************************************************
    !                  else                                                      !here I am treating the surface water as another layer of groundwater
     !                     kcell(i,j) = 1*(wtd(i,j) + 1.5+fdepth_sent(i,j))                       !this is kind of just a placefiller equation to make it do something above the surface. I'm sure we can do better. 
                      endif
                  endif
              end do
          end do

     !    write(6,*) 'kcell done'

          do j=2,nmax
!write(6,*)'lets see',j
              do i=2,n2-1

                  if(landmask(i,j) .gt. 0 ) then
                     

                     ! qnorth=0.
                     ! qsouth=0.
                     ! qeast=0.
                     ! qwest=0.
                      !north
                      qnorth  = (kcell(i,j+1)+kcell(i,j))*(head(i,j+1)-head(i,j)) * cos(xlat(j)+pi/(180.*delta_xy*2.))   
!it seems like we are getting the total which will be discharged from each cell
                      !south
                      qsouth  =(kcell(i,j-1)+kcell(i,j))*(head(i,j-1)-head(i,j)) * cos(xlat(j)-pi/(180.*delta_xy*2.))
                      !west
                      qwest  = (kcell(i-1,j)+kcell(i,j))*(head(i-1,j)-head(i,j)) / cos(xlat(j))
                      !east
                      qeast  = (kcell(i+1,j)+kcell(i,j))*(head(i+1,j)-head(i,j)) / cos(xlat(j))

                     

                      qlat_north(i,j) = alpha(j)*qnorth
                      qlat_south(i,j) = alpha(j)*qsouth
                      qlat_east(i,j) = alpha(j)*qeast
                      qlat_west(i,j) = alpha(j)*qwest
                      

                      wtdnew(i,j) = wtd(i,j) + (qlat_north(i,j)+qlat_south(i,j)+qlat_east(i,j)+qlat_west(i,j))    !check all of your signs! I'm not totally sure if it should be - or + here! 
                      
                      wtdnew(i,j+1) = wtd(i,j+1) - qlat_north(i,j)
                      wtdnew(i,j-1) = wtd(i,j-1) - qlat_south(i,j)
                      wtdnew(i-1,j) = wtd(i-1,j) - qlat_west(i,j)
                      wtdnew(i+1,j) = wtd(i+1,j) - qlat_east(i,j)


               !       if(wtdnew(i,j).gt.0)then
               !           wtdnew(i,j)= max((wtdnew(i,j)-Evaporation_sent(i,j)),0.)
               !       endif





                  endif
              end do
          end do
 
 !        write(6,*)'wtdnew updated',pid

      
  !        wtd = min(wtdnew,0.)
         



!sending and receiving the lines on either side of each section to allow for flow across those lines. 

          if(pid .eq. 1) then
!write(6,*)'if'
              call MPI_send(wtd(1,nend(1)-1),1,columntype,2,0,MPI_COMM_WORLD,ierr)
     
              call MPI_recv(wtd(1,nend(1)),1,columntype,2,0,MPI_COMM_WORLD,status,ierr)
     
     

          elseif (pid .eq. numtasks-1) then     !and, continue to do this for all of the different tasks. We separate odd and even to allow our sends and receives to work without creating a blockage. 
!write(6,*)'elseif'
              if (mod(pid,2) .eq.numtasks - 1) then
                  call MPI_recv(wtd(1,1),1,columntype,pid-1,0,MPI_COMM_WORLD,status,ierr)
              
                  call MPI_send(wtd(1,2),1,columntype,pid-1,0,MPI_COMM_WORLD,ierr)
     
    

              else
                  call MPI_send(wtd(1,2),1,columntype,pid-1,0,MPI_COMM_WORLD,ierr)
    
                  call MPI_recv(wtd(1,1),1,columntype,pid-1,0,MPI_COMM_WORLD,status,ierr)
    
    
              endif
          else
!write(6,*)'else'
              nmax = nend(pid) - nini(pid)+1
              if(mod(pid,2).eq.0) then
                  call MPI_recv(wtd(1,nmax),1,columntype,pid+1,0,MPI_COMM_WORLD,status,ierr)
                  call MPI_recv(wtd(1,1),1,columntype,pid-1,0,MPI_COMM_WORLD,status,ierr)
    
    
                  call MPI_send(wtd(1,nmax-1),1,columntype,pid+1,0,MPI_COMM_WORLD,ierr)
                  call MPI_send(wtd(1,2),1,columntype,pid-1,0,MPI_COMM_WORLD,ierr)
              else
                  call MPI_send(wtd(1,nmax-1),1,columntype,pid+1,0,MPI_COMM_WORLD,ierr)
                  call MPI_send(wtd(1,2),1,columntype,pid-1,0,MPI_COMM_WORLD,ierr)
    
                  call MPI_recv(wtd(1,nmax),1,columntype,pid+1,0,MPI_COMM_WORLD,status,ierr)
                  call MPI_recv(wtd(1,1),1,columntype,pid-1,0,MPI_COMM_WORLD,status,ierr)
    

              endif
        endif





      ENDIF 









  END DO EWTLOOP









  if(pid.eq.0) then
      write(6,*) 'done; iterations = ',iter
      write(15,*) 'done; iterations = ',iter

      allocate(wtdglob(n2,n3))
      wtdglob = 0.

      open(23,file = trim(time_end)//'_wtd_forward.dat',form='unformatted',access='direct',recl=n2*n3) !do the final write - create the file

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


if(pid .gt. 0) then
deallocate(topo_sent)
deallocate(rech_sent)
deallocate(fdepth_sent)
!write(6,*)'deallocated'
endif




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
 
 !deallocate(latent_heat_start_read,stat=error)
 ! if (error.ne.0) then
 !     print *, 'latent_heat_start_read error'
 ! endif
!deallocate(latent_heat_end_read,stat=error)
 ! if (error.ne.0) then
  !    print *, 'latent_heat_end_read error'
  !endif

  deallocate(landmask,stat=error)
  if (error.ne.0)then
      print *,'landmask error'
  endif
  deallocate(wtd,stat=error)
  if (error.ne.0)then
      print *,'wtd error'
  endif
  deallocate(wtdnew,stat=error)
  if (error.ne.0) then
      print *, 'wtdnew error'
  endif
 
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
  deallocate(domblockint,stat=error)
  if (error.ne.0) then
      print *,'domblockint error'
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


  write(6,*)'done'
  write(15,*)'done'


  call MPI_FINALIZE(ierr)
end program ewtd







!********************************************************************************************************


subroutine dividedomain(n2,n3,numtasks,nini,filemask,ntotal)
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
  call check_err(iret)
  write(6,*)'first call'

  iret = nf90_inq_varid(ncid,'value',varid) !get the ID of the value layer 
  call check_err(iret)
  write(6,*) 'second call'

  iret = nf90_get_var(ncid,varid,varread) !read the actual values into the array called varread
  call check_err(iret)
  write(6,*) 'third call'

  iret = nf90_close(ncid) !close the mask file
  call check_err(iret)
  write(6,*)'fourth call'

  ntotal = count(varread>0.5) !count the number of land cells. I changed this slightly since I am using mask here rather than topo; all cells with a value of 1 should be included. 



  allocate(ncells(n3))

 
  ncells = count(varread>0.5,1) !The number of cells which are defined in the 1 dimension

  ncount=0

  nini(1) = 1

  n=2 !counter
  
  do j=1,n3
      ncount=ncount+ncells(j) !add the number of cells defined in the current column
      if (ncount .ge. ntotal/(numtasks-1)) then !>= total number of defined cells/number of threads
          nini(n) = j-1 !Telling it which column it needs to start in, to divide up the work equally, so each task handles the same number of cells.
          ncount = ncells(j) !reset the ncount for the next task
          n = n+1
      endif
      if(n .eq. numtasks) exit
  end do

  deallocate(varread,ncells)

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
