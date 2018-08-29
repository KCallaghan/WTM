!This is for obtaining an equilibrium depth to water table under given climatic and topographic conditions
!You are assuming the water table is at equilibrium; for non-equilibrium water tables, you will need a starting
!depth to water table and the transient version of this model. 
!Note this version models groundwater ONLY. Any water above ground surface/exfiltrating from the ground is assumed to run off/evaporate and leave the land surface. 

!Required inputs are topography (metres), estimates of hydraulic conductivity, recharge (precipitation - evapotranspiration) (metres), 
!slope, and winter temperature, in netcdf (.nc) format. 
!Output is depth to water table in metres. 

!For compilation, this requires netcdf and mpi libraries. 


!declare variables ******************************************************************************************************************************

program ewtd

  use netcdf
  use mpi

  implicit none !ensures you don't create any variables without explicitly defining their type

!define all variables:

  REAL,PARAMETER :: UNDEF = -1.0E+7
  REAL(KIND=8),PARAMETER :: pi=3.141592653589793D0+0

  integer :: n2,n3,ntotal,iret,ncid,varid,n,j,nmax,i,numbertotal,numbercount,iter,error,numberold
  integer ierr,pid,numtasks,rc,tasktype,status(MPI_STATUS_SIZE),columntype,columntypeint,Wtd_input

  real :: dltxy = 120 !there are 120 30 arc-second pieces in one degree
  real thres,d0,d1,d2,d3,deltat,dx,dy,xn,xs,area,wtdmax,q,total

  character*60 :: surfdatadir,initdatadir,output_string
  character*10 :: region
  character*100 :: filetopo,filewtd,fileslope,fileksat,filerech,filemask,filetemp

  REAL(KIND=8) :: SEDGE !kind defines precision, number of bytes

  integer,allocatable :: domblock(:),domblocksmall(:),domblockint(:),nini(:),nend(:)
  integer*1,allocatable,dimension(:) :: maskline
  integer*1, allocatable,dimension(:,:) :: landmask,mask,maskold,watermask

  real, allocatable :: alpha(:),xlat(:),alphamonth(:)
  real, allocatable, dimension(:,:) :: topo,wtd,kcell,ksat,head,qlat,rechmean,fdepth,varread,wtd_read,fdepth_read,ksat_read,&
  rech_read,rech_month_read,rech_month,wtdglob,fdepth_start,temp_read,slope_read,temp_sent

   

 
  call MPI_INIT(ierr)       !Initialise MPI
  if (ierr .ne. MPI_SUCCESS) then       !Error catching - this part should never run.
      print *,'Error starting MPI program. Terminating.'
      call MPI_ABORT(MPI_COMM_WORLD,rc,ierr)
  end if


!************************************************************************************************************************************************

!All of your global variables need to be set here. Change these values as needed for your application:

  n3 = 17400 !number of cells in x and y directions
  n2 = 43200
  region = '000000' !first few characters in name for the 4 files. Limit 10 characters. E.g. time period '000500' or location 'N_America'
  surfdatadir = 'surfdata/' !folder names for data, limit 60 characters. All of your input data files go here.
  initdatadir = 'initdata/' !mask and wtd go here. wtd = water table depth. This is optional, and only used if you have starting estimates to make code faster. 
!ksat is not within a folder since the same one is being used for all times.


  Thres = 0.001 !When to stop calculating - can try with different values and see how it goes. This is the finest option (choose 0.01, 0.005 or 0.001). 
  SEDGE = -60 !southern most latitude in degrees, used for calculating changing cell areas

  Wtd_input = 0 !0 means no starting estimate, wtd will start as 0 everywhere. Change to 1 if you would like to use a starting estimated depth to water table file input. 

!input files:

  filetopo  = trim(surfdatadir)//trim(region)//'_topo.nc'     !topography file
  fileslope = trim(surfdatadir)//trim(region)//'_slope.nc'    !slope file
  fileksat  = trim(surfdatadir)//trim(region)//'_ksat.nc'     !estimated hydraulic conductivity
  filerech  = trim(surfdatadir)//trim(region)//'_rech.nc'     !recharge as precipitation - evapotranspiration
  filemask  = trim(surfdatadir)//trim(region)//'_mask.nc'     !mask file, 1 for land and 0 for ocean
  filewtd   = trim(surfdatadir)//trim(region)//'_wtd.nc'      !OPTIONAL depth to water table. 
  filetemp  = trim(surfdatadir)//trim(region)//'_temp.nc'     !WINTER temperature in degrees Celsius


  output_string = '_text_output.txt'                          !name of your output text file

!End of variable definitions ***********************************************************************************************************************************************

!output text file

open (15,file=trim(region)//trim(output_string)) !creating a text file to store the values of iterations etc

!different increments to use to move water table depth:

  d0 = 0.005 
  d1 = 0.02 
  d2 = 0.1
  d3 = 0.25

  IF (Thres .eq. 0.01) THEN
      WRITE(15,*) 'Doing the coarse adjustment based on the 10 mm threshold'
      d0=d0
  ELSEIF (Thres .eq. 0.005) THEN
      WRITE(15,*) 'Medium adjustment for the 5 mm threshold'
      d0=d0/5.
  ELSEIF (Thres .eq. 0.001) THEN
      WRITE(15,*) 'Medium adjustment for the 5 mm threshold'
      d0=d0/10.
  ELSE
      WRITE(15,*) 'The threshold input does not exist, the program will now quit.'
      STOP
  ENDIF

!Behind-the-scenes MPI work. No actual processing to see here! ****************************************************************************************************

 
!get the number of tasks available
  call MPI_COMM_RANK(MPI_COMM_WORLD,pid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,numtasks,ierr)
  print *,'Number of tasks=',numtasks,'My rank=',pid
 
!allocate array space, the size of the number of tasks available

  allocate(nini(1:numtasks-1))
  allocate(nend(1:numtasks-1))
  allocate(domblock(1:numtasks-1))
  allocate(domblockint(1:numtasks-1))
  allocate(domblocksmall(1:numtasks-1))
  allocate(maskline(n2))               

  call MPI_TYPE_CONTIGUOUS(numtasks-1,MPI_INTEGER,tasktype,ierr) !creates a contiguous datatype. The new datatype is called tasktype
  call MPI_Type_commit(tasktype,ierr)


!divide the domain among tasks:

  if(pid .eq. 0) then
      write(15,*) 'PID = 0'
      call dividedomain(n2,n3,numtasks,nini,filemask,ntotal) !Divides up the work equally among all of the ranks, by number of defined land cells.
      write (15,*) 'Dividedomain done'
      
      do n=1,numtasks-1
          call MPI_send(nini(1),1,tasktype,n,1,MPI_COMM_WORLD,ierr) !because only PID=0 has these values right now, so we have to send them out. 
          call MPI_send(ntotal,1,MPI_INTEGER,n,2001,MPI_COMM_WORLD,ierr)
      end do
  else
      call MPI_recv(nini(1),1,tasktype,0,1,MPI_COMM_WORLD,status,ierr) !receive what was sent above.
      call MPI_recv(ntotal,1,MPI_INTEGER,0,2001,MPI_COMM_WORLD,status,ierr)

  endif   
write(15,*)'sent here'

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

  !Done with MPI setup. Calculation for changing cell area with latitude: 


  deltat = (365.*24.*3600.) !Seconds in an annual timestep

  dy = 6370000.*pi/(180.*dltxy) !radius of the earth * pi / number of possible cells in the y-direction. This should equal the height of each cell in the N-S direction.
  dx=dy
  
  if(pid .gt. 0) then
      nmax = nend(pid) - nini(pid) +1
      allocate(xlat(nmax))
      allocate(alpha(nmax))
      allocate(alphamonth(nmax))

      do j=1,nmax  !changing area of cell depending on its latitude. 
          xlat(j) = (float(j+nini(pid)-2)/dltxy+SEDGE)*pi/180.
          xs = (float(2*(j+nini(pid)-2)-1)/(dltxy*2.)+SEDGE)*pi/180. !Latitude number 1
          xn = (float(2*(j+nini(pid)-2)+1)/(dltxy*2.)+SEDGE)*pi/180. !latitude number 2
          area = dy*6370000.*(sin(xn)-sin(xs)) !final cell area for that latitude

          alpha(j) = 0.5*deltat/area 
          alphamonth(j) = 0.5*(deltat/12)/area 
      end do
  end if



  wtdmax = 0 !Max water table depth is 0. This version of code is for groundwater ONLY. 


!Import data files ********************************************************************************************************************************************************  

  if(pid .eq. 0) then
      allocate(varread(n2,n3))
     
      iret = nf90_open(filetopo,0,ncid) !reading in the topography
      call check_err(iret)

      iret = nf90_inq_varid(ncid,'value',varid)
      call check_err(iret)

      iret = nf90_get_var(ncid,varid,varread)
      call check_err(iret)

      iret = nf90_close(ncid)
      call check_err(iret)


      where(varread .le. UNDEF) varread = 0. !change undefined cells to 0


      allocate(temp_read(n2,n3))
    
      iret = nf90_open(filetemp,0,ncid) !reading in the temperature
      call check_err(iret)

      iret = nf90_inq_varid(ncid,'value',varid)
      call check_err(iret)

      iret = nf90_get_var(ncid,varid,temp_read)
      call check_err(iret)

      iret = nf90_close(ncid)
      call check_err(iret)


      allocate(watermask(n2,n3)) 
     
      iret = nf90_open(filemask,0,ncid) !reading in the mask
      call check_err(iret)

      iret = nf90_inq_varid(ncid,'value',varid)
      call check_err(iret)

      iret = nf90_get_var(ncid,varid,watermask)
      call check_err(iret)

      iret = nf90_close(ncid)
      call check_err(iret)
      


      allocate(wtd_read(n2,n3))
      if(wtd_input .eq. 0)then
          wtd_read = 0           !no longer using an input water table! 
      else
          iret = nf90_open(filewtd,0,ncid) !reading in the wtd
          call check_err(iret)

          iret = nf90_inq_varid(ncid,'value',varid)
          call check_err(iret)
    
          iret = nf90_get_var(ncid,varid,wtd_read)
          call check_err(iret)
    
              iret = nf90_close(ncid)
          call check_err(iret)
      endif
     
    
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
      
     
      iret = nf90_open(filerech,0,ncid) !reading in the recharge file
      call check_err(iret)

      iret = nf90_inq_varid(ncid,'value',varid)
      call check_err(iret)

      iret = nf90_get_var(ncid,varid,rech_read)
      call check_err(iret)

      iret = nf90_close(ncid)
      call check_err(iret)


      rech_read = max(rech_read,0.) !setting it so recharge can only be positive

      where(rech_read .ge. 10000) rech_read = 0. 


      allocate(rech_month_read(n2,n3))

      rech_month_read = (rech_read/12.)  



!now send everything we have opened:

      do n=1,numtasks-1
          call MPI_send(varread(1,nini(n)),1,domblock(n),n,1,MPI_COMM_WORLD,ierr)
          call MPI_send(watermask(1,nini(n)),1,domblockint(n),n,2,MPI_COMM_WORLD,ierr)
          call MPI_send(wtd_read(1,nini(n)),1,domblock(n),n,3,MPI_COMM_WORLD,ierr)
          call MPI_send(ksat_read(1,nini(n)),1,domblock(n),n,5,MPI_COMM_WORLD,ierr)
          call MPI_send(rech_read(1,nini(n)),1,domblock(n),n,6,MPI_COMM_WORLD,ierr)
          call MPI_send(rech_month_read(1,nini(n)),1,domblock(n),n,7,MPI_COMM_WORLD,ierr)
      end do


      deallocate(varread)
      deallocate(ksat_read)
      deallocate(wtd_read)
      deallocate(rech_read)
      deallocate(watermask)
      deallocate(rech_month_read)


  else
  
      nmax = nend(pid) - nini(pid) +1
      allocate(topo(n2,nmax))
      allocate(landmask(n2,nmax))
      allocate(wtd(n2,nmax))
      allocate(ksat(n2,nmax))
      allocate(rechmean(n2,nmax))
      allocate(rech_month(n2,nmax))
      allocate(kcell(n2,nmax))
      allocate(head(n2,nmax))
      allocate(qlat(n2,nmax))
      allocate(mask(n2,nmax),maskold(n2,nmax))
      write(15,*) 'allocated all'


      call MPI_recv(topo(1,1),1,domblock(pid),0,1,MPI_COMM_WORLD,status,ierr) !receiving everthing that was sent above     
      call MPI_recv(landmask(1,1),1,domblockint(pid),0,2,MPI_COMM_WORLD,status,ierr)
      call MPI_recv(wtd(1,1),1,domblock(pid),0,3,MPI_COMM_WORLD,status,ierr)
      call MPI_recv(ksat(1,1),1,domblock(pid),0,5,MPI_COMM_WORLD,status,ierr)
      call MPI_recv(rechmean(1,1),1,domblock(pid),0,6,MPI_COMM_WORLD,status,ierr)
      call MPI_recv(rech_month(1,1),1,domblock(pid),0,7,MPI_COMM_WORLD,status,ierr)
 
      maskold = 1 !This is just being done ahead of the big loop, maskold is used to show which cells still need to be processed. 

  endif


  if (pid.eq.0) then
      allocate(slope_read(n2,n3))
      
      iret = nf90_open(fileslope,0,ncid) !reading in the slope
      call check_err(iret)
  
      iret = nf90_inq_varid(ncid,'value',varid)
      call check_err(iret)
  
      iret = nf90_get_var(ncid,varid,slope_read)
      call check_err(iret)
  
      iret = nf90_close(ncid)
      call check_err(iret)
  
      allocate(fdepth_start(n2,n3))
  
      write(15,*)'fdepth_Start allocated'
   
      do j=1,n3
          do i=1,n2
              if (temp_read(i,j) .gt. -5) then
                  fdepth_start(i,j) = 100/(1+150*slope_read(i,j))
              elseif (temp_read(i,j) .lt. -14) then
                  fdepth_start(i,j) =  (100/(1+150*slope_read(i,j))) * (0.17+0.005*temp_read(i,j))
              else
                  fdepth_start(i,j) = (100/(1+150*slope_read(i,j)))  * (1.5 + 0.1*temp_read(i,j))
              endif
          end do
      end do
      write(15,*) 'fdepth_start created'
  

      where(fdepth_start .le. 0.0001) fdepth_start = 0.0001 !change undefined cells to 0

    
      do n=1,numtasks-1
          write (15,*) 'trying to send'
             
          call MPI_send(fdepth_start(1,nini(n)),1,domblock(n),n,15,MPI_COMM_WORLD,ierr)
             
      end do
 
      write(15,*)'sent initial'

      deallocate(fdepth_start)

  else
  
      nmax = nend(pid) - nini(pid) +1

      allocate(fdepth(n2,nmax))
       
      call MPI_recv(fdepth(1,1),1,domblock(pid),0,15,MPI_COMM_WORLD,status,ierr)
 
      write(15,*)'received initial'

  endif





 
  iter = 0              !count the number of iterations - this is used to know when to switch from annual to monthly cycles. 
  numbercount = 0
  numbertotal = ntotal !we will use this to get numbertotal less than x % land cells as equilibrium condition
 

 
  EWTLOOP : DO while(numbertotal > ntotal/100 .and. iter<400000 ) !start of the main loop with condition - either number of iterations or % equilibrium. I am going for 99% equilibrium. 
      
 

      if (iter .eq. 30000) then           !I have switched to 30000 iterations before switching to monthly processing, because with 50000 it always seemed to be doing nothing for a long time. 
!here we automatically switch to monthly processing, 
          write(15,*) '30000 iterations,adjusting the values for monthly processing'
          thres = thres/12.
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

      iter = iter + 1 !count the iterations

      if (pid .eq. 0) then !This is not completely necessary, but it's nice to see how far the run is. 
          write (15,*) 'PID = 0. Numbertotal:',numbertotal,'ntotal',ntotal
          write (15,*) 'iter number',iter
      endif

      numbercount = 0


      if(mod(float(iter),10000.).eq.0.)then       !save the data
          if(pid.eq.0) then
              allocate(wtdglob(n2,n3))
              open(23,file=trim(region)//'dollewt1.dat',form='unformatted',access='direct',recl=n2*n3)

              do n=1,numtasks-1
                  write(15,*)'receiving wtd',n
                  call MPI_recv(wtdglob(1,nini(n)+1),1,domblocksmall(n),n,8,MPI_COMM_WORLD,status,ierr)

              end do
              
              write(23,rec=1)((wtdglob(i,j),i=1,n2),j=1,n3)

              deallocate(wtdglob)
          else
              call MPI_send(wtd(1,2),1,domblocksmall(pid),0,8,MPI_COMM_WORLD,ierr)
          endif
      endif


      if (pid .gt. 0) wtd = min(wtd,0.) !any water table above ground gets reset to 0


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
                      else
                          kcell(i,j) = ksat(i,j)*(wtd(i,j)+1.5+fdepth(i,j)) !equation S4 from the paper 
                      endif
                  endif
              end do
          end do



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
                     
                      total = rechmean(i,j) + qlat(i,j)
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
                      elseif (total .lt. -0.05) then
                          wtd(i,j) = wtd(i,j) -d1 
                          numbercount = numbercount + 1
                          
                      elseif (total .lt. -thres) then
                          wtd(i,j) = wtd(i,j) -d0
                          numbercount = numbercount + 1

                      elseif(total .gt. 1. .and.wtd(i,j).lt.wtdmax) then                 
                          wtd(i,j) = wtd(i,j) +d3
                          numbercount = numbercount + 1
                      elseif (total .gt. 0.25 .and.wtd(i,j).lt.wtdmax) then
                          wtd(i,j) = wtd(i,j) +d2
                          numbercount = numbercount + 1
                      elseif (total .gt. 0.05 .and.wtd(i,j).lt.wtdmax) then
                          wtd(i,j) = wtd(i,j) +d1 
                          numbercount = numbercount + 1
                      elseif (total .gt. thres .and.wtd(i,j).lt.wtdmax) then
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


      
         maskold = mask !to record which cells still need to be processed.

         numbertotal = numbercount !total number of cells still out of equilibrium this iteration
                             


!sending and receiving the lines on either side of each section to allow for flow across those lines. 

          if(pid .eq. 1) then
              call MPI_send(wtd(1,nend(1)-1),1,columntype,2,0,MPI_COMM_WORLD,ierr)
              call MPI_send(maskold(1,nend(1)),1,columntypeint,2,1,MPI_COMM_WORLD,ierr) !so we send out the border part of our wtd and our mask

              call MPI_recv(wtd(1,nend(1)),1,columntype,2,0,MPI_COMM_WORLD,status,ierr)
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

 

  END DO EWTLOOP !done with the big piece of processing! We now have 98% equilibrium! 


  if(pid.eq.0) then
      write(15,*) 'done; numbertotal = ',numbertotal,'iterations = ',iter

      allocate(wtdglob(n2,n3))
      wtdglob = 0.

      open(23,file = trim(region)//'_dollewt1.dat',form='unformatted',access='direct',recl=n2*n3) !do the final write - create the file

      do n=1,numtasks-1
          call MPI_recv(wtdglob(1,nini(n)+1),1,domblocksmall(n),n,8,MPI_COMM_WORLD,status,ierr) !receive the final wtd data from everyone
      end do

      write(23,rec=1)((wtdglob(i,j),i=1,n2),j=1,n3) !and write it to file

      close(23)

  else
      call MPI_send(wtd(1,2),1,domblocksmall(pid),0,8,MPI_COMM_WORLD,ierr) !everyone sends out the final result to pid 0! 

  endif


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
  deallocate(slope_read,stat=error)
  if (error.ne.0) then
      print *,'slope error'
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

 deallocate(temp_sent,stat=error)
  if (error.ne.0) then
      print *, 'temp error'
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

  write(15,*)'reading in the mask to divide the domain'

  iret = nf90_open(filemask,0,ncid)  !open the mask file
  call check_err(iret)
  write(15,*)'first call'

  iret = nf90_inq_varid(ncid,'value',varid) !get the ID of the value layer
  call check_err(iret)
  write(15,*) 'second call'

  iret = nf90_get_var(ncid,varid,varread) !read the actual values into the array called varread
  call check_err(iret)
  write(15,*) 'third call'

  iret = nf90_close(ncid) !close the mask file
  call check_err(iret)
  write(15,*)'fourth call'

  ntotal = count(varread>0.5) !count the number of land cells. I changed this slightly since I am using mask here rather than topo; all cells with a value of 1 should be included. 



  allocate(ncells(n3))

 
  ncells = count(varread>0.5,1) !The number of cells which are defined in the 1 dimension

  ncount=0

  nini(1) = 1

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




