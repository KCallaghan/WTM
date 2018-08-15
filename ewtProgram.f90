program ewtd
  
  use mpi
  use Globalpars
  !use mpi
  implicit none
  !include 'mpif.h'

  integer :: n2,n3
  real :: dltxy=120 !=1 degree X resolution (30 arc second) = 120
  integer :: i,j,iun,iter,numberold,n,nmax,ntotal,number
  real dx,dy,total,deltat,d0,d1,d2,d3,wtdmax,q,xn,xs,area,thres!pi,
  real, allocatable, dimension(:,:) ::  topo,wtd,kcell,ksat,head,qlat,rechmean,fdepth,varread,wtdglob
  integer*1, allocatable, dimension(:,:) :: veg,varreadint,mask,maskold
  integer*1,allocatable, dimension(:) :: maskline
  real, allocatable :: alpha(:),xlat(:)
  real :: walltime,wtime_start,wtime1,wtime2,t1,t2
  integer ierr,rc,pid,numtasks
  integer columntype,columntypeint,tasktype,status(MPI_STATUS_SIZE)
  !integer(KIND=8) tasktype1
  integer, allocatable :: domblock(:),domblocksmall(:),domblockint(:),nini(:),nend(:)
  character*60 :: surfdatadir
  character*60 :: initdatadir
  !character*100 :: surfile
  character*4 :: region
  real :: xlatsw
  !integer :: POWER
  !integer(KIND=8) :: numbertotal
  integer :: numbertotal
  REAL(KIND=8) :: SEDGE,ULYMAP,YDIM!,ntotal

  logical :: InitIs !Where there are initial ewtd estimates to use. True or False.

  call MPI_INIT(ierr)
  if (ierr .ne. MPI_SUCCESS) then
     print *,'Error starting MPI program. Terminating.'
     call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
  end if
  
  !Global parameters that need to be set
  n2=%nnx%
  n3=%nny%
  region='%regionLW%'
  surfdatadir='%surfdatadir%'
  initdatadir='%initdatadir%'
  Thres=%threshold%
  InitIs=%InitIs%
  ULYMAP=%ULYMAP%
  YDIM=%YDIM%

  SEDGE=ULYMAP-(n3-1)*YDIM
  !WRITE(6,*) SEDGE

  d0=0.01
  d1=0.05
  d2=0.1
  d3=0.25
  IF(Thres.eq.0.01) THEN
     !WRITE(6,*) 'Doing the coarse adjustment based on the 10 mm threshold'
     !   POWER=1   
     d0=d0
  ELSEIF(Thres.eq.0.005) THEN
     !WRITE(6,*) 'Medium adjustment for 5 mm threshold'
     !   POWER=2
     d0=d0/5.
  ELSEIF(Thres.eq.0.001) THEN
     !WRITE(6,*) 'Fine adjustment for 1 mm threshold'
     !   POWER=10
     d0=d0/10.
  ELSE
     WRITE(6,*) 'The threshold input does not exist, now quit the program ...'
     STOP
  ENDIF
   
  call MPI_COMM_RANK(MPI_COMM_WORLD, pid, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
  print *, 'Number of tasks=',numtasks,' My rank=',pid
  IF(pid.eq.0) WRITE(6,*) 'Threshold = ',Thres

  !gmmdivide the domain

  allocate(nini(1:numtasks-1))
  allocate(nend(1:numtasks-1))
  allocate(domblock(numtasks-1))
  allocate(domblockint(numtasks-1))
  allocate(domblocksmall(numtasks-1))
  allocate(maskline(n2))

  call MPI_TYPE_CONTIGUOUS(numtasks-1,MPI_INTEGER,tasktype,ierr)
  call MPI_Type_commit(tasktype,ierr)

  if(pid.eq.0)then
     call dividedomain(n2,n3,numtasks,nini,surfdatadir,region,ntotal)
     !  call MPI_BCAST(ntotal,1,tasktype,0,MPI_COMM_WORLD,ierr) 
     do n=1,numtasks-1
        call MPI_send(nini(1),1,tasktype,n,1,MPI_COMM_WORLD,ierr) 
        call MPI_send(ntotal,1,MPI_INTEGER,n,2001,MPI_COMM_WORLD,ierr)    
        !     call MPI_send(number,1,tasktype,n,2002,MPI_COMM_WORLD,ierr)
     end do

  else   
     call MPI_recv(nini(1),1,tasktype,0,1,MPI_COMM_WORLD,status,ierr)
     call MPI_recv(ntotal,1,MPI_INTEGER,0,2001,MPI_COMM_WORLD,status,ierr)
     !  call MPI_send(number,1,tasktype,0,2002,MPI_COMM_WORLD,ierr)
  endif

  if(pid.eq.numtasks-1) write(6,*) 'total land cells',ntotal

  nend(numtasks-1)=n3

  do n=2,numtasks-1
     nend(n-1)=nini(n)+1
  end do

  !gmmdeclare pieces to be send and received
  do n=1,numtasks-1
     if(pid.eq.0)write(6,*)nini(n),nend(n),n,pid
     nmax=nend(n)-nini(n)+1
     call MPI_TYPE_CONTIGUOUS(n2*nmax,MPI_REAL,domblock(n),ierr)
     call MPI_Type_commit(domblock(n),ierr)

     call MPI_TYPE_CONTIGUOUS(n2*nmax,MPI_INTEGER1,domblockint(n),ierr)
     call MPI_Type_commit(domblockint(n),ierr)

     nmax=nmax-2
     call MPI_TYPE_CONTIGUOUS(n2*nmax,MPI_REAL,domblocksmall(n),ierr)
     call MPI_Type_commit(domblocksmall(n),ierr)
  end do

  call MPI_TYPE_CONTIGUOUS(n2,MPI_REAL,columntype,ierr)
  call MPI_Type_commit(columntype,ierr)

  call MPI_TYPE_CONTIGUOUS(n2,MPI_INTEGER1,columntypeint,ierr)
  call MPI_Type_commit(columntypeint,ierr)


  deltat=365.*24.*3600.

  dy=6370000.*pi/(180.*dltxy)
  dx=dy

  if(pid.gt.0)then
     nmax=nend(pid)-nini(pid)+1
     allocate(xlat(nmax))
     allocate(alpha(nmax))

     do j=1,nmax
        xlat(j)=(float(j+nini(pid)-2)/dltxy+SEDGE)*pi/180.
        xs=(float(2*(j+nini(pid)-2)-1)/(dltxy*2.)+SEDGE)*pi/180.
        xn=(float(2*(j+nini(pid)-2)+1)/(dltxy*2.)+SEDGE)*pi/180.
        area=dy*6370000.*(sin(xn)-sin(xs))
        alpha(j)=0.5*deltat/area
     end do

  endif

  wtdmax=0.

  if(pid.eq.0)then

     allocate(varread(n2,n3))

     open(23,file=trim(surfdatadir)//trim(region)//'_topo1.dat',form='unformatted',access='direct',recl=n2*n3)
     read(23,rec=1)((varread(i,j),i=1,n2),j=1,n3/4)

     open(24,file=trim(surfdatadir)//trim(region)//'_topo2.dat',form='unformatted',access='direct',recl=n2*n3)
     read(24,rec=1)((varread(i,j),i=1,n2),j=n3/4+1,2*n3/4)

     open(25,file=trim(surfdatadir)//trim(region)//'_topo3.dat',form='unformatted',access='direct',recl=n2*n3)
     read(25,rec=1)((varread(i,j),i=1,n2),j=2*n3/4+1,3*n3/4)

     open(26,file=trim(surfdatadir)//trim(region)//'_topo4.dat',form='unformatted',access='direct',recl=n2*n3)
     read(26,rec=1)((varread(i,j),i=1,n2),j=3*n3/4+1,n3)
     !sea mask
     allocate(varreadint(n2,n3))
     varreadint=1
     where(varread.le.UNDEF) varreadint=0
     where(varread.le.UNDEF) varread=0.

     do i=1,n2
        do j=1,n3
           if(varread(i,j)<-1000) write(6,*) i,j,varread(i,j)
        end do
     end do

     write(6,*)'reading topo'

     do n=1,numtasks-1
        write(6,*)'sending topo to pid ',n
        call MPI_send(varread(1,nini(n)),1,domblock(n),n,1,MPI_COMM_WORLD,ierr)
     end do

     do n=1,numtasks-1
        write(6,*)'sending veg to pid ',n
        call MPI_send(varreadint(1,nini(n)),1,domblockint(n),n,2,MPI_COMM_WORLD,ierr)
     end do

     deallocate(varreadint)

     !gmmread here initial water table
     if(InitIs) then
        write(6,*)'reading initial wtd'
        open(27,file=trim(initdatadir)//trim(region)//'_init_ewt1.dat',form='unformatted',access='direct',recl=n2*n3)
        read(27,rec=1)((varread(i,j),i=1,n2),j=1,n3/4)
        close(27)

        open(27,file=trim(initdatadir)//trim(region)//'_init_ewt2.dat',form='unformatted',access='direct',recl=n2*n3)
        read(27,rec=1)((varread(i,j),i=1,n2),j=n3/4+1,2*n3/4)
        close(27)

        open(27,file=trim(initdatadir)//trim(region)//'_init_ewt3.dat',form='unformatted',access='direct',recl=n2*n3)
        read(27,rec=1)((varread(i,j),i=1,n2),j=2*n3/4+1,3*n3/4)
        close(27)

        open(27,file=trim(initdatadir)//trim(region)//'_init_ewt4.dat',form='unformatted',access='direct',recl=n2*n3)
        read(27,rec=1)((varread(i,j),i=1,n2),j=3*n3/4+1,n3)
        close(27)
     else
 varread=0.
 write(6,*) 'Initial wtd set to zero'
     endif


     do n=1,numtasks-1
        write(6,*)'sending wtd to pid ',n
        call MPI_send(varread(1,nini(n)),1,domblock(n),n,3,MPI_COMM_WORLD,ierr)
     end do

     !gmmfdepth
     read(23,rec=3)((varread(i,j),i=1,n2),j=1,n3/4)
     read(24,rec=3)((varread(i,j),i=1,n2),j=n3/4+1,2*n3/4)
     read(25,rec=3)((varread(i,j),i=1,n2),j=2*n3/4+1,3*n3/4)
     read(26,rec=3)((varread(i,j),i=1,n2),j=3*n3/4+1,n3)


     do n=1,numtasks-1
        write(6,*)'sending fdepth to pid ',n
        call MPI_send(varread(1,nini(n)),1,domblock(n),n,4,MPI_COMM_WORLD,ierr)
     end do

     !gmmksat
     read(23,rec=2)((varread(i,j),i=1,n2),j=1,n3/4)
     read(24,rec=2)((varread(i,j),i=1,n2),j=n3/4+1,2*n3/4)
     read(25,rec=2)((varread(i,j),i=1,n2),j=2*n3/4+1,3*n3/4)
     read(26,rec=2)((varread(i,j),i=1,n2),j=3*n3/4+1,n3)

     do n=1,numtasks-1
        write(6,*)'sending ksat to pid ',n
        call MPI_send(varread(1,nini(n)),1,domblock(n),n,5,MPI_COMM_WORLD,ierr)
     end do

     close(23)
     close(24)
     close(25)
     close(26)

     !gmmrechmean

     open(23,file=trim(surfdatadir)//trim(region)//'_dollmean1.dat',form='unformatted',access='direct',recl=n2*n3)
     read(23,rec=1)((varread(i,j),i=1,n2),j=1,n3/4)

     open(24,file=trim(surfdatadir)//trim(region)//'_dollmean2.dat',form='unformatted',access='direct',recl=n2*n3)
     read(24,rec=1)((varread(i,j),i=1,n2),j=n3/4+1,2*n3/4)

     open(25,file=trim(surfdatadir)//trim(region)//'_dollmean3.dat',form='unformatted',access='direct',recl=n2*n3)
     read(25,rec=1)((varread(i,j),i=1,n2),j=2*n3/4+1,3*n3/4)

     open(26,file=trim(surfdatadir)//trim(region)//'_dollmean4.dat',form='unformatted',access='direct',recl=n2*n3)
     read(26,rec=1)((varread(i,j),i=1,n2),j=3*n3/4+1,n3)
     varread=max(varread,0.)

     do n=1,numtasks-1
        write(6,*)'sending rechmean to pid ',n
        call MPI_send(varread(1,nini(n)),1,domblock(n),n,6,MPI_COMM_WORLD,ierr)
     end do

     close(23)
     close(24)
     close(25)
     close(26)

     deallocate(varread)

  else !rest of nodes, pid !=0

     !gmmdeclare variables
     nmax=nend(pid)-nini(pid)+1
     !    write(6,*)'x dimension node',nmax,pid
     allocate(topo(n2,nmax))
     allocate(veg(n2,nmax))
     allocate(wtd(n2,nmax))
     allocate(fdepth(n2,nmax))
     allocate(ksat(n2,nmax))
     allocate(rechmean(n2,nmax))
     allocate(kcell(n2,nmax))
     allocate(head(n2,nmax))
     allocate(qlat(n2,nmax))
     allocate(mask(n2,nmax),maskold(n2,nmax))

     !gmmnodes receive the pieces
     write(6,*)'receiving variables ',pid
     call MPI_recv(topo(1,1),1,domblock(pid),0,1,MPI_COMM_WORLD,status,ierr)
     write(6,*)'topo received',pid
     call MPI_recv(veg(1,1),1,domblockint(pid),0,2,MPI_COMM_WORLD,status,ierr)
     write(6,*)'veg received',pid
     call MPI_recv(wtd(1,1),1,domblock(pid),0,3,MPI_COMM_WORLD,status,ierr)
     write(6,*)'wtd received',pid
     call MPI_recv(fdepth(1,1),1,domblock(pid),0,4,MPI_COMM_WORLD,status,ierr)
     write(6,*)'fdepth received',pid
     call MPI_recv(ksat(1,1),1,domblock(pid),0,5,MPI_COMM_WORLD,status,ierr)
     write(6,*)'ksat received',pid
     call MPI_recv(rechmean(1,1),1,domblock(pid),0,6,MPI_COMM_WORLD,status,ierr)
     !write(6,*)'rechmean received',pid
     write(6,*)'variables received',pid

     maskold=1

  endif

  !Now do the calculation.
  !now calculate the lateral flow

  if(pid.eq.0)write(6,*)'now calculating the lateral flow'

  iter=0
  number=0
  numbertotal=0
  !To make numbertotal should be less than 1% of land cells as equilibrium condition
  numbertotal=ntotal!*100D+0
  !write(6,*) numbertotal,ntotal
  if(pid.eq.0)wtime_start=walltime(0.)

  EWTLOOP : DO while(numbertotal > ntotal/100.)

     if(pid.eq.0)then
        call timing(1,t1)
        wtime1=walltime(wtime_start)
     endif

     iter=iter+1
     if(pid.eq.0)write(6,*)'iter number ',iter,numbertotal, & 
                 numbertotal*1D+0/ntotal
     !STOP
     number=0

     if(pid.gt.0) wtd = min(wtd, 0.)

     !gmmlateral flow calculation

     IF(pid.gt.0)then
             
        nmax=nend(pid)-nini(pid)

        do j=1,nmax+1
           do i=1,n2
              if(fdepth(i,j).gt.0..and.(maskold(i,j).gt.0.or.j.eq.1.or.j.eq.nmax+1))then

                 head(i,j)=topo(i,j)+wtd(i,j)
                 if(wtd(i,j).lt.-1.5)then
                    kcell(i,j)=fdepth(i,j)*ksat(i,j)*exp((wtd(i,j)+1.5)/fdepth(i,j))
                 else
                    kcell(i,j)=ksat(i,j)*(wtd(i,j)+1.5+fdepth(i,j))
                 endif

              endif

           end do
        end do

        !write(6,*)'kcell calculated'

        !head=topo+wtd

        !nmax=nend(pid)-nini(pid)

        do j=2,nmax
           do i=2,n2-1

              IF(veg(i,j).gt.0.and.maskold(i,j).gt.0) then
                 q=0.
                 !north
                 q  = q + (kcell(i,j+1)+kcell(i,j))*(head(i,j+1)-head(i,j)) * cos(xlat(j)+pi/(180.*dltxy*2.))
                 !south
                 q  = q + (kcell(i,j-1)+kcell(i,j))*(head(i,j-1)-head(i,j)) * cos(xlat(j)-pi/(180.*dltxy*2.))
                 !west
                 q  = q + (kcell(i-1,j)+kcell(i,j))*(head(i-1,j)-head(i,j)) / cos(xlat(j))
                 !east
                 q  = q + (kcell(i+1,j)+kcell(i,j))*(head(i+1,j)-head(i,j)) / cos(xlat(j))
                 qlat(i,j)=alpha(j)*q
              ENDIF

           end do
        end do

     ENDIF

     if(mod(float(iter),500.).eq.0.)then

        if(pid.eq.0)then
           allocate(wtdglob(n2,n3))
           !     wtdglob=0.

           open(23,file=trim(region)//'_dollewt1.dat',form='unformatted',access='direct',recl=n2*n3)
           open(24,file=trim(region)//'_dollewt2.dat',form='unformatted',access='direct',recl=n2*n3)
           open(25,file=trim(region)//'_dollewt3.dat',form='unformatted',access='direct',recl=n2*n3)
           open(26,file=trim(region)//'_dollewt4.dat',form='unformatted',access='direct',recl=n2*n3)

           do n=1,numtasks-1
              write(6,*)'receiving wtd ',n
              call MPI_recv(wtdglob(1,nini(n)+1),1,domblocksmall(n),n,8,MPI_COMM_WORLD,status,ierr)
           end do

           write(23,rec=1)((wtdglob(i,j),i=1,n2),j=1,n3/4)
           write(24,rec=1)((wtdglob(i,j),i=1,n2),j=n3/4+1,2*n3/4)
           write(25,rec=1)((wtdglob(i,j),i=1,n2),j=2*n3/4+1,3*n3/4)
           write(26,rec=1)((wtdglob(i,j),i=1,n2),j=3*n3/4+1,n3)

           do n=1,numtasks-1
              write(6,*)'receiving qlat ',n
              call MPI_recv(wtdglob(1,nini(n)+1),1,domblocksmall(n),n,9,MPI_COMM_WORLD,status,ierr)
           end do

           write(23,rec=2)((wtdglob(i,j),i=1,n2),j=1,n3/4)
           write(24,rec=2)((wtdglob(i,j),i=1,n2),j=n3/4+1,2*n3/4)
           write(25,rec=2)((wtdglob(i,j),i=1,n2),j=2*n3/4+1,3*n3/4)
           write(26,rec=2)((wtdglob(i,j),i=1,n2),j=3*n3/4+1,n3)

           close(23)
           close(24)
           close(25)
           close(26)

           deallocate(wtdglob)

        else
           call MPI_send(wtd(1,2),1,domblocksmall(pid),0,8,MPI_COMM_WORLD,ierr)
           call MPI_send(qlat(1,2),1,domblocksmall(pid),0,9,MPI_COMM_WORLD,ierr)
        endif

     endif

     IF(pid.gt.0)then

        mask=0

        do j=2,nmax
           do i=2,n2-1

              if(veg(i,j).gt.0.and.maskold(i,j).gt.0)then
                 !mm -> m
                 total=rechmean(i,j)*1.e-3 + qlat(i,j)
                 numberold=number
                 !       As recharge is fixed, the following applies:
                 !       (a) if total <0, meaning too much lateral flows i.e., water table is too high.
                 !       (b) if total >0, meaning too little lateral flow, i.e., water table is too low.

                 if(total .lt.-1.)then
                    !wtd(i,j)=wtd(i,j)+total/2
                    wtd(i,j)=wtd(i,j)-d3
                    !              wtd(i,j)=wtd(i,j)+total/2
                    number=number+1
                 elseif(total .lt.-.25)then
                    wtd(i,j)=wtd(i,j)-d2
                    !wtd(i,j)=wtd(i,j)+total/4  
                    !               wtd(i,j)=wtd(i,j)+total/2
                    number=number+1
                 elseif(total .lt.-0.05) then
                    wtd(i,j)=wtd(i,j)-d1
                    number=number+1
                 elseif(total .lt.-thres)then            
                    wtd(i,j)=wtd(i,j)-d0
                    number=number+1
                 elseif(total.gt.1.and.wtd(i,j).lt.wtdmax)then
                    wtd(i,j)=wtd(i,j)+d3
                    number=number+1
                 elseif(total.gt.0.25.and.wtd(i,j).lt.wtdmax)then
                    wtd(i,j)=wtd(i,j)+d2
                    number=number+1
                 elseif(total.gt.0.05.and.wtd(i,j).lt.wtdmax)then
                    wtd(i,j)=wtd(i,j)+d1
                    number=number+1
                 elseif(total.gt.thres.and.wtd(i,j).lt.wtdmax)then
                    wtd(i,j)=wtd(i,j)+d0
                    number=number+1
                 endif

                 if(numberold.ne.number)then
                    mask(i+1,j)=1
                    mask(i-1,j)=1
                    mask(i,j+1)=1
                    mask(i,j-1)=1
                    mask(i,j)=1
                 endif

              endif
              !      WRITE(6,*) 'OK here'
           end do
        end do

        maskold=mask
        numbertotal=number
       
        !1111 continue
        !gmm send each other the boundaries and add up number

        !if(numtasks.eq.1)goto 2222
        !goto 2222
            if(pid.eq.1)then
           !write(6,*)'sending border'
           call MPI_send(wtd(1,nend(1)-1),1,columntype,2,0,MPI_COMM_WORLD,ierr)
           call MPI_send(maskold(1,nend(1)),1,columntypeint,2,1,MPI_COMM_WORLD,ierr)

           call MPI_recv(wtd(1,nend(1)),1,columntype,2,0,MPI_COMM_WORLD,status,ierr)
           !write(6,*)'border received'
           call MPI_recv(maskline(1),1,columntypeint,2,1,MPI_COMM_WORLD,status,ierr)
           do i=1,n2
              maskold(i,nend(1)-1)=maskold(i,nend(1)-1)+maskline(i)
           end do

        elseif(pid.eq.numtasks-1)then
	   if(mod(pid,2).eq.0) then
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
	   
	   nmax=nend(pid)-nini(pid)+1
	   if(mod(pid,2).eq.0) then

           call MPI_recv(wtd(1,nmax),1,columntype,pid+1,0,MPI_COMM_WORLD,status,ierr)
           call MPI_recv(wtd(1,1),1,columntype,pid-1,0,MPI_COMM_WORLD,status,ierr)
           call MPI_recv(maskline(1),1,columntypeint,pid+1,1,MPI_COMM_WORLD,status,ierr)
           do i=1,n2
              maskold(i,nmax-1)=maskold(i,nmax-1)+maskline(i)
           enddo
           call MPI_recv(maskline(1),1,columntypeint,pid-1,1,MPI_COMM_WORLD,status,ierr)
           do i=1,n2
              maskold(i,2)=maskold(i,2)+maskline(i)
           enddo

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
           enddo
           call MPI_recv(maskline(1),1,columntypeint,pid-1,1,MPI_COMM_WORLD,status,ierr)
           do i=1,n2
              maskold(i,2)=maskold(i,2)+maskline(i)
           enddo
	   endif
	   
        endif

        !2222 continue
     ENDIF


    call  MPI_ALLREDUCE (number,numbertotal,1,MPI_INTEGER,mpi_sum,MPI_COMM_WORLD,ierr)
    !write(6,*) pid,'update numbertotal is good'
     if(pid.eq.0)then
        wtime2=walltime(wtime_start)
        CALL TIMING(2,T2)
        write(6,'(a25,2f10.3)')'Time,CPU(sec),Wall(sec)',T2-T1,wtime2-wtime1
     endif

     !When less than 1% land cells are not in equilibrium, exit the program.
     !if(numbertotal.lt.0.01*ntotal)EXIT
     !WRITE(6,*) numbertotal
     !write(6,*) 'done a lxxp',iter 
  END DO EWTLOOP

  !1111 continue

  if(pid.eq.0)then

     write(6,*)'done',numbertotal

     allocate(wtdglob(n2,n3))
     wtdglob=0.

     open(23,file=trim(region)//'_dollewt1.dat',form='unformatted',access='direct',recl=n2*n3)
     open(24,file=trim(region)//'_dollewt2.dat',form='unformatted',access='direct',recl=n2*n3)
     open(25,file=trim(region)//'_dollewt3.dat',form='unformatted',access='direct',recl=n2*n3)
     open(26,file=trim(region)//'_dollewt4.dat',form='unformatted',access='direct',recl=n2*n3)


     do n=1,numtasks-1
        !write(6,*)'receiving wtd ',n
        call MPI_recv(wtdglob(1,nini(n)+1),1,domblocksmall(n),n,8,MPI_COMM_WORLD,status,ierr)
     end do

     write(23,rec=1)((wtdglob(i,j),i=1,n2),j=1,n3/4)
     write(24,rec=1)((wtdglob(i,j),i=1,n2),j=n3/4+1,2*n3/4)
     write(25,rec=1)((wtdglob(i,j),i=1,n2),j=2*n3/4+1,3*n3/4)
     write(26,rec=1)((wtdglob(i,j),i=1,n2),j=3*n3/4+1,n3)

     do n=1,numtasks-1
        write(6,*)'receiving qlat',n
        call MPI_recv(wtdglob(1,nini(n)+1),1,domblocksmall(n),n,9,MPI_COMM_WORLD,status,ierr)
        write(6,*)'received block',n
     end do

     write(23,rec=2)((wtdglob(i,j),i=1,n2),j=1,n3/4)
     write(24,rec=2)((wtdglob(i,j),i=1,n2),j=n3/4+1,2*n3/4)
     write(25,rec=2)((wtdglob(i,j),i=1,n2),j=2*n3/4+1,3*n3/4)
     write(26,rec=2)((wtdglob(i,j),i=1,n2),j=3*n3/4+1,n3)

     close(23)
     close(24)
     close(25)
     close(26)


  else
     call MPI_send(wtd(1,2),1,domblocksmall(pid),0,8,MPI_COMM_WORLD,ierr)
     call MPI_send(qlat(1,2),1,domblocksmall(pid),0,9,MPI_COMM_WORLD,ierr)
  endif
deallocate(topo)
deallocate(veg)
deallocate(wtd)
deallocate(fdepth)
deallocate(ksat)
deallocate(rechmean)
deallocate(kcell)
deallocate(head)
deallocate(qlat)
deallocate(mask,maskold)
deallocate(nini)
deallocate(nend)
deallocate(domblock)
deallocate(domblockint)
deallocate(domblocksmall)
deallocate(xlat)
deallocate(alpha)
  call MPI_FINALIZE(ierr)

end program ewtd

!**********************************************************************

subroutine dividedomain(n1,n2,numtasks,nini,surfdatadir,region,ntotal)
  USE GlobalPars
  implicit none
  integer :: n1,n2,numtasks
  integer :: nini(1:numtasks-1)
  real, allocatable, dimension(:,:) :: varread
  integer, allocatable, dimension(:) :: ncells
  integer :: ncount,i,j,n,ntotal
  !integer(KIND=8) :: ntotal
  character*4 :: region
  character*60 :: surfdatadir

  allocate(varread(n1,n2))

  write(6,*)'reading soil data'
  !read in topo data to use as mask

  open(21,file=trim(surfdatadir)//trim(region)//'_topo1.dat',form='unformatted',access='direct',recl=n1*n2)
  read(21,rec=1)((varread(i,j),i=1,n1),j=1,n2/4)
  close(21)

  open(24,file=trim(surfdatadir)//trim(region)//'_topo2.dat',form='unformatted',access='direct',recl=n1*n2)
  read(24,rec=1)((varread(i,j),i=1,n1),j=n2/4+1,2*n2/4)
  close(24)

  open(25,file=trim(surfdatadir)//trim(region)//'_topo3.dat',form='unformatted',access='direct',recl=n1*n2)
  read(25,rec=1)((varread(i,j),i=1,n1),j=2*n2/4+1,3*n2/4)
  close(25)

  open(26,file=trim(surfdatadir)//trim(region)//'_topo4.dat',form='unformatted',access='direct',recl=n1*n2)
  read(26,rec=1)((varread(i,j),i=1,n1),j=3*n2/4+1,n2)
  close(26)

  ntotal=count(varread>UNDEF)

  !write(6,*)'total number of land cells',ntotal

  allocate(ncells(n2))

  ncells=count(varread>UNDEF,1)
  ncount=0

  ! nini(0)=1
  ! n=1
  nini(1)=1
  n=2
  do j=1,n2
     ncount=ncount+ncells(j)
     !   if(ncount.ge.ntotal/numtasks)then
     if(ncount.ge.ntotal/(numtasks-1))then
        !        write(6,*)ncount,ncount-ncells(j-1),j,n
        nini(n)=j-1
        ncount=ncells(j)
        n=n+1
     endif
     if(n.eq.numtasks)exit
  end do

  !ntotal=ncount

  deallocate(varread,ncells)

  return
end subroutine dividedomain
!***************************************************************************

real function walltime(wstart)

  call system_clock(count=ii,count_rate=ir)
  walltime=float(ii)/float(ir) - wstart
  return
end function walltime

!***************************************************************************
SUBROUTINE TIMING(ICALL,T1)

  !Routine returns CPU time.  Called with ICALL=1 at beginning
  !of timestep, ICALL=2 at end of timestep.

  dimension et(2)

  IF(ICALL.EQ.1) THEN
     T1=etime(et)
  ELSEIF(ICALL.EQ.2) THEN
     T1=etime(et)
  ENDIF

  RETURN
END SUBROUTINE TIMING


