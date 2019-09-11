!correct and tested way of calculating area

program area_test

implicit none

integer :: i,j,n2

real dx,dy,xn,xs

real,allocatable,dimension(:) :: area

REAL(KIND=8) :: SEDGE
REAL(KIND=8),PARAMETER :: pi=3.141592653589793D0+0
real :: dltxy = 120 !there are 120 30 arc-second pieces in one degree

SEDGE=-60
n2=17400


dy = 6370000.*pi/(180.*dltxy) !radius of the earth * pi / number of possible cells in the y-direction. This should equal the height of each cell in the N-S direction.
dx=dy
  

    allocate(area(n2))
    
    do j=1,n2  !changing area of cell depending on its latitude. Not totally sure what each of the different variables here represents...
   
        xs = (float(2*(j-1))/(dltxy*2.)+SEDGE)*pi/180. !Latitude on southern cell edge
   
!print *,xs
     xn = (float(2*(j+1))/(dltxy*2.)+SEDGE)*pi/180. !latitude on northern cell edge
        area(j) = dy * 6370000.*(sin(xn)-sin(xs))/2. !final cell area for that latitude: trapezoid dy * dx
    end do



print *,area



end program area_test
