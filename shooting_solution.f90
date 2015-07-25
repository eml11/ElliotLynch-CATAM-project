subroutine readinputfile (unitinput,filename)

  use mod_shared

  integer :: unitinput
  character(len=256) :: filename
  character (len=256) :: buf, label
  integer :: ios = 0
  integer :: pos = 0

  open(unit,file=filename)

  do while (ios.eq.0)

    read(unitinput, '(A)', iostat=ios) buf
    
    if (ios.eq.0) then
    
      pos = scan(buf, '    ')
      label = buf(1:pos)
      buf = buf(pos+1:)

      select case (label)
        case ('RADIUSS')

        case ('GAMMA')

        case ('MASSM')

        case('MASST')

        case('PRESSUREC')

        

      end select

    end if
  end do
  
  close(unitinput)

end subroutine

subroutine rk4step (f,x,y,dx)

  double precision, external :: f(:)
  double precision :: y(:)
  double precision :: x, dx

  double precision, allocatable :: k1(:),k2(:),k3(:),k4(:) !these are arrays

  allocate (k1,SIZE(y)) !obviously not efficient to continually allocate/deallocate arrays
  allocate (k2,SIZE(y))
  allocate (k3,SIZE(y))
  allocate (k4,SIZE(y))

  k1 = f(x,y)
  k2 = f(x + (dx/2.0),y + (dx/2.0)*k1)
  k3 = f(x + (dx/2.0),y + (dx/2.0)*k2)
  k4 = f(x + dx,y + dx*k3)

  y = y + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4)
  x = x + dx

  deallocate(k1)
  deallocate(k2)
  deallocate(k3)
  deallocate(k4)

end subroutine 

subroutine shootingmethod(f,xar,yinner,youter,unitinner,unitouter)
  
  double precision, external :: f(:)
  double precision :: yinner(:),youter
  double precision :: xar(3)
  double precision :: xinner, xouter

  double precision :: dxinner, dxouter

  integer :: unitinner,unitouter

  xinner = xar(1)
  xouter = xar(3)

  dxinner = 1000.0
  dxouter = -1000.0

  while (xinner.lt.xar(2)).or.(xouter.gt.xar(2)) do

    if (xinner.lt.xar(2))
      rk4step (f,xinner,yinner,dxinner)
      write(*,unitinner) xinner, yinner
    end if 

    if (xouter.gt.xar(2))
      rk4step (f,xouter,youter,dxouter)
      write(*,unitouter) xinner, youter
    end if

  end do

end subroutine

program stellarstructure

  use mod_shared

  implicit none

  double precision :: xar(3)
  double precision :: yinner(2), youter(2)
  integer :: unitinner = 101,unitouter = 102
  integer :: unitinput = 103

  call readinputfile (unitinput,'star.inp')

  open(unitinner,'stellarstructure_inner.txt')
  open(unitouter,'stellarstructure_outer.txt')

  !set up problem
  xar(1) = 0
  xar(2) = massm
  xar(3) = masst

  yinner(1) = 0.0
  youter(1) = radiusstar**3.

  yinner(2) = DLOG(pressurec)
  youter(2) = 

  !run solver


  close(unitinner)
  close(unitouter)

  contains

  function question4rkfunc (x,y)

    double precision :: x
    double precision :: y(2)
    double precision :: question4rkfunc(2)

    question4rkfunc(1) = 12*PI*(y(1)**(4.0/3.0))*densitycenter*DEXP((y(2) - y2center)/ratiospecificheat)
    question4rkfunc(2) = - ((G*x)/(4*PI))*(y(1)**(4.0/3.0))*DEXP(y(2))

  end function

end program

