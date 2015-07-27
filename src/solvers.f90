subroutine rk4step (f,x,y,dx)

  implicit none

  interface
    function f(xin,yout)
      double precision :: xin
      double precision :: yout(:)
      double precision, dimension(size(yout)) :: f    
    end function
  end interface
  !double precision, external :: f(:)
  double precision :: y(:)
  double precision :: x, dx

  double precision, allocatable :: k1(:),k2(:),k3(:),k4(:) !these are arrays

  allocate (k1(SIZE(y))) !obviously not efficient to continually allocate/deallocate arrays
  allocate (k2(SIZE(y)))
  allocate (k3(SIZE(y)))
  allocate (k4(SIZE(y)))

  k1 = f(x,y)
  k2 = f(x + (dx/2.0),y + (dx/2.0)*k1)
  k3 = f(x + (dx/2.0),y + (dx/2.0)*k2) !appears to be failing here despite the others being secable
  k4 = f(x + dx,y + dx*k3)



  !print *, y

  !print *, (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4)

  y = y + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4)
  x = x + dx

  deallocate(k1)
  deallocate(k2)
  deallocate(k3)
  deallocate(k4)

end subroutine 

subroutine shootingmethod(f,xar,yinner,youter,unitinner,unitouter)
 
  implicit none
 
  interface
    function f(xin,yout)
      double precision :: xin
      double precision :: yout(:)
      double precision, dimension(size(yout)) :: f
    end function
    subroutine rk4step (usrf,xin,yout,dxin)
      double precision :: yout(:)
      double precision :: xin,dxin
      interface
        function usrf(xin,yout)
          double precision :: xin
          double precision :: yout(:)
          double precision, dimension(size(yout)) :: usrf
      end function
      end interface
    end subroutine
  end interface
  double precision :: yinner(:),youter(:)
  double precision :: xar(3)
  double precision :: xinner, xouter

  double precision :: dxinner, dxouter

  integer :: unitinner,unitouter

  xinner = xar(1)
  xouter = -xar(3)

  dxinner = 0.001
  dxouter = 0.001

  do while (xinner.lt.xar(2) .or. xouter.lt.-xar(2))

    if (xinner.lt.xar(2)) then
      call rk4step (f,xinner,yinner,dxinner)
      write(unitinner,*) xinner, yinner
      flush(unitinner)
    end if 

    if (xouter.lt.-xar(2)) then
      call rk4step (fouter,xouter,youter,dxouter)
      write(unitouter,*) -xouter, youter
      flush(unitouter)
    end if

    !stop

  end do

  contains

  function fouter(x,y)
     double precision :: x
     double precision :: y(:)
     double precision, dimension(size(y)) :: fouter

     fouter = -1.0*f(-x,y)

  end function

end subroutine
