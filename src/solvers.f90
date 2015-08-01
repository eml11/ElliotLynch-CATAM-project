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
  k3 = f(x + (dx/2.0),y + (dx/2.0)*k2) 
  k4 = f(x + dx,y + dx*k3)

  print *, k1
  print *, k2
  print *, k3
  print *, k4


  !print *, y

  print *, (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4)

  y = y + (dx/6.0)*(k1 + 2*k2 + 2*k3 + k4)
  x = x + dx

  deallocate(k1)
  deallocate(k2)
  deallocate(k3)
  deallocate(k4)

end subroutine 

subroutine shootingmethod(f,finnerboundary,xar,yinner,youter,unitinner,unitouter)
 
  implicit none
 
  interface
    function f(xin,yout)
      double precision :: xin
      double precision :: yout(:)
      double precision, dimension(size(yout)) :: f
    end function
    function finnerboundary(xin,yout)
      double precision :: xin
      double precision :: yout(:)
      double precision, dimension(size(yout)) :: finnerboundary
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
  double precision :: pyinner(SIZE(yinner)), pyouter(SIZE(youter))
  double precision :: xar(3)
  double precision :: xinner, xouter

  double precision :: dxinner, dxouter

  double precision :: maxchange = 1.0D-4

  integer :: unitinner,unitouter

  integer :: start = 1

  xinner = xar(1)
  xouter = -xar(3)

  dxinner = 1.0D-10
  dxouter = 1.0D-10

  write(unitinner,*) xinner, yinner
  write(unitouter,*) -xouter, youter

  do while (xinner.lt.xar(2) .or. xouter.lt.-xar(2))

    if (xinner.lt.xar(2)) then

      if (start.eq.1) then
  
        print *, 1
        xinner = xinner + dxinner
        print *, 2
        yinner = finnerboundary (xinner,yinner)
        print *, 3
        write(unitinner,*) xinner, yinner
        print *, 4
        start = 0
      end if
 
      pyinner = yinner

      call rk4step (f,xinner,yinner,dxinner)
      write(unitinner,*) xinner, yinner
      flush(unitinner)
! this adaptive stepping need a heck of alot of work
      !dxinner = maxchange*DABS((dxinner*pyinner(1))/((yinner(1) - pyinner(1))))
      dxinner = MINVAL(maxchange*DABS((dxinner*pyinner)/((yinner - pyinner))))

    end if 

    print *,
    print *,

    if (xouter.lt.-xar(2)) then

      pyouter = youter

      call rk4step (fouter,xouter,youter,dxouter)
      write(unitouter,*) -xouter, youter
      flush(unitouter)

      dxouter = MINVAL(maxchange*DABS(dxouter*pyouter/((youter - pyouter))))
      !dxouter = maxchange*DABS(dxouter*pyouter(1)/((youter(1) - pyouter(1))))
      !print *, dxouter

    end if

    !stop

  end do

  contains

  function fouter(x,y)
     double precision :: x
     double precision :: y(:)
     double precision, dimension(size(y)) :: fouter

     fouter = -1.0*f(-1.0*x,y)

  end function

end subroutine
