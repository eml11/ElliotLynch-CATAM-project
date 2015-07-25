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
