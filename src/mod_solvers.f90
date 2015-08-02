module mod_solvers

  implicit none

  contains

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

    !print *, (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4)

    y = y + (dx/6.0)*(k1 + 2*k2 + 2*k3 + k4)
    x = x + dx

    deallocate(k1)
    deallocate(k2)
    deallocate(k3)
    deallocate(k4)

  end subroutine 

  subroutine householderinversion (G)

    double precision :: G(:,:)
    double precision :: P(SHAPE(G)(1),SHAPE(G)(2))
    double precision :: u(SHAPE(G)(2)),eorthoganal(SHAPE(G)(2))
    double precision :: QT(SHAPE(G)(1),SHAPE(G)(2))
    !SHAPE(G)
    integer i,j

    QT = 0
    forall (j = 1,SHAPE(G)(1)) QT(I,I) = 1.0

    !alot of muddling of dimensions
    do i=1,SHAPE(G)(1)     
  
      eorthoganal = (/(j, j=1,SHAPE(G)(1))/)
      where (eorthoganal.eq.i)
        eorthoganal=-1.0*SIGN(1.0,G(i,1))NORM2(G(i,:))
      elsewhere
        eorthoganal=0.0
      end where

      u = G(i,:) - eorthoganal
      
      P = 0.0
      forall (j = 1,SHAPE(G)(1)) P(I,I) = 1.0

      P = P - 2.0*(spread(u,1,SHAPE(G)(1))*spread(TRANSPOSE(u),2,SHAPE(G)(2)))/(DOT_PRODUCT(u,u))

      QT = MATMUL(P,QT) !this is most likely currently wrong
      G(i,:) = eorthoganal 
    end do  

    !pivots!
    !watch out for zeros
    do i=SHAPE(G)(1),1,-1!this needs to be done in a correct order
      do j=SHAPE(G)(2),1,-1
       if (i=j) then
         G(i,j) = 1
         QT(i,:) = QT(i,:)/G(i,j) 
       else
         G(i,:) = G(i,:) - G(i,j)*G(j,:)
         QT(i,:) = QT(i,:) - G(i,j)*QT(j,:)
       end if 
      end do
    end do

    G = QT

  end subroutine

  subroutine shootingmethod(f,finnerboundary,xar,yinner,youter,unitinner,unitouter,stype_flag,write_flag)
 
    implicit none

    interface
      function f(xin,yout)
        double precision :: xin
        double precision :: yout(:)
        double precision, dimension(size(yout)) :: f
      end function
    end interface 
    integer, optional :: stype_flag,write_flag
    integer :: stype_var = 0
    integer :: write_var = 1
    double precision :: yinner(:),youter(:)
    double precision :: pyinner(SIZE(yinner)), pyouter(SIZE(youter))
    double precision :: xar(3)
    double precision :: xinner, xouter
    double precision :: xc

    double precision :: dxinner, dxouter

    double precision :: maxchange = 1.0D-4

    integer :: unitinner,unitouter

    integer :: start = 1

    if (present(stype_flag)) then !is this nesicary?
      stype_var = stype_flag
    end if

    if (present(write_flag)) then
      write_var = write_flag
    end if

    xinner = xar(1)
    xouter = -xar(3)

    dxinner = 1.0D-10
    dxouter = 1.0D-10

    if (write_var.eq.1) then
      write(unitinner,*) xinner, yinner
      write(unitouter,*) -xouter, youter
    end if

    if (stype_var.eq.0) then
      xc = xar(2)
    else if (stype_var.eq.1) then
      xc = xar(3)
    else if (stype_var.eq.2) then
      xc = xar(1)
    end if

    do while (xinner.lt.xc .or. xouter.lt.-xc)

      if (xinner.lt.xc .and. (stype_var.eq.0 .or. stype_var.eq.1)) then

        if (start.eq.1) then
 
          xinner = xinner + dxinner
          yinner = finnerboundary (xinner,yinner)
          if (write_var.eq.1) then
            write(unitinner,*) xinner, yinner
          end if 
          start = 0
        end if
 
        pyinner = yinner

        call rk4step (f,xinner,yinner,dxinner)
        if (write_var.eq.1) then
          write(unitinner,*) xinner, yinner
          flush(unitinner)
        end if
        ! this adaptive stepping need a heck of alot of work
        !dxinner = maxchange*DABS((dxinner*pyinner(1))/((yinner(1) - pyinner(1))))
        dxinner = MINVAL(maxchange*DABS((dxinner*pyinner)/((yinner - pyinner))))

      end if 

      if (xouter.lt.-xc .and. (stype_var.eq.0 .or. stype_var.eq.2)) then

        pyouter = youter

        call rk4step (fouter,xouter,youter,dxouter)
        if (write_var.eq.1) then
          write(unitouter,*) -xouter, youter
          flush(unitouter)
        end if

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

  subroutine jacobian_solver(f,finnerboundary,xar,yparam,unitinner,unitouter,stype_flag)

    implicit none

    interface
      function f(xin,yout)
        double precision :: xin
        double precision :: yout(:)
        double precision, dimension(size(yout)) :: f
      end function
    end interface

    integer, optional :: stype_flag
    integer :: stype_var = 0
    double precision :: yparam(:)
    double precision :: yinner(yparam),youter(yparam)
    double precision :: pyinner(SIZE(yparam)), pyouter(SIZE(yparam))
    double precision :: xar(3)

    double precision :: yinnervar(SIZE(yparam)),youtervar(SIZE(yparam))
    double precision :: delta_er(SIZE(yparam))

    double precision :: dxinner, dxouter

    integer :: unitinner,unitouter

    double precision :: jacobian(SIZE(yparam),SIZE(yparam))

    double precision :: er(SIZE(yparam))
    double precision :: dyparam(yparam)
    double precision :: targeter = 1.0D-4

    integer :: i

    if (present(stype_flag)) then
      stype_var = stype_flag
    end if

    do
      yinner(1) = 0.0 !cannot do it this way
      youter(1) = yparam(1)

      yinner(2) = yparam(2)
      youter(2) = ((2.0*G)/(3.0*0.034))*(xar(3)/(yparam(1)**2.0))*(MSUN/(RSUN**2.0)) !from eddington approx

      yinner(3) = yparam(3)
      youter(3) = ((yparam(4)/(4.0*PI*STBOLTZ*(yparam(1)**2)))*(LSUN/(RSUN**2.0)))**(1.0/4.0)

      yinner(4) = 0.0
      youter(4) = yparam(4)

      call shootingmethod(f,finnerboundary,xar,yinner,youter,unitinner,unitouter,stype_flag,0)

      er = (youter - yinner)/yinner
  
      !main solver run with correct initial conditions
      if (MAXVAL(DABS(er)) .le. targeter) then
        call shootingmethod(f,finnerboundary,xar,yinner,youter,unitinner,unitouter,stype_flag,1)
        return
      end if

      do i=1,SIZE(yparam)
        yinnervar = yinner
        youtervar = youter
        yinnervar(i) = (1 + targeter)*yinner(i)
        youtervar(i) = (1 + targeter)*youter(i)
        call shootingmethod(f,finnerboundary,xar,yinnervar,youtervar,unitinner,unitouter,stype_flag,0)

        jacobian(i) = (youtervar - yinnervar)/yinnervar
      
        yinnervar = yinner
        youtervar = youter
        yinnervar(i) = (1 - targeter)*yinner(i)
        youtervar(i) = (1 - targeter)*youter(i)
        call shootingmethod(f,finnerboundary,xar,yinnervar,youtervar,unitinner,unitouter,stype_flag,0)

        jacobian(i) = delta_er - (youtervar - yinnervar)/yinnervar
        jacobian(i) = delta_er/(2.0*targeter)!check the rows and columbs are the correct way round.

      end do
      
      call householderinversion(jacobian)

      dyparam = -1.0*MATMUL(jacobian,er)
      yparam = yparam + dyparam

    end do

  end subroutine

end module
