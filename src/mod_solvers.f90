module mod_solvers

  use mod_shared

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

    double precision, dimension(SIZE(y)) :: k1,k2,k3,k4 !these are arrays

    !allocate (k1(SIZE(y))) !obviously not efficient to continually allocate/deallocate arrays
    !allocate (k2(SIZE(y)))
    !allocate (k3(SIZE(y)))
    !allocate (k4(SIZE(y)))

    k1 = f(x,y)
    k2 = f(x + (dx/2.0),y + (dx/2.0)*k1)
    k3 = f(x + (dx/2.0),y + (dx/2.0)*k2) 
    k4 = f(x + dx,y + dx*k3)

    !print *, k1
    !print *, k2
    !print *, k3
    !print *, k4


    !print *, y

    !print *, (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4)

    y = y + (dx/6.0)*(k1 + 2*k2 + 2*k3 + k4)
    x = x + dx

 !   deallocate(k1)
 !   deallocate(k2)
 !   deallocate(k3)
 !   deallocate(k4)

  end subroutine 

  subroutine householderinversion (G)

    double precision :: G(:,:)
    double precision :: P(SIZE(G(:,1)),SIZE(G(1,:)))
    double precision :: u(SIZE(G(1,:))),eorthoganal(SIZE(G(1,:)))
    double precision :: QT(SIZE(G(:,1)),SIZE(G(1,:)))
    double precision :: bufar(SIZE(G(:,1)))
    double precision :: pivot(SIZE(G(:,1)),SIZE(G(:,1)))
    !SHAPE(G)
    integer i,j

    

    !pivot G
    !very simple to get to work

    pivot = 0
    QT = 0
    forall (j = 1:SIZE(G(:,1))) QT(j,j) = 1.0
    forall (j = 1:SIZE(G(:,1))) pivot(j,j) = 1.0

    !alot of muddling of dimensions
    !also requires square - actual 
    !householder does not
    do i=1,SIZE(G(:,1))     
  
      !preform simple pivot
      if (i.ge.SIZE(G(:,1))-1) then
        if (G(SIZE(G(:,1)),SIZE(G(:,1))).eq.0 .or. G(SIZE(G(:,1))-1,SIZE(G(:,1))-1).eq.0) then
          j = SIZE(G(:,1)) - 1

          bufar=G(:,SIZE(G(:,1)))
          G(:,SIZE(G(:,1)))=G(:,j)
          G(:,j)=bufar
          
          bufar=pivot(:,SIZE(pivot(:,1)))
          pivot(:,SIZE(pivot(:,1)))=pivot(:,j)
          pivot(:,j)=bufar
        end if
      else 
        if (G(i,i).eq.0) then 
          j=0
          do
            j=j+1
            if (G(j,i).ne.0) exit !should add error for not being possible to find pivot
          end do
          
          bufar=G(:,i)
          G(:,i)=G(:,i+j)
          G(:,i+j)=bufar
          
          bufar=pivot(:,i)
          pivot(:,i)=pivot(:,i+j)
          pivot(:,i+j)=bufar
        end if
      end if

      print *, 'testdim ',G(i,:)
      print *, G(i,i:)

      eorthoganal = (/(j, j=1,SIZE(G(:,1)))/)
      where (eorthoganal.eq.i)
        eorthoganal=-1.0*SIGN(1.0D0,G(i,i))*NORM2(G(i,i:))
      elsewhere
        eorthoganal=0.0
      end where

      u = G(i,:)
      u(i:) = G(i,i:) - eorthoganal(i:)
  
      print *, 'rho: ',u(i)

      P = 0.0
      forall (j = 1:SIZE(G(:,1))) P(j,j) = 1.0

      P(i:,i:) = P(i:,i:) - 2.0*(spread(u(i:),1,SIZE(G(:,1)) + 1 - i)*spread(u(i:),2,SIZE(G(1,:)) + 1 - i)) &
     &/(DOT_PRODUCT(u(i:),u(i:)))

      print *, 'P: ',P
      print *, 'QT',QT

      !P is correct

      !QT(i:,i:) = MATMUL(P(i:,i:),QT(i:,i:)) !this is most likely currently wrong
      !QT = MATMUL(P,QT)
       QT = MATMUL(QT,P)

      !print *, 'ident ', MATMUL(QT,TRANSPOSE(QT))
      !print *, 'apply householder: ', MATMUL(G,P)

      !G(i,i:) = eorthoganal(i:)!this is wrong
      G = MATMUL(G,P)
  
      print *, 'R ', G
    end do  

    print *, 'umat: ',G

    print *, 'QT ',QT!this is too small

    !appears to be preforming rotations but by incorrect value

    print *,
    print *, 'original G: ', MATMUL(G,TRANSPOSE(QT)) !returning original G
    print *,

    !print *, pivot

    !do i=1,SIZE(QT(:,1))!this isn't pivoting, check which is supposed to be pivoted
    !  if (pivot(i).eq.0) CYCLE
    !  print *, QT(:,i)
    !  print *, QT(:,pivot(i))
    !  print *, 'pivot',i
    !  bufar=QT(:,i)
    !  QT(:,i)=QT(:,pivot(i))
    !  QT(:,pivot(i))=bufar
    !end do

    !print *, 'QT pivot ',QT

    !print *, 'test'
    !watch out for zeros
    !only works for square
    do i=SIZE(G(:,1)),1,-1 !this inversion now correct
      !this is possibly wrong dim
      !print *,
      !print *, QT(:,i+1:)
      !print *, G(i,i)
      !print *, SUM(SPREAD(G(i+1:,i),1,SIZE(G(:,1)))*QT(:,i+1:),2)
      !print *, QT(:,i)
      !print *,
      QT(:,i) = (QT(:,i) - SUM(SPREAD(G(i+1:,i),1,SIZE(G(:,1)))*QT(:,i+1:),2))/G(i,i)
      !print *, 'inverse step: ',QT
    end do

    !want to pivot after multiplication

    !print *, 'test2'

    !print *, SHAPE(G)

    
    !G = QT


    !need to pivot back
    G = MATMUL(QT,pivot)!this should be swapted?
    !do i=1,SIZE(G(:,1))
    !  if (pivot(i).eq.0) CYCLE
    !  bufar=G(:,i)
    !  G(:,i)=G(:,pivot(i))
    !  G(:,pivot(i))=bufar
    !end do
    !print *, G
    !print *, QT

    print *, 'test3'

  end subroutine

  subroutine shootingmethod(f,finnerboundary,xar,yinner,youter,unitinner,unitouter,stype_flag,write_flag)
 
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

  subroutine jacobian_solver(f,finnerboundary,boundaryconditions,xar,yparam,unitinner,unitouter,stype_flag)

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
      subroutine boundaryconditions(xar,yinner,youter,yparam)
        double precision xar(3)
        double precision yinner(:),youter(:)
        double precision yparam(:)
      end subroutine
    end interface

    integer, optional :: stype_flag
    integer :: stype_var = 0
    double precision :: yparam(:)
    double precision :: yinner(SIZE(yparam)),youter(SIZE(yparam))
    double precision :: xar(3)
  
    double precision :: yparamtmp(SIZE(yparam))
    !double precision :: delta_er(SIZE(yparam))

    double precision :: dxinner, dxouter

    integer :: unitinner,unitouter

    double precision :: jacobian(SIZE(yparam),SIZE(yparam))

    double precision :: er(SIZE(yparam))
    double precision :: dyparam(SIZE(yparam))
    double precision :: targeter = 1.0D-4

    integer :: i

    if (present(stype_flag)) then
      stype_var = stype_flag
    end if

    call boundaryconditions(xar,yinner,youter,yparam)

    do !note currently only supports INOUT

      call shootingmethod(f,finnerboundary,xar,yinner,youter,unitinner,unitouter,stype_flag,0)

      er = (youter - yinner)/yinner
  
      !main solver run with correct initial conditions
      if (MAXVAL(DABS(er)) .le. targeter) return
        !call shootingmethod(f,finnerboundary,xar,yinner,youter,unitinner,unitouter,stype_flag,1)
        !return
      !end if

      do i=1,SIZE(yparam)
        yparamtmp = yparam
        yparamtmp(i) = (1.0 + targeter)*yparam(i)
        call boundaryconditions(xar,yinner,youter,yparamtmp)
        call shootingmethod(f,finnerboundary,xar,yinner,youter,unitinner,unitouter,stype_flag,0)

        jacobian(i,:) = (youter - yinner)/yinner
      
        yparamtmp = yparam
        yparamtmp(i) = (1.0 - targeter)*yparam(i)
        call boundaryconditions(xar,yinner,youter,yparamtmp)
        call shootingmethod(f,finnerboundary,xar,yinner,youter,unitinner,unitouter,stype_flag,0)

        jacobian(i,:) = jacobian(i,:) - (youter - yinner)/yinner
        jacobian(i,:) = jacobian(i,:)/(2.0*targeter)!check the rows and columbs are the correct way round.

      end do
      
      call householderinversion(jacobian)

      dyparam = -1.0*MATMUL(jacobian,er)
      yparam = yparam + dyparam

    end do

  end subroutine

end module

!unit tests
#ifdef unittest
program solvers_unit_tests

  use mod_solvers

  implicit none

  !test matrix inversion
  double precision GMAT(3,3),GINV(3,3)

  GMAT = reshape((/1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,0.0/), shape(GMAT)) 
  GINV = (1.0/9.0)*reshape((/-16.0,8.0,-1.0,14.0,-7.0,2.0,-1.0,2.0,-1.0/), shape(GMAT))

  print *, 'Testting: householderinversion'
  print *,
  print *, 'G: ',GMAT 
  call householderinversion(GMAT)   
  print *, 'G^-1: ',GMAT
  print *, 'ACTUAL: ', GINV
  print *, 'er: ',(GMAT - GINV)/GINV
  print *,

end program
#endif
