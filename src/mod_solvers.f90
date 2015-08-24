!
! MODULE: mod_solvers 
!
!> module containing the numerical  
!! integrators and matrix solvers
module mod_solvers

  use mod_shared

  implicit none

  contains

  !> Fourth order Runge-Kutta step
  !! @param f function to be integrated
  !! @param x independent variable
  !! @param y array of dependent variables
  !! @param dx step size
  subroutine rk4step (f,x,y,dx)

    implicit none

    interface 
      function f(xin,yout)
        double precision :: xin
        double precision :: yout(:)
        double precision, dimension(size(yout)) :: f    
      end function
    end interface
    
    double precision :: y(:)
    double precision :: x, dx

    double precision, dimension(SIZE(y)) :: k1,k2,k3,k4

    k1 = f(x,y)
    k2 = f(x + (dx/2.0),y + (dx/2.0)*k1)
    k3 = f(x + (dx/2.0),y + (dx/2.0)*k2) 
    k4 = f(x + dx,y + dx*k3)

    y = y + (dx/6.0)*(k1 + 2*k2 + 2*k3 + k4)
    x = x + dx

  end subroutine 

  !> Computes matrix determinant using householder
  !! algorithm, transforms to L Matrix and computes
  !! product of the diagonal
  !! @param G input matrix
  !! @return householderdeterminant determinant of G
  function householderdeterminant(G)

    double precision :: G(:,:)
    double precision :: P(SIZE(G(:,1)),SIZE(G(1,:)))
    double precision :: u(SIZE(G(1,:))),eorthoganal(SIZE(G(1,:)))
    double precision :: QT(SIZE(G(:,1)),SIZE(G(1,:)))
    double precision :: bufar(SIZE(G(:,1)))
    double precision :: pivot(SIZE(G(:,1)),SIZE(G(:,1)))
    double precision :: UMAT(SIZE(G(:,1)),SIZE(G(1,:)))

    double precision :: householderdeterminant    

    integer i,j

    pivot = 0
    forall (j = 1:SIZE(G(:,1))) pivot(j,j) = 1.0

    UMAT = G
    householderdeterminant = 1.0
    do i=1,SIZE(G(:,1))

      !preform simple pivot
      if (i.ge.SIZE(G(:,1))-1) then
        if (UMAT(SIZE(G(:,1)),SIZE(G(:,1))).eq.0 .or. UMAT(SIZE(G(:,1))-1,SIZE(G(:,1))-1).eq.0) then
          j = SIZE(G(:,1)) - 1

          bufar=UMAT(:,SIZE(G(:,1)))
          UMAT(:,SIZE(G(:,1)))=UMAT(:,j)
          UMAT(:,j)=bufar

          bufar=pivot(:,SIZE(pivot(:,1)))
          pivot(:,SIZE(pivot(:,1)))=pivot(:,j)
          pivot(:,j)=bufar
        end if
      else
        if (UMAT(i,i).eq.0) then
          j=0
          do
            j=j+1
            if (UMAT(j,i).ne.0) exit
          end do

          bufar=UMAT(:,i)
          UMAT(:,i)=UMAT(:,i+j)
          UMAT(:,i+j)=bufar

          bufar=pivot(:,i)
          pivot(:,i)=pivot(:,i+j)
          pivot(:,i+j)=bufar
        end if
      end if

      !compute orthoganal vector
      eorthoganal = (/(j, j=1,SIZE(G(:,1)))/)
      where (eorthoganal.eq.i)
        eorthoganal=-1.0*SIGN(1.0D0,UMAT(i,i))*NORM2(UMAT(i,i:))
      elsewhere
        eorthoganal=0.0
      end where

      u = UMAT(i,:)
      u(i:) = UMAT(i,i:) - eorthoganal(i:)

      !compute rotation matrix
      P = 0.0
      forall (j = 1:SIZE(G(:,1))) P(j,j) = 1.0

      P(i:,i:) = P(i:,i:) - 2.0*(spread(u(i:),1,SIZE(G(:,1)) + 1 - i)*spread(u(i:),2,SIZE(G(1,:)) + 1 - i)) &
     &/(DOT_PRODUCT(u(i:),u(i:)))

      UMAT = MATMUL(UMAT,P)

      !Multiply diagonal element to return value
      householderdeterminant = householderdeterminant*eorthoganal(i)
    end do

  end function
 
  !> Inverts matrix using a householder reflection
  !! algorithm, transforms to Lower Matrix and solves
  !! by back substitution.
  !! @param G matrix to be inverted, currently must be square
  subroutine householderinversion (G)

    double precision :: G(:,:)
    double precision :: P(SIZE(G(:,1)),SIZE(G(1,:)))
    double precision :: u(SIZE(G(1,:))),eorthoganal(SIZE(G(1,:)))
    double precision :: QT(SIZE(G(:,1)),SIZE(G(1,:)))
    double precision :: bufar(SIZE(G(:,1)))
    double precision :: pivot(SIZE(G(:,1)),SIZE(G(:,1)))
    integer i,j

    pivot = 0
    QT = 0
    forall (j = 1:SIZE(G(:,1))) QT(j,j) = 1.0
    forall (j = 1:SIZE(G(:,1))) pivot(j,j) = 1.0

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
            if (G(j,i).ne.0) exit
          end do
          
          bufar=G(:,i)
          G(:,i)=G(:,i+j)
          G(:,i+j)=bufar
          
          bufar=pivot(:,i)
          pivot(:,i)=pivot(:,i+j)
          pivot(:,i+j)=bufar
        end if
      end if
     
      !compute orthoganal vector
      eorthoganal = (/(j, j=1,SIZE(G(:,1)))/)
      where (eorthoganal.eq.i)
        eorthoganal=-1.0*SIGN(1.0D0,G(i,i))*NORM2(G(i,i:))
      elsewhere
        eorthoganal=0.0
      end where

      u = G(i,:)
      u(i:) = G(i,i:) - eorthoganal(i:)
  
      !compute rotation matrix
      P = 0.0
      forall (j = 1:SIZE(G(:,1))) P(j,j) = 1.0

      P(i:,i:) = P(i:,i:) - 2.0*(spread(u(i:),1,SIZE(G(:,1)) + 1 - i)*spread(u(i:),2,SIZE(G(1,:)) + 1 - i)) &
     &/(DOT_PRODUCT(u(i:),u(i:)))

      !compute transpose of Q matrix
      QT = MATMUL(QT,P)

      !update R matrix computation
      G = MATMUL(G,P)
  
    end do  

    !back substitution
    do i=SIZE(G(:,1)),1,-1
      QT(:,i) = (QT(:,i) - SUM(SPREAD(G(i+1:,i),1,SIZE(G(:,1)))*QT(:,i+1:),2))/G(i,i)
    end do

    !reverse pivot
    G = MATMUL(QT,pivot)

  end subroutine

  !> Shooting method integration using 4th order Runge-Kutta
  !! with an adaptive step size.
  !! @param f function to be integrated
  !! @param finnerboundary function for initial step for inner integral
  !! @param xar integration bounds
  !! @param yinner dependent variables for inner integral
  !! @param youter dependent variables for outer integral
  !! @param unitinner file unit for inner integral output
  !! @param unitouter file unit for outer integral output
  !! @param stype_flag optional integration direction; in,out,inout
  !! @param write_flag optional write integration to file
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
    integer :: stype_var
    integer :: write_var
    double precision :: yinner(:),youter(:)
    double precision :: pyinner(SIZE(yinner)), pyouter(SIZE(youter))
    double precision :: xar(3)
    double precision :: xinner, xouter
    double precision :: xc

    double precision :: dxinner, dxouter

    double precision :: maxchange

    integer :: unitinner,unitouter

    integer :: start

    integer :: iteration !test code

    start = 1
    maxchange = 1.0D-4
    stype_var = 0
    write_var = 1 

    !optional variables
    if (present(stype_flag)) then
      stype_var = stype_flag
    end if

    if (present(write_flag)) then
      write_var = write_flag
    end if

    !setup integration
    xinner = xar(1)
    xouter = -xar(3)

    dxinner = 1.0D-10
    dxouter = 1.0D-10

    !write initial conditions
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

      !inner integral
      if (xinner.lt.xc .and. (stype_var.eq.0 .or. stype_var.eq.1)) then

        !first step - avoid divide by zero issue
        if (start.eq.1) then
        
          xinner = xinner + dxinner
          yinner = finnerboundary (xinner,yinner)
          if (write_var.eq.1) then
            write(unitinner,*) xinner, yinner
          end if 
          start = 0

        end if

        pyinner = yinner
        
        !preform rk step
        call rk4step (f,xinner,yinner,dxinner)
        if (write_var.eq.1) then
          write(unitinner,*) xinner, yinner
          flush(unitinner)
        end if

        !compute stepsize
        dxinner = MINVAL(maxchange*DABS((dxinner*pyinner)/((yinner - pyinner))))

        if (ANY(yinner .ne. yinner)) then
          print *, 'convergence fail inner : ', yinner
          stop
        end if

      end if 

      !outer integral step
      if (xouter.lt.-xc .and. (stype_var.eq.0 .or. stype_var.eq.2)) then

        pyouter = youter

        !preform rk step
        call rk4step (fouter,xouter,youter,dxouter)
        if (write_var.eq.1) then
          write(unitouter,*) -xouter, youter
          flush(unitouter)
        end if

        !compute stepsize
        dxouter = MINVAL(maxchange*DABS(dxouter*pyouter/((youter - pyouter))))

        if (ANY(youter .ne. youter)) then
          print *, 'convergence fail outer : ', youter
          stop
        end if

      end if

    end do

    contains

    !> Helper function to deal with reversal of
    !! integration direction for outer integration
    !! @param x independent variable
    !! @param y array of dependent variables
    !! @return fouter array of derivatives of y wrt x
    function fouter(x,y)
      double precision :: x
      double precision :: y(:)
      double precision, dimension(size(y)) :: fouter

      fouter = -1.0*f(-1.0*x,y)

    end function

  end subroutine

  !> Uses shooting method and jacobian iteration to refine
  !! array of input parameters to yeild desired relative error
  !! @param f function to be integrated
  !! @param finnerboundary function for initial step for inner integral
  !! @param boundaryconditions function to compute initial conditions from
  !! yparam
  !! @param xar integration bounds
  !! @param yparam array of parameters to be refined
  !! @param unitinner file unit for inner integral output
  !! @param unitouter file unit for outer integral output
  !! @param stype_flag optional integration direction; in,out,inout
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

    double precision :: dxinner, dxouter

    integer :: unitinner,unitouter

    double precision :: jacobian(SIZE(yparam),SIZE(yparam))
    double precision :: jmax(SIZE(yparam))
    integer :: r_j_num

    double precision :: er(SIZE(yparam)),forward_er(SIZE(yparam)),backward_er(SIZE(yparam))
    double precision :: dyparam(SIZE(yparam))
    double precision :: targeter = 1.0D-5

    integer :: i

    !optional variable
    if (present(stype_flag)) then
      stype_var = stype_flag
      print *, "stype_flag is set"
    end if

    do !note currently only supports INOUT

      !call shooting method to compute error
      call boundaryconditions(xar,yinner,youter,yparam)
      call shootingmethod(f,finnerboundary,xar,yinner,youter,unitinner,unitouter,stype_var,0)

      if (ANY(yinner .ne. yinner)) then
        call boundaryconditions(xar,yinner,youter,yparam)
        call shootingmethod(f,finnerboundary,xar,yinner,youter,unitinner,unitouter,stype_var,1)
        print *, 'nan in yinner ', yinner
        stop
      end if
 
      er = yinner - youter

      !exit condition
      if (ALL(DABS(er) .le. targeter)) return

      !loop over parameters in jacobian
      do i=1,SIZE(yparam)

        !forward step
        yparamtmp = yparam
        yparamtmp(i) = yparam(i) + targeter

        call boundaryconditions(xar,yinner,youter,yparamtmp)
        call shootingmethod(f,finnerboundary,xar,yinner,youter,unitinner,unitouter,stype_var,0)

        forward_er = yinner - youter

        !backward step
        yparamtmp = yparam
        yparamtmp(i) = yparam(i) - targeter 

        call boundaryconditions(xar,yinner,youter,yparamtmp)
        call shootingmethod(f,finnerboundary,xar,yinner,youter,unitinner,unitouter,stype_var,0)

        backward_er = yinner - youter
        
        jacobian(:,i) = (forward_er - backward_er)/(2.0*targeter)
 
      end do 

      !invert jacobian
      call householderinversion(jacobian)
      
      !preform parameter update
      dyparam = -1.0*MATMUL(jacobian,er)
      yparam = DABS(yparam + dyparam) 

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
  double precision xar(3),yinner(1),youter(1)

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

  print *, 

  print *, 'Testting :: shootingsolution'
  print *,

  xar(1) = 0.0D0
  xar(2) = 1.0D0
  xar(3) = 2.0D0

  yinner(1) = 1.0D0

  !test in shoot
  call shootingmethod(decay,decayinner,xar,yinner,youter,6,6,1,0)
  print *, 'test inner'
  print *, 'y(2.0) : ',yinner,'ACTUAL: ', DEXP(-0.2D0)


  youter(1) = DEXP(-0.2D0)

  !test out shoot
  call shootingmethod(decay,decayinner,xar,yinner,youter,6,6,2,0)
  print *, 'y(0.0) : ',youter,'ACTUAL: ', 1.0
  print *, 

  print *, 
  

  contains

    function decay(x,y)

      double precision :: x
      double precision :: y(:)
      double precision, dimension(SIZE(y)) :: decay

      decay(1) = -0.1*y(1)

    end function

    function decayinner(x,y)

      double precision :: x
      double precision :: y(:)
      double precision, dimension(SIZE(y)) :: decayinner
    
      decayinner(1) = y(1)

    end function

end program
#endif
