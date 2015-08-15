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
    !print *, 'k1'
    k2 = f(x + (dx/2.0),y + (dx/2.0)*k1)
    !print *, 'k2'
    k3 = f(x + (dx/2.0),y + (dx/2.0)*k2) 
    !print *, 'k3'
    k4 = f(x + dx,y + dx*k3)
    !print *, 'k4'

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

      !print *, 'testdim ',G(i,:)
      !print *, G(i,i:)

      eorthoganal = (/(j, j=1,SIZE(G(:,1)))/)
      where (eorthoganal.eq.i)
        eorthoganal=-1.0*SIGN(1.0D0,G(i,i))*NORM2(G(i,i:))
      elsewhere
        eorthoganal=0.0
      end where

      u = G(i,:)
      u(i:) = G(i,i:) - eorthoganal(i:)
  
      !print *, 'rho: ',u(i)

      P = 0.0
      forall (j = 1:SIZE(G(:,1))) P(j,j) = 1.0

      P(i:,i:) = P(i:,i:) - 2.0*(spread(u(i:),1,SIZE(G(:,1)) + 1 - i)*spread(u(i:),2,SIZE(G(1,:)) + 1 - i)) &
     &/(DOT_PRODUCT(u(i:),u(i:)))

      !print *, 'P: ',P
      !print *, 'QT',QT

      !P is correct

      !QT(i:,i:) = MATMUL(P(i:,i:),QT(i:,i:)) !this is most likely currently wrong
      !QT = MATMUL(P,QT)
       QT = MATMUL(QT,P)

      !print *, 'ident ', MATMUL(QT,TRANSPOSE(QT))
      !print *, 'apply householder: ', MATMUL(G,P)

      !G(i,i:) = eorthoganal(i:)!this is wrong
      G = MATMUL(G,P)
  
      !print *, 'R ', G
    end do  

    !print *, 'umat: ',G

    !print *, 'QT ',QT!this is too small

    !appears to be preforming rotations but by incorrect value

    !print *,
    !print *, 'original G: ', MATMUL(G,TRANSPOSE(QT)) !returning original G
    !print *,

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

    !print *, 'test3'

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

    integer :: start != 1

    integer :: iteration !test code


    start = 1
    iteration = 1 !test code
 

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

    !print *, 'input ',xar,yinner

    do while (xinner.lt.xc .or. xouter.lt.-xc)



      if (xinner.lt.xc .and. (stype_var.eq.0 .or. stype_var.eq.1)) then

        !this is most likely unstable

        if (start.eq.1) then
 
          xinner = xinner + dxinner
          yinner = finnerboundary (xinner,yinner)
          if (write_var.eq.1) then
            write(unitinner,*) xinner, yinner
          end if 
   
          !random test fix
          !if (yinner(1) .ge. 0.1) then
          start = 0
          !stop
          !end if
          !print *, 'start ',start
        end if

        !print *, yinner

        pyinner = yinner

     !   print *, "iteration - inner loop",iteration 

        call rk4step (f,xinner,yinner,dxinner)
        if (write_var.eq.1) then
          write(unitinner,*) xinner, yinner
          flush(unitinner)
        end if
        ! this adaptive stepping need a heck of alot of work
        !dxinner = maxchange*DABS((dxinner*pyinner(1))/((yinner(1) - pyinner(1))))
        dxinner = MINVAL(maxchange*DABS((dxinner*pyinner)/((yinner - pyinner))))

      end if 

     ! stop

      if (xouter.lt.-xc .and. (stype_var.eq.0 .or. stype_var.eq.2)) then

        pyouter = youter

      !  print *, "iteration - outer loop",iteration
        call rk4step (fouter,xouter,youter,dxouter)
        if (write_var.eq.1) then
          write(unitouter,*) -xouter, youter
          flush(unitouter)
        end if

        dxouter = MINVAL(maxchange*DABS(dxouter*pyouter/((youter - pyouter))))
        !dxouter = maxchange*DABS(dxouter*pyouter(1)/((youter(1) - pyouter(1))))
        !print *, dxouter

      end if

      iteration = iteration + 1 !test code
      !print *,

      !call sleep(5)

    end do

    !print *, '\n\n\n#######################\n\n\n#################\n\n\n'


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
    double precision :: rescale(SIZE(yparam),SIZE(yparam))
    !double precision :: delta_er(SIZE(yparam))

    double precision :: dxinner, dxouter

    integer :: unitinner,unitouter

    double precision :: jacobian(SIZE(yparam),SIZE(yparam))
    double precision :: jmax(SIZE(yparam))
    double precision :: regulerised_jacobian(SIZE(yparam),SIZE(yparam))
    integer :: r_j_num

    double precision :: er(SIZE(yparam))
    double precision :: dyparam(SIZE(yparam))
    double precision :: targeter = 1.0D-5

    integer :: i

    if (present(stype_flag)) then
      stype_var = stype_flag
      print *, "stype_flag is set"
    end if

    !call boundaryconditions(xar,yinner,youter,yparam)

    do !note currently only supports INOUT

      print *, "setting boundary conditions"
      call boundaryconditions(xar,yinner,youter,yparam)

      
      print *, 'xar=',xar, " yinner=", yinner, " youter=", youter, " yparam=", yparam

      print *, 'calling shooting method'
      call shootingmethod(f,finnerboundary,xar,yinner,youter,unitinner,unitouter,stype_var,0)
      
      print *, "yinner=", yinner, " youter=", youter

     ! stop

      !potentially changing one of the read in parameters?
      !print *, xar
      !print *, yinner
      !print *, '-'
      
      !call boundaryconditions(xar,yinner,youter,yparam)
      !call shootingmethod(f,finnerboundary,xar,yinner,youter,unitinner,unitouter,stype_var,0)!this seems to break something for no sensible reason
      !potentially changing one of the read in parameters?
      !print *, xar
      !print *, yinner
      !print *, '-'

      er = yinner - youter
      print *, 'error =', er
      

      !main solver run with correct initial conditions
      if (ALL(DABS(er) .le. targeter*yinner)) return
        !call shootingmethod(f,finnerboundary,xar,yinner,youter,unitinner,unitouter,stype_flag,1)
        !return
      !end if

      !rescale=0.0

      !yparamtmp = yparam
      !yparamtmp(2) = 1.0*yparam(2) + targeter
      !call boundaryconditions(xar,yinner,youter,yparamtmp)
      !call shootingmethod(f,finnerboundary,xar,yinner,youter,unitinner,unitouter,stype_var,0)

      !er = - yinner + youter

      !yparamtmp = yparam
      !yparamtmp(2) = 1.0*yparam(2) - targeter
      !call boundaryconditions(xar,yinner,youter,yparamtmp)
      !call shootingmethod(f,finnerboundary,xar,yinner,youter,unitinner,unitouter,stype_var,0)

      !er = er + yinner - youter

      !print *, er/targeter

      !stop

      call sleep(10)
      !need to rescale variables so they are all the same magnitude, somehow.
      jacobian = 0.0
      print *, 'do'
      do i=1,SIZE(yparam)
        print *, 'parameter ',i

        yparamtmp = yparam
        yparamtmp(i) = yparam(i) + targeter
 
        print *, 'f delta yparam', yparamtmp
        print *, 'calling boundaryconditions'
        call boundaryconditions(xar,yinner,youter,yparamtmp)
        print *, 'yinner=',yinner,' youter=',youter
        print *, 'calling shootingmethod'
        call shootingmethod(f,finnerboundary,xar,yinner,youter,unitinner,unitouter,stype_var,0)
        print *, 'yinner=',yinner,' youter=',youter      

        !print *,
        !print *, youter

        jacobian(:,i) = -(- yinner + youter)      
        print *, 'jacobian 1',jacobian
        !jacobian(i,i) = - yinner(i) + youter(i)!/yinner

        yparamtmp = yparam
        yparamtmp(i) = yparam(i) - targeter
        print *, 'b delta yparam', yparamtmp
        print *, 'calling boundaryconditions'
        call boundaryconditions(xar,yinner,youter,yparamtmp)
        print *, 'yinner=',yinner,' youter=',youter
        print *, 'calling shootingmethod'
        call shootingmethod(f,finnerboundary,xar,yinner,youter,unitinner,unitouter,stype_var,0)
        print *, 'yinner=',yinner,' youter=',youter

        !print *, youter
        !print *,

        jacobian(:,i) = jacobian(:,i) - (yinner - youter)
        print *, 'jacobian 2',jacobian
        jacobian(:,i) = jacobian(:,i)/(2.0*targeter)!check the rows and columbs are the correct way round.
        print *, 'jacobian 3',jacobian

        !jacobian(i,i) = jacobian(i,i) + yinner(i) - youter(i)!/yinner
        !jacobian(i,i) = jacobian(i,i)/(2.0*targeter)

        !jacobian(i,:) = jacobian(i,:)
        call sleep(10)
      end do
 
      !issue with luminosity

      !print *, 'lum ', jacobian(:,4)

      !have to regulerise first
      !forall (i = 1:SIZE(yparam)) jacobian(i,i) = jacobian(i,i) + targeter

      print *, 'jacobian', jacobian

      print *,
      !do i=1,SIZE(yparam)
      !  print *, DABS(jacobian(i,:))
      !end do
      print *,

      !forall (i = 1:SIZE(yparam)) rescale(i,i) = 1.0D0/DABS(yparam(i) + targeter)
 
      forall (i = 1:SIZE(yparam)) jacobian(:,i) = jacobian(:,i)/MAXVAL(DABS(jacobian(:,i)) + targeter)

      if (ANY(jacobian .ne. jacobian)) then
         print *, '1/mx jacobian: ', jacobian
         stop
      end if

      !basic regulerisation sceme:

      !jacobian = MATMUL(jacobian,rescale)

      !jacobian = jacobian*1.0D6

      !print *, 'rescale', rescale
      print *, 'jacobian', jacobian

      if (ANY(jacobian .ne. jacobian)) then
         print *, 'rescale jacobian: ', jacobian
         stop
      end if


      !print *, 'jacobiani: ', jacobian
    
      !jmax = max(jacobian,1)

      !r_j_num = 0

      !regulerised_jacobian = MATMUL(jacobian,TRANSPOSE(jacobian)) !note squares the condition number so watch out for that
      !forall (i = 1:SIZE(yparam)) regulerised_jacobian(i,i) = regulerised_jacobian(i,i) + 1.0

      !forall (i = 1:SIZE(yparam)) jacobian(i,i) = jacobian(i,i) + 0.01

      regulerised_jacobian = jacobian

      print *, 'calling householderinversion'
      call householderinversion(jacobian)

      if (ANY(jacobian .ne. jacobian)) then
         print *, 'arcjacobian: ', jacobian
         stop
      end if

      !rint *, 'arcjacobian: ', jacobian
      
      !jacobian = MATMUL(jacobian,rescale)

      print *, 'ident ;', MATMUL(jacobian,regulerised_jacobian)

      print *, 'arcjacobian: ', jacobian

      print *, 'yparam: ', yparam
      print *, er
    
      !er = MATMUL(er,TRANSPOSE(jacobian))

      dyparam = -1.0*MATMUL(jacobian,er)

      print *, 'dyparam: ', dyparam
      yparam = yparam + dyparam

      !stop

      call sleep(10)

    end do

  end subroutine

  !function condition_number(G)

  !  double precision :: G(:,:)
  !  double precision :: condition_number

    

  !end function

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
