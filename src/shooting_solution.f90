
program stellarstructure

  use mod_shared

  implicit none

  interface
    subroutine shootingmethod(f,xar,yinner,youter,unitinner,unitouter)
      interface
        function f(xin,yout)
          double precision :: xin
          double precision :: yout(:)
          double precision, dimension(lbound(yout,dim=1)) :: f
        end function
      end interface
      double precision :: yinner(:),youter(:)
      double precision :: xar(3)
      double precision :: dxinner, dxouter
      integer :: unitinner,unitouter
    end subroutine
  end interface

  double precision :: xar(3)
  double precision :: yinner(2), youter(2)
  integer :: unitinner = 101,unitouter = 102
  integer :: unitinput = 103

  call readinputfile (unitinput,'star.inp')

  open(unit=unitinner,FILE='stellarstructure_inner.txt')
  open(unit=unitouter,FILE='stellarstructure_outer.txt')

  !set up problem
  xar(1) = 0
  xar(2) = massm
  xar(3) = masst

  yinner(1) = 0.0
  youter(1) = radiusstar

  yinner(2) = pressurec
  youter(2) = ((2.0*G)/(3*0.034))*(MSUN/(pressurec*RSUN**2))*(masst/(radiusstar**2)) !from eddington approx

  !unitinner header
  write(unitinner,*) "! inner stellar structure"
  write(unitinner,*) "mass radius^3 lnP"
  flush(unitinner) !not sure why this is suddenly necessary


  !unitouter header
  write(unitouter,*) "! outer stellar structure"
  write(unitouter,*) "mass radius^3 lnP"
  flush(unitouter)

  !run solver
  call shootingmethod(question4rkfunc,xar,yinner,youter,unitinner,unitouter) 

  close(unitinner)
  close(unitouter)

  contains

  function question4rkfunc (x,y)

    double precision :: x
    double precision :: y(:)  !again possibly bad practice 
    double precision, dimension(lbound(y,dim=1)) :: question4rkfunc

    !question4rkfunc(1) = (3.0/(4.0*PI))*((pressurec**(1/ratiospecificheat))/densityc)*DEXP(-y(2)/ratiospecificheat)
    !question4rkfunc(2) = - (G/(4.0*PI))*x*(y(1)**(-4.0/3.0))*DEXP(-y(2))

    question4rkfunc(1) = (1/(4.0*PI))*(1/(y(1)**2))*(1/densityc)*((y(2))**(-1/ratiospecificheat))
    question4rkfunc(2) = -(G/(4.0*PI))*(x/(y(1)**4))*((MSUN**2)/(pressurec*RSUN**4))

  end function

end program

