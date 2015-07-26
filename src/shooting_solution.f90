
program stellarstructure

  use mod_shared

  implicit none

  double precision :: xar(3)
  double precision :: yinner(2), youter(2)
  integer :: unitinner = 101,unitouter = 102
  integer :: unitinput = 103

  call readinputfile (unitinput,'star.inp')

  open(unitinner,FILE='stellarstructure_inner.txt')
  open(unitouter,FILE='stellarstructure_outer.txt')

  !set up problem
  xar(1) = 0
  xar(2) = massm
  xar(3) = masst

  yinner(1) = 0.0
  youter(1) = radiusstar**3.

  yinner(2) = DLOG(pressurec)
  youter(2) = -10.0

  !unitinner header
  write(unitinner,*) '! inner stellar structure' 
  write(unitinner,*) 'mass radius^3 lnP'

  !unitouter header
  write(unitinner,*) '! outer stellar structure'
  write(unitinner,*) 'mass radius^3 lnP'

  !run solver
  call shootingmethod(question4rkfunc,xar,yinner,youter,unitinner,unitouter) 

  close(unitinner)
  close(unitouter)

  contains

  function question4rkfunc (x,y)

    double precision :: x
    double precision :: y(2)
    double precision :: question4rkfunc(2)

    question4rkfunc(1) = 12*PI*(y(1)**(4.0/3.0))*densityc*DEXP((y(2) - y2c)/ratiospecificheat)
    question4rkfunc(2) = - ((G*x)/(4*PI))*(y(1)**(4.0/3.0))*DEXP(y(2))

  end function

end program

