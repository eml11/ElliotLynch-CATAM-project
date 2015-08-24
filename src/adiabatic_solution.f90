
program adiabatic_solver

  use mod_shared
  use mod_solvers

  implicit none
  
  double precision :: xar(3)
  double precision :: yinner(2), youter(2)
  integer :: unitinner = 101,unitouter = 102
  integer :: unitinput = 103

  double precision :: opacity = 0.034
  double precision :: mmolecweight = 1.0/1.625

  call readinputfile (unitinput,'star.inp')

  open(unit=unitinner,FILE='stellarstructure_inner.txt')
  open(unit=unitouter,FILE='stellarstructure_outer.txt')

  !composition dependent parameters
  mmolecweight = 1/(2.0*hmassfrac + (3.0/4.0)*(1 - hmassfrac))
  opacity = 0.02*(1 + hmassfrac)
  
  !set up problem
  xar(1) = 0
  xar(2) = massm
  xar(3) = masst

  yinner(1) = 0.0
  youter(1) = radiusstar

  yinner(2) = pressurec
  youter(2) = (1.0/Gbar)*(((2.0*G)/(3.0*opacity))*(masst/(radiusstar**2.0))*(MSUN/(RSUN**2.0))) !from eddington approx

  !unitinner header
  write(unitinner,*) "! inner stellar structure"
  write(unitinner,*) "Mass Radius P"
  flush(unitinner) !not sure why this is suddenly necessary


  !unitouter header
  write(unitouter,*) "! outer stellar structure"
  write(unitouter,*) "Mass Radius P"
  flush(unitouter)

  !run solver
  call shootingmethod(question4rkfunc,question4innerboundary,xar,yinner,youter,unitinner,unitouter,stype,1) 

  close(unitinner)
  close(unitouter)

  contains

  function question4rkfunc (x,y)

    double precision :: x
    double precision :: y(:)  !again possibly bad practice
    double precision, dimension(size(y)) :: question4rkfunc
    double precision :: density 

    !question4rkfunc(1) = (3.0/(4.0*PI))*((pressurec**(1/ratiospecificheat))/densityc)*DEXP(-y(2)/ratiospecificheat)
    !question4rkfunc(2) = - (G/(4.0*PI))*x*DEXP(-y(2))*(y(1)**(-4.0/3.0))

    density = ((mmolecweight*pressurec*Gbar)/(RGAS*tempc*MKELVIN))*((y(2)/pressurec)**(1.0/ratiospecificheat))

    question4rkfunc(1) = (1.0/(4.0*PI*density))*(1.0/(y(1)**2))*(MSUN/(RSUN**3))
    question4rkfunc(2) = - (G/(4.0*PI))*(x/(y(1)**4))*((MSUN**2)/(Gbar*RSUN**4))

  end function

  function question4innerboundary (x,y)

    double precision :: x
    double precision :: y(:)  !again possibly bad practice
    double precision, dimension(size(y)) :: question4innerboundary
    double precision :: density
   
    density = ((mmolecweight*pressurec*Gbar)/(RGAS*tempc*MKELVIN))

    question4innerboundary(1) = 1.0D-4

    x = (1.0/MSUN)*((4.0*PI/3.0)*density*(question4innerboundary(1)*RSUN)**3)

    !question4innerboundary(1) = (1.0/RSUN)*(((3.0/(4.0*PI*density))*x*MSUN)**(1.0/3.0))
    question4innerboundary(2) = pressurec - (1.0/Gbar)*(2.0*PI/3.0)*G*((question4innerboundary(1)*RSUN)**2)*density**2

  end function

end program

