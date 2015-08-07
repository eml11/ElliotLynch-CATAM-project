
program stellarstructure

  use mod_shared
  use mod_solvers

  implicit none
  
  double precision :: xar(3)
  double precision :: yinner(4), youter(4),  yparam(4)
  integer :: unitinner = 101,unitouter = 102
  integer :: unitinput = 103

  double precision :: opacity
  double precision :: mmolecweight

  call readinputfile (unitinput,'star.inp')

  !set up problem

  mmolecweight = 1/(2.0*hmassfrac + (3.0/4.0)*(1 - hmassfrac))
  opacity = 0.02*(1 + hmassfrac)

  xar(1) = 0
  xar(2) = massm
  xar(3) = masst

  !yinner(1) = 0.0
  !youter(1) = radiusstar

  !yinner(2) = pressurec
  !youter(2) = ((2.0*G)/(3.0*0.034))*(masst/(radiusstar**2.0))*(MSUN/(RSUN**2.0)) !from eddington approx

  !yinner(3) = tempc
  !youter(3) = ((luminosity0/(4.0*PI*STBOLTZ*(radiusstar**2)))*(LSUN/(RSUN**2.0)))**(1.0/4.0)

  !yinner(4) = 0.0
  !youter(4) = luminosity0

  yparam = (/radiusstar,pressurec,tempc,luminosity0/)

  if (stype.eq.0 .or. stype.eq.1) then
    !unitinner header
    open(unit=unitinner,FILE='stellarstructure_inner.txt')

    write(unitinner,*) "! inner stellar structure"
    write(unitinner,*) "mass radius P Temp L"
    flush(unitinner) !not sure why this is suddenly necessary
  end if

  if (stype.eq.0 .or. stype.eq.2) then
    !unitouter header
    open(unit=unitouter,FILE='stellarstructure_outer.txt')

    write(unitouter,*) "! outer stellar structure"
    write(unitouter,*) "mass radius P Temp L"
    flush(unitouter)
  end if

  !run jacobian solver to find paramters
  call jacobian_solver(question5rkfunc,question5innerboundary,question5boundcond,xar,yparam,unitinner,unitouter,stype) 
 
  ! write to file using shooting solver with correct parameters 
  call boundaryconditions(xar,yinner,youter,yparam)
  call shootingmethod(question5rkfunc,question5innerboundary,xar,yinner,youter,unitinner,unitouter,stype,1)

  close(unitinner)
  close(unitouter)

  contains

  subroutine question5boundcond(xar,yinner,youter,yparam)

    double precision :: xar(3)
    double precision :: yinner(:),youter(:)
    double precision :: yparam(:)
    double precision :: masst

    masst = xar(3)

    yinner(1) = 0.0
    youter(1) = yparam(1)

    yinner(2) = yparam(2)
    youter(2) = ((2.0*G)/(3.0*0.034))*(masst/(yparam(1)**2.0))*(MSUN/(RSUN**2.0)) !from eddington approx

    yinner(3) = yparam(3)
    youter(3) = ((yparam(4)/(4.0*PI*STBOLTZ*(yparam(1)**2)))*(LSUN/(RSUN**2.0)))**(1.0/4.0)

    yinner(4) = 0.0
    youter(4) = yparam(4)

  end subroutine

  function question5rkfunc (x,y)

    double precision :: x
    double precision :: y(:)  !again possibly bad practice
    double precision, dimension(size(y)) :: question5rkfunc
    double precision :: ppCNO_energypdc, density

    !question4rkfunc(1) = (3.0/(4.0*PI))*((pressurec**(1/ratiospecificheat))/densityc)*DEXP(-y(2)/ratiospecificheat)
    !question4rkfunc(2) = - (G/(4.0*PI))*x*DEXP(-y(2))*(y(1)**(-4.0/3.0))
    
    density = ((mmolecweight*y(2))/(RGAS*y(3)))

    ppCNO_energypdc = ((0.25*(hmassfrac**2)*DEXP(-33.8 - 33.8*(y(3)/1.0D6)**(-1.0/3.0)) + & 
   & 8.8D18*DEXP(-152.2*(y(3)/1.0D6)**(-1.0/3.0)))*((y(3)/1.0D6)**(-2.0/3.0)))*density

    question5rkfunc(1) = (1/(4.0*PI))*(1/(y(1)**2))*(1/density) * &
   & (MSUN/(RSUN**3))
    question5rkfunc(2) = -(G/(4.0*PI))*(x/(y(1)**4))*((MSUN**2)/(RSUN**4))
    question5rkfunc(3) = -((3*opacity*y(4))/(((16*PI)**2)*STBOLTZ*(y(1)**4)*(y(3)**3)))*((MSUN*LSUN)/(RSUN**4))
    question5rkfunc(4) = ppCNO_energypdc*(MSUN/LSUN)

  end function

  function question5innerboundary (x,y) !this is inconsistent

    double precision :: x
    double precision :: y(:)  !again possibly bad practice
    double precision, dimension(size(y)) :: question5innerboundary
    double precision :: ppCNO_energypdc, density

    density = (mmolecweight*y(2))/(RGAS*y(3))

    ppCNO_energypdc = ((0.25*(hmassfrac**2)*DEXP(-33.8 - 33.8*(y(3)/1.0D6)**(-1.0/3.0)) + & 
   & 8.8D18*DEXP(-152.2*(y(3)/1.0D6)**(-1.0/3.0)))*((y(3)/1.0D6)**(-2.0/3.0)))*density*((y(3)/1.0D6)**(-2.0/3.0))

    question5innerboundary(1) = ((3.0/(4.0*PI))*((RGAS*tempc)/(mmolecweight*pressurec))*x)**(1.0/3.0)
    question5innerboundary(2) = pressurec - (2.0*PI/3.0)*G*(y(1)**2)*((mmolecweight*pressurec)/(RGAS*tempc))**2
    question5innerboundary(3) = ((4.0*PI/3.0)*((mmolecweight*pressurec)/((RGAS*tempc)))**(1.0/3.0)) * &
    & ((((12.0*opacity)/(STBOLTZ*(16*PI)**2)) * &
    & (3.0*luminosity0*LSUN*((x*MSUN)**(-1.0/3.0)) - 1.5*ppCNO_energypdc*(2.0/3.0)))**(0.25)) * MSUN
    question5innerboundary(4) = ppCNO_energypdc*x*(MSUN/LSUN)

  end function

end program

