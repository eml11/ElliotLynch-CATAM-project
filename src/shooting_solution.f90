
program shootingsolution

  use mod_shared
  use mod_solvers

  implicit none

  double precision :: xar(3)
  double precision :: yinner(4), youter(4), yparam(4)
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

  yparam = (/radiusstar,pressurec,tempc,luminosity0/)

  call question5boundcond(xar,yinner,youter,yparam)

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

  !run solver
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
    youter(2) = (1.0/Gbar)*((2.0*G)/(3.0*opacity))*(masst/(yparam(1)**2.0))*(MSUN/(RSUN**2.0)) !from eddington approx

    yinner(3) = yparam(3)
    youter(3) = (1.0/MKELVIN)*(((yparam(4)/(4.0*PI*STBOLTZ*(yparam(1)**2)))*(LSUN/(RSUN**2.0)))**(1.0/4.0))

    yinner(4) = 0.0
    youter(4) = yparam(4)

  end subroutine

  function question5rkfunc (x,y)

    double precision :: x
    double precision :: y(:)  !again possibly bad practice
    double precision, dimension(size(y)) :: question5rkfunc
    double precision :: ppCNO_energypdc, density


    density = ((mmolecweight*y(2)*Gbar)/(RGAS*y(3)*MKELVIN))

    ppCNO_energypdc = (0.25*(hmassfrac**2)*DEXP(-33.8*(y(3)**(-1.0/3.0))) + &
   & 8.8D18*hmassfrac*DEXP(-152.28*(y(3)**(-1.0/3.0)))) * density * &
   & (y(3)**(-2.0/3.0)) 

    question5rkfunc(1) = (1/(4.0*PI))*(1/(y(1)**2))*(1/density) * &
   & (MSUN/(RSUN**3))
    question5rkfunc(2) = -(G/(4.0*PI))*(x/(y(1)**4))*((MSUN**2)/(Gbar*(RSUN**4)))
    question5rkfunc(3) = -((3*opacity*y(4))/(((16*PI)**2)*STBOLTZ*(y(1)**4)*((y(3)*MKELVIN)**3)))*((MSUN*LSUN)/(MKELVIN*(RSUN**4)))
    question5rkfunc(4) = ppCNO_energypdc*(MSUN/LSUN)

  end function

  function question5innerboundary (x,y) !this is inconsistent

    double precision :: x
    double precision :: y(:)  !again possibly bad practice
    double precision, dimension(size(y)) :: question5innerboundary
    double precision :: ppCNO_energypdc, density

    density = (mmolecweight*y(2))*Gbar/(RGAS*y(3)*MKELVIN)
   
    ppCNO_energypdc = (0.25*(hmassfrac**2)*DEXP(-33.8*(y(3)**(-1.0/3.0))) + &
   & 8.8D18*hmassfrac*DEXP(-152.28*(y(3)**(-1.0/3.0))) )* density * &
   & (y(3)**(-2.0/3.0))


    !hacky way of ensuring r is sufficiently large
    question5innerboundary(1) = 1.0D-4

    x = ((4.0*PI/3.0)*density*question5innerboundary(1)**3)

    question5innerboundary(2) = pressurec - (1.0/Gbar)*(2.0*PI/3.0)*G*(y(1)**2)*((mmolecweight*pressurec*Gbar)/(RGAS*tempc*MKELVIN))**2
   
    question5innerboundary(3) = (((tempc*MKELVIN)**4 - ((3.0*opacity)/(32.0*PI*STBOLTZ))*(1/density) * &
   & ppCNO_energypdc*((x*MSUN)**(2.0/3.0)))**(0.25))/MKELVIN

    question5innerboundary(4) = ppCNO_energypdc*x*(MSUN/LSUN)

   ! print *, question5innerboundary

  end function

end program

