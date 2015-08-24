!
! PROGRAM: stellarstructure 
!
!> Solves the full equations of stellar structure 
!! using the jacobian interation method
program stellarstructure

  use mod_shared
  use mod_solvers

  implicit none
  
  double precision :: xar(3)
  double precision :: yinner(4), youter(4),  yparam(4), er(4)
  integer :: unitinner = 101,unitouter = 102,unitparameters=103
  integer :: unitinput = 104

  double precision :: opacity
  double precision :: mmolecweight

  call readinputfile (unitinput,'star.inp')

  !compute composition dependent variables
  mmolecweight = 1/(2.0*hmassfrac + (3.0/4.0)*(1 - hmassfrac))
  opacity = 0.02*(1 + hmassfrac)

  !set up problem
  xar(1) = 0
  xar(2) = massm
  xar(3) = masst

  yparam = (/radiusstar,pressurec,tempc,luminosity0/)

  !write approriate ouput headers
  if (stype.eq.0 .or. stype.eq.1) then
    !unitinner header
    open(unit=unitinner,FILE='stellarstructure_inner.txt')

    write(unitinner,*) "! inner stellar structure"
    write(unitinner,*) "Mass Radius P Temp L"
    flush(unitinner) !not sure why this is suddenly necessary
  end if

  if (stype.eq.0 .or. stype.eq.2) then
    !unitouter header
    open(unit=unitouter,FILE='stellarstructure_outer.txt')

    write(unitouter,*) "! outer stellar structure"
    write(unitouter,*) "Mass Radius P Temp L"
    flush(unitouter)
  end if

  !run jacobian solver to find paramters
  call jacobian_solver(question5rkfunc,question5innerboundary,question5boundcond,xar,yparam,unitinner,unitouter,stype) 
 
  ! write to file using shooting solver with correct parameters 
  call question5boundcond(xar,yinner,youter,yparam)
  call shootingmethod(question5rkfunc,question5innerboundary,xar,yinner,youter,unitinner,unitouter,stype,1)

  er = 2.0*((yinner - youter)/(yinner + youter))

  !write parameter output
  open(unit=unitparameters,FILE='parameters.txt')

  write(unitparameters,*) '! stellar parameters and errors'
  write(unitparameters,*)

  write(unitparameters,*) 'Parameter Value IMidpoint OMidpoint Error'

  write(unitparameters,*) 'Radius', yparam(1), yinner(1), youter(1), er(1)
  write(unitparameters,*) 'Pressure', yparam(2), yinner(2), youter(2), er(2)
  write(unitparameters,*) 'Temperature', yparam(3), yinner(3), youter(3), er(3)
  write(unitparameters,*) 'Luminosity', yparam(4), yinner(4), youter(4), er(4)

  !close output files
  close(unitinner)
  close(unitouter)
  close(unitparameters)

  contains

  !> Sets up the initial conditions for the stellar variables
  !! calculated from the parameters R*,Pc,Tc,L*
  !! @param xar array of mass bounds for the integration 
  !! @param yinner initial values for the inner integral
  !! @param youter initial values for the outer integral
  !! @param yparam stellar parameters R*,Pc,Tc,L*
  subroutine question5boundcond(xar,yinner,youter,yparam)

    double precision :: xar(3)
    double precision :: yinner(:),youter(:)
    double precision :: yparam(:)
    double precision :: masst

    masst = xar(3)

    yinner(1) = 0.0
    youter(1) = yparam(1)

    yinner(2) = yparam(2)
    youter(2) = (1.0/Gbar)*((2.0*G)/(3.0*opacity))*(masst/(yparam(1)**2))*(MSUN/(RSUN**2)) !from eddington approx

    yinner(3) = yparam(3)
    youter(3) = (1.0/MKELVIN)*(((yparam(4)/(4.0*PI*STBOLTZ*(yparam(1)**2)))*(LSUN/(RSUN**2)))**(1.0/4.0))

    yinner(4) = 0.0
    youter(4) = yparam(4)

  end subroutine

  !> Derivatives of the stellar variables with respect to mass.
  !! Used in Runge-Kutta integration
  !! @param x mass indepenent variable
  !! @parma y array of dependent variables
  !! @return question5rkfunc array of derivatives of y wrt x
  function question5rkfunc (x,y)

    double precision :: x
    double precision :: y(:)
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

  !> Analytic solutions to the linearised equation for the stellar variables
  !! near the core. Gives stellar variables at R=1.0D-5 
  !! @param x mass indepenent variable
  !! @param y array of dependent variables
  !! @return question4innerboundary array of dependent variables at R=1.0D-5
  function question5innerboundary (x,y) 

    double precision :: x
    double precision :: y(:)
    double precision, dimension(size(y)) :: question5innerboundary
    double precision :: ppCNO_energypdc, density

    density = (mmolecweight*y(2)*Gbar)/(RGAS*y(3)*MKELVIN)

    ppCNO_energypdc = (0.25*(hmassfrac**2)*DEXP(-33.8*(y(3)**(-1.0/3.0))) + &
   & 8.8D18*hmassfrac*DEXP(-152.28*(y(3)**(-1.0/3.0))) )* density * &
   & (y(3)**(-2.0/3.0))

    question5innerboundary(1) = 1.0D-5

    x = (1.0/MSUN)*((4.0*PI/3.0)*density*(question5innerboundary(1)*RSUN)**3)

    question5innerboundary(2) = y(2) - &
   & (1.0/Gbar)*(2.0*PI/3.0)*G*((y(1)*RSUN)**2)*(density**2)
   
    question5innerboundary(3) = (((y(3)*MKELVIN)**4 - ((3.0*opacity)/(32.0*PI*STBOLTZ))*(1/density) * &
   & ppCNO_energypdc*((x*MSUN)**(2.0/3.0)))**(0.25))/MKELVIN
    question5innerboundary(4) = ppCNO_energypdc*x*(MSUN/LSUN)

  end function

end program

