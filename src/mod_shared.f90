!
! MODULE: mod_shared 
!
!> Contains common varibles. Physical Constants, 
!! Scaling parameters and variables for file input
module mod_shared

  double precision :: massm
  double precision :: masst
  double precision :: radiusstar
  double precision :: pressurec
  double precision :: ratiospecificheat
  double precision :: y2c   
  double precision :: densityc
  double precision :: tempc
  double precision :: hmassfrac
  double precision :: luminosity0
  integer :: stype = 0

  double precision, parameter :: G = 6.674D-11
  double precision, parameter :: PI = 4*DATAN(1.D0)
  double precision, parameter :: STBOLTZ = 5.67D-8 
  double precision, parameter :: RGAS = 8.3145D3
  double precision, parameter :: MSUN = 1.9891D30
  double precision, parameter :: RSUN = 6.9598D8
  double precision, parameter :: LSUN = 3.8515D26
  double precision, parameter :: PSUN = ((2*G)/(3*0.034))*(MSUN/(RSUN**2))
  double precision, parameter :: DENSUN = (0.61538*PSUN/RGAS)*((LSUN/(4*PI*RSUN*RSUN*STBOLTZ))**(1.0/4.0))
  double precision, parameter :: Gbar = 1.0D14
  double precision, parameter :: MKELVIN = 1.0D6

end module
