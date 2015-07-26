module mod_shared

  double precision :: massm
  double precision :: masst
  double precision :: radiusstar
  double precision :: pressurec
  double precision :: ratiospecificheat
  double precision :: y2c   
  double precision :: densityc

  double precision, parameter :: G = 6.674D-11
  double precision, parameter :: PI = 4*DATAN(1.D0) 
  double precision, parameter :: MSUN = 1.9891D30
  double precision, parameter :: RSUN = 6.9598D8

end module
